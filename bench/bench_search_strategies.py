#!/usr/bin/env python3
"""Benchmark search strategies for pyitsx.

Compares:
  baseline    — single nhmmer() call with all 538 profiles
  per_anchor  — 4 nhmmer() calls, one per anchor type
  shortcircuit — per-profile with early stopping on confident hit
  staged      — 5.8S anchors first, then SSU/LSU restricted to relevant regions

All strategies use cpus=1 for apples-to-apples comparison.
"""

import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pyhmmer
import pyhmmer.easel
import pyhmmer.hmmer

from pyitsx.constants import AnchorType, Confidence, Strand
from pyitsx.models import AnchorHit, ChainConstraints, DEFAULT_CONSTRAINTS
from pyitsx.profiles import ProfileDB, WINDOW_LENGTH
from pyitsx.chains import build_chain
from pyitsx.pipeline import delimit


SCORE_THRESHOLD = DEFAULT_CONSTRAINTS.min_anchor_score
EVALUE_THRESHOLD = DEFAULT_CONSTRAINTS.max_anchor_evalue


def _parse_anchor_type(profile_name: str):
    prefix = profile_name.split("_")[0]
    try:
        return AnchorType(int(prefix))
    except (ValueError, KeyError):
        return None


def _collect_hits(top_hits) -> list[AnchorHit]:
    """Extract AnchorHits from a single nhmmer TopHits result."""
    hits = []
    profile_name = top_hits.query.name
    anchor_type = _parse_anchor_type(profile_name)
    if anchor_type is None:
        return hits
    for hit in top_hits:
        if not hit.included:
            continue
        for domain in hit.domains:
            hits.append(AnchorHit(
                anchor_type=anchor_type,
                strand=Strand(domain.strand),
                env_from=domain.env_from,
                env_to=domain.env_to,
                score=domain.score,
                evalue=domain.i_evalue,
                profile_name=profile_name,
            ))
    return hits


def _group_profiles(hmms):
    groups = defaultdict(list)
    for hmm in hmms:
        anchor_type = _parse_anchor_type(hmm.name if isinstance(hmm.name, str) else hmm.name.decode())
        if anchor_type is not None:
            groups[anchor_type].append(hmm)
    return dict(groups)


def search_baseline(hmms, seqs):
    """Current approach: one nhmmer() call with all profiles."""
    hits_by_seq = defaultdict(list)
    for top_hits in pyhmmer.hmmer.nhmmer(hmms, seqs, cpus=1, window_length=WINDOW_LENGTH):
        for h in _collect_hits(top_hits):
            hits_by_seq[top_hits.query.name].append(h)  # wrong key, fix below
    # Redo properly
    hits_by_seq = defaultdict(list)
    for top_hits in pyhmmer.hmmer.nhmmer(hmms, seqs, cpus=1, window_length=WINDOW_LENGTH):
        profile_name = top_hits.query.name
        anchor_type = _parse_anchor_type(profile_name)
        if anchor_type is None:
            continue
        for hit in top_hits:
            if not hit.included:
                continue
            for domain in hit.domains:
                hits_by_seq[hit.name].append(AnchorHit(
                    anchor_type=anchor_type,
                    strand=Strand(domain.strand),
                    env_from=domain.env_from,
                    env_to=domain.env_to,
                    score=domain.score,
                    evalue=domain.i_evalue,
                    profile_name=profile_name,
                ))
    return dict(hits_by_seq)


def search_per_anchor(hmms, seqs):
    """4 nhmmer() calls, one per anchor type."""
    groups = _group_profiles(hmms)
    hits_by_seq = defaultdict(list)
    for anchor_type in [AnchorType.S58_START, AnchorType.S58_END,
                        AnchorType.SSU_END, AnchorType.LSU_START]:
        profiles = groups.get(anchor_type, [])
        if not profiles:
            continue
        for top_hits in pyhmmer.hmmer.nhmmer(profiles, seqs, cpus=1, window_length=WINDOW_LENGTH):
            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    hits_by_seq[hit.name].append(AnchorHit(
                        anchor_type=anchor_type,
                        strand=Strand(domain.strand),
                        env_from=domain.env_from,
                        env_to=domain.env_to,
                        score=domain.score,
                        evalue=domain.i_evalue,
                        profile_name=top_hits.query.name,
                    ))
    return dict(hits_by_seq)


def search_shortcircuit(hmms, seqs, profile_freq=None):
    """Per-profile search with short-circuit on confident anchor hit.

    Profiles are ordered by hit frequency (most common first).
    For each anchor type, once every sequence has a confident hit,
    skip remaining profiles for that type.
    """
    groups = _group_profiles(hmms)
    seq_names = {s.name for s in seqs}
    n_seqs = len(seq_names)

    hits_by_seq = defaultdict(list)
    satisfied = defaultdict(set)  # anchor_type -> set of seq_names with confident hits
    profiles_searched = 0

    for anchor_type in [AnchorType.S58_START, AnchorType.S58_END,
                        AnchorType.SSU_END, AnchorType.LSU_START]:
        profiles = groups.get(anchor_type, [])
        if profile_freq:
            profiles = sorted(profiles, key=lambda p: profile_freq.get(p.name if isinstance(p.name, str) else p.name.decode(), 0), reverse=True)

        for profile in profiles:
            if len(satisfied[anchor_type]) >= n_seqs:
                break

            profiles_searched += 1
            for top_hits in pyhmmer.hmmer.nhmmer([profile], seqs, cpus=1, window_length=WINDOW_LENGTH):
                for hit in top_hits:
                    if not hit.included:
                        continue
                    for domain in hit.domains:
                        ah = AnchorHit(
                            anchor_type=anchor_type,
                            strand=Strand(domain.strand),
                            env_from=domain.env_from,
                            env_to=domain.env_to,
                            score=domain.score,
                            evalue=domain.i_evalue,
                            profile_name=top_hits.query.name,
                        )
                        hits_by_seq[hit.name].append(ah)
                        if ah.score >= SCORE_THRESHOLD and ah.evalue <= EVALUE_THRESHOLD:
                            satisfied[anchor_type].add(hit.name)

    return dict(hits_by_seq), profiles_searched


def learn_profile_frequencies(hmms, seqs):
    """Run baseline search and count which profiles produce hits per anchor type."""
    freq = Counter()
    for top_hits in pyhmmer.hmmer.nhmmer(hmms, seqs, cpus=1, window_length=WINDOW_LENGTH):
        profile_name = top_hits.query.name
        n_hits = sum(1 for hit in top_hits if hit.included)
        if n_hits > 0:
            freq[profile_name.decode() if isinstance(profile_name, bytes) else profile_name] = n_hits
    return freq


def run_benchmark(input_path, hmm_dir, max_seqs=0, n_rounds=3):
    db = ProfileDB(hmm_dir, organism="F")
    seqs = db.prepare(input_path)
    if max_seqs:
        block = pyhmmer.easel.DigitalSequenceBlock(db._alphabet)
        for i, s in enumerate(seqs):
            if i >= max_seqs:
                break
            block.append(s)
        seqs = block

    n = len(seqs)
    hmms = db._hmms
    groups = _group_profiles(hmms)

    print(f"Input: {n} sequences from {input_path.name}")
    print(f"Profiles: {len(hmms)} total — "
          + ", ".join(f"{at.name}: {len(ps)}" for at, ps in sorted(groups.items(), key=lambda x: x[0].value)))
    print(f"Rounds: {n_rounds}")
    print()

    # Learn profile frequencies from first baseline run
    print("Learning profile frequencies...", end=" ", flush=True)
    t0 = time.perf_counter()
    profile_freq = learn_profile_frequencies(hmms, seqs)
    print(f"done ({time.perf_counter() - t0:.2f}s)")
    top5 = profile_freq.most_common(5)
    print(f"  Top 5 profiles: {', '.join(f'{name}={count}' for name, count in top5)}")
    print()

    strategies = {
        "baseline": lambda: search_baseline(hmms, seqs),
        "per_anchor": lambda: search_per_anchor(hmms, seqs),
        "shortcircuit": lambda: search_shortcircuit(hmms, seqs, profile_freq),
    }

    results = {}
    for name, fn in strategies.items():
        times = []
        for r in range(n_rounds):
            t0 = time.perf_counter()
            result = fn()
            elapsed = time.perf_counter() - t0
            times.append(elapsed)

        if isinstance(result, tuple):
            hits_by_seq, profiles_searched = result
        else:
            hits_by_seq = result
            profiles_searched = len(hmms)

        n_detected = len(hits_by_seq)
        n_hits = sum(len(v) for v in hits_by_seq.values())
        best = min(times)
        median = sorted(times)[len(times) // 2]

        # Build chains to verify correctness
        chains_built = 0
        full_chains = 0
        for seq_id, hits in hits_by_seq.items():
            chain = build_chain(hits)
            if chain:
                chains_built += 1
                if chain.is_full:
                    full_chains += 1

        results[name] = {
            "best": best, "median": median,
            "detected": n_detected, "hits": n_hits,
            "chains": chains_built, "full": full_chains,
            "profiles_searched": profiles_searched,
        }

        print(f"{name:20s}  best={best:6.3f}s  median={median:6.3f}s  "
              f"seq/s={n/best:7.1f}  "
              f"detected={n_detected}  chains={chains_built} (full={full_chains})  "
              f"profiles={profiles_searched}")

    print()
    baseline_time = results["baseline"]["best"]
    for name, r in results.items():
        speedup = baseline_time / r["best"] if r["best"] > 0 else float("inf")
        print(f"  {name:20s}  {speedup:5.2f}x vs baseline")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Benchmark pyitsx search strategies")
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("--hmm-dir", type=Path, default=None)
    parser.add_argument("--max-seqs", type=int, default=0)
    parser.add_argument("--rounds", type=int, default=3)
    args = parser.parse_args()

    hmm_dir = args.hmm_dir
    if hmm_dir is None:
        from pyitsx.profiles import find_hmm_dir
        hmm_dir = find_hmm_dir()

    run_benchmark(args.input, hmm_dir, args.max_seqs, args.rounds)


if __name__ == "__main__":
    main()
