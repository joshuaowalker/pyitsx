#!/usr/bin/env python3
"""Empirical analysis of HMM profile hit distributions.

Runs all profiles against the full dataset (no short-circuit) and analyzes:
1. Profile frequency — which profiles dominate per anchor type?
2. Score/E-value distributions — what's a good short-circuit threshold?
3. Position agreement — do top profiles agree on boundary positions?
4. Recommended static profile ordering for cold-start optimization.
"""

import json
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import pyhmmer.easel
import pyhmmer.hmmer
import pyhmmer.plan7

from pyitsx.constants import AnchorType
from pyitsx.profiles import WINDOW_LENGTH, _NHMMER_Z, _parse_anchor_type, find_hmm_dir


@dataclass
class RawHit:
    seq_id: str
    anchor_type: AnchorType
    profile_name: str
    score: float
    evalue: float
    env_from: int
    env_to: int
    strand: str


def collect_all_hits(hmm_path, fasta_path, cpus=4):
    """Run all profiles against all sequences, return every included hit."""
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as f:
        hmms = list(f)

    alphabet = pyhmmer.easel.Alphabet.dna()
    with pyhmmer.easel.SequenceFile(str(fasta_path), digital=True, alphabet=alphabet) as sf:
        seqs = sf.read_block()

    n_seqs = len(seqs)
    n_profiles = len(hmms)
    print(f"Scanning {n_seqs:,} sequences × {n_profiles} profiles (cpus={cpus})...",
          file=sys.stderr, flush=True)

    hits = []
    t0 = time.perf_counter()
    for i, top_hits in enumerate(pyhmmer.hmmer.nhmmer(
        hmms, seqs, cpus=cpus, window_length=WINDOW_LENGTH, Z=_NHMMER_Z
    )):
        profile_name = top_hits.query.name
        anchor_type = _parse_anchor_type(profile_name)
        if anchor_type is None:
            continue
        for hit in top_hits:
            if not hit.included:
                continue
            for domain in hit.domains:
                hits.append(RawHit(
                    seq_id=hit.name,
                    anchor_type=anchor_type,
                    profile_name=profile_name,
                    score=domain.score,
                    evalue=domain.i_evalue,
                    env_from=domain.env_from,
                    env_to=domain.env_to,
                    strand="+" if domain.strand == "watson" else "-",
                ))
        if (i + 1) % 100 == 0:
            elapsed = time.perf_counter() - t0
            print(f"  {i+1}/{n_profiles} profiles done ({elapsed:.0f}s)",
                  file=sys.stderr, flush=True)

    elapsed = time.perf_counter() - t0
    print(f"Done: {len(hits):,} hits in {elapsed:.1f}s", file=sys.stderr, flush=True)
    return hits, n_seqs


def analyze_profile_frequency(hits):
    """Which profiles produce hits, and how often?"""
    print("\n" + "=" * 80)
    print("1. PROFILE FREQUENCY (how many sequences each profile matches)")
    print("=" * 80)

    for anchor_type in AnchorType:
        anchor_hits = [h for h in hits if h.anchor_type == anchor_type]
        if not anchor_hits:
            continue

        # Count unique sequences per profile
        profile_seqs = defaultdict(set)
        for h in anchor_hits:
            profile_seqs[h.profile_name].add(h.seq_id)

        total_seqs_with_hits = len({h.seq_id for h in anchor_hits})
        print(f"\n--- {anchor_type.name} ({total_seqs_with_hits:,} seqs with any hit) ---")
        print(f"  {'Profile':<55s}  {'Seqs':>7s}  {'%':>6s}  {'MeanScore':>9s}")

        ranked = sorted(profile_seqs.items(), key=lambda x: len(x[1]), reverse=True)
        for profile_name, seq_ids in ranked[:20]:
            scores = [h.score for h in anchor_hits if h.profile_name == profile_name]
            mean_score = sum(scores) / len(scores)
            pct = 100 * len(seq_ids) / total_seqs_with_hits
            print(f"  {profile_name:<55s}  {len(seq_ids):>7,}  {pct:>5.1f}%  {mean_score:>9.1f}")

        n_with_hits = len([p for p, s in profile_seqs.items() if len(s) > 0])
        n_dominant = len([p for p, s in profile_seqs.items() if len(s) >= total_seqs_with_hits * 0.01])
        print(f"  Total profiles with any hit: {n_with_hits}/{len(ranked)}")
        print(f"  Profiles matching >=1% of seqs: {n_dominant}")


def analyze_score_distributions(hits):
    """Score distributions to inform short-circuit threshold."""
    print("\n" + "=" * 80)
    print("2. SCORE DISTRIBUTIONS (per anchor type)")
    print("=" * 80)

    for anchor_type in AnchorType:
        anchor_hits = [h for h in hits if h.anchor_type == anchor_type]
        if not anchor_hits:
            continue

        scores = sorted([h.score for h in anchor_hits])
        n = len(scores)
        print(f"\n--- {anchor_type.name} ({n:,} hits) ---")
        print(f"  min={scores[0]:.1f}  p5={scores[int(n*0.05)]:.1f}  "
              f"p25={scores[int(n*0.25)]:.1f}  median={scores[int(n*0.5)]:.1f}  "
              f"p75={scores[int(n*0.75)]:.1f}  p95={scores[int(n*0.95)]:.1f}  "
              f"max={scores[-1]:.1f}")

        # Per-sequence best score distribution
        best_by_seq = {}
        for h in anchor_hits:
            key = h.seq_id
            if key not in best_by_seq or h.score > best_by_seq[key]:
                best_by_seq[key] = h.score
        best_scores = sorted(best_by_seq.values())
        n_best = len(best_scores)
        print(f"  Best-per-seq ({n_best:,} seqs): "
              f"min={best_scores[0]:.1f}  p5={best_scores[int(n_best*0.05)]:.1f}  "
              f"median={best_scores[int(n_best*0.5)]:.1f}  "
              f"p95={best_scores[int(n_best*0.95)]:.1f}  max={best_scores[-1]:.1f}")

        # How many seqs have a hit >= various thresholds?
        for threshold in [10, 15, 20, 25, 30]:
            n_above = sum(1 for s in best_scores if s >= threshold)
            print(f"  Seqs with best score >= {threshold}: {n_above:,} ({100*n_above/n_best:.1f}%)")


def analyze_position_agreement(hits):
    """Do different profiles for the same anchor agree on boundary position?"""
    print("\n" + "=" * 80)
    print("3. POSITION AGREEMENT (boundary consistency across profiles)")
    print("=" * 80)

    for anchor_type in AnchorType:
        anchor_hits = [h for h in hits if h.anchor_type == anchor_type and h.strand == "+"]
        if not anchor_hits:
            continue

        # Group by sequence, find best and second-best profile
        by_seq = defaultdict(list)
        for h in anchor_hits:
            by_seq[h.seq_id].append(h)

        diffs = []
        n_multi = 0
        for seq_id, seq_hits in by_seq.items():
            if len(seq_hits) < 2:
                continue
            n_multi += 1
            seq_hits.sort(key=lambda h: h.score, reverse=True)
            best = seq_hits[0]
            second = seq_hits[1]
            # Compare env_from and env_to positions
            from_diff = abs(best.env_from - second.env_from)
            to_diff = abs(best.env_to - second.env_to)
            diffs.append((max(from_diff, to_diff), from_diff, to_diff,
                         best.score, second.score, best.profile_name, second.profile_name))

        if not diffs:
            print(f"\n--- {anchor_type.name}: no sequences with multiple profile hits ---")
            continue

        diffs.sort()
        max_diffs = [d[0] for d in diffs]
        n = len(max_diffs)
        print(f"\n--- {anchor_type.name} ({n_multi:,} seqs with 2+ profile hits) ---")
        print(f"  Position diff (best vs 2nd-best profile):")
        print(f"    exact={sum(1 for d in max_diffs if d == 0):,}  "
              f"within1={sum(1 for d in max_diffs if d <= 1):,}  "
              f"within5={sum(1 for d in max_diffs if d <= 5):,}  "
              f"within10={sum(1 for d in max_diffs if d <= 10):,}  "
              f">10={sum(1 for d in max_diffs if d > 10):,}")

        # Score gap between best and second-best
        score_gaps = [d[3] - d[4] for d in diffs]
        score_gaps.sort()
        print(f"  Score gap (best - 2nd): "
              f"min={score_gaps[0]:.1f}  median={score_gaps[len(score_gaps)//2]:.1f}  "
              f"max={score_gaps[-1]:.1f}")

        # Show worst position disagreements
        worst = sorted(diffs, reverse=True)[:5]
        if worst[0][0] > 5:
            print(f"  Worst disagreements:")
            for max_d, from_d, to_d, s1, s2, p1, p2 in worst:
                if max_d <= 5:
                    break
                print(f"    {max_d}bp: {p1} (score={s1:.1f}) vs {p2} (score={s2:.1f})")


def analyze_shortcircuit_safety(hits):
    """Would short-circuiting at threshold T miss the best profile?"""
    print("\n" + "=" * 80)
    print("4. SHORT-CIRCUIT SAFETY (does first-above-threshold match best profile?)")
    print("=" * 80)

    for anchor_type in AnchorType:
        anchor_hits = [h for h in hits if h.anchor_type == anchor_type and h.strand == "+"]
        if not anchor_hits:
            continue

        by_seq = defaultdict(list)
        for h in anchor_hits:
            by_seq[h.seq_id].append(h)

        for threshold in [15, 20, 25]:
            n_seqs_above = 0
            n_first_is_best = 0
            position_diffs = []

            for seq_id, seq_hits in by_seq.items():
                above = [h for h in seq_hits if h.score >= threshold]
                if not above:
                    continue
                n_seqs_above += 1

                best_overall = max(seq_hits, key=lambda h: h.score)
                # "first above threshold" = any profile that clears; in practice
                # we'd stop at the first one found in profile order, but for this
                # analysis we check if ANY above-threshold profile agrees with best
                best_above = max(above, key=lambda h: h.score)

                if best_above.profile_name == best_overall.profile_name:
                    n_first_is_best += 1

                pos_diff = max(
                    abs(best_above.env_from - best_overall.env_from),
                    abs(best_above.env_to - best_overall.env_to)
                )
                position_diffs.append(pos_diff)

            if not position_diffs:
                continue

            n = len(position_diffs)
            position_diffs.sort()
            exact = sum(1 for d in position_diffs if d == 0)
            within5 = sum(1 for d in position_diffs if d <= 5)
            print(f"\n  {anchor_type.name} @ threshold={threshold}: "
                  f"{n_seqs_above:,} seqs above threshold")
            print(f"    Best-above == overall-best profile: "
                  f"{n_first_is_best:,}/{n_seqs_above:,} ({100*n_first_is_best/n_seqs_above:.1f}%)")
            print(f"    Position agreement (best-above vs overall-best): "
                  f"exact={exact:,} within5={within5:,} "
                  f"median={position_diffs[n//2]}bp max={position_diffs[-1]}bp")


def generate_default_ordering(hits, n_seqs):
    """Produce a recommended static profile ordering per anchor type."""
    print("\n" + "=" * 80)
    print("5. RECOMMENDED DEFAULT PROFILE ORDERING")
    print("=" * 80)

    ordering = {}
    for anchor_type in AnchorType:
        anchor_hits = [h for h in hits if h.anchor_type == anchor_type]
        if not anchor_hits:
            continue

        # Rank by: number of sequences with a confident hit (score >= 20)
        confident_seqs = defaultdict(set)
        all_seqs = defaultdict(set)
        for h in anchor_hits:
            all_seqs[h.profile_name].add(h.seq_id)
            if h.score >= 20.0:
                confident_seqs[h.profile_name].add(h.seq_id)

        ranked = sorted(confident_seqs.items(), key=lambda x: len(x[1]), reverse=True)
        total = len({h.seq_id for h in anchor_hits})

        print(f"\n--- {anchor_type.name} ---")
        cumulative = set()
        ordering[anchor_type.name] = []
        for profile_name, seq_ids in ranked:
            if not seq_ids:
                break
            prev_coverage = len(cumulative)
            cumulative |= seq_ids
            marginal = len(cumulative) - prev_coverage
            ordering[anchor_type.name].append(profile_name)
            if marginal > 0:
                print(f"  {profile_name:<55s}  confident={len(seq_ids):>7,}  "
                      f"cumulative={len(cumulative):>7,}/{total:,}  "
                      f"marginal=+{marginal:,}")

        print(f"  Coverage with top profiles: {len(cumulative):,}/{total:,} "
              f"({100*len(cumulative)/total:.1f}%)")

    return ordering


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--organism", default="F")
    args = parser.parse_args()

    hmm_dir = find_hmm_dir()
    hmm_path = hmm_dir / f"{args.organism}.hmm"

    hits, n_seqs = collect_all_hits(hmm_path, args.input, cpus=args.cpus)

    analyze_profile_frequency(hits)
    analyze_score_distributions(hits)
    analyze_position_agreement(hits)
    analyze_shortcircuit_safety(hits)
    ordering = generate_default_ordering(hits, n_seqs)

    # Write ordering to JSON for potential inclusion in pyitsx
    out_path = Path(__file__).parent / "profile_ordering_F.json"
    with open(out_path, "w") as f:
        json.dump(ordering, f, indent=2)
    print(f"\nOrdering written to {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
