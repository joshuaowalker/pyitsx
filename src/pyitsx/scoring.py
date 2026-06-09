"""Multi-organism scoring for ITS sequences.

Scores sequences against multiple organism HMM profile sets to determine
the best-matching taxonomic group. This is expensive — one full HMM search
per organism — and is separated from the main pipeline to make that cost
explicit.
"""

import logging
from collections import defaultdict
from pathlib import Path
from typing import Optional, Sequence

from pyitsx.constants import AnchorType, Organism, SearchMode
from pyitsx.models import AnchorHit, OrganismResult, OrganismScore
from pyitsx.profiles import DEFAULT_BATCH_SIZE, ProfileDB, SequenceInput, find_hmm_dir

logger = logging.getLogger(__name__)


def score_organisms(
    sequences: SequenceInput,
    hmm_dir: Optional[Path] = None,
    organisms: Optional[Sequence[Organism]] = None,
    cpus: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
    mode: SearchMode = SearchMode.FAST,
) -> list[OrganismResult]:
    if hmm_dir is None:
        hmm_dir = find_hmm_dir()
    if organisms is None:
        organisms = list(Organism)

    dbs: list[tuple[Organism, ProfileDB]] = []
    for org in organisms:
        try:
            db = ProfileDB(hmm_dir, organism=org.name)
        except (FileNotFoundError, EOFError):
            logger.debug("Skipping organism %s: no usable HMM profiles", org.name)
            continue
        if db.n_profiles == 0:
            continue
        dbs.append((org, db))

    if not dbs:
        raise FileNotFoundError("No usable organism HMM profiles found")

    seqs = dbs[0][1].prepare(sequences)
    all_seq_ids = [s.name for s in seqs]

    scores_by_seq: dict[str, list[OrganismScore]] = defaultdict(list)

    for org, db in dbs:
        hits_by_seq = db.search(seqs, cpus=cpus, batch_size=batch_size, mode=mode)

        for seq_id in all_seq_ids:
            hits = hits_by_seq.get(seq_id, [])
            if not hits:
                continue
            score = _score_hits(org, hits)
            scores_by_seq[seq_id].append(score)

    results = []
    for seq_id in all_seq_ids:
        org_scores = sorted(
            scores_by_seq.get(seq_id, []),
            key=lambda s: s.total_score,
            reverse=True,
        )
        results.append(
            OrganismResult(
                seq_id=seq_id,
                best=org_scores[0] if org_scores else None,
                scores=tuple(org_scores),
            )
        )

    return results


def _score_hits(org: Organism, hits: list[AnchorHit]) -> OrganismScore:
    best_by_anchor: dict[AnchorType, AnchorHit] = {}
    for h in hits:
        prev = best_by_anchor.get(h.anchor_type)
        if prev is None or h.score > prev.score:
            best_by_anchor[h.anchor_type] = h
    anchors = list(best_by_anchor.values())
    return OrganismScore(
        organism=org,
        n_anchors=len(anchors),
        total_score=sum(h.score for h in anchors),
        best_evalue=min(h.evalue for h in anchors),
    )
