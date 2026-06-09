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

from pyitsx.constants import Organism, SearchMode
from pyitsx.models import OrganismResult, OrganismScore
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
            scores_by_seq[seq_id].append(
                OrganismScore(
                    organism=org,
                    n_anchors=len(hits),
                    total_score=sum(h.score for h in hits),
                    best_evalue=min(h.evalue for h in hits),
                )
            )

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
