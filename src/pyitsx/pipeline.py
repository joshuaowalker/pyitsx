from collections import defaultdict
from typing import Iterable, Union

import pyhmmer.easel

from pyitsx.chains import build_chain
from pyitsx.constants import AnchorType, Confidence, Strand
from pyitsx.models import (
    AnchorHit,
    ChainConstraints,
    ClassifyResult,
    DEFAULT_CONSTRAINTS,
    DelimitResult,
    OrientResult,
)
from pyitsx.profiles import ProfileDB
from pyitsx.regions import extract_regions


def orient(
    sequences: Union[
        pyhmmer.easel.DigitalSequenceBlock,
        Iterable[pyhmmer.easel.DigitalSequence],
    ],
    db: ProfileDB,
    cpus: int = 0,
) -> list[OrientResult]:
    hits_by_seq = db.search(sequences, cpus=cpus)
    results = []
    for seq_id, hits in hits_by_seq.items():
        result = _orient_from_hits(seq_id, hits)
        if result is not None:
            results.append(result)
    return results


def _orient_from_hits(
    seq_id: str, hits: list[AnchorHit]
) -> OrientResult | None:
    if not hits:
        return None

    score_by_strand: dict[Strand, float] = defaultdict(float)
    count_by_strand: dict[Strand, int] = defaultdict(int)
    top_score_by_strand: dict[Strand, float] = defaultdict(float)

    for hit in hits:
        score_by_strand[hit.strand] += hit.score
        count_by_strand[hit.strand] += 1
        if hit.score > top_score_by_strand[hit.strand]:
            top_score_by_strand[hit.strand] = hit.score

    best_strand = max(score_by_strand, key=lambda s: score_by_strand[s])
    return OrientResult(
        seq_id=seq_id,
        strand=best_strand,
        top_score=top_score_by_strand[best_strand],
        n_anchors=count_by_strand[best_strand],
    )


def classify(
    sequences: Union[
        pyhmmer.easel.DigitalSequenceBlock,
        Iterable[pyhmmer.easel.DigitalSequence],
    ],
    db: ProfileDB,
    cpus: int = 0,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
) -> list[ClassifyResult]:
    hits_by_seq = db.search(sequences, cpus=cpus)
    results = []
    for seq_id, hits in hits_by_seq.items():
        chain = build_chain(hits, constraints)
        if chain is None:
            continue
        has_its1 = (
            chain.anchors[AnchorType.SSU_END.value - 1] is not None
            and chain.anchors[AnchorType.S58_START.value - 1] is not None
        )
        has_its2 = (
            chain.anchors[AnchorType.S58_END.value - 1] is not None
            and chain.anchors[AnchorType.LSU_START.value - 1] is not None
        )
        results.append(
            ClassifyResult(
                seq_id=seq_id,
                strand=chain.strand,
                has_its1=has_its1,
                has_its2=has_its2,
                confidence=chain.confidence,
                chain=chain,
            )
        )
    return results


def delimit(
    sequences: Union[
        pyhmmer.easel.DigitalSequenceBlock,
        Iterable[pyhmmer.easel.DigitalSequence],
    ],
    db: ProfileDB,
    cpus: int = 0,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
) -> list[DelimitResult]:
    seq_lengths = _collect_seq_lengths(sequences)
    hits_by_seq = db.search(sequences, cpus=cpus)
    results = []
    for seq_id, hits in hits_by_seq.items():
        chain = build_chain(hits, constraints)
        if chain is None:
            continue
        seq_length = seq_lengths.get(seq_id, 0)
        bounds = extract_regions(chain, seq_length)
        results.append(
            DelimitResult(
                seq_id=seq_id,
                seq_length=seq_length,
                strand=chain.strand,
                chain=chain,
                bounds=bounds,
                confidence=chain.confidence,
            )
        )
    return results


def _collect_seq_lengths(
    sequences: Union[
        pyhmmer.easel.DigitalSequenceBlock,
        Iterable[pyhmmer.easel.DigitalSequence],
    ],
) -> dict[str, int]:
    lengths = {}
    try:
        for seq in sequences:
            lengths[seq.name] = len(seq)
    except TypeError:
        pass
    return lengths
