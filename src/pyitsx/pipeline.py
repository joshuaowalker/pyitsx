from collections import defaultdict

from pyitsx.chains import build_chain
from pyitsx.constants import Confidence, Region, Strand
from pyitsx.models import (
    AnchorHit,
    ChainConstraints,
    ClassifyResult,
    DEFAULT_CONSTRAINTS,
    DelimitResult,
    OrientResult,
)
from pyitsx.profiles import ProfileDB, SequenceInput
from pyitsx.regions import extract_regions


def orient(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 0,
) -> list[OrientResult]:
    seqs = db.prepare(sequences)
    hits_by_seq = db.search(seqs, cpus=cpus)
    results = []
    for seq_id, hits in hits_by_seq.items():
        result = _orient_from_hits(seq_id, hits)
        if result is not None:
            results.append(result)
    return results


def classify(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 0,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
) -> list[ClassifyResult]:
    seqs = db.prepare(sequences)
    seq_lengths = {s.name: len(s) for s in seqs}
    hits_by_seq = db.search(seqs, cpus=cpus)
    results = []
    for seq_id, hits in hits_by_seq.items():
        chain = build_chain(hits, constraints)
        if chain is None:
            continue
        seq_length = seq_lengths.get(seq_id, 0)
        bounds = extract_regions(chain, seq_length)
        region_set = {b.region for b in bounds}
        results.append(
            ClassifyResult(
                seq_id=seq_id,
                strand=chain.strand,
                has_its1=Region.ITS1 in region_set,
                has_its2=Region.ITS2 in region_set,
                confidence=chain.confidence,
                chain=chain,
            )
        )
    return results


def delimit(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 0,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
) -> list[DelimitResult]:
    seqs = db.prepare(sequences)
    seq_lengths = {s.name: len(s) for s in seqs}
    hits_by_seq = db.search(seqs, cpus=cpus)
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
