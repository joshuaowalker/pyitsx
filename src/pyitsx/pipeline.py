from collections import defaultdict
from typing import Optional, Sequence

from pyitsx.chains import build_chain, detect_chimera
from pyitsx.constants import Confidence, Region, SearchMode, Strand
from pyitsx.models import (
    AnchorHit,
    ChainConstraints,
    ClassifyResult,
    DEFAULT_CONSTRAINTS,
    DelimitResult,
    ExtractionResult,
    OrientResult,
)
from pyitsx.profiles import DEFAULT_BATCH_SIZE, ProfileDB, SequenceInput
from pyitsx.regions import extract_regions

_COMPLEMENT = str.maketrans(
    "ACGTacgtNnRYSWKMBDHVryswkmbdhv",
    "TGCAtgcaNnYRSWMKVHDByrswmkvhdb",
)


def orient(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
    mode: SearchMode = SearchMode.FAST,
) -> list[OrientResult]:
    seqs = db.prepare(sequences)
    all_seq_ids = [s.name for s in seqs]
    hits_by_seq = db.search(seqs, cpus=cpus, batch_size=batch_size, mode=mode)
    results = []
    detected_ids: set[str] = set()
    for seq_id, hits in hits_by_seq.items():
        result = _orient_from_hits(seq_id, hits)
        if result is not None:
            detected_ids.add(seq_id)
            results.append(result)
    for seq_id in all_seq_ids:
        if seq_id not in detected_ids:
            results.append(
                OrientResult(
                    seq_id=seq_id,
                    strand=None,
                    top_score=0.0,
                    n_anchors=0,
                )
            )
    return results


def classify(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
    mode: SearchMode = SearchMode.FAST,
) -> list[ClassifyResult]:
    seqs = db.prepare(sequences)
    seq_lengths = {s.name: len(s) for s in seqs}
    all_seq_ids = [s.name for s in seqs]
    hits_by_seq = db.search(seqs, cpus=cpus, batch_size=batch_size, mode=mode)
    results = []
    detected_ids: set[str] = set()
    for seq_id, hits in hits_by_seq.items():
        chain = build_chain(hits, constraints)
        if chain is None:
            continue
        detected_ids.add(seq_id)
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
                chimeric=detect_chimera(hits),
            )
        )
    for seq_id in all_seq_ids:
        if seq_id not in detected_ids:
            results.append(
                ClassifyResult(
                    seq_id=seq_id,
                    strand=None,
                    has_its1=False,
                    has_its2=False,
                    confidence=Confidence.NONE,
                    chain=None,
                )
            )
    return results


def delimit(
    sequences: SequenceInput,
    db: ProfileDB,
    cpus: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
    mode: SearchMode = SearchMode.FAST,
) -> list[DelimitResult]:
    seqs = db.prepare(sequences)
    seq_lengths = {s.name: len(s) for s in seqs}
    all_seq_ids = [s.name for s in seqs]
    hits_by_seq = db.search(seqs, cpus=cpus, batch_size=batch_size, mode=mode)
    results = []
    detected_ids: set[str] = set()
    for seq_id, hits in hits_by_seq.items():
        chain = build_chain(hits, constraints)
        if chain is None:
            continue
        detected_ids.add(seq_id)
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
                chimeric=detect_chimera(hits),
            )
        )
    for seq_id in all_seq_ids:
        if seq_id not in detected_ids:
            results.append(
                DelimitResult(
                    seq_id=seq_id,
                    seq_length=seq_lengths.get(seq_id, 0),
                    strand=None,
                    chain=None,
                    bounds=(),
                    confidence=Confidence.NONE,
                )
            )
    return results


def extract(
    sequences: SequenceInput,
    db: ProfileDB,
    regions: Optional[Sequence[Region]] = None,
    cpus: int = 1,
    batch_size: int = DEFAULT_BATCH_SIZE,
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
    mode: SearchMode = SearchMode.FAST,
) -> list[ExtractionResult]:
    seqs = db.prepare(sequences)
    text_by_id: dict[str, str] = {}
    for ds in seqs:
        ts = ds.textize()
        text_by_id[ds.name] = ts.sequence

    delimit_results = delimit(seqs, db, cpus=cpus, batch_size=batch_size, constraints=constraints, mode=mode)
    target_regions = set(regions) if regions is not None else None
    wants_full_its = target_regions is not None and Region.FULL_ITS in target_regions
    normal_regions = (target_regions - {Region.FULL_ITS}) if target_regions is not None else None

    extracted: list[ExtractionResult] = []
    for r in delimit_results:
        if r.confidence == Confidence.NONE:
            continue
        seq_text = text_by_id.get(r.seq_id)
        if seq_text is None:
            continue
        if r.strand == Strand.MINUS:
            seq_text = seq_text.translate(_COMPLEMENT)[::-1]
        for b in r.bounds:
            if normal_regions is not None and b.region not in normal_regions:
                continue
            extracted.append(
                ExtractionResult(
                    seq_id=r.seq_id,
                    region=b.region,
                    start=b.start,
                    end=b.end,
                    sequence=seq_text[b.start - 1 : b.end],
                )
            )
        if wants_full_its:
            full = r.full_its
            if full is not None:
                extracted.append(
                    ExtractionResult(
                        seq_id=r.seq_id,
                        region=Region.FULL_ITS,
                        start=full.start,
                        end=full.end,
                        sequence=seq_text[full.start - 1 : full.end],
                    )
                )
    return extracted


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
        chimeric=detect_chimera(hits),
    )
