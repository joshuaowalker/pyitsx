from collections import defaultdict
from itertools import product
from typing import Optional

from pyitsx.constants import AnchorType, Confidence, Strand
from pyitsx.models import AnchorChain, AnchorHit, ChainConstraints, DEFAULT_CONSTRAINTS


def build_chain(
    hits: list[AnchorHit],
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
    max_per_anchor: int = 8,
) -> Optional[AnchorChain]:
    if not hits:
        return None

    grouped = _group_hits(hits, max_per_anchor)

    best_full: Optional[AnchorChain] = None
    for strand in (Strand.PLUS, Strand.MINUS):
        chain = _best_full_chain(grouped, strand, constraints)
        if chain and (best_full is None or chain.total_score > best_full.total_score):
            best_full = chain

    if best_full is not None:
        return best_full

    best_partial: Optional[AnchorChain] = None
    for strand in (Strand.PLUS, Strand.MINUS):
        chain = _best_partial_chain(grouped, strand, constraints)
        if chain and (best_partial is None or chain.total_score > best_partial.total_score):
            best_partial = chain

    return best_partial


def _group_hits(
    hits: list[AnchorHit],
    max_per_anchor: int,
) -> dict[tuple[Strand, AnchorType], list[AnchorHit]]:
    groups: dict[tuple[Strand, AnchorType], list[AnchorHit]] = defaultdict(list)
    for h in hits:
        groups[(h.strand, h.anchor_type)].append(h)

    for key in groups:
        groups[key].sort(key=lambda h: h.score, reverse=True)
        groups[key] = groups[key][:max_per_anchor]

    return dict(groups)


def _best_full_chain(
    grouped: dict[tuple[Strand, AnchorType], list[AnchorHit]],
    strand: Strand,
    constraints: ChainConstraints,
) -> Optional[AnchorChain]:
    a1s = grouped.get((strand, AnchorType.SSU_END), [])
    a2s = grouped.get((strand, AnchorType.S58_START), [])
    a3s = grouped.get((strand, AnchorType.S58_END), [])
    a4s = grouped.get((strand, AnchorType.LSU_START), [])

    if not (a1s and a2s and a3s and a4s):
        return None

    best: Optional[AnchorChain] = None

    for a1, a2, a3, a4 in product(a1s, a2s, a3s, a4s):
        if not (a1.env_to < a2.env_from and a2.env_to < a3.env_from and a3.env_to < a4.env_from):
            continue

        its1_len = a2.env_from - a1.env_to - 1
        if not (constraints.min_its1 <= its1_len <= constraints.max_its1):
            continue

        its2_len = a4.env_from - a3.env_to - 1
        if not (constraints.min_its2 <= its2_len <= constraints.max_its2):
            continue

        full_len = a4.env_from - a1.env_to - 1
        if not (constraints.min_full <= full_len <= constraints.max_full):
            continue

        total = a1.score + a2.score + a3.score + a4.score
        if best is None or total > best.total_score:
            confidence = _assess_confidence(
                [a1, a2, a3, a4], constraints
            )
            best = AnchorChain(
                anchors=(a1, a2, a3, a4),
                strand=strand,
                total_score=total,
                confidence=confidence,
            )

    return best


def _best_partial_chain(
    grouped: dict[tuple[Strand, AnchorType], list[AnchorHit]],
    strand: Strand,
    constraints: ChainConstraints,
) -> Optional[AnchorChain]:
    pairs = [
        (AnchorType.SSU_END, AnchorType.S58_START, constraints.min_its1, constraints.max_its1),
        (AnchorType.S58_END, AnchorType.LSU_START, constraints.min_its2, constraints.max_its2),
        (AnchorType.SSU_END, AnchorType.LSU_START, constraints.min_full, constraints.max_full),
    ]

    best: Optional[AnchorChain] = None

    for left_type, right_type, min_len, max_len in pairs:
        lefts = grouped.get((strand, left_type), [])
        rights = grouped.get((strand, right_type), [])

        for left, right in product(lefts, rights):
            if left.env_to >= right.env_from:
                continue

            region_len = right.env_from - left.env_to - 1
            if not (min_len <= region_len <= max_len):
                continue

            total = left.score + right.score
            if best is None or total > best.total_score:
                anchors = [None, None, None, None]
                anchors[left_type.value - 1] = left
                anchors[right_type.value - 1] = right
                best = AnchorChain(
                    anchors=tuple(anchors),
                    strand=strand,
                    total_score=total,
                    confidence=Confidence.PARTIAL,
                )

    return best


def _assess_confidence(
    anchors: list[AnchorHit],
    constraints: ChainConstraints,
) -> Confidence:
    for a in anchors:
        if a.score < constraints.min_anchor_score or a.evalue > constraints.max_anchor_evalue:
            return Confidence.AMBIGUOUS
    return Confidence.CONFIDENT
