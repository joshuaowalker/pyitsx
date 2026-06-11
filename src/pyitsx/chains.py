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

    if best_partial is not None:
        return best_partial

    best_single: Optional[AnchorChain] = None
    for strand in (Strand.PLUS, Strand.MINUS):
        chain = _best_single_anchor_chain(grouped, strand, constraints)
        if chain and (best_single is None or chain.total_score > best_single.total_score):
            best_single = chain

    return best_single


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
        # 5.8S-only: anchors 2+3 span the interior of 5.8S (~67bp gap between ~45bp models)
        (AnchorType.S58_START, AnchorType.S58_END, 20, 200),
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


def _best_single_anchor_chain(
    grouped: dict[tuple[Strand, AnchorType], list[AnchorHit]],
    strand: Strand,
    constraints: ChainConstraints,
) -> Optional[AnchorChain]:
    best: Optional[AnchorChain] = None
    for anchor_type in AnchorType:
        candidates = grouped.get((strand, anchor_type), [])
        if not candidates:
            continue
        hit = candidates[0]
        if (hit.score < constraints.min_anchor_score
                or hit.evalue > constraints.max_anchor_evalue):
            continue
        if best is not None and hit.score <= best.total_score:
            continue
        anchors = [None, None, None, None]
        anchors[anchor_type.value - 1] = hit
        best = AnchorChain(
            anchors=tuple(anchors),
            strand=strand,
            total_score=hit.score,
            confidence=Confidence.PARTIAL,
        )
    return best


def detect_chimera(
    hits: list[AnchorHit],
    score_threshold: float = 20.0,
) -> bool:
    if not hits:
        return False

    plus_score = sum(h.score for h in hits if h.strand == Strand.PLUS and h.score >= score_threshold)
    minus_score = sum(h.score for h in hits if h.strand == Strand.MINUS and h.score >= score_threshold)
    if plus_score >= score_threshold and minus_score >= score_threshold:
        return True

    dominant = Strand.PLUS if plus_score >= minus_score else Strand.MINUS
    best_by_anchor: dict[AnchorType, AnchorHit] = {}
    for h in hits:
        if h.strand != dominant or h.score < score_threshold:
            continue
        if h.anchor_type not in best_by_anchor or h.score > best_by_anchor[h.anchor_type].score:
            best_by_anchor[h.anchor_type] = h

    ordered = sorted(best_by_anchor.values(), key=lambda h: h.anchor_type.value)
    for i in range(len(ordered) - 1):
        if dominant == Strand.PLUS:
            if max(ordered[i].env_from, ordered[i].env_to) >= min(ordered[i + 1].env_from, ordered[i + 1].env_to):
                return True
        else:
            if min(ordered[i].env_from, ordered[i].env_to) <= max(ordered[i + 1].env_from, ordered[i + 1].env_to):
                return True

    return False


def _assess_confidence(
    anchors: list[AnchorHit],
    constraints: ChainConstraints,
) -> Confidence:
    for a in anchors:
        if a.score < constraints.min_anchor_score or a.evalue > constraints.max_anchor_evalue:
            return Confidence.AMBIGUOUS
    return Confidence.CONFIDENT
