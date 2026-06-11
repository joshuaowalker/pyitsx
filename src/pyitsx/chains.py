from collections import defaultdict
from itertools import product
from typing import Optional

from pyitsx.constants import AnchorType, Confidence, Strand
from pyitsx.models import AnchorChain, AnchorHit, ChainConstraints, DEFAULT_CONSTRAINTS


def _normalize_minus_hits(
    hits: list[AnchorHit], seq_length: int
) -> list[AnchorHit]:
    out = []
    for h in hits:
        if h.strand == Strand.MINUS:
            h = AnchorHit(
                anchor_type=h.anchor_type,
                strand=h.strand,
                env_from=seq_length - h.env_from + 1,
                env_to=seq_length - h.env_to + 1,
                score=h.score,
                evalue=h.evalue,
                profile_name=h.profile_name,
            )
        out.append(h)
    return out


def build_chain(
    hits: list[AnchorHit],
    constraints: ChainConstraints = DEFAULT_CONSTRAINTS,
    max_per_anchor: int = 8,
    seq_length: int = 0,
) -> Optional[AnchorChain]:
    if not hits:
        return None

    has_minus = any(h.strand == Strand.MINUS for h in hits)
    if has_minus and seq_length <= 0:
        raise ValueError(
            "seq_length is required when hits contain minus-strand entries"
        )
    if has_minus:
        hits = _normalize_minus_hits(hits, seq_length)

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


def _hit_left(h: AnchorHit) -> int:
    return min(h.env_from, h.env_to)


def _hit_right(h: AnchorHit) -> int:
    return max(h.env_from, h.env_to)


def detect_chimera(
    hits: list[AnchorHit],
    score_threshold: float = 20.0,
    min_cross_strand_hits: int = 1,
) -> bool:
    if not hits:
        return False

    significant = [h for h in hits if h.score >= score_threshold]
    if not significant:
        return False

    plus_hits = [h for h in significant if h.strand == Strand.PLUS]
    minus_hits = [h for h in significant if h.strand == Strand.MINUS]
    plus_score = sum(h.score for h in plus_hits)
    minus_score = sum(h.score for h in minus_hits)

    if plus_hits and minus_hits:
        if plus_score >= minus_score:
            dominant_hits, minor_hits = plus_hits, minus_hits
        else:
            dominant_hits, minor_hits = minus_hits, plus_hits

        dom_left = min(_hit_left(h) for h in dominant_hits)
        dom_right = max(_hit_right(h) for h in dominant_hits)

        n_outside = sum(
            1 for h in minor_hits
            if _hit_right(h) < dom_left or _hit_left(h) > dom_right
        )
        if n_outside >= min_cross_strand_hits:
            return True

    dominant = Strand.PLUS if plus_score >= minus_score else Strand.MINUS
    best_by_anchor: dict[AnchorType, AnchorHit] = {}
    for h in significant:
        if h.strand != dominant:
            continue
        if h.anchor_type not in best_by_anchor or h.score > best_by_anchor[h.anchor_type].score:
            best_by_anchor[h.anchor_type] = h

    ordered = sorted(best_by_anchor.values(), key=lambda h: h.anchor_type.value)
    for i in range(len(ordered) - 1):
        if dominant == Strand.PLUS:
            if _hit_right(ordered[i]) >= _hit_left(ordered[i + 1]):
                return True
        else:
            if _hit_left(ordered[i]) <= _hit_right(ordered[i + 1]):
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
