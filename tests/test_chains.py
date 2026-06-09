"""Unit tests for chain building — uses synthetic AnchorHit objects, no HMMs needed."""

import pytest

from pyitsx.chains import build_chain, _group_hits
from pyitsx.constants import AnchorType, Confidence, Strand
from pyitsx.models import AnchorChain, AnchorHit, ChainConstraints


def _hit(
    anchor: AnchorType,
    env_from: int,
    env_to: int,
    score: float = 50.0,
    evalue: float = 1e-10,
    strand: Strand = Strand.PLUS,
    profile: str = "test",
) -> AnchorHit:
    return AnchorHit(
        anchor_type=anchor,
        strand=strand,
        env_from=env_from,
        env_to=env_to,
        score=score,
        evalue=evalue,
        profile_name=profile,
    )


CONSTRAINTS = ChainConstraints(
    min_its1=50,
    max_its1=1500,
    min_its2=50,
    max_its2=2000,
    min_full=150,
    max_full=4000,
)


class TestGroupHits:

    def test_groups_by_strand_and_type(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, score=30),
            _hit(AnchorType.SSU_END, 10, 50, score=40),
            _hit(AnchorType.S58_START, 200, 240, score=50),
        ]
        grouped = _group_hits(hits, max_per_anchor=8)
        assert len(grouped[(Strand.PLUS, AnchorType.SSU_END)]) == 2
        assert len(grouped[(Strand.PLUS, AnchorType.S58_START)]) == 1

    def test_keeps_top_k(self):
        hits = [_hit(AnchorType.SSU_END, 10, 50, score=i) for i in range(20)]
        grouped = _group_hits(hits, max_per_anchor=3)
        ssu_hits = grouped[(Strand.PLUS, AnchorType.SSU_END)]
        assert len(ssu_hits) == 3
        assert ssu_hits[0].score == 19

    def test_separates_strands(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.MINUS),
        ]
        grouped = _group_hits(hits, max_per_anchor=8)
        assert len(grouped[(Strand.PLUS, AnchorType.SSU_END)]) == 1
        assert len(grouped[(Strand.MINUS, AnchorType.SSU_END)]) == 1


class TestBuildChainFull:

    def _four_anchors(self, strand=Strand.PLUS, scores=(50, 50, 50, 50)):
        return [
            _hit(AnchorType.SSU_END, 10, 54, score=scores[0], strand=strand),
            _hit(AnchorType.S58_START, 317, 361, score=scores[1], strand=strand),
            _hit(AnchorType.S58_END, 430, 474, score=scores[2], strand=strand),
            _hit(AnchorType.LSU_START, 672, 724, score=scores[3], strand=strand),
        ]

    def test_valid_four_anchor_chain(self):
        hits = self._four_anchors()
        chain = build_chain(hits, CONSTRAINTS)

        assert chain is not None
        assert chain.is_full
        assert chain.n_anchors == 4
        assert chain.strand == Strand.PLUS
        assert chain.total_score == 200.0
        assert chain.confidence == Confidence.CONFIDENT

    def test_selects_higher_scoring_strand(self):
        plus_hits = self._four_anchors(Strand.PLUS, scores=(40, 40, 40, 40))
        minus_hits = self._four_anchors(Strand.MINUS, scores=(60, 60, 60, 60))
        chain = build_chain(plus_hits + minus_hits, CONSTRAINTS)

        assert chain is not None
        assert chain.strand == Strand.MINUS
        assert chain.total_score == 240.0

    def test_rejects_wrong_anchor_order(self):
        hits = [
            _hit(AnchorType.SSU_END, 400, 450),
            _hit(AnchorType.S58_START, 100, 140),
            _hit(AnchorType.S58_END, 200, 240),
            _hit(AnchorType.LSU_START, 300, 350),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is None or not chain.is_full

    def test_rejects_its1_too_short(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50),
            _hit(AnchorType.S58_START, 60, 100),  # ITS1 = 60-50-1 = 9bp < 50
            _hit(AnchorType.S58_END, 170, 210),
            _hit(AnchorType.LSU_START, 400, 450),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is None or not chain.is_full

    def test_rejects_its2_too_short(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50),
            _hit(AnchorType.S58_START, 200, 240),
            _hit(AnchorType.S58_END, 310, 350),
            _hit(AnchorType.LSU_START, 360, 410),  # ITS2 = 360-350-1 = 9bp < 50
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is None or not chain.is_full

    def test_rejects_full_too_long(self):
        constraints = ChainConstraints(max_full=500)
        hits = [
            _hit(AnchorType.SSU_END, 10, 50),
            _hit(AnchorType.S58_START, 200, 240),
            _hit(AnchorType.S58_END, 310, 350),
            _hit(AnchorType.LSU_START, 4200, 4250),
        ]
        chain = build_chain(hits, constraints)
        assert chain is None or not chain.is_full

    def test_picks_best_among_multiple_candidates(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 54, score=50),
            _hit(AnchorType.S58_START, 317, 361, score=50),
            _hit(AnchorType.S58_START, 320, 364, score=80),
            _hit(AnchorType.S58_END, 430, 474, score=50),
            _hit(AnchorType.LSU_START, 672, 724, score=50),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is not None
        assert chain.anchors[1].score == 80


class TestBuildChainPartial:

    def test_partial_its1(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, score=40),
            _hit(AnchorType.S58_START, 300, 340, score=40),
        ]
        chain = build_chain(hits, CONSTRAINTS)

        assert chain is not None
        assert chain.confidence == Confidence.PARTIAL
        assert chain.n_anchors == 2
        assert chain.anchors[0] is not None  # SSU_END
        assert chain.anchors[1] is not None  # S58_START
        assert chain.anchors[2] is None
        assert chain.anchors[3] is None

    def test_partial_its2(self):
        hits = [
            _hit(AnchorType.S58_END, 100, 140, score=40),
            _hit(AnchorType.LSU_START, 350, 400, score=40),
        ]
        chain = build_chain(hits, CONSTRAINTS)

        assert chain is not None
        assert chain.confidence == Confidence.PARTIAL
        assert chain.anchors[0] is None
        assert chain.anchors[1] is None
        assert chain.anchors[2] is not None  # S58_END
        assert chain.anchors[3] is not None  # LSU_START

    def test_full_chain_preferred_over_partial(self):
        full_hits = [
            _hit(AnchorType.SSU_END, 10, 54, score=30),
            _hit(AnchorType.S58_START, 317, 361, score=30),
            _hit(AnchorType.S58_END, 430, 474, score=30),
            _hit(AnchorType.LSU_START, 672, 724, score=30),
        ]
        chain = build_chain(full_hits, CONSTRAINTS)
        assert chain is not None
        assert chain.is_full
        assert chain.confidence != Confidence.PARTIAL

    def test_no_chain_from_empty(self):
        assert build_chain([], CONSTRAINTS) is None

    def test_single_ssu_end_builds_partial_chain(self):
        hits = [_hit(AnchorType.SSU_END, 10, 50, score=100)]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is not None
        assert chain.confidence == Confidence.PARTIAL
        assert chain.n_anchors == 1
        assert chain.anchors[0] is not None
        assert chain.anchors[1] is None

    def test_single_weak_anchor_no_chain(self):
        hits = [_hit(AnchorType.SSU_END, 10, 50, score=5, evalue=0.1)]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is None

    def test_single_s58_anchor_no_chain(self):
        hits = [_hit(AnchorType.S58_START, 100, 150, score=100)]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain is None


class TestConfidence:

    def test_confident_when_all_above_thresholds(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 54, score=50, evalue=1e-10),
            _hit(AnchorType.S58_START, 317, 361, score=50, evalue=1e-10),
            _hit(AnchorType.S58_END, 430, 474, score=50, evalue=1e-10),
            _hit(AnchorType.LSU_START, 672, 724, score=50, evalue=1e-10),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain.confidence == Confidence.CONFIDENT

    def test_ambiguous_when_score_below_threshold(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 54, score=5, evalue=1e-10),
            _hit(AnchorType.S58_START, 317, 361, score=50, evalue=1e-10),
            _hit(AnchorType.S58_END, 430, 474, score=50, evalue=1e-10),
            _hit(AnchorType.LSU_START, 672, 724, score=50, evalue=1e-10),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain.confidence == Confidence.AMBIGUOUS

    def test_ambiguous_when_evalue_above_threshold(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 54, score=50, evalue=1e-10),
            _hit(AnchorType.S58_START, 317, 361, score=50, evalue=0.1),
            _hit(AnchorType.S58_END, 430, 474, score=50, evalue=1e-10),
            _hit(AnchorType.LSU_START, 672, 724, score=50, evalue=1e-10),
        ]
        chain = build_chain(hits, CONSTRAINTS)
        assert chain.confidence == Confidence.AMBIGUOUS
