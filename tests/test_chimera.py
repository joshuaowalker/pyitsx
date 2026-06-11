"""Unit tests for chimera detection — synthetic AnchorHit objects, no HMMs."""

from pyitsx.chains import detect_chimera
from pyitsx.constants import AnchorType, Strand
from pyitsx.models import AnchorHit


def _hit(
    anchor: AnchorType,
    env_from: int,
    env_to: int,
    score: float = 50.0,
    strand: Strand = Strand.PLUS,
) -> AnchorHit:
    return AnchorHit(
        anchor_type=anchor,
        strand=strand,
        env_from=env_from,
        env_to=env_to,
        score=score,
        evalue=1e-10,
        profile_name="test",
    )


class TestDetectChimera:

    def test_empty_hits(self):
        assert detect_chimera([]) is False

    def test_single_anchor(self):
        hits = [_hit(AnchorType.SSU_END, 10, 50)]
        assert detect_chimera(hits) is False

    def test_normal_chain_not_chimeric(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50),
            _hit(AnchorType.S58_START, 200, 240),
            _hit(AnchorType.S58_END, 300, 340),
            _hit(AnchorType.LSU_START, 500, 550),
        ]
        assert detect_chimera(hits) is False

    def test_cross_strand_chimera(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 200, 240, strand=Strand.PLUS),
            _hit(AnchorType.S58_END, 300, 340, strand=Strand.MINUS),
            _hit(AnchorType.LSU_START, 500, 550, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is True

    def test_out_of_order_chimera(self):
        hits = [
            _hit(AnchorType.SSU_END, 400, 450),
            _hit(AnchorType.S58_START, 100, 140),
            _hit(AnchorType.S58_END, 200, 240),
            _hit(AnchorType.LSU_START, 500, 550),
        ]
        assert detect_chimera(hits) is True

    def test_weak_cross_strand_not_chimeric(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, score=50.0, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 200, 240, score=50.0, strand=Strand.PLUS),
            _hit(AnchorType.S58_END, 300, 340, score=5.0, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is False

    def test_normal_minus_strand_not_chimeric(self):
        hits = [
            _hit(AnchorType.SSU_END, 788, 744, strand=Strand.MINUS),
            _hit(AnchorType.S58_START, 486, 443, strand=Strand.MINUS),
            _hit(AnchorType.S58_END, 373, 329, strand=Strand.MINUS),
            _hit(AnchorType.LSU_START, 109, 65, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is False

    def test_out_of_order_minus_strand_chimera(self):
        hits = [
            _hit(AnchorType.SSU_END, 200, 150, strand=Strand.MINUS),
            _hit(AnchorType.S58_START, 600, 550, strand=Strand.MINUS),
            _hit(AnchorType.S58_END, 373, 329, strand=Strand.MINUS),
            _hit(AnchorType.LSU_START, 109, 65, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is True

    def test_single_minor_hit_outside_is_chimeric(self):
        """One strong minor-strand hit outside dominant span — chimeric."""
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 200, 240, strand=Strand.PLUS),
            _hit(AnchorType.S58_END, 300, 340, strand=Strand.PLUS),
            _hit(AnchorType.LSU_START, 500, 550, strand=Strand.PLUS),
            _hit(AnchorType.SSU_END, 600, 640, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is True

    def test_single_minor_hit_inside_not_chimeric(self):
        """One strong minor-strand hit inside dominant span — background complementarity."""
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 200, 240, strand=Strand.PLUS),
            _hit(AnchorType.S58_END, 300, 340, strand=Strand.PLUS),
            _hit(AnchorType.LSU_START, 500, 550, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 250, 290, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is False

    def test_two_minor_hits_inside_not_chimeric(self):
        """Two minor-strand hits inside dominant span — not disjoint chimera."""
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 200, 240, strand=Strand.PLUS),
            _hit(AnchorType.S58_END, 300, 340, strand=Strand.PLUS),
            _hit(AnchorType.LSU_START, 500, 550, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 100, 140, strand=Strand.MINUS),
            _hit(AnchorType.S58_END, 400, 440, strand=Strand.MINUS),
        ]
        assert detect_chimera(hits) is False

    def test_minus_dominant_cross_strand_chimera(self):
        """Symmetry: minus dominant with two plus hits outside its span."""
        hits = [
            _hit(AnchorType.SSU_END, 788, 744, strand=Strand.MINUS),
            _hit(AnchorType.S58_START, 486, 443, strand=Strand.MINUS),
            _hit(AnchorType.S58_END, 373, 329, strand=Strand.MINUS),
            _hit(AnchorType.LSU_START, 109, 65, strand=Strand.MINUS),
            _hit(AnchorType.SSU_END, 10, 50, strand=Strand.PLUS),
            _hit(AnchorType.S58_START, 850, 890, strand=Strand.PLUS),
        ]
        assert detect_chimera(hits) is True

    def test_all_below_threshold_not_chimeric(self):
        hits = [
            _hit(AnchorType.SSU_END, 10, 50, score=5.0),
            _hit(AnchorType.S58_START, 200, 240, score=5.0),
        ]
        assert detect_chimera(hits) is False
