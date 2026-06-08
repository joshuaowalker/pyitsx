"""Unit tests for region extraction — uses synthetic AnchorChain objects."""

from pyitsx.constants import AnchorType, Confidence, Region, Strand
from pyitsx.models import AnchorChain, AnchorHit, RegionBounds
from pyitsx.regions import extract_regions


def _hit(anchor: AnchorType, env_from: int, env_to: int) -> AnchorHit:
    return AnchorHit(
        anchor_type=anchor,
        strand=Strand.PLUS,
        env_from=env_from,
        env_to=env_to,
        score=50.0,
        evalue=1e-10,
        profile_name="test",
    )


def _chain(*anchors_or_none) -> AnchorChain:
    return AnchorChain(
        anchors=tuple(anchors_or_none),
        strand=Strand.PLUS,
        total_score=sum(a.score for a in anchors_or_none if a is not None),
        confidence=Confidence.CONFIDENT,
    )


class TestExtractRegionsFull:

    def test_full_chain_five_regions(self):
        chain = _chain(
            _hit(AnchorType.SSU_END, 10, 54),
            _hit(AnchorType.S58_START, 317, 361),
            _hit(AnchorType.S58_END, 430, 474),
            _hit(AnchorType.LSU_START, 672, 724),
        )
        bounds = extract_regions(chain, seq_length=800)

        assert len(bounds) == 5
        regions = {b.region: b for b in bounds}

        assert regions[Region.SSU] == RegionBounds(Region.SSU, 1, 54)
        assert regions[Region.ITS1] == RegionBounds(Region.ITS1, 55, 316)
        assert regions[Region.S58] == RegionBounds(Region.S58, 317, 474)
        assert regions[Region.ITS2] == RegionBounds(Region.ITS2, 475, 671)
        assert regions[Region.LSU] == RegionBounds(Region.LSU, 672, 800)

    def test_its1_length(self):
        chain = _chain(
            _hit(AnchorType.SSU_END, 10, 54),
            _hit(AnchorType.S58_START, 317, 361),
            _hit(AnchorType.S58_END, 430, 474),
            _hit(AnchorType.LSU_START, 672, 724),
        )
        bounds = extract_regions(chain, seq_length=800)
        regions = {b.region: b for b in bounds}
        assert regions[Region.ITS1].length == 262  # 316 - 55 + 1

    def test_its2_length(self):
        chain = _chain(
            _hit(AnchorType.SSU_END, 10, 54),
            _hit(AnchorType.S58_START, 317, 361),
            _hit(AnchorType.S58_END, 430, 474),
            _hit(AnchorType.LSU_START, 672, 724),
        )
        bounds = extract_regions(chain, seq_length=800)
        regions = {b.region: b for b in bounds}
        assert regions[Region.ITS2].length == 197  # 671 - 475 + 1


class TestExtractRegionsPartial:

    def test_its1_only(self):
        chain = _chain(
            _hit(AnchorType.SSU_END, 10, 54),
            _hit(AnchorType.S58_START, 317, 361),
            None,
            None,
        )
        bounds = extract_regions(chain, seq_length=400)

        regions = {b.region: b for b in bounds}
        assert Region.SSU in regions
        assert Region.ITS1 in regions
        assert Region.ITS2 not in regions
        assert Region.LSU not in regions

    def test_its2_only(self):
        chain = _chain(
            None,
            None,
            _hit(AnchorType.S58_END, 100, 140),
            _hit(AnchorType.LSU_START, 350, 400),
        )
        bounds = extract_regions(chain, seq_length=500)

        regions = {b.region: b for b in bounds}
        assert Region.SSU not in regions
        assert Region.ITS1 not in regions
        assert Region.ITS2 in regions
        assert Region.LSU in regions

    def test_no_58s_without_both_anchors(self):
        chain = _chain(
            _hit(AnchorType.SSU_END, 10, 54),
            None,
            None,
            _hit(AnchorType.LSU_START, 672, 724),
        )
        bounds = extract_regions(chain, seq_length=800)

        regions = {b.region: b for b in bounds}
        assert Region.S58 not in regions
        assert Region.SSU in regions
        assert Region.LSU in regions


class TestRegionBounds:

    def test_length_property(self):
        b = RegionBounds(Region.ITS1, 55, 316)
        assert b.length == 262

    def test_single_base_region(self):
        b = RegionBounds(Region.ITS1, 100, 100)
        assert b.length == 1
