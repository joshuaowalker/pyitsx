from typing import Optional

from pyitsx.constants import AnchorType, Region
from pyitsx.models import AnchorChain, RegionBounds


def extract_regions(
    chain: AnchorChain, seq_length: int
) -> tuple[RegionBounds, ...]:
    a1 = chain.anchors[AnchorType.SSU_END.value - 1]
    a2 = chain.anchors[AnchorType.S58_START.value - 1]
    a3 = chain.anchors[AnchorType.S58_END.value - 1]
    a4 = chain.anchors[AnchorType.LSU_START.value - 1]

    bounds: list[RegionBounds] = []

    if a1 is not None:
        bounds.append(RegionBounds(Region.SSU, 1, a1.env_to))

    if a1 is not None and a2 is not None:
        start = a1.env_to + 1
        end = a2.env_from - 1
        if start <= end:
            bounds.append(RegionBounds(Region.ITS1, start, end))

    if a2 is not None and a3 is not None:
        bounds.append(RegionBounds(Region.S58, a2.env_from, a3.env_to))

    if a3 is not None and a4 is not None:
        start = a3.env_to + 1
        end = a4.env_from - 1
        if start <= end:
            bounds.append(RegionBounds(Region.ITS2, start, end))

    if a4 is not None:
        bounds.append(RegionBounds(Region.LSU, a4.env_from, seq_length))

    return tuple(bounds)
