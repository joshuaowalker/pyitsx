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

    # ITS1: prefer anchor1..anchor2 boundary; infer boundaries from seq edges
    if a1 is not None and a2 is not None:
        start = a1.env_to + 1
        end = a2.env_from - 1
        if start <= end:
            bounds.append(RegionBounds(Region.ITS1, start, end))
    elif a1 is None and a2 is not None and a2.env_from > 1:
        bounds.append(RegionBounds(Region.ITS1, 1, a2.env_from - 1))
    elif a1 is not None and a2 is None and a3 is None and a4 is None and a1.env_to < seq_length:
        bounds.append(RegionBounds(Region.ITS1, a1.env_to + 1, seq_length))

    if a2 is not None and a3 is not None:
        bounds.append(RegionBounds(Region.S58, a2.env_from, a3.env_to))

    # ITS2: prefer anchor3..anchor4 boundary; infer boundaries from seq edges
    if a3 is not None and a4 is not None:
        start = a3.env_to + 1
        end = a4.env_from - 1
        if start <= end:
            bounds.append(RegionBounds(Region.ITS2, start, end))
    elif a3 is not None and a4 is None and a3.env_to < seq_length:
        bounds.append(RegionBounds(Region.ITS2, a3.env_to + 1, seq_length))
    elif a4 is not None and a3 is None and a2 is None and a1 is None and a4.env_from > 1:
        bounds.append(RegionBounds(Region.ITS2, 1, a4.env_from - 1))

    if a4 is not None:
        bounds.append(RegionBounds(Region.LSU, a4.env_from, seq_length))

    return tuple(bounds)
