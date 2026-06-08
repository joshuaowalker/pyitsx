from dataclasses import dataclass
from typing import Optional

from pyitsx.constants import AnchorType, Confidence, Region, Strand


@dataclass(frozen=True)
class AnchorHit:
    anchor_type: AnchorType
    strand: Strand
    env_from: int
    env_to: int
    score: float
    evalue: float
    profile_name: str


@dataclass(frozen=True)
class OrientResult:
    seq_id: str
    strand: Strand
    top_score: float
    n_anchors: int


@dataclass(frozen=True)
class ChainConstraints:
    min_its1: int = 50
    max_its1: int = 1500
    min_its2: int = 50
    max_its2: int = 2000
    min_full: int = 150
    max_full: int = 4000
    min_anchor_score: float = 20.0
    max_anchor_evalue: float = 1e-4


DEFAULT_CONSTRAINTS = ChainConstraints()


@dataclass(frozen=True)
class AnchorChain:
    anchors: tuple[Optional[AnchorHit], ...]
    strand: Strand
    total_score: float
    confidence: Confidence

    @property
    def is_full(self) -> bool:
        return all(a is not None for a in self.anchors)

    @property
    def n_anchors(self) -> int:
        return sum(1 for a in self.anchors if a is not None)


@dataclass(frozen=True)
class RegionBounds:
    region: Region
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class ClassifyResult:
    seq_id: str
    strand: Strand
    has_its1: bool
    has_its2: bool
    confidence: Confidence
    chain: AnchorChain


@dataclass(frozen=True)
class DelimitResult:
    seq_id: str
    seq_length: int
    strand: Strand
    chain: AnchorChain
    bounds: tuple[RegionBounds, ...]
    confidence: Confidence

    def get_region(self, region: Region) -> Optional[RegionBounds]:
        for b in self.bounds:
            if b.region == region:
                return b
        return None
