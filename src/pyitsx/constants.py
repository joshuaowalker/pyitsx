from enum import Enum


class AnchorType(Enum):
    SSU_END = 1
    S58_START = 2
    S58_END = 3
    LSU_START = 4

    @classmethod
    def from_profile_prefix(cls, prefix: str) -> "AnchorType":
        return cls(int(prefix))


class Strand(Enum):
    PLUS = "+"
    MINUS = "-"


class Confidence(Enum):
    CONFIDENT = "confident"
    AMBIGUOUS = "ambiguous"
    PARTIAL = "partial"


class Region(Enum):
    SSU = "SSU"
    ITS1 = "ITS1"
    S58 = "5.8S"
    ITS2 = "ITS2"
    LSU = "LSU"
