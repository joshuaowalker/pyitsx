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
    NONE = "none"


class SearchMode(Enum):
    FAST = "fast"
    BEST = "best"


class Region(Enum):
    SSU = "SSU"
    ITS1 = "ITS1"
    S58 = "5.8S"
    ITS2 = "ITS2"
    LSU = "LSU"
    FULL_ITS = "full_ITS"


class Organism(Enum):
    """ITSx organism group codes."""

    A = "Alveolata"
    B = "Bryophyta"
    C = "Bacillariophyta"
    D = "Amoebozoa"
    E = "Euglenozoa"
    F = "Fungi"
    G = "Chlorophyta"
    H = "Rhodophyta"
    I = "Phaeophyceae"
    L = "Marchantiophyta"
    M = "Metazoa"
    N = "Microsporidia"
    O = "Oomycota"
    P = "Haptophyceae"
    Q = "Raphidophyceae"
    R = "Rhizaria"
    S = "Synurophyceae"
    T = "Tracheophyta"
    U = "Eustigmatophyceae"

    @classmethod
    def from_code(cls, code: str) -> "Organism":
        try:
            return cls[code.upper()]
        except KeyError:
            valid = ", ".join(f"{m.name} ({m.value})" for m in cls)
            raise ValueError(f"Unknown organism code {code!r}. Valid codes: {valid}")
