import logging
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Union

import pyhmmer
import pyhmmer.easel
import pyhmmer.hmmer
import pyhmmer.plan7

from pyitsx.constants import AnchorType, Strand
from pyitsx.models import AnchorHit

logger = logging.getLogger(__name__)

WINDOW_LENGTH = 200


class ProfileDB:
    """Loaded HMM profiles for a single organism group."""

    def __init__(self, hmm_dir: Path, organism: str = "F"):
        hmm_path = hmm_dir / f"{organism}.hmm"
        if not hmm_path.exists():
            raise FileNotFoundError(f"HMM profile not found: {hmm_path}")
        self._hmms = _load_hmms(hmm_path)
        self._alphabet = pyhmmer.easel.Alphabet.dna()
        self.organism = organism
        logger.info(
            "Loaded %d profiles for organism %s", len(self._hmms), organism
        )

    @classmethod
    def from_hmm_file(cls, hmm_path: Path) -> "ProfileDB":
        db = cls.__new__(cls)
        db._hmms = _load_hmms(hmm_path)
        db._alphabet = pyhmmer.easel.Alphabet.dna()
        db.organism = hmm_path.stem
        return db

    @property
    def n_profiles(self) -> int:
        return len(self._hmms)

    def search(
        self,
        sequences: Union[
            pyhmmer.easel.DigitalSequenceBlock,
            Iterable[pyhmmer.easel.DigitalSequence],
        ],
        cpus: int = 0,
    ) -> dict[str, list[AnchorHit]]:
        hits_by_seq: dict[str, list[AnchorHit]] = defaultdict(list)

        for top_hits in pyhmmer.hmmer.nhmmer(
            self._hmms, sequences, cpus=cpus, window_length=WINDOW_LENGTH
        ):
            profile_name = top_hits.query.name
            anchor_prefix = profile_name.split("_")[0]
            try:
                anchor_type = AnchorType.from_profile_prefix(anchor_prefix)
            except (ValueError, KeyError):
                continue

            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    seq_name = hit.name
                    anchor_hit = AnchorHit(
                        anchor_type=anchor_type,
                        strand=Strand(domain.strand),
                        env_from=domain.env_from,
                        env_to=domain.env_to,
                        score=domain.score,
                        evalue=domain.i_evalue,
                        profile_name=profile_name,
                    )
                    hits_by_seq[seq_name].append(anchor_hit)

        return dict(hits_by_seq)

    def load_sequences(
        self, path: Path
    ) -> pyhmmer.easel.DigitalSequenceBlock:
        with pyhmmer.easel.SequenceFile(
            str(path), digital=True, alphabet=self._alphabet
        ) as sf:
            return sf.read_block()


def _load_hmms(hmm_path: Path) -> list[pyhmmer.plan7.HMM]:
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as f:
        return list(f)
