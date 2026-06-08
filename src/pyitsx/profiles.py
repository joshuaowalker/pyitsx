import logging
import os
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional, Sequence, Union

import pyhmmer
import pyhmmer.easel
import pyhmmer.hmmer
import pyhmmer.plan7

from pyitsx.constants import AnchorType, Organism, Strand
from pyitsx.models import AnchorHit

logger = logging.getLogger(__name__)

WINDOW_LENGTH = 200

SequenceInput = Union[
    Path,
    pyhmmer.easel.DigitalSequenceBlock,
    Iterable[pyhmmer.easel.DigitalSequence],
    Sequence[str],
    Sequence[tuple[str, str]],
]


class ProfileDB:
    """Loaded HMM profiles for a single organism group."""

    def __init__(
        self,
        hmm_dir: Optional[Path] = None,
        organism: Union[Organism, str] = Organism.F,
    ):
        if isinstance(organism, str):
            organism = Organism.from_code(organism)
        if hmm_dir is None:
            hmm_dir = find_hmm_dir()
        hmm_path = Path(hmm_dir) / f"{organism.name}.hmm"
        if not hmm_path.exists():
            raise FileNotFoundError(f"HMM profile not found: {hmm_path}")
        self._hmms = _load_hmms(hmm_path)
        self._alphabet = pyhmmer.easel.Alphabet.dna()
        self.organism = organism
        logger.info(
            "Loaded %d profiles for organism %s from %s",
            len(self._hmms), organism, hmm_dir,
        )

    @classmethod
    def from_hmm_file(cls, hmm_path: Path) -> "ProfileDB":
        db = cls.__new__(cls)
        db._hmms = _load_hmms(hmm_path)
        db._alphabet = pyhmmer.easel.Alphabet.dna()
        try:
            db.organism = Organism.from_code(hmm_path.stem)
        except ValueError:
            db.organism = hmm_path.stem
        return db

    @property
    def n_profiles(self) -> int:
        return len(self._hmms)

    def prepare(self, sequences: SequenceInput) -> pyhmmer.easel.DigitalSequenceBlock:
        """Normalize input sequences to a DigitalSequenceBlock.

        Accepts:
          - Path to a FASTA/FASTQ file
          - pyhmmer DigitalSequenceBlock (returned as-is)
          - Iterable of pyhmmer DigitalSequence
          - List of Bio.SeqRecord.SeqRecord
          - List of (name, sequence) string tuples
          - List of bare sequence strings (auto-named seq_0, seq_1, ...)
        """
        if isinstance(sequences, pyhmmer.easel.DigitalSequenceBlock):
            return sequences

        if isinstance(sequences, Path):
            return self._load_from_file(sequences)

        items = list(sequences)
        if not items:
            return pyhmmer.easel.DigitalSequenceBlock(self._alphabet)

        first = items[0]

        if isinstance(first, pyhmmer.easel.DigitalSequence):
            block = pyhmmer.easel.DigitalSequenceBlock(self._alphabet)
            for s in items:
                block.append(s)
            return block

        if isinstance(first, str):
            return self._from_strings(items)

        if isinstance(first, tuple):
            return self._from_tuples(items)

        # Try BioPython SeqRecord (duck-type to avoid hard import)
        if hasattr(first, "seq") and hasattr(first, "id"):
            return self._from_seqrecords(items)

        raise TypeError(
            f"Cannot prepare sequences from {type(first).__name__}. "
            "Expected Path, str, (name, seq) tuple, Bio.SeqRecord, "
            "or pyhmmer DigitalSequence."
        )

    def search(
        self,
        sequences: pyhmmer.easel.DigitalSequenceBlock,
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
                    anchor_hit = AnchorHit(
                        anchor_type=anchor_type,
                        strand=Strand(domain.strand),
                        env_from=domain.env_from,
                        env_to=domain.env_to,
                        score=domain.score,
                        evalue=domain.i_evalue,
                        profile_name=profile_name,
                    )
                    hits_by_seq[hit.name].append(anchor_hit)

        return dict(hits_by_seq)

    def load_sequences(
        self, path: Path
    ) -> pyhmmer.easel.DigitalSequenceBlock:
        return self.prepare(path)

    def _load_from_file(self, path: Path) -> pyhmmer.easel.DigitalSequenceBlock:
        with pyhmmer.easel.SequenceFile(
            str(path), digital=True, alphabet=self._alphabet
        ) as sf:
            return sf.read_block()

    def _from_strings(
        self, sequences: list[str]
    ) -> pyhmmer.easel.DigitalSequenceBlock:
        block = pyhmmer.easel.DigitalSequenceBlock(self._alphabet)
        for i, seq in enumerate(sequences):
            ts = pyhmmer.easel.TextSequence(
                name=f"seq_{i}".encode(), sequence=seq
            )
            block.append(ts.digitize(self._alphabet))
        return block

    def _from_tuples(
        self, pairs: list[tuple[str, str]]
    ) -> pyhmmer.easel.DigitalSequenceBlock:
        block = pyhmmer.easel.DigitalSequenceBlock(self._alphabet)
        for name, seq in pairs:
            ts = pyhmmer.easel.TextSequence(
                name=name.encode(), sequence=seq
            )
            block.append(ts.digitize(self._alphabet))
        return block

    def _from_seqrecords(self, records: list) -> pyhmmer.easel.DigitalSequenceBlock:
        block = pyhmmer.easel.DigitalSequenceBlock(self._alphabet)
        for rec in records:
            ts = pyhmmer.easel.TextSequence(
                name=str(rec.id).encode(), sequence=str(rec.seq)
            )
            block.append(ts.digitize(self._alphabet))
        return block


def find_hmm_dir() -> Path:
    """Auto-detect the ITSx HMM profile directory.

    Search order:
      1. PYITSX_HMM_DIR environment variable
      2. ITSx_db/HMMs/ alongside the ITSx executable on PATH
    """
    env_dir = os.environ.get("PYITSX_HMM_DIR")
    if env_dir:
        p = Path(env_dir)
        if p.is_dir():
            return p
        raise FileNotFoundError(
            f"PYITSX_HMM_DIR is set to {env_dir} but directory does not exist"
        )

    itsx_bin = shutil.which("ITSx")
    if itsx_bin:
        candidate = Path(itsx_bin).resolve().parent / "ITSx_db" / "HMMs"
        if candidate.is_dir():
            return candidate

    raise FileNotFoundError(
        "Cannot find ITSx HMM profiles. Either pass hmm_dir explicitly, "
        "set PYITSX_HMM_DIR, or install ITSx so it is on PATH."
    )


def _load_hmms(hmm_path: Path) -> list[pyhmmer.plan7.HMM]:
    with pyhmmer.plan7.HMMFile(str(hmm_path)) as f:
        return list(f)
