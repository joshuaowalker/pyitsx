import logging
import os
import shutil
from collections import Counter, defaultdict
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
DEFAULT_BATCH_SIZE = 1
SHORTCIRCUIT_SCORE = 20.0
SHORTCIRCUIT_EVALUE = 1e-4
# Fixed database size (Mb) for nhmmer E-value calculation.
# Normalizes E-values to per-sequence scale (~1 Kb) so detection
# sensitivity is independent of how many sequences are in the input.
_NHMMER_Z = 0.001

# Empirical profile ordering for fungi, derived from 103K macrofungi dataset.
# Profiles listed first produce confident hits (score >= 20) on the most
# sequences. Used to seed _profile_freq so the first batch gets near-optimal
# short-circuit behavior. Profiles absent from the loaded HMMs are ignored;
# profiles present in HMMs but not listed here start at priority 0.
_DEFAULT_PROFILE_ORDER: dict[str, dict[str, list[str]]] = {
    "F": {
        "SSU_END": [
            "1_SSU_Fungi_Ge_1_SSU_end_onea_long",
            "1_SSU_Fungi_Cr_1_SSU_end_long",
            "1_SSU_Fungi_Tu_1_SSU_end_four_long",
            "1_SSU_Fungi_Ge_1_SSU_end_oneb_long",
            "1_SSU_Fungi_Ca_1_SSU_end_long",
            "1_SSU_Fungi_Ge_1_SSU_end_threeh_long",
            "1_SSU_Fungi_Ge_1_SSU_end_five_long",
            "1_SSU_Fungi_Ge_1_SSU_end_zero_long",
            "1_SSU_Fungi_Ge_1_SSU_end_threej_long",
            "1_SSU_Fungi_Ge_1_SSU_end_nine_long",
        ],
        "S58_START": [
            "2_5.8_Fungi_Ge_58S_start_zero_long",
            "2_5.8_Fungi_Ge_58S_start_threed_long",
            "2_5.8_Fungi_Ca_58S_start_long",
            "2_5.8_Fungi_Ge_58S_start_two_long",
            "2_5.8_Fungi_Ge_58S_start_oned_long",
            "2_5.8_Tedersoo_58S_start_long_one",
        ],
        "S58_END": [
            "3_End_Fungi_Ge_58S_end_fourteenb_long",
            "3_End_Fungi_Ge_58S_end_thirtyfived_long",
            "3_End_Fungi_Ge_58S_end_oneg_long",
            "3_End_Fungi_Ge_58S_end_twentyc_long",
            "3_End_Fungi_Ge_58S_end_thirtythree_long",
            "3_End_Fungi_Ge_58S_end_thirtyfivea_long",
            "3_End_Fungi_Ge_58S_end_thirtyoneb_long",
            "3_End_Fungi_Ge_58S_end_foura_long",
            "3_End_Fungi_Ge_58S_end_thirtyfour_long",
            "3_End_Fungi_Ge_58S_end_thirtytwo_long",
            "3_End_Tedersoo_58S_end_long_ten",
            "3_End_Fungi_Ge_58S_end_threea_long",
            "3_End_Fungi_Ge_58S_end_oned_long",
            "3_End_Fungi_Ge_58S_end_thirtyfivec_long",
            "3_End_Fungi_Ge_58S_end_ninetythree_long",
            "3_End_Fungi_Ge_58S_end_twentysix_long",
            "3_End_Fungi_Ge_58S_end_seventyonea_long",
            "3_End_Fungi_Ge_58S_end_seventyoneb_long",
            "3_end_Calocera_58S_end_long_three",
            "3_End_Fungi_Ge_58S_end_sixtyfive_long",
            "3_End_Fungi_Ge_58S_end_fiftythree_long",
            "3_End_Fungi_Ge_58S_end_nineteen_long",
            "3_End_Fungi_Ge_58S_end_fiftyone_long",
            "3_End_Fungi_Ca_58S_end_one_long",
            "3_End_Fungi_Ge_58S_end_sixtyc_long",
            "3_end_Calocera_58S_end_long_two",
            "3_End_Fungi_Cr_58S_end_one_long",
        ],
        "LSU_START": [
            "4_LSU_Fungi_Ge_4_LSU_start_zeroe_long",
            "4_LSU_Fungi_Ge_4_LSU_start_fortythree_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zerof_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zeroq_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zerok_long",
            "4_LSU_Fungi_Ge_4_LSU_start_sixteen_long",
            "4_LSU_Fungi_Ge_4_LSU_start_fortyfive_long",
            "4_LSU_Fungi_Ge_4_LSU_start_seventeen_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zerog_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zeroj_long",
            "4_LSU_Fungi_Ge_4_LSU_start_fortyfour_long",
            "4_LSU_Fungi_Ge_4_LSU_start_tena_long",
            "4_LSU_Fungi_Ge_4_LSU_start_nine_long",
            "4_LSU_Fungi_Ge_4_LSU_start_elevenc_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zeroc_long",
            "4_LSU_Fungi_Ge_4_LSU_start_twentythreea_long",
            "4_LSU_Fungi_Ge_4_LSU_start_six_long",
            "4_LSU_Fungi_Ge_4_LSU_start_zerob_long",
            "4_LSU_Fungi_Ge_4_LSU_start_seven_long",
            "4_LSU_Fungi_Ge_4_LSU_start_two_long",
            "4_LSU_Fungi_Ge_4_LSU_start_twentyfive_long",
            "4_LSU_Fungi_Cr_4_LSU_start_two_long",
            "4_LSU_Fungi_Ca_4_LSU_start_three_long",
            "4_LSU_Fungi_Tu_4_LSU_start_fourteen_long",
        ],
    },
}

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
        self._profiles_by_anchor = _group_profiles(self._hmms)
        self._profile_freq: Counter = Counter()
        _seed_profile_freq(self._profile_freq, organism.name, self._hmms)
        logger.info(
            "Loaded %d profiles for organism %s from %s",
            len(self._hmms), organism, hmm_dir,
        )

    @classmethod
    def from_hmm_file(cls, hmm_path: Path) -> "ProfileDB":
        db = cls.__new__(cls)
        db._hmms = _load_hmms(hmm_path)
        db._alphabet = pyhmmer.easel.Alphabet.dna()
        db._profiles_by_anchor = _group_profiles(db._hmms)
        db._profile_freq = Counter()
        try:
            _seed_profile_freq(db._profile_freq, hmm_path.stem, db._hmms)
        except (ValueError, KeyError):
            pass
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
        cpus: int = 1,
        batch_size: int = DEFAULT_BATCH_SIZE,
    ) -> dict[str, list[AnchorHit]]:
        if batch_size <= 0:
            return self._search_bulk(sequences, cpus)
        return self._search_batched(sequences, cpus, batch_size)

    def _search_bulk(
        self,
        sequences: pyhmmer.easel.DigitalSequenceBlock,
        cpus: int,
    ) -> dict[str, list[AnchorHit]]:
        hits_by_seq: dict[str, list[AnchorHit]] = defaultdict(list)
        for top_hits in pyhmmer.hmmer.nhmmer(
            self._hmms, sequences, cpus=cpus, window_length=WINDOW_LENGTH,
            Z=_NHMMER_Z,
        ):
            profile_name = top_hits.query.name
            anchor_type = _parse_anchor_type(profile_name)
            if anchor_type is None:
                continue
            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    hits_by_seq[hit.name].append(AnchorHit(
                        anchor_type=anchor_type,
                        strand=Strand(domain.strand),
                        env_from=domain.env_from,
                        env_to=domain.env_to,
                        score=domain.score,
                        evalue=domain.i_evalue,
                        profile_name=profile_name,
                    ))
        return dict(hits_by_seq)

    def _search_batched(
        self,
        sequences: pyhmmer.easel.DigitalSequenceBlock,
        cpus: int,
        batch_size: int,
    ) -> dict[str, list[AnchorHit]]:
        hits_by_seq: dict[str, list[AnchorHit]] = defaultdict(list)
        batches = _make_batches(sequences, batch_size, self._alphabet)

        for anchor_type in AnchorType:
            profiles = self._profiles_by_anchor.get(anchor_type, [])
            if not profiles:
                continue
            profiles = sorted(
                profiles,
                key=lambda p: self._profile_freq.get(p.name, 0),
                reverse=True,
            )
            for block in batches:
                satisfied: set[str] = set()
                n_seqs = len(block)
                for profile in profiles:
                    if len(satisfied) >= n_seqs:
                        break
                    profile_name = profile.name
                    for top_hits in pyhmmer.hmmer.nhmmer(
                        [profile], block, cpus=cpus, window_length=WINDOW_LENGTH,
                        Z=_NHMMER_Z,
                    ):
                        for hit in top_hits:
                            if not hit.included:
                                continue
                            for domain in hit.domains:
                                ah = AnchorHit(
                                    anchor_type=anchor_type,
                                    strand=Strand(domain.strand),
                                    env_from=domain.env_from,
                                    env_to=domain.env_to,
                                    score=domain.score,
                                    evalue=domain.i_evalue,
                                    profile_name=profile_name,
                                )
                                hits_by_seq[hit.name].append(ah)
                                if (ah.score >= SHORTCIRCUIT_SCORE
                                        and ah.evalue <= SHORTCIRCUIT_EVALUE):
                                    satisfied.add(hit.name)
                                    self._profile_freq[profile_name] += 1

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


def _parse_anchor_type(profile_name: str) -> Optional[AnchorType]:
    prefix = profile_name.split("_")[0]
    try:
        return AnchorType(int(prefix))
    except (ValueError, KeyError):
        return None


def _group_profiles(
    hmms: list[pyhmmer.plan7.HMM],
) -> dict[AnchorType, list[pyhmmer.plan7.HMM]]:
    groups: dict[AnchorType, list[pyhmmer.plan7.HMM]] = defaultdict(list)
    for hmm in hmms:
        anchor_type = _parse_anchor_type(hmm.name)
        if anchor_type is not None:
            groups[anchor_type].append(hmm)
    return dict(groups)


def _make_batches(
    sequences: pyhmmer.easel.DigitalSequenceBlock,
    batch_size: int,
    alphabet: pyhmmer.easel.Alphabet,
) -> list[pyhmmer.easel.DigitalSequenceBlock]:
    batches = []
    block = pyhmmer.easel.DigitalSequenceBlock(alphabet)
    for s in sequences:
        block.append(s)
        if len(block) >= batch_size:
            batches.append(block)
            block = pyhmmer.easel.DigitalSequenceBlock(alphabet)
    if len(block) > 0:
        batches.append(block)
    return batches


def _seed_profile_freq(
    freq: Counter,
    organism_code: str,
    hmms: list[pyhmmer.plan7.HMM],
) -> None:
    """Seed profile frequency counter from empirical ordering.

    Assigns decreasing synthetic counts so profiles sort in the
    recommended order for short-circuit. Only profiles actually present
    in the loaded HMMs are seeded; unknown names are silently skipped.
    """
    ordering = _DEFAULT_PROFILE_ORDER.get(organism_code)
    if ordering is None:
        return
    loaded_names = {hmm.name for hmm in hmms}
    for anchor_names in ordering.values():
        n = len(anchor_names)
        for i, name in enumerate(anchor_names):
            if name in loaded_names:
                freq[name] = n - i


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
