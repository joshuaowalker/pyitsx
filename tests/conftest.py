"""Shared test fixtures for pyitsx tests."""

import os
import pytest
from pathlib import Path


DATA_DIR = Path(os.environ.get("PYITSX_TEST_DATA", Path.home() / "mm" / "data"))
ITSX_DB = Path(os.environ.get(
    "PYITSX_HMM_DB",
    Path(__file__).parent.parent / "ITSx_1.1.3" / "ITSx_db" / "HMMs",
))


def _has_test_data():
    return DATA_DIR.exists()


def _has_hmm_db():
    return ITSX_DB.exists() and (ITSX_DB / "F.hmm").exists()


requires_test_data = pytest.mark.skipif(
    not _has_test_data(), reason=f"Test data not found at {DATA_DIR}"
)
requires_hmm_db = pytest.mark.skipif(
    not _has_hmm_db(), reason=f"HMM database not found at {ITSX_DB}"
)


@pytest.fixture(scope="session")
def hmm_db_path() -> Path:
    """Path to ITSx HMM profile directory."""
    if not _has_hmm_db():
        pytest.skip(f"HMM database not found at {ITSX_DB}")
    return ITSX_DB


@pytest.fixture(scope="session")
def fungi_hmm_path(hmm_db_path) -> Path:
    """Path to the fungi (F) HMM profile file."""
    return hmm_db_path / "F.hmm"


@pytest.fixture(scope="session")
def consensus_fasta() -> Path:
    """Path to consensus ITS sequences (verified.fasta from ont98)."""
    path = DATA_DIR / "ont98" / "data" / "verified.fasta"
    if not path.exists():
        pytest.skip(f"Consensus test data not found at {path}")
    return path


@pytest.fixture(scope="session")
def nanopore_fastq() -> Path:
    """Path to raw nanopore ITS reads (all25k.fastq from ont98)."""
    path = DATA_DIR / "ont98" / "scale-test" / "all25k.fastq"
    if not path.exists():
        pytest.skip(f"Nanopore test data not found at {path}")
    return path


@pytest.fixture(scope="session")
def non_its_data_dir() -> Path:
    """Path to non-ITS data directory (ont_chants RPB2/TEF1 reads)."""
    path = DATA_DIR / "ont_chants" / "demux0130"
    if not path.exists():
        pytest.skip(f"Non-ITS test data not found at {path}")
    return path


@pytest.fixture(scope="session")
def itsx_test_fasta() -> Path:
    """Path to the small test.fasta bundled with ITSx."""
    path = Path(__file__).parent.parent / "ITSx_1.1.3" / "test.fasta"
    if not path.exists():
        pytest.skip(f"ITSx test.fasta not found at {path}")
    return path
