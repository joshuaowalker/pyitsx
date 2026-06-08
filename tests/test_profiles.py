import pytest
from collections import Counter

from pyitsx.constants import AnchorType, Organism, Strand
from pyitsx.profiles import ProfileDB, find_hmm_dir
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestProfileDB:

    def test_load_profiles(self):
        db = ProfileDB(ITSX_DB, organism="F")
        assert db.n_profiles == 538
        assert db.organism == Organism.F

    def test_search_returns_hits(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        hits_by_seq = db.search(seqs, cpus=1)
        assert len(hits_by_seq) > 40

    def test_all_anchor_types_found(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        hits_by_seq = db.search(seqs, cpus=1)

        all_types = set()
        for hits in hits_by_seq.values():
            for h in hits:
                all_types.add(h.anchor_type)

        assert all_types == {
            AnchorType.SSU_END,
            AnchorType.S58_START,
            AnchorType.S58_END,
            AnchorType.LSU_START,
        }

    def test_hit_fields_populated(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        hits_by_seq = db.search(seqs, cpus=1)

        first_hits = next(iter(hits_by_seq.values()))
        h = first_hits[0]
        assert isinstance(h.anchor_type, AnchorType)
        assert isinstance(h.strand, Strand)
        assert h.env_from >= 1
        assert h.env_to > h.env_from
        assert h.score > 0
        assert h.evalue >= 0
        assert len(h.profile_name) > 0

    def test_strand_detection(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        hits_by_seq = db.search(seqs, cpus=1)

        strand_counts = Counter()
        for hits in hits_by_seq.values():
            for h in hits:
                strand_counts[h.strand] += 1

        assert strand_counts[Strand.PLUS] > 0


class TestFindHmmDir:

    def test_env_var_override(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PYITSX_HMM_DIR", str(tmp_path))
        assert find_hmm_dir() == tmp_path

    def test_env_var_missing_dir_raises(self, monkeypatch):
        monkeypatch.setenv("PYITSX_HMM_DIR", "/nonexistent/path")
        with pytest.raises(FileNotFoundError, match="PYITSX_HMM_DIR"):
            find_hmm_dir()

    def test_itsx_on_path(self, tmp_path, monkeypatch):
        monkeypatch.delenv("PYITSX_HMM_DIR", raising=False)
        hmm_dir = tmp_path / "ITSx_db" / "HMMs"
        hmm_dir.mkdir(parents=True)
        fake_itsx = tmp_path / "ITSx"
        fake_itsx.write_text("#!/bin/sh\n")
        fake_itsx.chmod(0o755)
        monkeypatch.setenv("PATH", str(tmp_path))
        assert find_hmm_dir() == hmm_dir

    def test_nothing_found_raises(self, monkeypatch):
        monkeypatch.delenv("PYITSX_HMM_DIR", raising=False)
        monkeypatch.setenv("PATH", "/nonexistent")
        with pytest.raises(FileNotFoundError, match="Cannot find ITSx HMM profiles"):
            find_hmm_dir()

    def test_auto_detect_from_real_itsx(self):
        """If ITSx is actually installed, find_hmm_dir should locate it."""
        import shutil
        if not shutil.which("ITSx"):
            pytest.skip("ITSx not installed")
        result = find_hmm_dir()
        assert result.is_dir()
        assert (result / "F.hmm").exists()

    def test_profiledb_auto_detect(self):
        """ProfileDB() with no hmm_dir should auto-detect when ITSx is installed."""
        import shutil
        if not shutil.which("ITSx"):
            pytest.skip("ITSx not installed")
        db = ProfileDB(organism="F")
        assert db.n_profiles > 0
