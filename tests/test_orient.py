import pytest
import pyhmmer.easel

from pyitsx.constants import Strand
from pyitsx.models import OrientResult
from pyitsx.pipeline import orient
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestOrient:

    def test_orient_test_fasta(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = orient(itsx_test_fasta, db, cpus=1)

        detected = [r for r in results if r.strand is not None]
        assert len(detected) > 40
        for r in detected:
            assert isinstance(r, OrientResult)
            assert isinstance(r.strand, Strand)
            assert r.top_score > 0
            assert r.n_hits > 0

    def test_orient_mostly_plus_strand(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = orient(itsx_test_fasta, db, cpus=1)

        detected = [r for r in results if r.strand is not None]
        plus_count = sum(1 for r in detected if r.strand == Strand.PLUS)
        assert plus_count > len(detected) * 0.8

    def test_orient_includes_undetected(self):
        db = ProfileDB(ITSX_DB, organism="F")
        results = orient([("random_seq", "ACGTACGTACGT" * 10)], db, cpus=1)
        assert len(results) == 1
        assert results[0].strand is None
        assert results[0].top_score == 0.0
        assert results[0].n_hits == 0


@requires_hmm_db
class TestOrientConsensus:

    @pytest.fixture(scope="class")
    def orient_results(self, consensus_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        return orient(consensus_fasta, db, cpus=4)

    def test_high_detection_rate(self, orient_results, consensus_fasta):
        from Bio import SeqIO

        n_seqs = sum(1 for _ in SeqIO.parse(consensus_fasta, "fasta"))
        detected = [r for r in orient_results if r.strand is not None]
        assert len(detected) > n_seqs * 0.9

    def test_all_detected_have_strand(self, orient_results):
        detected = [r for r in orient_results if r.strand is not None]
        for r in detected:
            assert r.strand in (Strand.PLUS, Strand.MINUS)
