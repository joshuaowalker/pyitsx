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
        seqs = db.load_sequences(itsx_test_fasta)
        results = orient(seqs, db, cpus=1)

        assert len(results) > 40
        for r in results:
            assert isinstance(r, OrientResult)
            assert isinstance(r.strand, Strand)
            assert r.top_score > 0
            assert r.n_anchors > 0

    def test_orient_mostly_plus_strand(self, itsx_test_fasta):
        """ITSx test.fasta sequences are mostly in forward orientation."""
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        results = orient(seqs, db, cpus=1)

        plus_count = sum(1 for r in results if r.strand == Strand.PLUS)
        assert plus_count > len(results) * 0.8

    def test_orient_returns_nothing_for_no_hits(self):
        """Sequences with no anchor hits should not appear in results."""
        db = ProfileDB(ITSX_DB, organism="F")
        alphabet = db._alphabet
        seq = pyhmmer.easel.TextSequence(
            name=b"random_seq", sequence="ACGTACGTACGT" * 10
        ).digitize(alphabet)

        results = orient([seq], db, cpus=1)
        assert len(results) == 0


@requires_hmm_db
class TestOrientConsensus:

    @pytest.fixture(scope="class")
    def orient_results(self, consensus_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(consensus_fasta)
        return orient(seqs, db, cpus=4)

    def test_high_detection_rate(self, orient_results, consensus_fasta):
        from Bio import SeqIO

        n_seqs = sum(1 for _ in SeqIO.parse(consensus_fasta, "fasta"))
        assert len(orient_results) > n_seqs * 0.9

    def test_all_have_strand(self, orient_results):
        for r in orient_results:
            assert r.strand in (Strand.PLUS, Strand.MINUS)
