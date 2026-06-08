import pytest

from pyitsx.constants import Confidence, Strand
from pyitsx.models import ClassifyResult
from pyitsx.pipeline import classify
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestClassify:

    def test_classify_test_fasta(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        results = classify(seqs, db, cpus=1)

        assert len(results) > 0
        for r in results:
            assert isinstance(r, ClassifyResult)
            assert isinstance(r.strand, Strand)
            assert isinstance(r.confidence, Confidence)

    def test_most_test_seqs_have_both_its_regions(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        results = classify(seqs, db, cpus=1)

        both = [r for r in results if r.has_its1 and r.has_its2]
        assert len(both) > 15

    def test_confident_results_exist(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.load_sequences(itsx_test_fasta)
        results = classify(seqs, db, cpus=1)

        confident = [r for r in results if r.confidence == Confidence.CONFIDENT]
        assert len(confident) > 10
