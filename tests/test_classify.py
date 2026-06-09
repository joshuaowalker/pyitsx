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
        results = classify(itsx_test_fasta, db, cpus=1)

        assert len(results) > 0
        detected = [r for r in results if r.confidence != Confidence.NONE]
        assert len(detected) > 0
        for r in detected:
            assert isinstance(r, ClassifyResult)
            assert isinstance(r.strand, Strand)
            assert isinstance(r.confidence, Confidence)

    def test_most_test_seqs_have_both_its_regions(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = classify(itsx_test_fasta, db, cpus=1)

        both = [r for r in results if r.has_its1 and r.has_its2]
        assert len(both) > 15

    def test_confident_results_exist(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = classify(itsx_test_fasta, db, cpus=1)

        confident = [r for r in results if r.confidence == Confidence.CONFIDENT]
        assert len(confident) > 10

    def test_58s_only_seqs_have_inferred_its(self, itsx_test_fasta):
        """Sequences with only 5.8S anchors should still report has_its1/has_its2
        via boundary inference."""
        db = ProfileDB(ITSX_DB, organism="F")
        results = classify(itsx_test_fasta, db, cpus=1)

        by_id = {r.seq_id: r for r in results}
        # Pyricularia_sp_Br38 has only anchors 2+3 (5.8S)
        pyricularia = [r for r in results if "Pyricularia_sp_Br38" in r.seq_id]
        assert len(pyricularia) == 1
        r = pyricularia[0]
        assert r.has_its1
        assert r.has_its2
        assert r.confidence == Confidence.PARTIAL
