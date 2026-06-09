import pytest

from pyitsx.constants import Organism
from pyitsx.models import OrganismResult, OrganismScore
from pyitsx.scoring import score_organisms
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestScoreOrganisms:

    def test_returns_results_for_all_sequences(self, itsx_test_fasta):
        from pyitsx.profiles import ProfileDB
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.prepare(itsx_test_fasta)
        n_input = len(seqs)

        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.F], cpus=1,
        )
        assert len(results) == n_input
        for r in results:
            assert isinstance(r, OrganismResult)

    def test_best_organism_is_fungi_for_fungal_seqs(self, itsx_test_fasta):
        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.F, Organism.T], cpus=1,
        )
        detected = [r for r in results if r.best is not None]
        assert len(detected) > 0
        for r in detected:
            assert r.best.organism == Organism.F

    def test_scores_sorted_descending(self, itsx_test_fasta):
        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.F, Organism.T, Organism.M], cpus=1,
        )
        for r in results:
            if len(r.scores) > 1:
                for i in range(len(r.scores) - 1):
                    assert r.scores[i].total_score >= r.scores[i + 1].total_score

    def test_organism_filter(self, itsx_test_fasta):
        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.F], cpus=1,
        )
        for r in results:
            for s in r.scores:
                assert s.organism == Organism.F

    def test_organism_score_fields(self, itsx_test_fasta):
        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.F], cpus=1,
        )
        detected = [r for r in results if r.best is not None]
        assert len(detected) > 0
        for r in detected:
            assert isinstance(r.best, OrganismScore)
            assert 1 <= r.best.n_anchors <= 4
            assert r.best.total_score > 0
            assert r.best.best_evalue > 0

    def test_empty_hmm_skipped(self, itsx_test_fasta):
        results = score_organisms(
            itsx_test_fasta, hmm_dir=ITSX_DB,
            organisms=[Organism.N, Organism.F], cpus=1,
        )
        for r in results:
            for s in r.scores:
                assert s.organism != Organism.N
