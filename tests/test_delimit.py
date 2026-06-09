import pytest

from pyitsx.constants import Confidence, Region, Strand
from pyitsx.models import DelimitResult, RegionBounds
from pyitsx.pipeline import delimit
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestDelimit:

    def test_delimit_test_fasta(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        assert len(results) > 0
        detected = [r for r in results if r.confidence != Confidence.NONE]
        assert len(detected) > 0
        for r in detected:
            assert isinstance(r, DelimitResult)
            assert r.seq_length > 0
            assert len(r.bounds) > 0

    def test_full_chains_have_five_regions(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        full_results = [r for r in results if r.chain and r.chain.is_full]
        assert len(full_results) > 10

        for r in full_results:
            regions = {b.region for b in r.bounds}
            assert Region.SSU in regions
            assert Region.ITS1 in regions
            assert Region.S58 in regions
            assert Region.ITS2 in regions
            assert Region.LSU in regions

    def test_region_coordinates_are_ordered(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        for r in results:
            if r.chain and r.chain.is_full:
                sorted_bounds = sorted(r.bounds, key=lambda b: b.start)
                for i in range(len(sorted_bounds) - 1):
                    assert sorted_bounds[i].end < sorted_bounds[i + 1].start

    def test_its_regions_have_reasonable_lengths(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        for r in results:
            its1 = r.get_region(Region.ITS1)
            its2 = r.get_region(Region.ITS2)
            if its1:
                assert 50 <= its1.length <= 1500
            if its2:
                assert 50 <= its2.length <= 2000

    def test_get_region_returns_none_for_missing(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        partial_results = [r for r in results if r.chain and not r.chain.is_full]
        if partial_results:
            r = partial_results[0]
            missing = [
                reg for reg in Region
                if r.get_region(reg) is None
            ]
            assert len(missing) > 0

    def test_regions_dict_property(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        full_results = [r for r in results if r.chain and r.chain.is_full]
        assert len(full_results) > 0
        r = full_results[0]

        regions = r.regions
        assert isinstance(regions, dict)
        assert Region.ITS1 in regions
        assert Region.ITS2 in regions
        assert isinstance(regions[Region.ITS1], RegionBounds)
        assert regions[Region.ITS1] == r.get_region(Region.ITS1)

    def test_full_its_property(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        full_results = [r for r in results if r.chain and r.chain.is_full]
        assert len(full_results) > 0
        for r in full_results:
            full = r.full_its
            assert full is not None
            assert full.region == Region.FULL_ITS
            its1 = r.get_region(Region.ITS1)
            its2 = r.get_region(Region.ITS2)
            assert full.start == its1.start
            assert full.end == its2.end

    def test_full_its_none_when_missing_region(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = delimit(itsx_test_fasta, db, cpus=1)

        undetected = [r for r in results if r.confidence == Confidence.NONE]
        for r in undetected:
            assert r.full_its is None

    def test_undetected_sequences_included(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        seqs = db.prepare(itsx_test_fasta)
        n_input = len(seqs)

        results = delimit(seqs, db, cpus=1)
        assert len(results) == n_input

        detected = [r for r in results if r.confidence != Confidence.NONE]
        undetected = [r for r in results if r.confidence == Confidence.NONE]

        assert len(detected) + len(undetected) == n_input
        for r in undetected:
            assert r.strand is None
            assert r.chain is None
            assert r.bounds == ()
            assert len(r.regions) == 0
