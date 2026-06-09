import pytest

from pyitsx.constants import Confidence, Region
from pyitsx.models import ExtractionResult
from pyitsx.pipeline import extract
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestExtract:

    def test_extract_returns_results(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = extract(itsx_test_fasta, db, cpus=1)

        assert len(results) > 0
        for r in results:
            assert isinstance(r, ExtractionResult)
            assert len(r.sequence) > 0
            assert r.end - r.start + 1 == len(r.sequence)

    def test_extract_specific_region(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = extract(itsx_test_fasta, db, regions=[Region.ITS2], cpus=1)

        assert len(results) > 0
        for r in results:
            assert r.region == Region.ITS2

    def test_extract_multiple_regions(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = extract(
            itsx_test_fasta, db,
            regions=[Region.ITS1, Region.ITS2], cpus=1,
        )

        regions_found = {r.region for r in results}
        assert Region.ITS1 in regions_found
        assert Region.ITS2 in regions_found

    def test_extracted_sequences_are_dna(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = extract(itsx_test_fasta, db, regions=[Region.ITS2], cpus=1)

        for r in results:
            assert all(c in "ACGTacgtNn" for c in r.sequence), (
                f"Non-DNA character in {r.seq_id}: {r.sequence[:50]}"
            )

    def test_its_region_lengths_are_reasonable(self, itsx_test_fasta):
        db = ProfileDB(ITSX_DB, organism="F")
        results = extract(itsx_test_fasta, db, cpus=1)

        for r in results:
            if r.region == Region.ITS1:
                assert 50 <= len(r.sequence) <= 1500
            elif r.region == Region.ITS2:
                assert 50 <= len(r.sequence) <= 2000
