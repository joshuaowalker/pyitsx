"""Tests for flexible sequence input via ProfileDB.prepare()."""

import pytest
from pathlib import Path
from Bio import SeqIO

from pyitsx.pipeline import orient, delimit
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestPrepareInputTypes:

    @pytest.fixture(scope="class")
    def db(self):
        return ProfileDB(ITSX_DB, organism="F")

    @pytest.fixture(scope="class")
    def reference_seqs(self, db, itsx_test_fasta):
        """Load a few test sequences as (name, sequence) pairs."""
        records = list(SeqIO.parse(itsx_test_fasta, "fasta"))[:5]
        return [(rec.id, str(rec.seq)) for rec in records]

    @pytest.fixture(scope="class")
    def reference_names(self, reference_seqs):
        return {name for name, _ in reference_seqs}

    def test_from_path(self, db, itsx_test_fasta):
        results = orient(itsx_test_fasta, db, cpus=1)
        assert len(results) == 50

    def test_from_digital_block(self, db, itsx_test_fasta):
        block = db.prepare(itsx_test_fasta)
        results = orient(block, db, cpus=1)
        assert len(results) == 50

    def test_from_tuples(self, db, reference_seqs, reference_names):
        results = orient(reference_seqs, db, cpus=1)
        assert len(results) == len(reference_seqs)
        assert {r.seq_id for r in results} == reference_names

    def test_from_bare_strings(self, db, reference_seqs):
        strings = [seq for _, seq in reference_seqs]
        results = orient(strings, db, cpus=1)
        assert len(results) == len(reference_seqs)
        assert all(r.seq_id.startswith("seq_") for r in results)

    def test_from_biopython_records(self, db, itsx_test_fasta, reference_names):
        records = [r for r in SeqIO.parse(itsx_test_fasta, "fasta")
                    if r.id in reference_names]
        results = orient(records, db, cpus=1)
        assert len(results) == len(reference_names)
        assert {r.seq_id for r in results} == reference_names

    def test_empty_input(self, db):
        results = orient([], db, cpus=1)
        assert results == []

    def test_bad_type_raises(self, db):
        with pytest.raises(TypeError, match="Cannot prepare"):
            db.prepare([42])
