import pytest
from collections import Counter

from pyitsx.constants import AnchorType, Strand
from pyitsx.profiles import ProfileDB
from tests.conftest import requires_hmm_db, ITSX_DB


@requires_hmm_db
class TestProfileDB:

    def test_load_profiles(self):
        db = ProfileDB(ITSX_DB, organism="F")
        assert db.n_profiles == 538
        assert db.organism == "F"

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
