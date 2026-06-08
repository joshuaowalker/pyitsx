"""Smoke tests to verify pyhmmer can load ITSx profiles and search sequences."""

import pytest
import pyhmmer
from collections import defaultdict
from tests.conftest import requires_hmm_db


class TestPyhmmerSmoke:
    """Verify pyhmmer integration with ITSx HMM profiles."""

    @requires_hmm_db
    def test_load_fungi_profiles(self, fungi_hmm_path):
        with pyhmmer.plan7.HMMFile(str(fungi_hmm_path)) as f:
            hmms = list(f)

        assert len(hmms) == 538
        assert hmms[0].alphabet.type == "DNA"

    @requires_hmm_db
    def test_anchor_types_present(self, fungi_hmm_path):
        with pyhmmer.plan7.HMMFile(str(fungi_hmm_path)) as f:
            hmms = list(f)

        prefixes = set()
        for h in hmms:
            parts = h.name.split("_")
            prefixes.add(parts[0])

        assert "1" in prefixes, "Missing SSU (1_SSU) profiles"
        assert "2" in prefixes, "Missing 5.8S start (2_5.8) profiles"
        assert "3" in prefixes, "Missing 5.8S end (3_End) profiles"
        assert "4" in prefixes, "Missing LSU (4_LSU) profiles"

    @requires_hmm_db
    def test_hmmsearch_on_test_fasta(self, fungi_hmm_path, itsx_test_fasta):
        with pyhmmer.plan7.HMMFile(str(fungi_hmm_path)) as f:
            hmms = list(f)

        alphabet = pyhmmer.easel.Alphabet.dna()
        with pyhmmer.easel.SequenceFile(
            str(itsx_test_fasta), digital=True, alphabet=alphabet
        ) as sf:
            seqs = sf.read_block()

        assert len(seqs) == 50

        hit_count = 0
        seqs_with_hits = set()
        for top_hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=1):
            for hit in top_hits:
                if hit.included:
                    hit_count += 1
                    seqs_with_hits.add(hit.name)

        assert hit_count > 0, "No hits found"
        assert len(seqs_with_hits) > 40, f"Expected most sequences to have hits, got {len(seqs_with_hits)}/50"

    @requires_hmm_db
    def test_all_four_anchors_found(self, fungi_hmm_path, itsx_test_fasta):
        """Verify that a known ITS sequence produces hits for all 4 anchor types."""
        with pyhmmer.plan7.HMMFile(str(fungi_hmm_path)) as f:
            hmms = list(f)

        alphabet = pyhmmer.easel.Alphabet.dna()
        with pyhmmer.easel.SequenceFile(
            str(itsx_test_fasta), digital=True, alphabet=alphabet
        ) as sf:
            seqs = sf.read_block()

        # Collect anchor types hit per sequence
        anchors_by_seq = defaultdict(set)
        for top_hits in pyhmmer.hmmer.hmmsearch(hmms, seqs, cpus=1):
            query_name = top_hits.query.name
            anchor_type = query_name.split("_")[0]
            for hit in top_hits:
                if hit.included:
                    anchors_by_seq[hit.name].add(anchor_type)

        # At least some sequences should have all 4 anchor types
        full_chain_seqs = [
            name for name, anchors in anchors_by_seq.items()
            if {"1", "2", "3", "4"}.issubset(anchors)
        ]
        assert len(full_chain_seqs) > 0, "No sequences found with all 4 anchor types"
