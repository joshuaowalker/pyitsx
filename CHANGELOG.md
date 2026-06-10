# Changelog

## 0.1.0 (2026-06-09)

Initial release.

- Orient, classify, delimit, and extract pipeline functions with typed result objects
- FAST mode with short-circuit search and adaptive profile ordering (~300x faster than ITSx)
- BEST mode for exhaustive search matching ITSx behavior
- Four-anchor chain model with partial-chain (2-anchor and single-anchor) fallback
- Chimera detection (cross-strand and out-of-order anchors)
- Full ITS extraction (ITS1 + 5.8S + ITS2 as a single region)
- Multi-organism scoring via `score_organisms()`
- Flexible input: file paths, strings, tuples, BioPython SeqRecords, pyhmmer blocks
- CLI with TSV and JSONL output formats
- Multiprocessing parallelism in CLI (`--cpus`)
- Fixed E-value normalization independent of input batch size
- Empirical profile ordering for fungi (organism group F)
