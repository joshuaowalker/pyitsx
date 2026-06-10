# pyitsx

Fast ITS region detection and extraction using in-process HMM search.

pyitsx is a Python replacement for [ITSx](https://microbiology.se/software/itsx/) that runs HMMER searches in-process via [pyhmmer](https://pyhmmer.readthedocs.io/), eliminating file I/O overhead and subprocess coordination. It uses the same HMM profiles as ITSx but does not bundle them — an ITSx installation is required.

## Installation

```bash
pip install pyitsx
```

For BioPython SeqRecord support:

```bash
pip install pyitsx[bio]
```

### Requirements

- Python 3.10+
- [ITSx](https://microbiology.se/software/itsx/) installed and on PATH (for HMM profiles), or set `PYITSX_HMM_DIR` to the `HMMs/` directory

## Python API

```python
import pyitsx

db = pyitsx.ProfileDB(organism="F")

# Orient sequences (determine 5'->3' strand)
results = pyitsx.orient("reads.fasta", db)

# Classify ITS region presence
results = pyitsx.classify("reads.fasta", db)

# Delimit all regions with coordinates
results = pyitsx.delimit("reads.fasta", db)

# Extract region sequences
results = pyitsx.extract("reads.fasta", db, regions=[pyitsx.Region.ITS2])

# Extract full ITS (ITS1 + 5.8S + ITS2)
results = pyitsx.extract("reads.fasta", db, regions=[pyitsx.Region.FULL_ITS])

# Score against multiple organism groups
results = pyitsx.score_organisms("reads.fasta")
```

All pipeline functions accept flexible input: file paths, `(name, sequence)` tuples, bare sequence strings, BioPython SeqRecords, or pyhmmer DigitalSequenceBlocks.

## CLI

```bash
pyitsx orient   -i reads.fasta
pyitsx classify -i reads.fasta
pyitsx delimit  -i reads.fasta
pyitsx extract  -i reads.fasta --region ITS2
pyitsx score    -i reads.fasta --organisms F T M
```

Use `--format jsonl` for machine-readable output, `--mode best` for exhaustive (non-short-circuit) search.

## How it works

pyitsx uses the same biological model as ITSx — profile-HMM search against conserved ribosomal flanks — but differs in implementation in a few ways that may be useful to users and other implementers.

**Short-circuit search with adaptive profile ordering.** ITSx organism groups contain many HMM profiles per anchor type (e.g., 538 for fungi). ITSx and ITSxRust search all of them. pyitsx searches profiles one at a time in frequency-ranked order and stops as soon as a confident hit is found (score >= 20, E-value <= 1e-4). An empirical ordering for fungi is built in; the ordering also adapts at runtime as hits accumulate, so even cold-start runs on unfamiliar datasets converge quickly. This is the primary source of the ~300x speedup over ITSx in FAST mode. BEST mode (`--mode best`) disables short-circuit and searches all profiles exhaustively, matching ITSx behavior.

**In-process HMM search.** ITSx invokes `hmmscan` as a subprocess and parses text output files. ITSxRust invokes `nhmmer` as a subprocess and streams its tblout. pyitsx calls HMMER's C library directly via pyhmmer, eliminating process startup, file I/O, and text parsing overhead. This also enables the per-profile short-circuit strategy above, which would be difficult to implement over subprocess boundaries.

**Fixed E-value normalization.** pyitsx passes a fixed database size (`Z=0.001` Mb) to nhmmer so that E-values reflect per-sequence significance and are stable regardless of how many sequences are in the input batch. Without this, E-values shift with total input size, which can cause detection to depend on batch composition.

**Single-anchor chain fallback.** ITSxRust introduced partial chains from two-anchor pairs when a full four-anchor chain cannot be formed. pyitsx extends this to single-anchor fallback: a sequence with one confident anchor (any of the four types) produces a PARTIAL chain with inferred region boundaries, recovering fragments that would otherwise be discarded.

**Use-case oriented API.** The Python API exposes `orient`, `classify`, `delimit`, and `extract` as separate functions with typed result objects, rather than a single monolithic entry point. Pipeline functions accept file paths, bare strings, `(name, sequence)` tuples, BioPython SeqRecords, or pyhmmer DigitalSequenceBlocks.

**Process-level parallelism.** The CLI splits input sequences into chunks and processes them in independent worker processes (`--cpus`). Each worker loads its own HMM profiles and runs single-threaded short-circuit search. This avoids GIL contention and scales near-linearly to 4+ cores. The library API is intentionally single-process — parallelization is left to user code, where workload and resource constraints are understood.

## Caveats

- **Tested on fungi only.** pyitsx supports all ITSx organism groups via `--organism`, but development and testing have focused exclusively on fungal ITS (organism group `F`). Other groups should work — the underlying search and chain logic is organism-agnostic — but have not been validated.
- **FAST mode profile ordering is tuned for fungi.** The short-circuit optimization in FAST mode searches the most frequently matched profiles first. The built-in profile ordering is derived from empirical fungal datasets. For other organism groups, the first search run will be slower as profiles are tried in an unoptimized order; subsequent searches within the same process benefit from learned ordering.
- **HMM profiles are not bundled.** An ITSx installation (or a copy of its `HMMs/` directory) is required. Set `--hmm-dir` or `PYITSX_HMM_DIR` to the profile directory if ITSx is not on PATH.

## Acknowledgments

pyitsx builds on the work of several projects:

**ITSx** established the HMM profile-based approach to ITS boundary detection that pyitsx (and all tools in this space) depends on. pyitsx uses the ITSx HMM profiles directly and would not exist without this foundational work.

> Bengtsson-Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A., De Wit, P., Sánchez-García, M., Ebersberger, I., de Sousa, F., Amend, A., Jumpponen, A., Unterseher, M., Kristiansson, E., Abarenkov, K., & Nilsson, R. H. (2013). Improved software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other eukaryotes for analysis of environmental sequencing data. *Methods in Ecology and Evolution*, 4(10), 914–919. https://doi.org/10.1111/2041-210X.12073

**ITSxRust** introduced the four-anchor chain model with partial-chain fallback and confidence labels that pyitsx's chain-building and region inference logic is based on.

> O'Brien, A., Lagos, C., Fernández, K., Ojeda, B., & Parada, P. (2026). ITSxRust: ITS region extraction with partial-chain recovery and structured diagnostics for long-read amplicon sequencing. *bioRxiv*. https://doi.org/10.64898/2026.02.25.707950

**PyHMMER** provides the in-process Cython bindings to HMMER3 that make pyitsx's performance possible — replacing subprocess calls and file I/O with direct C-level HMM search.

> Larralde, M. & Zeller, G. (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics*, 39(5), btad214. https://doi.org/10.1093/bioinformatics/btad214

## License

MIT
