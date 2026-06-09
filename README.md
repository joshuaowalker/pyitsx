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

## License

MIT
