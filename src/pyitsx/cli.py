"""Command-line interface for pyitsx."""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from pyitsx import __version__
from pyitsx.constants import Region
from pyitsx.models import ChainConstraints, DEFAULT_CONSTRAINTS
from pyitsx.pipeline import classify, delimit, orient
from pyitsx.profiles import ProfileDB

logger = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="pyitsx",
        description="Fast ITS region detection and extraction using in-process HMM search",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("-i", "--input", required=True, type=Path, help="Input FASTA file")
    shared.add_argument("-o", "--output", type=Path, help="Output file (default: stdout)")
    shared.add_argument("--hmm-dir", type=Path, default=None, help="Path to ITSx HMM profile directory (auto-detected if omitted)")
    shared.add_argument("--organism", default="F", help="Organism group (default: F for fungi)")
    shared.add_argument("--cpus", type=int, default=0, help="CPU threads (default: all available)")
    shared.add_argument("--format", choices=["tsv", "jsonl"], default="tsv", help="Output format")

    subparsers.add_parser("orient", parents=[shared], help="Determine 5'->3' orientation")
    subparsers.add_parser("classify", parents=[shared], help="Classify ITS region type")

    delimit_parser = subparsers.add_parser("delimit", parents=[shared], help="Full region delimitation")
    delimit_parser.add_argument(
        "--region", choices=["all", "its1", "its2", "full"],
        default="all", help="Which region(s) to report",
    )

    return parser.parse_args(argv)


def setup_logging(debug: bool = False) -> None:
    level = logging.DEBUG if debug else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
        stream=sys.stderr,
    )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    setup_logging()

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    logger.info("pyitsx %s — %s", __version__, args.command)
    t0 = time.perf_counter()

    db = ProfileDB(args.hmm_dir, organism=args.organism)
    seqs = db.load_sequences(args.input)
    logger.info("Loaded %d sequences", len(seqs))

    out = open(args.output, "w") if args.output else sys.stdout

    try:
        if args.command == "orient":
            results = orient(seqs, db, cpus=args.cpus)
            _write_orient(results, out, args.format)
        elif args.command == "classify":
            results = classify(seqs, db, cpus=args.cpus)
            _write_classify(results, out, args.format)
        elif args.command == "delimit":
            results = delimit(seqs, db, cpus=args.cpus)
            _write_delimit(results, out, args.format)
    finally:
        if args.output:
            out.close()

    elapsed = time.perf_counter() - t0
    logger.info(
        "Done: %d results in %.1fs (%.0f seq/s)",
        len(results), elapsed, len(seqs) / elapsed,
    )


def _write_orient(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            json.dump({"seq_id": r.seq_id, "strand": r.strand.value, "top_score": r.top_score, "n_anchors": r.n_anchors}, out)
            out.write("\n")
    else:
        out.write("seq_id\tstrand\ttop_score\tn_anchors\n")
        for r in results:
            out.write(f"{r.seq_id}\t{r.strand.value}\t{r.top_score:.1f}\t{r.n_anchors}\n")


def _write_classify(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            json.dump({
                "seq_id": r.seq_id, "strand": r.strand.value,
                "has_its1": r.has_its1, "has_its2": r.has_its2,
                "confidence": r.confidence.value,
            }, out)
            out.write("\n")
    else:
        out.write("seq_id\tstrand\thas_its1\thas_its2\tconfidence\n")
        for r in results:
            out.write(f"{r.seq_id}\t{r.strand.value}\t{r.has_its1}\t{r.has_its2}\t{r.confidence.value}\n")


def _write_delimit(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            rec = {
                "seq_id": r.seq_id, "seq_length": r.seq_length,
                "strand": r.strand.value, "confidence": r.confidence.value,
            }
            for b in r.bounds:
                rec[b.region.value] = f"{b.start}-{b.end}"
            json.dump(rec, out)
            out.write("\n")
    else:
        out.write("seq_id\tseq_length\tstrand\tconfidence\tSSU\tITS1\t5.8S\tITS2\tLSU\n")
        for r in results:
            regions = {b.region: b for b in r.bounds}
            cols = [r.seq_id, str(r.seq_length), r.strand.value, r.confidence.value]
            for reg in [Region.SSU, Region.ITS1, Region.S58, Region.ITS2, Region.LSU]:
                b = regions.get(reg)
                cols.append(f"{b.start}-{b.end}" if b else "-")
            out.write("\t".join(cols) + "\n")
