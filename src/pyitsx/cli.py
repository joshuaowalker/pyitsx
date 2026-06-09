"""Command-line interface for pyitsx."""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

from pyitsx import __version__
from pyitsx.constants import Organism, Region, SearchMode
from pyitsx.models import ChainConstraints, DEFAULT_CONSTRAINTS
from pyitsx.pipeline import classify, delimit, extract, orient
from pyitsx.profiles import ProfileDB
from pyitsx.scoring import score_organisms

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
    shared.add_argument(
        "--organism", default="F",
        choices=[o.name for o in Organism],
        help="Organism group code (default: F for fungi)",
    )
    shared.add_argument("--cpus", type=int, default=0, help="CPU threads for HMM search (default: all available)")
    shared.add_argument("--batch-size", type=int, default=1, help="Sequences per search batch (default: 1)")
    shared.add_argument("--format", choices=["tsv", "jsonl"], default="tsv", help="Output format")
    shared.add_argument("--mode", choices=["fast", "best"], default="fast", help="Search mode: fast (short-circuit) or best (exhaustive)")

    subparsers.add_parser("orient", parents=[shared], help="Determine 5'->3' orientation")
    subparsers.add_parser("classify", parents=[shared], help="Classify ITS region type")

    subparsers.add_parser("delimit", parents=[shared], help="Full region delimitation")

    extract_parser = subparsers.add_parser("extract", parents=[shared], help="Extract ITS region sequences")
    extract_parser.add_argument(
        "--region", choices=["SSU", "ITS1", "5.8S", "ITS2", "LSU", "full_ITS"],
        action="append", dest="regions", default=None,
        help="Region(s) to extract (default: all detected). Can be repeated.",
    )

    score_parser = subparsers.add_parser("score", parents=[shared], help="Score sequences against multiple organism groups")
    score_parser.add_argument(
        "--organisms", nargs="+", default=None,
        choices=[o.name for o in Organism],
        help="Organism group codes to score against (default: all)",
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

    mode = SearchMode(args.mode)
    out = open(args.output, "w") if args.output else sys.stdout

    try:
        if args.command == "score":
            org_enums = [Organism[o] for o in args.organisms] if args.organisms else None
            results = score_organisms(
                args.input, hmm_dir=args.hmm_dir, organisms=org_enums,
                cpus=args.cpus, batch_size=args.batch_size, mode=mode,
            )
            _write_score(results, out, args.format)
        else:
            db = ProfileDB(args.hmm_dir, organism=args.organism)
            seqs = db.load_sequences(args.input)
            logger.info("Loaded %d sequences", len(seqs))

            if args.command == "orient":
                results = orient(seqs, db, cpus=args.cpus, batch_size=args.batch_size, mode=mode)
                _write_orient(results, out, args.format)
            elif args.command == "classify":
                results = classify(seqs, db, cpus=args.cpus, batch_size=args.batch_size, mode=mode)
                _write_classify(results, out, args.format)
            elif args.command == "delimit":
                results = delimit(seqs, db, cpus=args.cpus, batch_size=args.batch_size, mode=mode)
                _write_delimit(results, out, args.format)
            elif args.command == "extract":
                region_enums = [Region(r) for r in args.regions] if args.regions else None
                results = extract(
                    seqs, db, regions=region_enums,
                    cpus=args.cpus, batch_size=args.batch_size, mode=mode,
                )
                _write_extract(results, out, args.format)
    finally:
        if args.output:
            out.close()

    elapsed = time.perf_counter() - t0
    logger.info(
        "Done: %d results in %.1fs (%.0f seq/s)",
        len(results), elapsed, len(results) / elapsed if elapsed > 0 else 0,
    )


def _write_orient(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            json.dump({"seq_id": r.seq_id, "strand": r.strand.value if r.strand else "-", "top_score": r.top_score, "n_anchors": r.n_anchors, "chimeric": r.chimeric}, out)
            out.write("\n")
    else:
        out.write("seq_id\tstrand\ttop_score\tn_anchors\tchimeric\n")
        for r in results:
            strand = r.strand.value if r.strand else "-"
            out.write(f"{r.seq_id}\t{strand}\t{r.top_score:.1f}\t{r.n_anchors}\t{r.chimeric}\n")


def _write_classify(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            json.dump({
                "seq_id": r.seq_id,
                "strand": r.strand.value if r.strand else "-",
                "has_its1": r.has_its1, "has_its2": r.has_its2,
                "confidence": r.confidence.value,
                "chimeric": r.chimeric,
            }, out)
            out.write("\n")
    else:
        out.write("seq_id\tstrand\thas_its1\thas_its2\tconfidence\tchimeric\n")
        for r in results:
            strand = r.strand.value if r.strand else "-"
            out.write(f"{r.seq_id}\t{strand}\t{r.has_its1}\t{r.has_its2}\t{r.confidence.value}\t{r.chimeric}\n")


def _write_delimit(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            rec = {
                "seq_id": r.seq_id, "seq_length": r.seq_length,
                "strand": r.strand.value if r.strand else "-",
                "confidence": r.confidence.value,
                "chimeric": r.chimeric,
            }
            for b in r.bounds:
                rec[b.region.value] = f"{b.start}-{b.end}"
            full = r.full_its
            rec["full_ITS"] = f"{full.start}-{full.end}" if full else "-"
            json.dump(rec, out)
            out.write("\n")
    else:
        out.write("seq_id\tseq_length\tstrand\tconfidence\tchimeric\tSSU\tITS1\t5.8S\tITS2\tLSU\tfull_ITS\n")
        for r in results:
            strand = r.strand.value if r.strand else "-"
            cols = [r.seq_id, str(r.seq_length), strand, r.confidence.value, str(r.chimeric)]
            for reg in [Region.SSU, Region.ITS1, Region.S58, Region.ITS2, Region.LSU]:
                b = r.regions.get(reg)
                cols.append(f"{b.start}-{b.end}" if b else "-")
            full = r.full_its
            cols.append(f"{full.start}-{full.end}" if full else "-")
            out.write("\t".join(cols) + "\n")


def _write_extract(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            json.dump({
                "seq_id": r.seq_id, "region": r.region.value,
                "start": r.start, "end": r.end,
                "sequence": r.sequence,
            }, out)
            out.write("\n")
    else:
        for r in results:
            out.write(f">{r.seq_id}|{r.region.value}\n{r.sequence}\n")


def _write_score(results, out, fmt):
    if fmt == "jsonl":
        for r in results:
            rec = {"seq_id": r.seq_id}
            rec["best_organism"] = r.best.organism.name if r.best else "-"
            rec["best_score"] = round(r.best.total_score, 1) if r.best else 0.0
            rec["scores"] = [
                {"organism": s.organism.name, "total_score": round(s.total_score, 1),
                 "n_anchors": s.n_anchors, "best_evalue": s.best_evalue}
                for s in r.scores
            ]
            json.dump(rec, out)
            out.write("\n")
    else:
        out.write("seq_id\tbest_organism\tbest_score\tbest_n_anchors\tscores\n")
        for r in results:
            if r.best:
                best_org = r.best.organism.name
                best_score = f"{r.best.total_score:.1f}"
                best_n = str(r.best.n_anchors)
            else:
                best_org = "-"
                best_score = "0.0"
                best_n = "0"
            score_parts = [
                f"{s.organism.name}:{s.total_score:.1f}:{s.n_anchors}"
                for s in r.scores
            ]
            out.write(f"{r.seq_id}\t{best_org}\t{best_score}\t{best_n}\t{';'.join(score_parts)}\n")
