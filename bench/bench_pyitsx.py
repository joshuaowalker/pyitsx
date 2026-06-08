#!/usr/bin/env python3
"""Benchmark pyitsx against ITSx and ITSxRust on real data."""

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO


def bench_pyitsx(input_path: Path, hmm_dir: Path, cpus: int, max_seqs: int = 0):
    from pyitsx.pipeline import delimit
    from pyitsx.profiles import ProfileDB

    db = ProfileDB(hmm_dir, organism="F")
    seqs = db.load_sequences(input_path)
    if max_seqs:
        seqs = list(seqs)[:max_seqs]
    n = len(seqs)

    t0 = time.perf_counter()
    results = delimit(seqs, db, cpus=cpus)
    elapsed = time.perf_counter() - t0

    full = sum(1 for r in results if r.chain.is_full)
    return {
        "tool": "pyitsx",
        "n_sequences": n,
        "elapsed_seconds": round(elapsed, 3),
        "seqs_per_second": round(n / elapsed, 1),
        "detected": len(results),
        "full_chains": full,
        "partial_chains": len(results) - full,
        "cpus": cpus,
    }


def bench_itsx(input_path: Path, cpus: int, max_seqs: int = 0):
    with tempfile.TemporaryDirectory(prefix="bench_itsx_") as tmpdir:
        tmpdir = Path(tmpdir)
        fasta_path = input_path
        fmt = "fasta" if input_path.suffix in (".fasta", ".fa", ".fna") else "fastq"

        if fmt == "fastq" or max_seqs:
            fasta_path = tmpdir / "input.fasta"
            parser = fmt
            n = 0
            with open(fasta_path, "w") as out:
                for rec in SeqIO.parse(input_path, parser):
                    out.write(f">{rec.id}\n{rec.seq}\n")
                    n += 1
                    if max_seqs and n >= max_seqs:
                        break
        else:
            n = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))

        cmd = [
            "ITSx", "-i", str(fasta_path), "-o", str(tmpdir / "out"),
            "-t", "F", "--cpu", str(cpus), "--graphical", "F",
            "--detailed_results", "T", "--positions", "T",
        ]

        t0 = time.perf_counter()
        subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        elapsed = time.perf_counter() - t0

        detected = 0
        pos_file = tmpdir / "out.positions.txt"
        if pos_file.exists():
            detected = sum(1 for line in open(pos_file) if line.strip())

        return {
            "tool": "ITSx",
            "n_sequences": n,
            "elapsed_seconds": round(elapsed, 3),
            "seqs_per_second": round(n / elapsed, 1),
            "detected": detected,
            "cpus": cpus,
        }


def bench_itsxrust(input_path: Path, hmm_path: Path, cpus: int, max_seqs: int = 0):
    with tempfile.TemporaryDirectory(prefix="bench_itsxrust_") as tmpdir:
        tmpdir = Path(tmpdir)

        if max_seqs:
            fmt = "fasta" if input_path.suffix in (".fasta", ".fa", ".fna") else "fastq"
            limited = tmpdir / "input.fasta"
            n = 0
            with open(limited, "w") as out:
                for rec in SeqIO.parse(input_path, fmt):
                    out.write(f">{rec.id}\n{rec.seq}\n")
                    n += 1
                    if n >= max_seqs:
                        break
            input_path = limited
        else:
            fmt = "fasta" if input_path.suffix in (".fasta", ".fa", ".fna") else "fastq"
            n = sum(1 for _ in SeqIO.parse(input_path, fmt))

        cmd = [
            "itsxrust", "extract",
            "-i", str(input_path), "--hmm", str(hmm_path),
            "-o", str(tmpdir / "out"), "--region", "all",
            "--preset", "ont", "--hmmer-cpu", str(cpus),
        ]

        t0 = time.perf_counter()
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        elapsed = time.perf_counter() - t0

        detected = 0
        for line in result.stderr.splitlines():
            if "Reads with computed bounds" in line:
                detected = int(line.split(":")[1].strip().split()[0])

        return {
            "tool": "ITSxRust",
            "n_sequences": n,
            "elapsed_seconds": round(elapsed, 3),
            "seqs_per_second": round(n / elapsed, 1),
            "detected": detected,
            "cpus": cpus,
        }


def main():
    parser = argparse.ArgumentParser(description="Benchmark pyitsx vs ITSx vs ITSxRust")
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("--hmm-dir", required=True, type=Path)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--max-seqs", type=int, default=0)
    parser.add_argument("-o", "--output", type=Path, help="JSON output file")
    parser.add_argument(
        "--tools", nargs="+", default=["pyitsx", "itsx", "itsxrust"],
        choices=["pyitsx", "itsx", "itsxrust"],
    )
    args = parser.parse_args()

    results = {}
    hmm_file = args.hmm_dir / "F.hmm"

    for tool in args.tools:
        print(f"\n{'='*60}", file=sys.stderr)
        print(f"Benchmarking {tool}...", file=sys.stderr)
        try:
            if tool == "pyitsx":
                results[tool] = bench_pyitsx(args.input, args.hmm_dir, args.cpus, args.max_seqs)
            elif tool == "itsx":
                results[tool] = bench_itsx(args.input, args.cpus, args.max_seqs)
            elif tool == "itsxrust":
                results[tool] = bench_itsxrust(args.input, hmm_file, args.cpus, args.max_seqs)
            r = results[tool]
            print(f"  {r['detected']}/{r['n_sequences']} detected in {r['elapsed_seconds']}s ({r['seqs_per_second']} seq/s)", file=sys.stderr)
        except Exception as e:
            print(f"  FAILED: {e}", file=sys.stderr)
            results[tool] = {"tool": tool, "error": str(e)}

    print(f"\n{'='*60}", file=sys.stderr)
    print("Summary:", file=sys.stderr)
    for tool, r in results.items():
        if "error" in r:
            print(f"  {tool}: FAILED — {r['error']}", file=sys.stderr)
        else:
            print(f"  {tool}: {r['seqs_per_second']} seq/s ({r['detected']}/{r['n_sequences']} detected)", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
