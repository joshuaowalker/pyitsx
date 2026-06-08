#!/usr/bin/env python3
"""
Baseline benchmark: run ITSx on test data and capture results + timing.

Produces a JSON results file that other benchmarks can compare against
for both correctness and performance.

Usage:
    python bench/bench_itsx_baseline.py --input DATA.fasta --output results.json [--cpu 8]
    python bench/bench_itsx_baseline.py --input DATA.fastq --output results.json --format fastq
"""

import argparse
import json
import subprocess
import sys
import tempfile
import shutil
import time
from pathlib import Path
from Bio import SeqIO


def fasta_from_fastq(fastq_path: Path, fasta_path: Path, max_seqs: int = 0) -> int:
    """Convert FASTQ to FASTA, optionally limiting sequence count."""
    count = 0
    with open(fasta_path, "w") as out:
        for rec in SeqIO.parse(fastq_path, "fastq"):
            out.write(f">{rec.id}\n{rec.seq}\n")
            count += 1
            if max_seqs and count >= max_seqs:
                break
    return count


def parse_positions_file(positions_path: Path) -> dict:
    """Parse ITSx positions.txt into a dict keyed by sequence ID."""
    results = {}
    with open(positions_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            if len(fields) < 7:
                continue
            seq_id = fields[0]
            results[seq_id] = {
                "length": fields[1],
                "ssu": fields[2],
                "its1": fields[3],
                "five_eight_s": fields[4],
                "its2": fields[5],
                "lsu": fields[6],
                "notes": fields[7] if len(fields) > 7 else "",
            }
    return results


def run_itsx_benchmark(
    input_path: Path,
    fmt: str = "fasta",
    cpu: int = 8,
    max_seqs: int = 0,
    organism: str = "F",
) -> dict:
    """Run ITSx and return structured results with timing."""
    temp_dir = Path(tempfile.mkdtemp(prefix="bench_itsx_"))

    try:
        if fmt == "fastq":
            fasta_path = temp_dir / "input.fasta"
            n_seqs = fasta_from_fastq(input_path, fasta_path, max_seqs)
            print(f"Converted {n_seqs} sequences from FASTQ to FASTA", file=sys.stderr)
        else:
            fasta_path = input_path
            n_seqs = sum(1 for _ in SeqIO.parse(fasta_path, "fasta"))
            if max_seqs and n_seqs > max_seqs:
                limited_path = temp_dir / "input.fasta"
                n_seqs = 0
                with open(limited_path, "w") as out:
                    for rec in SeqIO.parse(fasta_path, "fasta"):
                        out.write(f">{rec.id}\n{rec.seq}\n")
                        n_seqs += 1
                        if n_seqs >= max_seqs:
                            break
                fasta_path = limited_path

        output_prefix = temp_dir / "itsx_out"

        cmd = [
            "ITSx",
            "-i", str(fasta_path),
            "-o", str(output_prefix),
            "-t", organism,
            "--cpu", str(cpu),
            "--save_regions", "all",
            "--graphical", "F",
            "--detailed_results", "T",
            "--positions", "T",
            "--not_found", "T",
        ]

        print(f"Running ITSx on {n_seqs} sequences with {cpu} CPUs...", file=sys.stderr)
        t0 = time.perf_counter()

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)

        elapsed = time.perf_counter() - t0
        print(f"ITSx completed in {elapsed:.1f}s ({n_seqs / elapsed:.1f} seq/s)", file=sys.stderr)

        positions = {}
        pos_file = Path(f"{output_prefix}.positions.txt")
        if pos_file.exists():
            positions = parse_positions_file(pos_file)

        no_detections = []
        nodet_file = Path(f"{output_prefix}_no_detections.txt")
        if nodet_file.exists():
            no_detections = [line.strip() for line in open(nodet_file) if line.strip()]

        region_counts = {}
        for region in ["full", "ITS1", "ITS2", "5_8S", "SSU", "LSU"]:
            region_file = Path(f"{output_prefix}.{region}.fasta")
            if region_file.exists():
                region_counts[region] = sum(1 for _ in SeqIO.parse(region_file, "fasta"))
            else:
                region_counts[region] = 0

        complement_count = 0
        if pos_file.exists():
            with open(pos_file) as f:
                for line in f:
                    if "complementary" in line.lower():
                        complement_count += 1

        return {
            "tool": "ITSx",
            "version": "1.1.3",
            "input_file": str(input_path),
            "input_format": fmt,
            "organism": organism,
            "n_sequences": n_seqs,
            "cpu": cpu,
            "elapsed_seconds": round(elapsed, 3),
            "seqs_per_second": round(n_seqs / elapsed, 1),
            "detected": len(positions),
            "not_detected": len(no_detections),
            "complement_strand": complement_count,
            "region_counts": region_counts,
            "positions": positions,
            "no_detections": no_detections,
        }

    finally:
        shutil.rmtree(temp_dir)


def main():
    parser = argparse.ArgumentParser(description="Benchmark ITSx baseline performance")
    parser.add_argument("-i", "--input", required=True, type=Path, help="Input sequence file")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Output JSON results file")
    parser.add_argument("--format", choices=["fasta", "fastq"], default="fasta")
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--max-seqs", type=int, default=0, help="Limit input sequences (0=all)")
    parser.add_argument("-t", "--organism", default="F", help="Organism type (default: F for fungi)")

    args = parser.parse_args()

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    results = run_itsx_benchmark(
        args.input, fmt=args.format, cpu=args.cpu,
        max_seqs=args.max_seqs, organism=args.organism,
    )

    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults written to {args.output}", file=sys.stderr)
    print(f"  Detected: {results['detected']}/{results['n_sequences']}", file=sys.stderr)
    print(f"  Regions: {results['region_counts']}", file=sys.stderr)


if __name__ == "__main__":
    main()
