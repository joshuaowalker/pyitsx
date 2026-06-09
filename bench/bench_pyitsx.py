#!/usr/bin/env python3
"""Benchmark pyitsx against ITSx and ITSxRust on real data."""

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from Bio import SeqIO


def _worker_delimit(args):
    """Run delimit in a worker process on a chunk FASTA file."""
    chunk_path, hmm_dir, mode_str = args
    from pyitsx.constants import Confidence, SearchMode
    from pyitsx.pipeline import delimit
    from pyitsx.profiles import ProfileDB

    db = ProfileDB(Path(hmm_dir), organism="F")
    seqs = db.load_sequences(Path(chunk_path))
    results = delimit(seqs, db, cpus=1, mode=SearchMode(mode_str))
    detected = [r for r in results if r.confidence != Confidence.NONE]
    full = sum(1 for r in detected if r.chain and r.chain.is_full)
    return len(seqs), len(detected), full


def bench_pyitsx(input_path: Path, hmm_dir: Path, cpus: int, max_seqs: int = 0, mode: str = "fast"):
    from pyitsx.constants import Confidence, SearchMode
    from pyitsx.pipeline import delimit
    from pyitsx.profiles import ProfileDB

    db = ProfileDB(hmm_dir, organism="F")
    seqs = db.load_sequences(input_path)
    if max_seqs:
        seqs = list(seqs)[:max_seqs]
    n = len(seqs)

    if cpus <= 1:
        t0 = time.perf_counter()
        results = delimit(seqs, db, cpus=1, mode=SearchMode(mode))
        elapsed = time.perf_counter() - t0
        detected = [r for r in results if r.confidence != Confidence.NONE]
        full = sum(1 for r in detected if r.chain and r.chain.is_full)
        total_detected, total_full = len(detected), full
    else:
        from multiprocessing import Pool

        with tempfile.TemporaryDirectory(prefix="bench_mp_") as tmpdir:
            tmpdir = Path(tmpdir)
            chunk_paths = []
            chunk_size = max(1, n // cpus)
            idx = 0
            for i in range(cpus):
                end = min(idx + chunk_size, n) if i < cpus - 1 else n
                if idx >= n:
                    break
                chunk_path = tmpdir / f"chunk_{i}.fasta"
                with open(chunk_path, "w") as f:
                    for seq in list(seqs)[idx:end]:
                        ts = seq.textize()
                        f.write(f">{ts.name}\n{ts.sequence}\n")
                chunk_paths.append(chunk_path)
                idx = end

            worker_args = [(str(p), str(hmm_dir), mode) for p in chunk_paths]

            t0 = time.perf_counter()
            with Pool(len(chunk_paths)) as pool:
                chunk_results = pool.map(_worker_delimit, worker_args)
            elapsed = time.perf_counter() - t0

        total_detected = sum(r[1] for r in chunk_results)
        total_full = sum(r[2] for r in chunk_results)

    return {
        "tool": f"pyitsx ({mode})",
        "n_sequences": n,
        "elapsed_seconds": round(elapsed, 3),
        "seqs_per_second": round(n / elapsed, 1),
        "detected": total_detected,
        "full_chains": total_full,
        "partial_chains": total_detected - total_full,
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

        itsxrust_bin = shutil.which("itsxrust") or str(Path.home() / ".local" / "bin" / "itsxrust")
        cmd = [
            itsxrust_bin, "extract",
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

    all_results = {}
    hmm_file = args.hmm_dir / "F.hmm"

    for tool in args.tools:
        if tool == "pyitsx":
            for mode in ("fast", "best"):
                key = f"pyitsx_{mode}"
                print(f"\n{'='*60}", file=sys.stderr)
                print(f"Benchmarking pyitsx ({mode})...", file=sys.stderr)
                try:
                    all_results[key] = bench_pyitsx(args.input, args.hmm_dir, args.cpus, args.max_seqs, mode=mode)
                    r = all_results[key]
                    print(f"  {r['detected']}/{r['n_sequences']} detected in {r['elapsed_seconds']}s ({r['seqs_per_second']} seq/s)", file=sys.stderr)
                except Exception as e:
                    print(f"  FAILED: {e}", file=sys.stderr)
                    all_results[key] = {"tool": f"pyitsx ({mode})", "error": str(e)}
        else:
            print(f"\n{'='*60}", file=sys.stderr)
            print(f"Benchmarking {tool}...", file=sys.stderr)
            try:
                if tool == "itsx":
                    all_results[tool] = bench_itsx(args.input, args.cpus, args.max_seqs)
                elif tool == "itsxrust":
                    all_results[tool] = bench_itsxrust(args.input, hmm_file, args.cpus, args.max_seqs)
                r = all_results[tool]
                print(f"  {r['detected']}/{r['n_sequences']} detected in {r['elapsed_seconds']}s ({r['seqs_per_second']} seq/s)", file=sys.stderr)
            except Exception as e:
                print(f"  FAILED: {e}", file=sys.stderr)
                all_results[tool] = {"tool": tool, "error": str(e)}

    print(f"\n{'='*60}", file=sys.stderr)
    print("Summary:", file=sys.stderr)
    for key, r in all_results.items():
        if "error" in r:
            print(f"  {r['tool']}: FAILED — {r['error']}", file=sys.stderr)
        else:
            print(f"  {r['tool']}: {r['seqs_per_second']} seq/s ({r['detected']}/{r['n_sequences']} detected)", file=sys.stderr)

    results = all_results

    if args.output:
        with open(args.output, "w") as f:
            json.dump(results, f, indent=2)
        print(f"\nResults written to {args.output}", file=sys.stderr)
    else:
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
