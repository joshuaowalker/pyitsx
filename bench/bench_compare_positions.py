#!/usr/bin/env python3
"""Compare pyitsx vs ITSx: detection rate, anchor positions, and performance.

Runs ITSx to get ground-truth positions, then runs pyitsx in both batched
(short-circuit) and bulk modes, comparing detection and boundary accuracy.
"""

import subprocess
import sys
import tempfile
import time
from collections import Counter
from pathlib import Path

from pyitsx.constants import Region
from pyitsx.pipeline import delimit
from pyitsx.profiles import ProfileDB


def parse_itsx_positions(positions_file: Path) -> dict[str, dict[str, tuple]]:
    """Parse ITSx positions.txt into {seq_id: {region: (start, end)}}."""
    results = {}
    with open(positions_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            seq_id = parts[0]
            # Format: seq_id\tlength\tSSU\tITS1\t5.8S\tITS2\tLSU
            regions = {}
            region_names = ["SSU", "ITS1", "5.8S", "ITS2", "LSU"]
            for i, name in enumerate(region_names):
                field = parts[i + 2] if i + 2 < len(parts) else "Not found"
                if field and field != "Not found" and field != "No start" and field != "No end":
                    try:
                        start, end = field.split(": ")[1].split("-")
                        regions[name] = (int(start), int(end))
                    except (ValueError, IndexError):
                        pass
            if regions:
                results[seq_id] = regions
    return results


def run_itsx(input_path: Path, cpus: int) -> tuple[dict, float]:
    """Run ITSx and return parsed positions + elapsed time."""
    with tempfile.TemporaryDirectory(prefix="bench_itsx_") as tmpdir:
        tmpdir = Path(tmpdir)
        cmd = [
            "ITSx", "-i", str(input_path), "-o", str(tmpdir / "out"),
            "-t", "F", "--cpu", str(cpus), "--graphical", "F",
            "--save_regions", "none", "--fasta", "F",
            "--positions", "T",
        ]
        print(f"  Running: {' '.join(cmd[:6])}...", file=sys.stderr)
        t0 = time.perf_counter()
        subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        elapsed = time.perf_counter() - t0

        positions = parse_itsx_positions(tmpdir / "out.positions.txt")
        return positions, elapsed


def run_pyitsx(input_path: Path, db: ProfileDB, cpus: int, batch_size: int) -> tuple[dict, float]:
    """Run pyitsx delimit and return positions + elapsed time."""
    seqs = db.prepare(input_path)
    n = len(seqs)
    print(f"  Running pyitsx (batch_size={batch_size}, cpus={cpus}) on {n} seqs...",
          file=sys.stderr, flush=True)

    t0 = time.perf_counter()
    results = delimit(seqs, db, cpus=cpus, batch_size=batch_size)
    elapsed = time.perf_counter() - t0

    positions = {}
    for r in results:
        regions = {}
        for b in r.bounds:
            regions[b.region.value] = (b.start, b.end)
        positions[r.seq_id] = regions
    return positions, elapsed


def compare_positions(itsx_pos, pyitsx_pos, label):
    """Compare detected positions between ITSx and pyitsx."""
    itsx_ids = set(itsx_pos.keys())
    py_ids = set(pyitsx_pos.keys())

    both = itsx_ids & py_ids
    itsx_only = itsx_ids - py_ids
    py_only = py_ids - itsx_ids

    print(f"\n=== {label} ===")
    print(f"ITSx detected:   {len(itsx_ids):,}")
    print(f"pyitsx detected: {len(py_ids):,}")
    print(f"Both detected:   {len(both):,}")
    print(f"ITSx only:       {len(itsx_only):,}")
    print(f"pyitsx only:     {len(py_only):,}")

    if not both:
        return

    # Compare positions for sequences detected by both
    region_names = ["SSU", "ITS1", "5.8S", "ITS2", "LSU"]
    for region in region_names:
        diffs = []
        start_diffs = []
        end_diffs = []
        n_both_have = 0
        n_itsx_has = 0
        n_py_has = 0

        for seq_id in both:
            itsx_r = itsx_pos[seq_id].get(region)
            py_r = pyitsx_pos[seq_id].get(region)

            if itsx_r and py_r:
                n_both_have += 1
                sd = abs(py_r[0] - itsx_r[0])
                ed = abs(py_r[1] - itsx_r[1])
                start_diffs.append(sd)
                end_diffs.append(ed)
                diffs.append(max(sd, ed))
            elif itsx_r:
                n_itsx_has += 1
            elif py_r:
                n_py_has += 1

        if not diffs:
            print(f"  {region:5s}: no overlapping detections")
            continue

        exact = sum(1 for d in diffs if d == 0)
        within1 = sum(1 for d in diffs if d <= 1)
        within5 = sum(1 for d in diffs if d <= 5)
        max_diff = max(diffs)
        mean_start = sum(start_diffs) / len(start_diffs)
        mean_end = sum(end_diffs) / len(end_diffs)

        print(f"  {region:5s}: {n_both_have:,} shared | "
              f"exact={exact:,} within1={within1:,} within5={within5:,} | "
              f"mean_start_diff={mean_start:.1f} mean_end_diff={mean_end:.1f} max={max_diff} | "
              f"itsx_only={n_itsx_has:,} pyitsx_only={n_py_has:,}")

    # Show a few examples of large discrepancies
    big_diffs = []
    for seq_id in both:
        for region in region_names:
            itsx_r = itsx_pos[seq_id].get(region)
            py_r = pyitsx_pos[seq_id].get(region)
            if itsx_r and py_r:
                d = max(abs(py_r[0] - itsx_r[0]), abs(py_r[1] - itsx_r[1]))
                if d > 10:
                    big_diffs.append((d, seq_id, region, itsx_r, py_r))

    if big_diffs:
        big_diffs.sort(reverse=True)
        print(f"\n  Top discrepancies (>{10}bp):")
        for d, seq_id, region, itsx_r, py_r in big_diffs[:10]:
            print(f"    {seq_id[:60]:60s} {region:5s} ITSx={itsx_r[0]}-{itsx_r[1]} pyitsx={py_r[0]}-{py_r[1]} diff={d}")


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=Path)
    parser.add_argument("--hmm-dir", type=Path, default=None)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--max-seqs", type=int, default=0,
                        help="Limit sequences for faster testing")
    args = parser.parse_args()

    input_path = args.input
    if args.max_seqs:
        # Write a limited FASTA
        from Bio import SeqIO
        tmpf = Path(tempfile.mktemp(suffix=".fasta", prefix="bench_"))
        n = 0
        with open(tmpf, "w") as out:
            for rec in SeqIO.parse(input_path, "fasta"):
                out.write(f">{rec.id}\n{rec.seq}\n")
                n += 1
                if n >= args.max_seqs:
                    break
        input_path = tmpf
        print(f"Limited to {n} sequences in {tmpf}", file=sys.stderr)

    hmm_dir = args.hmm_dir
    if hmm_dir is None:
        from pyitsx.profiles import find_hmm_dir
        hmm_dir = find_hmm_dir()

    db = ProfileDB(hmm_dir, organism="F")

    # Run ITSx
    print("\n--- ITSx ---", file=sys.stderr)
    itsx_pos, itsx_time = run_itsx(input_path, args.cpus)
    n_seqs = sum(1 for _ in open(input_path) if _.startswith(">"))
    print(f"  {len(itsx_pos):,}/{n_seqs:,} detected in {itsx_time:.1f}s "
          f"({n_seqs/itsx_time:.1f} seq/s)", file=sys.stderr)

    # Run pyitsx batched (short-circuit, batch_size=1)
    print("\n--- pyitsx batched (batch_size=1) ---", file=sys.stderr)
    py_batched_pos, py_batched_time = run_pyitsx(input_path, db, cpus=1, batch_size=1)
    print(f"  {len(py_batched_pos):,}/{n_seqs:,} detected in {py_batched_time:.1f}s "
          f"({n_seqs/py_batched_time:.1f} seq/s)", file=sys.stderr)

    # Run pyitsx bulk (no short-circuit)
    print("\n--- pyitsx bulk (batch_size=0, cpus=1) ---", file=sys.stderr)
    py_bulk_pos, py_bulk_time = run_pyitsx(input_path, db, cpus=1, batch_size=0)
    print(f"  {len(py_bulk_pos):,}/{n_seqs:,} detected in {py_bulk_time:.1f}s "
          f"({n_seqs/py_bulk_time:.1f} seq/s)", file=sys.stderr)

    # Performance summary
    print("\n=== Performance ===")
    print(f"{'Tool':<30s}  {'Time':>8s}  {'seq/s':>8s}  {'Detected':>10s}")
    print("-" * 65)
    print(f"{'ITSx':<30s}  {itsx_time:7.1f}s  {n_seqs/itsx_time:8.1f}  {len(itsx_pos):>10,}")
    print(f"{'pyitsx bulk (batch=0, cpu=1)':<30s}  {py_bulk_time:7.1f}s  {n_seqs/py_bulk_time:8.1f}  {len(py_bulk_pos):>10,}")
    print(f"{'pyitsx batched (batch=1)':<30s}  {py_batched_time:7.1f}s  {n_seqs/py_batched_time:8.1f}  {len(py_batched_pos):>10,}")

    # Compare batched vs ITSx
    compare_positions(itsx_pos, py_batched_pos, "ITSx vs pyitsx batched (batch_size=1)")

    # Compare bulk vs ITSx
    compare_positions(itsx_pos, py_bulk_pos, "ITSx vs pyitsx bulk (batch_size=0)")

    # Compare batched vs bulk (internal consistency)
    compare_positions(py_bulk_pos, py_batched_pos, "pyitsx bulk vs batched (short-circuit accuracy)")

    if args.max_seqs:
        tmpf.unlink()


if __name__ == "__main__":
    main()
