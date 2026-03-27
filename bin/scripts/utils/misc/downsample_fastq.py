#!/usr/bin/env python3
"""
Standalone FASTQ downsample script.
Same logic as nipt_pipeline.downsample_fastq - uses seqtk sample.

Usage:
  python3 downsample_fastq.py R1.fastq.gz R2.fastq.gz [options]
  python3 downsample_fastq.py R1.fastq.gz R2.fastq.gz -o out_R1.fq.gz out_R2.fq.gz
"""

import argparse
import subprocess
import os
import sys


def count_fastq_reads(fastq_path):
    """Count number of reads in FASTQ file (4 lines per read)"""
    if not os.path.exists(fastq_path):
        print(f"Error: File not found: {fastq_path}", file=sys.stderr)
        return 0

    try:
        if fastq_path.endswith(".gz"):
            cmd = f"zcat '{fastq_path}' | wc -l"
        else:
            cmd = f"wc -l '{fastq_path}'"

        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            line_count = int(result.stdout.strip().split()[0])
            return line_count // 4
        return 0
    except Exception as e:
        print(f"Error counting reads: {e}", file=sys.stderr)
        return 0


def run_cmd(cmd, description):
    """Run shell command, return success bool"""
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            print(f"Failed: {description}", file=sys.stderr)
            if result.stderr:
                print(result.stderr, file=sys.stderr)
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"Timeout: {description}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return False


def downsample_fastq(r1_path, r2_path, output_r1, output_r2, target_reads, seed, backup):
    """
    Downsample paired FASTQ to target_reads.
    Uses seqtk sample -s{seed} {file} {fraction}
    """
    # Count reads
    r1_count = count_fastq_reads(r1_path)
    if r1_count == 0:
        return False

    print(f"R1 read count: {r1_count:,}")
    print(f"Target reads: {target_reads:,}")

    if r1_count <= target_reads:
        print("Read count already <= target. Copying to output...")
        if output_r1 != r1_path or output_r2 != r2_path:
            run_cmd(f"cp '{r1_path}' '{output_r1}'", "Copy R1")
            run_cmd(f"cp '{r2_path}' '{output_r2}'", "Copy R2")
        return True

    # Downsample
    sampling_fraction = target_reads / r1_count
    print(f"Downsampling to {target_reads:,} reads (fraction: {sampling_fraction:.4f})")
    print(f"Output: {output_r1}")

    # Backup if requested
    if backup:
        backup_dir = os.path.join(os.path.dirname(r1_path), "backup_original")
        os.makedirs(backup_dir, exist_ok=True)
        r1_backup = os.path.join(backup_dir, os.path.basename(r1_path))
        r2_backup = os.path.join(backup_dir, os.path.basename(r2_path))
        print(f"Backing up to {backup_dir}")
        if not run_cmd(f"cp '{r1_path}' '{r1_backup}'", "Backup R1"):
            return False
        if not run_cmd(f"cp '{r2_path}' '{r2_backup}'", "Backup R2"):
            return False

    # Ensure output dir exists
    for out in (output_r1, output_r2):
        out_dir = os.path.dirname(out)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

    # seqtk sample -s{seed} {file} {fraction} | gzip > output
    r1_cmd = f"seqtk sample -s{seed} '{r1_path}' {sampling_fraction} | gzip > '{output_r1}'"
    r2_cmd = f"seqtk sample -s{seed} '{r2_path}' {sampling_fraction} | gzip > '{output_r2}'"

    if not run_cmd(r1_cmd, "Downsample R1"):
        return False
    if not run_cmd(r2_cmd, "Downsample R2"):
        return False

    new_count = count_fastq_reads(output_r1)
    print(f"Downsampled R1 reads: {new_count:,}")
    print("Done.")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Downsample paired FASTQ files using seqtk (same logic as nipt_pipeline)"
    )
    parser.add_argument("r1", help="R1 FASTQ path (.gz or plain)")
    parser.add_argument("r2", help="R2 FASTQ path (.gz or plain)")
    parser.add_argument(
        "-o", "--output",
        nargs=2,
        metavar=("OUT_R1", "OUT_R2"),
        help="Output paths for downsampled R1 and R2 (default: input_dir/input_name_downsampled.fastq.gz)",
    )
    parser.add_argument(
        "--target-reads",
        type=int,
        default=7_500_000,
        help="Target read count after downsampling (default: 7500000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=100,
        help="seqtk sample seed (default: 100)",
    )
    parser.add_argument(
        "--backup",
        action="store_true",
        help="Backup original files to backup_original/ before downsampling",
    )
    parser.add_argument(
        "--in-place",
        action="store_true",
        help="Replace original files with downsampled (implies --backup)",
    )

    args = parser.parse_args()

    r1 = os.path.abspath(args.r1)
    r2 = os.path.abspath(args.r2)

    if not os.path.exists(r1):
        print(f"Error: R1 not found: {r1}", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(r2):
        print(f"Error: R2 not found: {r2}", file=sys.stderr)
        sys.exit(1)

    if args.in_place:
        output_r1, output_r2 = r1, r2
        backup = True
    elif args.output:
        output_r1 = os.path.abspath(args.output[0])
        output_r2 = os.path.abspath(args.output[1])
        backup = args.backup
    else:
        base1 = os.path.basename(r1).replace(".gz", "").replace(".fastq", "").replace(".fq", "").rstrip(".")
        base2 = os.path.basename(r2).replace(".gz", "").replace(".fastq", "").replace(".fq", "").rstrip(".")
        out_dir = os.path.dirname(r1)
        output_r1 = os.path.join(out_dir, f"{base1}_downsampled.fastq.gz")
        output_r2 = os.path.join(out_dir, f"{base2}_downsampled.fastq.gz")
        backup = args.backup

    ok = downsample_fastq(
        r1, r2,
        output_r1, output_r2,
        target_reads=args.target_reads,
        seed=args.seed,
        backup=backup,
    )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
