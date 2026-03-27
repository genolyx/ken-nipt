#!/usr/bin/env python3
"""
Wait until all male (LF_M_) artificial BAMs are generated, then run pipeline
for male samples only while female BAM generation continues.

This is designed to be run in parallel with BAM generation.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import subprocess
import time
from pathlib import Path
from typing import Optional


def ts() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def load_samples_from_manifest(manifest: Path) -> list[str]:
    with manifest.open() as f:
        r = csv.DictReader(f, delimiter="\t")
        out = []
        for row in r:
            s = (row.get("sample_name") or "").strip()
            if s:
                out.append(s)
    return out


def bam_ready(root: Path, work: str, sample: str) -> bool:
    bam = root / "analysis" / work / "artificial" / sample / f"{sample}.proper_paired.bam"
    bai = bam.with_suffix(bam.suffix + ".bai")
    return bam.exists() and bai.exists()


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--male-manifest", required=True, type=Path)
    ap.add_argument("--root", default=str(Path.cwd()))
    ap.add_argument("--work", default="lowff_test")
    ap.add_argument("--poll-seconds", type=int, default=60)
    ap.add_argument("--max-wait-minutes", type=int, default=0, help="0 means wait forever")
    ap.add_argument(
        "--check-reads-first",
        action="store_true",
        default=True,
        help="Check BAM read counts (QC) for LF_M_ samples before running pipeline (default: enabled).",
    )
    ap.add_argument(
        "--no-check-reads-first",
        dest="check_reads_first",
        action="store_false",
        help="Disable BAM read-count check before pipeline.",
    )
    ap.add_argument(
        "--check-only",
        action="store_true",
        help="Only wait + check BAM reads; do not run pipeline.",
    )
    ap.add_argument("--min-reads", type=int, default=10_000_000)
    ap.add_argument("--min-mapped-reads", type=int, default=9_500_000)

    # pipeline options
    ap.add_argument("--labcode", default="cordlife")
    ap.add_argument("--age", default="30")
    ap.add_argument("--max-workers", type=int, default=10)
    ap.add_argument("--config-dir", default=None)
    ap.add_argument("--force-pipeline", action="store_true", default=True)
    ap.add_argument("--no-force-pipeline", dest="force_pipeline", action="store_false")
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    male_manifest = args.male_manifest
    if not male_manifest.is_absolute():
        male_manifest = (root / male_manifest).resolve()

    samples = load_samples_from_manifest(male_manifest)
    male_samples = [s for s in samples if s.startswith("LF_M_")]
    if len(male_samples) == 0:
        raise SystemExit(f"No LF_M_ samples found in manifest: {male_manifest}")

    print(f"[{ts()}] male_manifest={male_manifest}", flush=True)
    print(f"[{ts()}] male_samples={len(male_samples)} poll={args.poll_seconds}s", flush=True)

    start = time.time()
    while True:
        ready = sum(1 for s in male_samples if bam_ready(root, args.work, s))
        print(f"[{ts()}] ready {ready}/{len(male_samples)}", flush=True)

        if ready == len(male_samples):
            break

        if args.max_wait_minutes and (time.time() - start) > args.max_wait_minutes * 60:
            raise SystemExit("Timed out waiting for male BAMs to be ready.")

        time.sleep(args.poll_seconds)

    print(f"[{ts()}] All male BAMs ready. Launching pipeline...", flush=True)

    if args.check_reads_first:
        checker = Path(__file__).resolve().parent / "lowff_check_bam_reads.py"
        out_tsv = root / "analysis" / args.work / "bam_read_check_male.tsv"
        check_cmd = [
            "PYTHONUNBUFFERED=1",
            "python3",
            str(checker),
            "--manifest",
            str(male_manifest),
            "--root",
            str(root),
            "--work",
            args.work,
            "--min-reads",
            str(args.min_reads),
            "--min-mapped-reads",
            str(args.min_mapped_reads),
            "--max-workers",
            str(args.max_workers),
            "--out",
            str(out_tsv),
        ]
        shell_check = " ".join(check_cmd)
        print(f"[{ts()}] CMD(check): {shell_check}", flush=True)
        rc = subprocess.call(shell_check, shell=True)
        if rc != 0:
            raise SystemExit(f"Read-count check failed with exit code {rc}. See: {out_tsv}")
        if args.check_only:
            print(f"[{ts()}] Check-only mode complete. TSV: {out_tsv}", flush=True)
            return 0

    runner = Path(__file__).resolve().parent / "lowff_run_parallel.py"
    cmd = [
        "PYTHONUNBUFFERED=1",
        "python3",
        str(runner),
        "--manifest",
        str(male_manifest),
        "--bg-bam-dir",
        str(root / "analysis" / "2510"),
        "--run-pipeline",
        "--labcode",
        args.labcode,
        "--age",
        str(args.age),
        "--root",
        str(root),
        "--work",
        args.work,
        "--max-workers",
        str(args.max_workers),
    ]
    if args.config_dir:
        cmd += ["--config-dir", str(Path(args.config_dir).resolve())]
    if args.force_pipeline:
        cmd.append("--force-pipeline")

    # Execute via shell so "PYTHONUNBUFFERED=1" works without env dict.
    shell_cmd = " ".join(cmd)
    print(f"[{ts()}] CMD: {shell_cmd}", flush=True)
    return subprocess.call(shell_cmd, shell=True)


if __name__ == "__main__":
    raise SystemExit(main())

