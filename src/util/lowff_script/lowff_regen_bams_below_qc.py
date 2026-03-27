#!/usr/bin/env python3
"""
Regenerate only "too short" artificial BAMs (below QC read threshold).

Why:
  Some previously generated artificial BAMs ended up with far fewer effective reads
  (e.g. ~1.5M) and therefore fail QC (number_of_reads / number_of_mapped_reads).

What this does:
  - Reads manifest TSV(s)
  - For each sample, counts effective reads in existing artificial BAM:
        samtools view -c -f 2 -F 256 -F 2048
  - If reads >= --min-reads: SKIP
  - Else: re-run lowff_make_artificial_bam.py to regenerate that sample BAM

Notes:
  - Does NOT run the NIPT pipeline. After regeneration, rerun pipeline for those samples.
  - Requires samtools on host (uses SAMTools env var if set).
"""

from __future__ import annotations

import argparse
import datetime as dt
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def ts() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def count_effective_reads(bam: Path, samtools: str, threads: int) -> int:
    cmd = [samtools, "view", "-c", "-f", "2", "-F", "256", "-F", "2048", "-@", str(threads), str(bam)]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return int(res.stdout.strip())


def regen_one(
    row: dict[str, Any],
    *,
    root: Path,
    work: str,
    samtools: str,
    threads: int,
    min_reads: int,
    dry_run: bool,
) -> tuple[str, str, Optional[int], Optional[int]]:
    sample = str(row["sample_name"]).strip()
    if not sample:
        return "", "SKIP_EMPTY", None, None

    # existing artificial BAM location (standard)
    bam = root / "analysis" / work / "artificial" / sample / f"{sample}.proper_paired.bam"
    if not bam.exists():
        return sample, "MISSING_BAM", None, None

    try:
        cur = count_effective_reads(bam, samtools, threads=max(1, min(threads, 4)))
    except Exception as e:
        return sample, f"COUNT_FAIL({e})", None, None

    if cur >= min_reads:
        return sample, "SKIP_OK", cur, None

    # Regenerate using manifest parameters
    make_py = Path(__file__).resolve().parent / "lowff_make_artificial_bam.py"
    analysis_dir_root = root / "analysis" / work / "artificial"

    cmd = [
        "python3",
        str(make_py),
        "--sample_name",
        sample,
        "--preg_bam",
        str(row["preg_bam"]),
        "--ff0",
        str(row["ff0"]),
        "--ff_target",
        str(row["ff_target"]),
        "--pairs",
        str(row["pairs"]),
        "--bg_bams",
        str(row["bg_bams"]),
        "--analysis_dir",
        str(analysis_dir_root),
        "--seed",
        str(row["seed"]),
        "--min_reads",
        str(min_reads),
        "--write_metadata",
    ]

    if dry_run:
        return sample, "WOULD_REGEN", cur, None

    p = subprocess.run(cmd, text=True)
    if p.returncode != 0:
        return sample, f"REGEN_FAIL({p.returncode})", cur, None

    # recount
    try:
        newc = count_effective_reads(bam, samtools, threads=max(1, min(threads, 4)))
    except Exception:
        newc = None
    return sample, "REGEN_OK", cur, newc


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, nargs="+")
    ap.add_argument("--root", default=str(Path.cwd()))
    ap.add_argument("--work", default="lowff_test")
    ap.add_argument("--min-reads", type=int, default=10_000_000)
    ap.add_argument("--max-workers", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8, help="Threads passed to generator; counting uses up to 4.")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    manifests = [(m if m.is_absolute() else (root / m)).resolve() for m in args.manifest]

    df_list = []
    for m in manifests:
        dfx = pd.read_csv(m, sep="\t")
        df_list.append(dfx)
    df = pd.concat(df_list, ignore_index=True)
    if args.limit is not None:
        df = df.head(args.limit)

    # Normalize bg_bams: manifest has bg_donors (ids). Convert to bam paths based on standard location.
    # For regen we prefer using the exact donor BAMs under /analysis/2510/<id>/<id>.proper_paired.bam,
    # since that's what was used during the batch.
    def donors_to_bams(s: str) -> str:
        s = str(s).strip()
        if not s:
            return ""
        ids = [x for x in s.split(",") if x]
        # hardcode 2510 as the donor month used in manifests (consistent with previous runs)
        bams = [str(root / "analysis" / "2510" / i / f"{i}.proper_paired.bam") for i in ids]
        return ",".join(bams)

    df = df.copy()
    if "bg_bams" not in df.columns:
        # manifest_highrisk_male.tsv uses bg_donors; convert.
        df["bg_bams"] = df["bg_donors"].apply(donors_to_bams)

    samtools = str(Path(subprocess.getoutput("echo ${SAMTools:-samtools}")).name)
    # If SAMTools env var is set to full path, prefer it
    samtools = str(Path(subprocess.getoutput("bash -lc 'echo ${SAMTools:-samtools}'")).as_posix()).strip() or "samtools"

    rows = df.to_dict(orient="records")
    print(f"[{ts()}] rows={len(rows)} min_reads={args.min_reads} max_workers={args.max_workers} dry_run={args.dry_run}", flush=True)

    regen = 0
    skipped = 0
    missing = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                regen_one,
                r,
                root=root,
                work=args.work,
                samtools=samtools,
                threads=args.threads,
                min_reads=args.min_reads,
                dry_run=args.dry_run,
            )
            for r in rows
        ]
        done = 0
        for f in as_completed(futs):
            sample, status, cur, newc = f.result()
            done += 1
            if status in ("SKIP_OK",):
                skipped += 1
            elif status in ("MISSING_BAM",):
                missing += 1
            elif status in ("REGEN_OK", "WOULD_REGEN"):
                regen += 1
            cur_s = "" if cur is None else str(cur)
            new_s = "" if newc is None else str(newc)
            print(f"[{ts()}] {done}/{len(futs)} {sample}\t{status}\tcur={cur_s}\tnew={new_s}", flush=True)

    print(f"[{ts()}] done skipped_ok={skipped} regen={regen} missing_bam={missing}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

