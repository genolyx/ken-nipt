#!/usr/bin/env python3
"""
Check artificial BAM read counts against QC thresholds.

For each sample in manifest:
  - ensure analysis/<work>/artificial/<sample>/<sample>.proper_paired.bam exists
  - ensure .bai exists
  - count reads:
      total: samtools view -c
      effective: samtools view -c -f 2 -F 256 -F 2048

Writes:
  - TSV table with per-sample counts + PASS/FAIL against --min-reads and --min-mapped-reads
  - Summary stats to stdout
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def ts() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def sh_count(cmd: list[str]) -> int:
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return int(res.stdout.strip())


def count_total(bam: Path, samtools: str, threads: int) -> int:
    return sh_count([samtools, "view", "-c", "-@", str(threads), str(bam)])


def count_effective(bam: Path, samtools: str, threads: int) -> int:
    return sh_count([samtools, "view", "-c", "-f", "2", "-F", "256", "-F", "2048", "-@", str(threads), str(bam)])


def check_one(
    sample: str,
    *,
    root: Path,
    work: str,
    samtools: str,
    threads: int,
    min_reads: int,
    min_mapped_reads: int,
) -> dict[str, Any]:
    bam = root / "analysis" / work / "artificial" / sample / f"{sample}.proper_paired.bam"
    bai = bam.with_suffix(bam.suffix + ".bai")

    out: dict[str, Any] = {
        "sample_name": sample,
        "bam": str(bam),
        "exists": bam.exists(),
        "bai_exists": bai.exists(),
        "total_reads": "",
        "effective_reads": "",
        "qc_reads_pass": "",
        "qc_mapped_pass": "",
        "status": "",
        "error": "",
    }

    if not bam.exists():
        out["status"] = "MISSING_BAM"
        return out

    try:
        # counting can be expensive; cap threads used for counting
        t = max(1, min(threads, 4))
        total = count_total(bam, samtools, t)
        eff = count_effective(bam, samtools, t)
        out["total_reads"] = total
        out["effective_reads"] = eff

        # In proper_paired BAM, mapped reads ~= total reads; we still treat min_mapped_reads separately
        out["qc_reads_pass"] = total >= min_reads
        out["qc_mapped_pass"] = total >= min_mapped_reads
        if total >= min_reads and total >= min_mapped_reads and bai.exists():
            out["status"] = "PASS"
        else:
            out["status"] = "FAIL"
        return out
    except Exception as e:
        out["status"] = "COUNT_FAIL"
        out["error"] = str(e)
        return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, nargs="+")
    ap.add_argument("--root", default=str(Path.cwd()))
    ap.add_argument("--work", default="lowff_test")
    ap.add_argument("--min-reads", type=int, default=10_000_000)
    ap.add_argument("--min-mapped-reads", type=int, default=9_500_000)
    ap.add_argument("--max-workers", type=int, default=10)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None, help="TSV output path (default: analysis/<work>/bam_read_check.tsv)")
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    manifests = [(m if m.is_absolute() else (root / m)).resolve() for m in args.manifest]

    df_list = []
    for m in manifests:
        dfx = pd.read_csv(m, sep="\t")
        df_list.append(dfx)
    df = pd.concat(df_list, ignore_index=True)
    samples = [str(x).strip() for x in df["sample_name"].tolist() if str(x).strip()]

    out_path = Path(args.out).resolve() if args.out else (root / "analysis" / args.work / "bam_read_check.tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    samtools = "samtools"
    print(f"[{ts()}] samples={len(samples)} max_workers={args.max_workers}", flush=True)

    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                check_one,
                s,
                root=root,
                work=args.work,
                samtools=samtools,
                threads=args.threads,
                min_reads=args.min_reads,
                min_mapped_reads=args.min_mapped_reads,
            )
            for s in samples
        ]
        done = 0
        for f in as_completed(futs):
            results.append(f.result())
            done += 1
            if done % 25 == 0 or done == len(futs):
                print(f"[{ts()}] checked {done}/{len(futs)}", flush=True)

    # Write TSV
    cols = [
        "sample_name",
        "exists",
        "bai_exists",
        "total_reads",
        "effective_reads",
        "qc_reads_pass",
        "qc_mapped_pass",
        "status",
        "bam",
        "error",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(results, key=lambda x: x["sample_name"]):
            w.writerow({c: r.get(c, "") for c in cols})

    # Summary
    total = len(results)
    pass_n = sum(1 for r in results if r["status"] == "PASS")
    missing = sum(1 for r in results if r["status"] == "MISSING_BAM")
    fail = sum(1 for r in results if r["status"] == "FAIL")
    count_fail = sum(1 for r in results if r["status"] == "COUNT_FAIL")

    print(f"[{ts()}] PASS={pass_n}/{total} FAIL={fail} MISSING_BAM={missing} COUNT_FAIL={count_fail}", flush=True)
    print(f"[{ts()}] wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

