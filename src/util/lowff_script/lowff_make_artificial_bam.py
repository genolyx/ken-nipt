#!/usr/bin/env python3
"""
Make low-FF artificial proper_paired.bam by diluting a pregnant proper_paired.bam
with non-pregnant proper_paired.bam donors (FF=0 background).

This creates:
  analysis/<work>/<sample>/<sample>.proper_paired.bam
and a metadata json if requested.

Method
  - Determine required pregnant fraction: preg_fraction = FF_target / FF0
  - Subsample pregnant BAM and background BAM(s) with samtools view -s seed.frac
  - Merge, sort, then final downsample to exact target read count (2*pairs)

Notes
  - Uses properly paired primary alignments only: -f 2 -F 256 -F 2048
  - No MAPQ filter by default (keep closer to pipeline input); if you want MAPQ>=30, set --mapq 30.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


def sh(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def count_reads(
    bam: Path,
    samtools: str,
    mapq: Optional[int],
    threads: int,
) -> int:
    cmd = [samtools, "view", "-c", "-f", "2", "-F", "256", "-F", "2048"]
    if mapq is not None:
        cmd += ["-q", str(mapq)]
    cmd += ["-@", str(threads), str(bam)]
    res = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return int(res.stdout.strip())


def frac_to_int(p: float) -> int:
    p = max(0.0, min(1.0, p))
    return int(p * 1_000_000)


def downsample_bam(
    in_bam: Path,
    out_bam: Path,
    samtools: str,
    seed: int,
    frac: float,
    mapq: Optional[int],
    threads: int,
) -> None:
    """
    Subsample with samtools view -s seed.frac.
    IMPORTANT: keep output option (-o) before input BAM for compatibility.
    """
    frac = max(0.0, min(1.0, float(frac)))
    cmd = [samtools, "view", "-b", "-f", "2", "-F", "256", "-F", "2048", "-@", str(threads)]
    if mapq is not None:
        cmd += ["-q", str(mapq)]

    # If frac==1.0, skip -s to avoid any parsing edge cases
    if frac < 1.0:
        frac_int = frac_to_int(frac)
        cmd += ["-s", f"{seed}.{frac_int:06d}"]

    cmd += ["-o", str(out_bam), str(in_bam)]
    sh(cmd)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Make low-FF artificial proper_paired.bam by dilution.")
    ap.add_argument("--sample_name", required=True)
    ap.add_argument("--preg_bam", required=True, type=Path)
    ap.add_argument("--ff0", required=True, type=float, help="FF0 percent for pregnant source (YFF_2 for male, M-SeqFF for female)")
    ap.add_argument("--ff_target", required=True, type=float)
    ap.add_argument("--pairs", required=True, type=int, help="Target read pairs in final BAM")
    ap.add_argument("--bg_bams", required=True, help="Comma-separated background BAM paths (non-pregnant donors)")
    ap.add_argument(
        "--analysis_dir",
        default="analysis/lowff_test/artificial",
        help="Output analysis dir for artificial samples (default: analysis/lowff_test/artificial)",
    )
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--mapq", type=int, default=None, help="Optional MAPQ cutoff (e.g. 30). Default: no MAPQ filter.")
    ap.add_argument(
        "--min_reads",
        type=int,
        default=10_000_000,
        help="Minimum effective reads required in final BAM. "
        "This should satisfy QC number_of_reads threshold (default: 10,000,000).",
    )
    ap.add_argument("--write_metadata", action="store_true")
    ap.add_argument("--dry_run", action="store_true")
    args = ap.parse_args(argv)

    if args.ff0 <= 0:
        raise SystemExit("ff0 must be > 0")
    if args.ff_target <= 0:
        raise SystemExit("ff_target must be > 0")
    if args.ff_target > args.ff0:
        raise SystemExit(f"ff_target ({args.ff_target}) > ff0 ({args.ff0}) not achievable by dilution with FF=0 background")

    samtools = os.environ.get("SAMTools", "samtools")

    bg_list = [Path(p) for p in args.bg_bams.split(",") if p.strip()]
    if not bg_list:
        raise SystemExit("No bg bams provided")

    analysis_dir = Path(args.analysis_dir) / args.sample_name
    out_bam = analysis_dir / f"{args.sample_name}.proper_paired.bam"
    meta_path = analysis_dir / f"{args.sample_name}.lowff_metadata.json"

    total_reads = args.pairs * 2
    if total_reads < args.min_reads:
        raise SystemExit(
            f"pairs*2={total_reads} < --min_reads={args.min_reads}. "
            "Increase --pairs or lower --min_reads."
        )
    preg_fraction = args.ff_target / args.ff0
    preg_reads_need = int(round(total_reads * preg_fraction))
    preg_reads_need = max(0, min(total_reads, preg_reads_need))
    bg_reads_need = total_reads - preg_reads_need

    if args.dry_run:
        print(f"[DRY-RUN] out_bam={out_bam}")
        print(f"[DRY-RUN] total_reads={total_reads}, preg_fraction={preg_fraction:.6f}, preg_reads={preg_reads_need}, bg_reads={bg_reads_need}")
        return 0

    analysis_dir.mkdir(parents=True, exist_ok=True)

    tmpd = Path(tempfile.mkdtemp(prefix=f"{args.sample_name}.", dir=str(analysis_dir)))
    try:
        # Subsample pregnant
        preg_tot = count_reads(args.preg_bam, samtools, args.mapq, args.threads)
        if preg_tot <= 0:
            raise SystemExit(f"preg effective reads=0: {args.preg_bam}")
        preg_frac = preg_reads_need / preg_tot
        preg_sub = tmpd / "preg.sub.bam"
        downsample_bam(args.preg_bam, preg_sub, samtools, args.seed, preg_frac, args.mapq, args.threads)
        preg_sub_tot = count_reads(preg_sub, samtools, args.mapq, args.threads)
        print(f"[INFO] preg_tot={preg_tot} preg_need={preg_reads_need} preg_frac={preg_frac:.6f} preg_sub={preg_sub_tot}")

        # Subsample background donors with equal quotas (+ remainder distribution)
        nbg = len(bg_list)
        base_quota = bg_reads_need // nbg
        rem = bg_reads_need % nbg
        bg_subs = []
        for i, b in enumerate(bg_list):
            quota = base_quota + (1 if i < rem else 0)
            bt = count_reads(b, samtools, args.mapq, args.threads)
            if bt <= 0:
                raise SystemExit(f"bg effective reads=0: {b}")
            # If quota exceeds available reads, just take all reads from that donor.
            frac = min(1.0, quota / bt) if bt > 0 else 0.0
            out = tmpd / f"bg{i+1}.sub.bam"
            downsample_bam(b, out, samtools, args.seed + 100 + i, frac, args.mapq, args.threads)
            out_tot = count_reads(out, samtools, args.mapq, args.threads)
            print(f"[INFO] bg{i+1}_tot={bt} quota={quota} frac={frac:.6f} bg_sub={out_tot} bam={b}")
            bg_subs.append(out)

        # Merge and sort
        mix_pre = tmpd / "mix.pre.bam"
        sh([samtools, "merge", "-@", str(args.threads), "-f", str(mix_pre), str(preg_sub), *map(str, bg_subs)])

        mix_sorted = tmpd / "mix.sorted.bam"
        sh([samtools, "sort", "-@", str(args.threads), "-o", str(mix_sorted), str(mix_pre)])

        # Final downsample to exact total_reads
        cur = count_reads(mix_sorted, samtools, args.mapq, args.threads)
        if cur <= 0:
            raise SystemExit("merged effective reads=0")
        if cur < args.min_reads:
            raise SystemExit(
                f"Not enough effective reads after merge for QC: cur={cur} < min_reads={args.min_reads}. "
                "This would fail number_of_reads QC."
            )

        # If we have enough reads, we may downsample to the desired target size.
        # If cur is slightly below total_reads, keep all reads (do not fail).
        final_input = mix_sorted
        final_frac: Optional[float] = None
        final_tmp = tmpd / "final.tmp.bam"
        if cur >= total_reads:
            final_frac = total_reads / cur
            downsample_bam(mix_sorted, final_tmp, samtools, args.seed + 999, final_frac, args.mapq, args.threads)
            final_tot = count_reads(final_tmp, samtools, args.mapq, args.threads)
            final_input = final_tmp
            print(
                f"[INFO] merged_cur={cur} target={total_reads} final_frac={final_frac:.6f} final_tmp={final_tot}"
            )
        else:
            # cur < total_reads but still passes QC minimums; keep all reads.
            print(f"[INFO] merged_cur={cur} < target={total_reads} (keeping all reads; QC-focused mode)")

        # Sort + index final output
        sh([samtools, "sort", "-@", str(args.threads), "-o", str(out_bam), str(final_input)])
        sh([samtools, "index", "-@", str(args.threads), str(out_bam)])

        out_tot = count_reads(out_bam, samtools, args.mapq, args.threads)
        if out_tot < args.min_reads:
            raise SystemExit(
                f"Final BAM fails QC minimum reads: out_reads={out_tot} < min_reads={args.min_reads}."
            )

        if args.write_metadata:
            meta = {
                "sample_name": args.sample_name,
                "analysis_dir": str(analysis_dir),
                "out_bam": str(out_bam),
                "preg_bam": str(args.preg_bam),
                "ff0": args.ff0,
                "ff_target": args.ff_target,
                "pairs": args.pairs,
                "total_reads": total_reads,
                "min_reads": args.min_reads,
                "merged_effective_reads": cur,
                "final_downsample_frac": final_frac,
                "preg_fraction": preg_fraction,
                "preg_reads_need": preg_reads_need,
                "bg_reads_need": bg_reads_need,
                "bg_bams": [str(p) for p in bg_list],
                "seed": args.seed,
                "mapq": args.mapq,
                "threads": args.threads,
            }
            with meta_path.open("w", encoding="utf-8") as f:
                json.dump(meta, f, indent=2)

        print(f"[OK] wrote {out_bam}")
        return 0
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)


if __name__ == "__main__":
    raise SystemExit(main())

