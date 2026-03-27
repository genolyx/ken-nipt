#!/usr/bin/env python3
"""
Batch runner for low-FF experiment:
  - reads manifest TSV from `lowff_make_manifest.py`
  - generates artificial proper_paired.bam for each row
  - optionally runs full NIPT pipeline from BAM using sibling `run_nipt_from_bam.sh`
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
from typing import Optional

import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument(
        "--bg-bam-dir",
        required=True,
        type=Path,
        help="Directory containing non-pregnant donor BAMs as analysis/<month>/<id>/<id>.proper_paired.bam or symlinks",
    )
    ap.add_argument(
        "--bg-donors",
        default="GNCI25100169,GNCI25100170,GNCI25100171,GNCI25100173,GNCI25100174",
    )
    ap.add_argument("--make-bams", action="store_true")
    ap.add_argument("--run-pipeline", action="store_true")
    ap.add_argument("--labcode", default="cordlife")
    ap.add_argument("--age", default="30")
    ap.add_argument("--root", default=str(Path.cwd()))
    ap.add_argument("--work", default="lowff_test")
    ap.add_argument(
        "--base-dir",
        default="analysis/lowff_test",
        help="LowFF workspace base directory (default: analysis/lowff_test)",
    )
    ap.add_argument("--detached", action="store_true")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--limit", type=int, default=None, help="Optional limit for quick testing")
    args = ap.parse_args(argv)

    df = pd.read_csv(args.manifest, sep="\t")
    if args.limit is not None:
        df = df.head(args.limit)

    donors = [d for d in args.bg_donors.replace(" ", "").split(",") if d]

    # Build donor bam paths
    donor_bams = {}
    for d in donors:
        # Expect: <bg-bam-dir>/<d>.proper_paired.bam OR <bg-bam-dir>/<d>/<d>.proper_paired.bam
        cand1 = args.bg_bam_dir / f"{d}.proper_paired.bam"
        cand2 = args.bg_bam_dir / d / f"{d}.proper_paired.bam"
        if cand1.exists():
            donor_bams[d] = cand1
        elif cand2.exists():
            donor_bams[d] = cand2
        else:
            # allow manifest to proceed; BAM creation will fail clearly if used
            donor_bams[d] = cand2

    make_bam_py = _SCRIPT_DIR / "lowff_make_artificial_bam.py"
    run_from_bam_sh = _SCRIPT_DIR / "run_nipt_from_bam.sh"
    base_dir = Path(args.base_dir)
    artificial_root = base_dir / "artificial"

    for _, r in df.iterrows():
        sample = str(r["sample_name"])
        preg_bam = str(r["preg_bam"])
        ff0 = float(r["ff0"])
        ff_target = float(r["ff_target"])
        pairs = int(r["pairs"])
        seed = int(r["seed"])
        bg_donors = str(r["bg_donors"]).split(",") if str(r["bg_donors"]).strip() else donors
        bg_bams = [str(donor_bams[d]) for d in bg_donors]

        if args.make_bams:
            cmd = [
                "python3",
                str(make_bam_py),
                "--sample_name",
                sample,
                "--preg_bam",
                preg_bam,
                "--ff0",
                str(ff0),
                "--ff_target",
                str(ff_target),
                "--pairs",
                str(pairs),
                "--bg_bams",
                ",".join(bg_bams),
                "--analysis_dir",
                str(artificial_root),
                "--seed",
                str(seed),
                "--write_metadata",
            ]
            run(cmd)

        if args.run_pipeline:
            # input BAM is in <base-dir>/artificial/<sample>/<sample>.proper_paired.bam
            in_bam = artificial_root / sample / f"{sample}.proper_paired.bam"
            cmd = [
                "bash",
                str(run_from_bam_sh),
                "-s",
                sample,
                "-l",
                args.labcode,
                "-a",
                str(args.age),
                "-root",
                str(args.root),
                "-work",
                args.work,
                "-b",
                str(in_bam),
            ]
            if args.detached:
                cmd.append("--detached")
            if args.force:
                cmd.append("-f")
            run(cmd)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

