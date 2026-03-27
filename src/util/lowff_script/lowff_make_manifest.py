#!/usr/bin/env python3
"""
Create a manifest for low-FF artificial samples.

Inputs
  - Pregnant source list TSV from `lowff_select_samples.py` (High Risk or Normal)
  - Non-pregnant donor IDs (5 donors) used as background (FF=0)

Outputs
  - Manifest TSV
  - Per-sample metadata JSON (optional)

Naming convention (sample_name == synthetic_id)
  LF_<SEX>_<PREGID>_FF<target>_BGPOOL5_S<seed>_P<pairs>M
  LF_<SEX>_<PREGID>_FF<target>_BG<NONPREGID>_S<seed>_P<pairs>M

where <SEX> is:
  - M for fetal_gender(gd_2)=XY
  - F for fetal_gender(gd_2)=XX
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import pandas as pd


def parse_targets(s: str) -> list[float]:
    return [float(x) for x in s.replace(" ", "").split(",") if x]


def ff_tag(x: float) -> str:
    s = f"{x}".rstrip("0").rstrip(".")
    return s.replace(".", "p")


def pairs_tag(pairs: int) -> str:
    # 7500000 -> 7p5M
    m = pairs / 1_000_000
    s = f"{m}".rstrip("0").rstrip(".")
    return "P" + s.replace(".", "p") + "M"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Create manifest for low-FF artificial BAM generation.")
    ap.add_argument("--preg-list", required=True, type=Path, help="TSV from lowff_select_samples.py")
    ap.add_argument("--targets", default="3.0,3.5,4.0,4.5,5.0")
    ap.add_argument("--pairs", type=int, default=7_500_000, help="Target read pairs in final proper_paired.bam")
    ap.add_argument("--seed", type=int, default=42, help="Base seed")
    ap.add_argument("--replicates", type=int, default=1, help="Replicate count per combination (different seeds)")
    ap.add_argument("--bg-mode", default="pool", choices=["pool", "per_donor", "pool_plus_per_donor"])
    ap.add_argument(
        "--bg-donors",
        default="GNCI25100169,GNCI25100170,GNCI25100171,GNCI25100173,GNCI25100174",
    )
    ap.add_argument("--out", required=True, type=Path, help="Output manifest TSV")
    ap.add_argument("--write-metadata", action="store_true", help="Write per-sample metadata JSONs next to output analysis dir")
    ap.add_argument(
        "--base-dir",
        default="analysis/lowff_test",
        help="LowFF workspace root directory (default: analysis/lowff_test). "
        "Artificial samples will be placed under <base-dir>/artificial/<sample_name>/...",
    )
    ap.add_argument(
        "--include-sex-tag",
        action="store_true",
        default=True,
        help="Include sex tag in sample_name as LF_<M/F>_<preg_id>_... (default: enabled).",
    )
    ap.add_argument(
        "--no-sex-tag",
        dest="include_sex_tag",
        action="store_false",
        help="Disable sex tag in sample_name (legacy naming).",
    )
    args = ap.parse_args(argv)

    preg_df = pd.read_csv(args.preg_list, sep="\t")
    targets = parse_targets(args.targets)
    donors = [d for d in args.bg_donors.replace(" ", "").split(",") if d]
    if len(donors) != 5:
        print(f"[WARN] bg donors count is {len(donors)} (expected 5)")

    bg_variants: list[tuple[str, list[str], str]] = []
    if args.bg_mode == "pool":
        bg_variants = [("pool", donors, "BGPOOL5")]
    elif args.bg_mode == "per_donor":
        bg_variants = [("donor", [d], f"BG{d}") for d in donors]
    else:
        bg_variants = [("pool", donors, "BGPOOL5")] + [("donor", [d], f"BG{d}") for d in donors]

    out_rows = []
    base_dir = Path(args.base_dir)
    artificial_root = base_dir / "artificial"
    pt = pairs_tag(args.pairs)

    for _, r in preg_df.iterrows():
        preg_id = str(r["sample_id"])
        ff0 = float(r["FF0"])
        ff0_method = str(r.get("FF0_method", ""))
        gender = str(r.get("fetal_gender(gd_2)", ""))
        preg_bam = str(r.get("bam_path", ""))
        sex_tag = ""
        if args.include_sex_tag:
            g = gender.strip().upper()
            if g == "XY":
                sex_tag = "M"
            elif g == "XX":
                sex_tag = "F"
            else:
                sex_tag = "U"
        for t in targets:
            if t > ff0:
                continue
            for rep in range(args.replicates):
                seed = args.seed + rep
                for bg_kind, bg_list, bg_tag in bg_variants:
                    if args.include_sex_tag:
                        sample_name = f"LF_{sex_tag}_{preg_id}_FF{ff_tag(t)}_{bg_tag}_S{seed}_{pt}"
                    else:
                        sample_name = f"LF_{preg_id}_FF{ff_tag(t)}_{bg_tag}_S{seed}_{pt}"
                    out_dir = artificial_root / sample_name
                    out_bam = out_dir / f"{sample_name}.proper_paired.bam"

                    meta = {
                        "sample_name": sample_name,
                        "analysis_dir": str(out_dir),
                        "output_bam": str(out_bam),
                        "preg_source": {
                            "sample_id": preg_id,
                            "bam_path": preg_bam,
                            "gender_gd2": gender,
                            "ff0_method": ff0_method,
                            "ff0": ff0,
                            "result": r.get("Result", ""),
                            "disease": r.get("Disease", ""),
                        },
                        "background": {
                            "mode": bg_kind,
                            "donors": bg_list,
                        },
                        "targets": {
                            "ff_target": t,
                            "pairs": args.pairs,
                            "seed": seed,
                        },
                        "note": "BAM dilution: FF_new ≈ FF0 * (preg_reads/total_reads), bg FF=0",
                    }

                    if args.write_metadata:
                        out_dir.mkdir(parents=True, exist_ok=True)
                        with (out_dir / f"{sample_name}.lowff_metadata.json").open("w", encoding="utf-8") as f:
                            json.dump(meta, f, indent=2)

                    out_rows.append(
                        {
                            "sample_name": sample_name,
                            "analysis_dir": str(out_dir),
                            "out_bam": str(out_bam),
                            "preg_id": preg_id,
                            "preg_bam": preg_bam,
                            "gender_gd2": gender,
                            "ff0_method": ff0_method,
                            "ff0": ff0,
                            "ff_target": t,
                            "pairs": args.pairs,
                            "seed": seed,
                            "bg_mode": bg_kind,
                            "bg_donors": ",".join(bg_list),
                        }
                    )

    out_df = pd.DataFrame(out_rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, sep="\t", index=False)
    print(f"[OK] wrote {args.out} (n={len(out_df)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

