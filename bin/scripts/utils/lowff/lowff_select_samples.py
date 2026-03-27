#!/usr/bin/env python3
"""
Select candidate pregnant samples (High Risk / Normal) for low-FF artificial dilution.

Primary input:
  data/refs/cordlife/reference_make/reference_sample_list_Cordlife_all.tsv

This TSV already aggregates values from output JSON + QC/FF files.

FF0 definition (aligned with NIPT pipeline QC logic):
  - male fetus (gd_2 == XY): FF0 = YFF_2
  - female fetus (gd_2 == XX): FF0 = M-SeqFF
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import pandas as pd


def norm(x: object) -> str:
    if x is None:
        return ""
    if isinstance(x, float) and pd.isna(x):
        return ""
    return str(x).strip()


def to_float(x: object) -> Optional[float]:
    s = norm(x)
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def compute_ff0(row: pd.Series) -> tuple[str, Optional[float]]:
    g = norm(row.get("fetal_gender(gd_2)")).upper()
    if g == "XY":
        return ("YFF_2", to_float(row.get("YFF_2")))
    return ("M-SeqFF", to_float(row.get("M-SeqFF")))


def build_bam_path(row: pd.Series) -> Path:
    sample_dir = Path(norm(row.get("sample_dir")))
    sample_id = norm(row.get("sample_id"))
    return sample_dir / f"{sample_id}.proper_paired.bam"


def parse_csv_list(s: str) -> list[str]:
    return [p for p in s.replace(" ", "").split(",") if p]


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Select High Risk / Normal samples for low-FF dilution experiments.")
    ap.add_argument("--sample-list", required=True, type=Path, help="Input TSV (reference_sample_list_*.tsv)")
    ap.add_argument("--mode", required=True, choices=["high_risk", "normal"])
    ap.add_argument(
        "--diseases",
        default="",
        help="Optional comma-separated disease filter for high_risk mode (substring match in Disease). "
        "Default: empty (no disease filter; include SCA/T16/Other, etc.).",
    )
    ap.add_argument("--gender", default="any", choices=["any", "male", "female"])
    ap.add_argument("--min-ff0", type=float, default=5.0)
    ap.add_argument("--min-mapping-rate", type=float, default=98.0)
    ap.add_argument("--max-dup-rate", type=float, default=None)
    ap.add_argument("--min-coverage", type=float, default=None)
    ap.add_argument(
        "--exclude-md-high-risk",
        action="store_true",
        default=True,
        help="Exclude samples with MDResult == High Risk or No call (default: enabled).",
    )
    ap.add_argument(
        "--include-md-high-risk",
        action="store_true",
        default=False,
        help="Include MDResult == High Risk (overrides --exclude-md-high-risk).",
    )
    ap.add_argument("--require-bam", action="store_true")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args(argv)

    df = pd.read_csv(args.sample_list, sep="\t")

    # Option resolution
    if args.include_md_high_risk:
        args.exclude_md_high_risk = False

    # Mode filters
    if args.mode == "high_risk":
        df = df[df["Result"].astype(str).str.strip().str.lower() == "high risk"]
        diseases = parse_csv_list(args.diseases) if args.diseases else []
        if diseases:
            df = df[df["Disease"].fillna("").astype(str).apply(lambda x: any(d in x for d in diseases))]
    else:
        df = df[df["Result"].astype(str).str.strip().str.lower() == "low risk"]
        if "MDResult" in df.columns:
            df = df[df["MDResult"].fillna("").astype(str).str.strip().str.lower().isin(["low risk", "none", ""])]

    # Gender filter
    if args.gender == "male":
        df = df[df["fetal_gender(gd_2)"].astype(str).str.upper() == "XY"]
    elif args.gender == "female":
        df = df[df["fetal_gender(gd_2)"].astype(str).str.upper() == "XX"]

    # QC metric filters
    if "mapping_rate(%)" in df.columns:
        df = df[pd.to_numeric(df["mapping_rate(%)"], errors="coerce") >= args.min_mapping_rate]
    if args.max_dup_rate is not None and "duplication_rate(%)" in df.columns:
        df = df[pd.to_numeric(df["duplication_rate(%)"], errors="coerce") <= args.max_dup_rate]
    if args.min_coverage is not None and "mean_coverageData(X)" in df.columns:
        df = df[pd.to_numeric(df["mean_coverageData(X)"], errors="coerce") >= args.min_coverage]

    # Exclude MD High Risk / No call (for both modes, but especially for high_risk as requested)
    if args.exclude_md_high_risk and "MDResult" in df.columns:
        md = df["MDResult"].fillna("").astype(str).str.strip().str.lower()
        df = df[~md.isin(["high risk", "no call"])]

    # Compute FF0 + BAM existence
    df = df.copy()
    methods, ff0s, bams, exists = [], [], [], []
    for _, r in df.iterrows():
        m, ff0 = compute_ff0(r)
        bp = build_bam_path(r)
        methods.append(m)
        ff0s.append(ff0)
        bams.append(str(bp))
        exists.append(bp.exists())

    df["FF0_method"] = methods
    df["FF0"] = ff0s
    df["bam_path"] = bams
    df["bam_exists"] = exists

    df = df[pd.to_numeric(df["FF0"], errors="coerce") >= args.min_ff0]
    if args.require_bam:
        df = df[df["bam_exists"] == True]  # noqa: E712

    keep = [
        "month",
        "sample_id",
        "Result",
        "Disease",
        "MDResult",
        "fetal_gender(gd_2)",
        "FF0_method",
        "FF0",
        "YFF_2",
        "M-SeqFF",
        "mapping_rate(%)",
        "duplication_rate(%)",
        "mean_coverageData(X)",
        "GC_content(%)",
        "bam_path",
        "bam_exists",
    ]
    keep = [c for c in keep if c in df.columns]

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values(["month", "sample_id"], inplace=True)
    df.to_csv(args.out, sep="\t", index=False, columns=keep)
    print(f"[OK] wrote {args.out} (n={len(df)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

