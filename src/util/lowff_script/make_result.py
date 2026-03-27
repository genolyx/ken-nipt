#!/usr/bin/env python3
"""
make_result.py

Aggregate low-FF artificial sample results and compare to original sample results.

Inputs
  - manifest TSV(s) (from lowff_make_manifest.py): contain sample_name, preg_id, preg_bam, ff_target, etc.
  - artificial output JSON: output/<work>/<sample>/<sample>.json
  - original output JSON: output/<month>/<preg_id>/<preg_id>.json (month inferred from preg_bam path)

Outputs
  - comparison table TSV
  - final_report.txt: summary (TP/FN/FP/TN/NoCall/MISSING), plus High-Risk disease-retention stats.
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def _norm_call(x: Optional[str]) -> Optional[str]:
    if x is None:
        return None
    s = str(x).strip().lower()
    if s in {"high risk", "highrisk"}:
        return "High Risk"
    if s in {"low risk", "lowrisk"}:
        return "Low Risk"
    if s in {"no call", "nocall", "no-call", "fail"}:
        return "No call"
    if s in {"", "none", "null", "nan"}:
        return None
    # keep unknown as-is (title-ish)
    return str(x).strip()


def _to_float(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        s = str(x).strip()
        if s == "" or s.lower() in {"nan", "none", "null"}:
            return None
        return float(s)
    except Exception:
        return None


def _as_list(x: Any) -> list[str]:
    if x is None:
        return []
    if isinstance(x, list):
        return [str(i) for i in x if str(i).strip() != ""]
    s = str(x).strip()
    if not s:
        return []
    # handle comma-separated strings defensively
    return [p.strip() for p in s.split(",") if p.strip()]


@dataclass
class NiptSummary:
    json_path: Path
    order_id: str = ""
    qc_result: Optional[str] = None
    fetal_gender: Optional[str] = None
    ff_yff: Optional[float] = None
    ff_seqff: Optional[float] = None
    reviewer_trisomy: Optional[str] = None
    reviewer_md: Optional[str] = None
    trisomy_list: list[str] = None  # type: ignore[assignment]
    md_list: list[str] = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.trisomy_list is None:
            self.trisomy_list = []
        if self.md_list is None:
            self.md_list = []

    @property
    def trisomy_call(self) -> str:
        # Prefer reviewer1 if present (this reflects UI/report behavior)
        if self.reviewer_trisomy is not None:
            return _norm_call(self.reviewer_trisomy) or "No call"
        # Fallback: QC gate
        if (self.qc_result or "").strip().upper() != "PASS":
            return "No call"
        return "High Risk" if len(self.trisomy_list) > 0 else "Low Risk"


def load_nipt_summary(path: Path) -> NiptSummary:
    raw = json.loads(path.read_text(encoding="utf-8"))
    nipt = raw.get("NIPT", {}) if isinstance(raw, dict) else {}

    final = nipt.get("final_results", {}) if isinstance(nipt, dict) else {}
    review = nipt.get("review", {}) if isinstance(nipt, dict) else {}
    reviewer1 = review.get("reviewer1", {}) if isinstance(review, dict) else {}

    order_id = str(final.get("order_id", "")) if isinstance(final, dict) else ""
    qc = str(final.get("QC_result", "")).strip() if isinstance(final, dict) else ""
    gender = str(final.get("fetal_gender", "")).strip() if isinstance(final, dict) else ""

    ff_yff = _to_float(final.get("fetal_fraction_yff")) if isinstance(final, dict) else None
    ff_seqff = _to_float(final.get("fetal_fraction_seqff")) if isinstance(final, dict) else None

    tri_list = _as_list(final.get("trisomy_result")) if isinstance(final, dict) else []
    md_list = _as_list(final.get("md_result")) if isinstance(final, dict) else []

    rev_tri = reviewer1.get("Trisomy_result") if isinstance(reviewer1, dict) else None
    rev_md = reviewer1.get("MD_result") if isinstance(reviewer1, dict) else None

    return NiptSummary(
        json_path=path,
        order_id=order_id,
        qc_result=qc or None,
        fetal_gender=gender or None,
        ff_yff=ff_yff,
        ff_seqff=ff_seqff,
        reviewer_trisomy=rev_tri,
        reviewer_md=rev_md,
        trisomy_list=tri_list,
        md_list=md_list,
    )


_MONTH_RE = re.compile(r"/analysis/(\d{4})/")


def infer_month_from_preg_bam(preg_bam: str) -> Optional[str]:
    m = _MONTH_RE.search(str(preg_bam))
    if m:
        return m.group(1)
    return None


def find_original_json(
    *,
    output_root: Path,
    preg_id: str,
    month_hint: Optional[str],
) -> Optional[Path]:
    candidates: list[Path] = []
    if month_hint:
        candidates.append(output_root / month_hint / preg_id / f"{preg_id}.json")
        candidates.append(output_root / month_hint / preg_id / f"{preg_id}_report.json")
        candidates.append(output_root / month_hint / f"{preg_id}.json")

    # fallback: scan one level of months
    for p in output_root.glob(f"*/{preg_id}/{preg_id}.json"):
        candidates.append(p)
    for p in output_root.glob(f"*/{preg_id}/{preg_id}_report.json"):
        candidates.append(p)

    for c in candidates:
        if c.exists():
            return c
    return None


def confusion_label(truth: str, pred: str) -> str:
    truth_n = _norm_call(truth) or "No call"
    pred_n = _norm_call(pred) or "No call"

    if pred_n == "No call":
        return "NoCall"

    if truth_n == "High Risk":
        if pred_n == "High Risk":
            return "TP"
        if pred_n == "Low Risk":
            return "FN"
        return f"UNK(truth={truth_n},pred={pred_n})"

    if truth_n == "Low Risk":
        if pred_n == "Low Risk":
            return "TN"
        if pred_n == "High Risk":
            return "FP"
        return f"UNK(truth={truth_n},pred={pred_n})"

    return f"UNK(truth={truth_n},pred={pred_n})"

def highrisk_outcome(
    *,
    truth_call: str,
    pred_call: str,
    truth_targets: list[str],
    pred_targets: list[str],
) -> str:
    """
    Outcome definition (per-sample), intended for Low-FF positive-retention experiments:
      - If original is High Risk:
          TP: artificial is High Risk AND contains ALL original targets
          HR_MISMATCH: artificial is High Risk but does NOT retain all original targets
          FN: artificial is Low Risk
          NoCall: artificial No call
      - If original is Low Risk:
          TN / FP / NoCall (based on artificial call)
      - Otherwise: UNKNOWN / MISSING
    """
    t = _norm_call(truth_call) or "No call"
    p = _norm_call(pred_call) or "No call"

    if p == "No call":
        return "NoCall"

    if t == "High Risk":
        if p == "Low Risk":
            return "FN"
        if p == "High Risk":
            ts = set(truth_targets)
            ps = set(pred_targets)
            if ts and ts.issubset(ps):
                return "TP"
            return "HR_MISMATCH"
        return f"UNK(truth={t},pred={p})"

    if t == "Low Risk":
        if p == "Low Risk":
            return "TN"
        if p == "High Risk":
            return "FP"
        return f"UNK(truth={t},pred={p})"

    return f"UNK(truth={t},pred={p})"


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, nargs="+", help="One or more manifest TSV paths")
    ap.add_argument("--root", default=str(Path.cwd()), help="Repo root (default: cwd)")
    ap.add_argument("--work", default="lowff_test", help="Artificial output work dir under output/<work>/")
    ap.add_argument(
        "--orig-output-root",
        default=None,
        help="Original output root directory (default: <root>/output).",
    )
    ap.add_argument(
        "--art-output-root",
        default=None,
        help="Artificial output root directory (default: <root>/output/<work>).",
    )
    ap.add_argument(
        "--out-table",
        default=None,
        help="Output TSV path (default: analysis/<work>/result_table.tsv).",
    )
    ap.add_argument(
        "--final-report",
        default=None,
        help="Summary report path (default: analysis/<work>/final_report.txt).",
    )
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)

    root = Path(args.root).resolve()
    manifests: list[Path] = []
    for m in args.manifest:
        mp = m if m.is_absolute() else (root / m)
        manifests.append(mp.resolve())

    orig_output_root = Path(args.orig_output_root).resolve() if args.orig_output_root else (root / "output")
    art_output_root = Path(args.art_output_root).resolve() if args.art_output_root else (root / "output" / args.work)

    out_table = Path(args.out_table).resolve() if args.out_table else (root / "analysis" / args.work / "result_table.tsv")
    final_report = (
        Path(args.final_report).resolve()
        if args.final_report
        else (root / "analysis" / args.work / "final_report.txt")
    )

    df_list = []
    for m in manifests:
        dfx = pd.read_csv(m, sep="\t")
        dfx["manifest_file"] = str(m)
        df_list.append(dfx)
    df = pd.concat(df_list, ignore_index=True)
    if args.limit is not None:
        df = df.head(args.limit)

    rows_out: list[dict[str, Any]] = []

    # cache per original sample (many artificial rows share same preg_id)
    orig_cache: dict[tuple[str, Optional[str]], tuple[Optional[Path], Optional[NiptSummary], Optional[str]]] = {}

    for _, r in df.iterrows():
        sample_name = str(r.get("sample_name", "")).strip()
        preg_id = str(r.get("preg_id", "")).strip()
        preg_bam = str(r.get("preg_bam", "")).strip()
        ff_target = _to_float(r.get("ff_target"))
        manifest_file = str(r.get("manifest_file", "")).strip()
        bg_mode = str(r.get("bg_mode", "")).strip()
        bg_donors = str(r.get("bg_donors", "")).strip()
        seed = str(r.get("seed", "")).strip()
        pairs = str(r.get("pairs", "")).strip()

        month_hint = infer_month_from_preg_bam(preg_bam)

        # original
        cache_key = (preg_id, month_hint)
        if cache_key not in orig_cache:
            ojson = find_original_json(output_root=orig_output_root, preg_id=preg_id, month_hint=month_hint)
            if ojson is None:
                orig_cache[cache_key] = (None, None, "missing_original_json")
            else:
                try:
                    osum = load_nipt_summary(ojson)
                    orig_cache[cache_key] = (ojson, osum, None)
                except Exception as e:
                    orig_cache[cache_key] = (ojson, None, f"failed_to_parse_original_json: {e}")

        ojson, osum, oerr = orig_cache[cache_key]

        # artificial
        ajson = art_output_root / sample_name / f"{sample_name}.json"
        asum: Optional[NiptSummary] = None
        aerr: Optional[str] = None
        if not ajson.exists():
            aerr = "missing_artificial_json"
        else:
            try:
                asum = load_nipt_summary(ajson)
            except Exception as e:
                aerr = f"failed_to_parse_artificial_json: {e}"

        # compare / metrics
        truth_call = osum.trisomy_call if osum else "UNKNOWN"
        pred_call = asum.trisomy_call if asum else "UNKNOWN"

        tri_truth = osum.trisomy_list if osum else []
        tri_pred = asum.trisomy_list if asum else []

        retained_any = bool(set(tri_truth) & set(tri_pred)) if tri_truth and tri_pred else False
        retained_all = set(tri_truth).issubset(set(tri_pred)) if tri_truth else False

        confusion = ""
        outcome = ""
        if osum and asum:
            confusion = confusion_label(truth_call, pred_call)
            outcome = highrisk_outcome(
                truth_call=truth_call,
                pred_call=pred_call,
                truth_targets=tri_truth,
                pred_targets=tri_pred,
            )
        elif oerr or aerr:
            confusion = "MISSING"
            outcome = "MISSING"
        else:
            confusion = "UNKNOWN"
            outcome = "UNKNOWN"

        rows_out.append(
            {
                "sample_name": sample_name,
                "manifest_file": manifest_file,
                "preg_id": preg_id,
                "month_hint": month_hint or "",
                "ff_target": ff_target if ff_target is not None else "",
                "bg_mode": bg_mode,
                "bg_donors": bg_donors,
                "seed": seed,
                "pairs": pairs,
                "orig_json": str(ojson) if ojson else "",
                "art_json": str(ajson),
                "orig_qc": (osum.qc_result if osum else "") or "",
                "art_qc": (asum.qc_result if asum else "") or "",
                "orig_gender": (osum.fetal_gender if osum else "") or "",
                "art_gender": (asum.fetal_gender if asum else "") or "",
                "orig_ff_yff": osum.ff_yff if osum else "",
                "art_ff_yff": asum.ff_yff if asum else "",
                "orig_ff_seqff": osum.ff_seqff if osum else "",
                "art_ff_seqff": asum.ff_seqff if asum else "",
                "orig_call_trisomy": truth_call,
                "art_call_trisomy": pred_call,
                "orig_trisomy_list": ",".join(tri_truth),
                "art_trisomy_list": ",".join(tri_pred),
                "retained_any": retained_any,
                "retained_all": retained_all,
                "confusion": confusion,
                "outcome": outcome,
                "orig_err": oerr or "",
                "art_err": aerr or "",
            }
        )

    out_df = pd.DataFrame(rows_out)
    out_table.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_table, sep="\t", index=False)

    # Summary report
    final_report.parent.mkdir(parents=True, exist_ok=True)

    def _count_conf(sub: pd.DataFrame) -> dict[str, int]:
        vc = sub["confusion"].fillna("MISSING").value_counts().to_dict()
        return {k: int(v) for k, v in vc.items()}

    with final_report.open("w", encoding="utf-8") as f:
        for m in manifests:
            f.write(f"manifest\t{m}\n")
        f.write(f"art_output_root\t{art_output_root}\n")
        f.write(f"orig_output_root\t{orig_output_root}\n")
        f.write(f"rows\t{len(out_df)}\n")
        f.write(f"table\t{out_table}\n")
        f.write("\n")

        # outcome (preferred for this experiment)
        f.write("## Outcome counts (disease-retention aware)\n")
        vc = out_df["outcome"].fillna("MISSING").value_counts().to_dict()
        for k, v in sorted(vc.items(), key=lambda x: (-x[1], x[0])):
            f.write(f"{k}\t{int(v)}\n")
        f.write("\n")

        # overall
        f.write("## Overall confusion counts (based on original call as truth)\n")
        for k, v in sorted(_count_conf(out_df).items(), key=lambda x: (-x[1], x[0])):
            f.write(f"{k}\t{v}\n")
        f.write("\n")

        # stratify by ff_target for High Risk truth
        truth_hr = out_df[out_df["orig_call_trisomy"].astype(str).str.lower() == "high risk"].copy()
        if len(truth_hr) > 0:
            f.write("## High Risk truth: TP/FN/HR_MISMATCH/NoCall/MISSING by ff_target\n")
            for ff, g in truth_hr.groupby("ff_target", dropna=False):
                oc = g["outcome"].fillna("MISSING").value_counts().to_dict()
                f.write(
                    "ff_target\t{ff}\tTP\t{tp}\tFN\t{fn}\tHR_MISMATCH\t{hm}\tNoCall\t{nc}\tMISSING\t{ms}\n".format(
                        ff=ff,
                        tp=int(oc.get("TP", 0)),
                        fn=int(oc.get("FN", 0)),
                        hm=int(oc.get("HR_MISMATCH", 0)),
                        nc=int(oc.get("NoCall", 0)),
                        ms=int(oc.get("MISSING", 0)),
                    )
                )
            f.write("\n")

        # stratify by ff_target for Low Risk truth
        truth_lr = out_df[out_df["orig_call_trisomy"].astype(str).str.lower() == "low risk"].copy()
        if len(truth_lr) > 0:
            f.write("## Low Risk truth: TN/FP/NoCall/MISSING by ff_target\n")
            for ff, g in truth_lr.groupby("ff_target", dropna=False):
                counts = _count_conf(g)
                f.write(f"ff_target\t{ff}\tTN\t{counts.get('TN',0)}\tFP\t{counts.get('FP',0)}\tNoCall\t{counts.get('NoCall',0)}\tMISSING\t{counts.get('MISSING',0)}\n")
            f.write("\n")

        # list false negatives (if any)
        fn = out_df[out_df["outcome"].isin(["FN", "HR_MISMATCH"])]
        if len(fn) > 0:
            f.write("## Failures to retain original High Risk targets (FN + HR_MISMATCH)\n")
            cols = [
                "sample_name",
                "preg_id",
                "ff_target",
                "orig_trisomy_list",
                "art_trisomy_list",
                "outcome",
                "orig_ff_yff",
                "art_ff_yff",
                "orig_ff_seqff",
                "art_ff_seqff",
                "orig_qc",
                "art_qc",
            ]
            f.write("\t".join(cols) + "\n")
            for _, rr in fn[cols].iterrows():
                f.write("\t".join(str(rr.get(c, "")) for c in cols) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

