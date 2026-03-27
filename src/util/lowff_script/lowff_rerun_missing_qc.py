#!/usr/bin/env python3
"""
Rerun only samples missing Output_QC/<sample>.qc.filter.txt (or JSON shows QC missing),
without recreating BAMs.

Strategy:
  - Read manifest TSV(s) (same format as manifest_highrisk_male.tsv)
  - For each sample:
      - If analysis/<work>/<sample>/Output_QC/<sample>.qc.filter.txt exists AND JSON reviewer is not "QC file not found",
        skip.
      - Else run run_nipt_from_bam.sh -f to regenerate QC + JSON, while pipeline skips expensive steps if files exist.
  - Parallelize with --max-workers (default 10)

This is intended after adding QC generation inside --from_proper_paired mode.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def ts() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def read_json(path: Path) -> Optional[dict[str, Any]]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def json_has_qc_missing(j: dict[str, Any]) -> bool:
    try:
        review = j.get("NIPT", {}).get("review", {})
        r1 = review.get("reviewer1", {})
        tri_c = str(r1.get("Trisomy_comment", "") or "")
        md_c = str(r1.get("MD_comment", "") or "")
        return ("QC file not found" in tri_c) or ("QC file not found" in md_c)
    except Exception:
        return False


def need_rerun(*, root: Path, work: str, sample: str) -> bool:
    qc_filter = root / "analysis" / work / sample / "Output_QC" / f"{sample}.qc.filter.txt"
    if not qc_filter.exists():
        return True

    out_json = root / "output" / work / sample / f"{sample}.json"
    if out_json.exists():
        j = read_json(out_json)
        if j and json_has_qc_missing(j):
            return True
    return False


def run_one(
    *,
    root: Path,
    work: str,
    sample: str,
    labcode: str,
    age: str,
    config_dir: Optional[Path],
    force: bool,
    log_dir: Path,
) -> tuple[str, str]:
    """
    Returns (sample, status)
    """
    in_bam = root / "analysis" / work / "artificial" / sample / f"{sample}.proper_paired.bam"
    if not in_bam.exists():
        return sample, "MISSING_BAM"

    sh = Path(__file__).resolve().parent / "run_nipt_from_bam.sh"
    if not sh.exists():
        return sample, "MISSING_RUNNER"

    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{sample}.log"

    cmd = [
        "bash",
        str(sh),
        "-s",
        sample,
        "-l",
        labcode,
        "-a",
        str(age),
        "-root",
        str(root),
        "-work",
        work,
        "-b",
        str(in_bam),
    ]
    if config_dir is not None:
        cmd += ["--config-dir", str(config_dir)]
    if force:
        cmd.append("-f")

    with log_path.open("a", encoding="utf-8") as lf:
        lf.write(f"\n[{ts()}] CMD: {' '.join(cmd)}\n")
        lf.flush()
        p = subprocess.run(cmd, stdout=lf, stderr=lf, text=True)
        lf.write(f"[{ts()}] EXIT: {p.returncode}\n")
        lf.flush()

    return sample, ("OK" if p.returncode == 0 else f"FAIL({p.returncode})")


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, nargs="+")
    ap.add_argument("--root", default=str(Path.cwd()))
    ap.add_argument("--work", default="lowff_test")
    ap.add_argument("--labcode", default="cordlife")
    ap.add_argument("--age", default="30")
    ap.add_argument("--max-workers", type=int, default=10)
    ap.add_argument("--config-dir", default=None, help="Optional config dir override (safe test config)")
    ap.add_argument("--force", action="store_true", default=True, help="Force rerun to refresh JSON (default: on)")
    ap.add_argument("--no-force", dest="force", action="store_false")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--dry-run", action="store_true")
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

    config_dir = Path(args.config_dir).resolve() if args.config_dir else None
    if config_dir is not None and not config_dir.exists():
        raise SystemExit(f"--config-dir not found: {config_dir}")

    samples = [str(x).strip() for x in df["sample_name"].tolist()]
    samples = [s for s in samples if s]

    targets = [s for s in samples if need_rerun(root=root, work=args.work, sample=s)]

    print(f"[{ts()}] total_samples={len(samples)} need_rerun={len(targets)} max_workers={args.max_workers}", flush=True)
    if args.dry_run:
        for s in targets[:20]:
            print(f"[DRY-RUN] {s}", flush=True)
        return 0

    log_dir = root / "analysis" / args.work / "rerun_logs"
    ok = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                run_one,
                root=root,
                work=args.work,
                sample=s,
                labcode=args.labcode,
                age=str(args.age),
                config_dir=config_dir,
                force=args.force,
                log_dir=log_dir,
            )
            for s in targets
        ]
        done = 0
        for f in as_completed(futs):
            sample, status = f.result()
            done += 1
            if status == "OK":
                ok += 1
            print(f"[{ts()}] {done}/{len(futs)} {sample} {status}", flush=True)

    print(f"[{ts()}] rerun_ok={ok}/{len(targets)} logs={log_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

