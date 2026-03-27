#!/usr/bin/env python3
"""
Parallel runner for Low-FF artificial sample experiment.

Reads a manifest TSV (from lowff_make_manifest.py) and for each row:
  - generates an artificial proper_paired.bam (optional)
  - runs full NIPT pipeline starting from proper_paired.bam (optional)

Runs up to N samples concurrently (default: 10).
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Optional

import pandas as pd


def _ts() -> str:
    return dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _run(cmd: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as lf:
        lf.write(f"\n[{_ts()}] CMD: {' '.join(cmd)}\n")
        lf.flush()
        p = subprocess.run(cmd, stdout=lf, stderr=lf, text=True)
        lf.write(f"[{_ts()}] EXIT: {p.returncode}\n")
        lf.flush()
        return p.returncode


def _build_donor_bams(bg_bam_dir: Path, donors: list[str]) -> dict[str, Path]:
    """
    Allow either:
      - <bg-bam-dir>/<ID>.proper_paired.bam
      - <bg-bam-dir>/<ID>/<ID>.proper_paired.bam
    """
    out: dict[str, Path] = {}
    for d in donors:
        cand1 = bg_bam_dir / f"{d}.proper_paired.bam"
        cand2 = bg_bam_dir / d / f"{d}.proper_paired.bam"
        if cand1.exists():
            out[d] = cand1
        else:
            out[d] = cand2
    return out


def _sample_done(marker: Path) -> bool:
    return marker.exists()


def _bam_done(out_bam: Path) -> bool:
    return out_bam.exists() and (out_bam.with_suffix(out_bam.suffix + ".bai")).exists()


def _process_one(
    row: dict[str, Any],
    *,
    repo_root: Path,
    work: str,
    labcode: str,
    age: str,
    config_dir: Optional[Path],
    make_bams: bool,
    run_pipeline: bool,
    force_bam: bool,
    force_pipeline: bool,
    donor_bams: dict[str, Path],
    donors_fallback: list[str],
    logs_dir: Path,
) -> dict[str, Any]:
    sample = str(row["sample_name"])
    log_path = logs_dir / f"{sample}.log"

    out_bam = Path(str(row.get("out_bam", "")))
    if not out_bam.is_absolute():
        out_bam = repo_root / out_bam

    analysis_sample_dir = repo_root / "analysis" / work / sample
    done_marker = analysis_sample_dir / f"{sample}.pipeline_completed.marker"

    result: dict[str, Any] = {
        "sample_name": sample,
        "out_bam": str(out_bam),
        "bam_status": "SKIP",
        "pipeline_status": "SKIP",
        "bam_exit_code": "",
        "pipeline_exit_code": "",
        "log": str(log_path),
        "started_at": _ts(),
        "finished_at": "",
    }

    # ---- BAM generation ----
    if make_bams:
        if (not force_bam) and _bam_done(out_bam):
            result["bam_status"] = "SKIP_EXISTS"
        else:
            make_bam_py = Path(__file__).resolve().parent / "lowff_make_artificial_bam.py"
            preg_bam = Path(str(row["preg_bam"]))
            ff0 = float(row["ff0"])
            ff_target = float(row["ff_target"])
            pairs = int(row["pairs"])
            seed = int(row["seed"])

            bg_donors = str(row.get("bg_donors", "")).strip()
            bg_ids = [x for x in bg_donors.split(",") if x] if bg_donors else donors_fallback
            bg_paths = [str(donor_bams[d]) for d in bg_ids]

            # analysis_dir arg expects root dir that will contain <sample_name>/
            analysis_dir_root = (repo_root / "analysis" / work / "artificial").as_posix()

            cmd = [
                "python3",
                str(make_bam_py),
                "--sample_name",
                sample,
                "--preg_bam",
                str(preg_bam),
                "--ff0",
                str(ff0),
                "--ff_target",
                str(ff_target),
                "--pairs",
                str(pairs),
                "--bg_bams",
                ",".join(bg_paths),
                "--analysis_dir",
                analysis_dir_root,
                "--seed",
                str(seed),
                "--write_metadata",
            ]
            rc = _run(cmd, log_path)
            result["bam_exit_code"] = rc
            result["bam_status"] = "OK" if rc == 0 else "FAIL"

            if rc != 0:
                result["finished_at"] = _ts()
                return result

    # ---- Pipeline run ----
    if run_pipeline:
        if (not force_pipeline) and _sample_done(done_marker):
            result["pipeline_status"] = "SKIP_DONE"
        else:
            run_from_bam_sh = Path(__file__).resolve().parent / "run_nipt_from_bam.sh"
            in_bam = repo_root / "analysis" / work / "artificial" / sample / f"{sample}.proper_paired.bam"

            cmd = [
                "bash",
                str(run_from_bam_sh),
                "-s",
                sample,
                "-l",
                labcode,
                "-a",
                str(age),
                "-root",
                str(repo_root),
                "-work",
                work,
                "-b",
                str(in_bam),
            ]
            if config_dir is not None:
                cmd += ["--config-dir", str(config_dir)]
            if force_pipeline:
                cmd.append("-f")

            rc = _run(cmd, log_path)
            result["pipeline_exit_code"] = rc
            result["pipeline_status"] = "OK" if rc == 0 else "FAIL"

    result["finished_at"] = _ts()
    return result


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, type=Path, nargs="+", help="One or more manifest TSVs")
    ap.add_argument(
        "--bg-bam-dir",
        required=True,
        type=Path,
        help="Directory containing non-pregnant donor BAMs (either <ID>.proper_paired.bam or <ID>/<ID>.proper_paired.bam).",
    )
    ap.add_argument(
        "--bg-donors",
        default="GNCI25100169,GNCI25100170,GNCI25100171,GNCI25100173,GNCI25100174",
        help="Comma-separated donor IDs (used as fallback if manifest bg_donors is empty).",
    )
    ap.add_argument("--labcode", default="cordlife")
    ap.add_argument("--age", default="30")
    ap.add_argument("--root", default=str(Path.cwd()), help="Repo root (default: cwd)")
    ap.add_argument("--work", default="lowff_test", help="analysis/<work>/..., output/<work>/...")
    ap.add_argument(
        "--config-dir",
        default=None,
        help="Optional config directory override passed to run_nipt_from_bam.sh (safe test config).",
    )
    ap.add_argument("--max-workers", type=int, default=10)
    ap.add_argument("--make-bams", action="store_true")
    ap.add_argument("--run-pipeline", action="store_true")
    ap.add_argument("--force-bam", action="store_true")
    ap.add_argument("--force-pipeline", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    args = ap.parse_args(argv)

    repo_root = Path(args.root).resolve()
    config_dir = Path(args.config_dir).resolve() if args.config_dir else None
    if config_dir is not None and not config_dir.exists():
        raise SystemExit(f"--config-dir not found: {config_dir}")
    manifests = [(m if m.is_absolute() else (repo_root / m)).resolve() for m in args.manifest]

    donors_fallback = [d for d in args.bg_donors.replace(" ", "").split(",") if d]
    donor_bams = _build_donor_bams(args.bg_bam_dir, donors_fallback)

    dfs = [pd.read_csv(m, sep="\t") for m in manifests]
    df = pd.concat(dfs, ignore_index=True)
    if args.limit is not None:
        df = df.head(args.limit)

    logs_dir = repo_root / "analysis" / args.work / "runner_logs"
    summary_dir = repo_root / "analysis" / args.work / "runner_summary"
    summary_dir.mkdir(parents=True, exist_ok=True)
    summary_path = summary_dir / f"lowff_parallel_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.tsv"

    rows = df.to_dict(orient="records")
    results: list[dict[str, Any]] = []

    print(
        f"[{_ts()}] manifests={len(manifests)} rows={len(rows)} max_workers={args.max_workers}",
        flush=True,
    )
    print(
        f"[{_ts()}] make_bams={args.make_bams} run_pipeline={args.run_pipeline} work={args.work}",
        flush=True,
    )

    with ThreadPoolExecutor(max_workers=args.max_workers) as ex:
        futs = [
            ex.submit(
                _process_one,
                r,
                repo_root=repo_root,
                work=args.work,
                labcode=args.labcode,
                age=str(args.age),
                config_dir=config_dir,
                make_bams=args.make_bams,
                run_pipeline=args.run_pipeline,
                force_bam=args.force_bam,
                force_pipeline=args.force_pipeline,
                donor_bams=donor_bams,
                donors_fallback=donors_fallback,
                logs_dir=logs_dir,
            )
            for r in rows
        ]

        done = 0
        for f in as_completed(futs):
            res = f.result()
            results.append(res)
            done += 1
            print(
                f"[{_ts()}] {done}/{len(futs)} {res['sample_name']} "
                f"bam={res['bam_status']} pipeline={res['pipeline_status']}",
                flush=True,
            )

    # Write summary TSV
    fieldnames = [
        "sample_name",
        "out_bam",
        "bam_status",
        "pipeline_status",
        "bam_exit_code",
        "pipeline_exit_code",
        "started_at",
        "finished_at",
        "log",
    ]
    with summary_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in sorted(results, key=lambda x: x["sample_name"]):
            w.writerow({k: r.get(k, "") for k in fieldnames})

    ok_bam = sum(1 for r in results if r["bam_status"] in ("OK", "SKIP_EXISTS", "SKIP"))
    ok_pipe = sum(1 for r in results if r["pipeline_status"] in ("OK", "SKIP_DONE", "SKIP"))
    print(f"[{_ts()}] wrote summary: {summary_path}", flush=True)
    print(
        f"[{_ts()}] bam_ok_or_skipped={ok_bam}/{len(results)} pipeline_ok_or_skipped={ok_pipe}/{len(results)}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

