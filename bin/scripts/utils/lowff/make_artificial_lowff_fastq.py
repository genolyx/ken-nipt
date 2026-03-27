#!/usr/bin/env python3
"""
Create artificial low-FF synthetic FASTQ pairs by diluting Positive samples
with a pool of non-pregnant cfDNA FASTQs.

Goal
  For each Positive sample, generate synthetic samples at target fetal fractions
  (e.g. 3%, 3.5%, 4%, 5%) by mixing:
    - Positive FASTQ (contains fetal+maternal mixture with FF0)
    - Non-pregnant FASTQ pool (FF=0 background)

Rationale
  If non-pregnant background has FF=0, then after dilution:
    FF_new ≈ FF0 * (reads_from_positive / total_reads)
  => reads_from_positive / total_reads = FF_target / FF0

FF0 source
  - Prefer reading existing pipeline outputs:
      {analysis_dir}/{sample}/Output_FF/{sample}.fetal_fraction.txt
      {analysis_dir}/{sample}/Output_FF/{sample}.gender.txt
    and use:
      - Male   : YFF_2 (percent)
      - Female : M-SeqFF (percent)
  - Or provide an explicit TSV/CSV mapping file via --ff-map.

Output layout (compatible with NIPT docker wrapper expectation)
  out_fastq_dir/
    <synthetic_sample_id>/
      <synthetic_sample_id>_R1.fastq.gz
      <synthetic_sample_id>_R2.fastq.gz

Notes / assumptions
  - Requires `seqtk` and `gzip` on PATH (NIPT docker image already includes seqtk).
  - Pairing assumption: R1 and R2 are in the same order (standard Illumina paired FASTQ),
    so using the same `seqtk sample -s <seed>` on both yields consistent paired subsets.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


@dataclass(frozen=True)
class FastqPair:
    sample_id: str
    r1: Path
    r2: Path


@dataclass(frozen=True)
class FfInfo:
    gender: str  # "Male" / "Female" / "Unknown"
    yff2: float  # percent
    mseqff: float  # percent

    def ff0_for_dilution(self) -> Tuple[str, float]:
        """
        Returns (method, ff0_percent) used for dilution planning.
        Male   -> YFF_2
        Female -> M-SeqFF
        Unknown -> M-SeqFF (fallback)
        """
        g = (self.gender or "Unknown").strip().lower()
        if g.startswith("m"):
            return ("YFF_2", float(self.yff2))
        if g.startswith("f"):
            return ("M-SeqFF", float(self.mseqff))
        return ("M-SeqFF", float(self.mseqff))


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def run_cmd(cmd: str, dry_run: bool = False) -> None:
    if dry_run:
        print(f"[DRY-RUN] {cmd}")
        return
    subprocess.run(cmd, shell=True, check=True)


def check_dep(tool: str) -> None:
    if shutil.which(tool) is None:
        raise SystemExit(
            f"[FATAL] Required tool not found on PATH: {tool}\n"
            f"        If you run inside NIPT docker, this should exist.\n"
            f"        Otherwise install it (e.g. apt install seqtk / conda / etc.)."
        )


def safe_mkdir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def parse_targets(s: str) -> List[float]:
    out: List[float] = []
    for part in re.split(r"[,\s]+", s.strip()):
        if not part:
            continue
        out.append(float(part))
    if not out:
        raise ValueError("No targets provided")
    return out


def ff_tag(ff_percent: float) -> str:
    # 3.5 -> 3p5, 4 -> 4, 5.25 -> 5p25
    x = float(ff_percent)
    if x.is_integer():
        return str(int(x))
    s = f"{x}".rstrip("0").rstrip(".")
    return s.replace(".", "p")


def count_fastq_reads(fq: Path) -> int:
    """
    Count reads in FASTQ (R1 only) using wc -l/4.
    This is expensive for huge FASTQs; use sparingly (verification / availability).
    """
    fq_str = str(fq)
    if fq_str.endswith(".gz"):
        cmd = f"zcat {sh_quote(fq_str)} | wc -l"
    else:
        cmd = f"wc -l {sh_quote(fq_str)}"
    res = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    line_count = int(res.stdout.strip().split()[0])
    return line_count // 4


def sh_quote(s: str) -> str:
    # Minimal safe shell quoting
    return "'" + s.replace("'", "'\"'\"'") + "'"


def discover_fastq_pairs(root: Path) -> Dict[str, FastqPair]:
    """
    Discover FASTQ pairs under `root`.

    Supported layouts:
      A) root/<sample_id>/*_R1*.fastq.gz and *_R2*.fastq.gz
      B) root/*_R1*.fastq.gz and root/*_R2*.fastq.gz (grouped by prefix before _R[12])
    """
    if not root.exists():
        raise FileNotFoundError(str(root))

    # If root has subdirectories with FASTQs, treat each subdir as one sample container
    subdirs = [p for p in root.iterdir() if p.is_dir()]
    pairs: Dict[str, FastqPair] = {}

    def find_pair_in_dir(d: Path, sample_id: str) -> Optional[FastqPair]:
        fqs = [p for p in d.iterdir() if p.is_file() and p.name.endswith(FASTQ_EXTS)]
        r1s = [p for p in fqs if re.search(r"(^|[_\.])R1([_\.]|$)", p.name)]
        r2s = [p for p in fqs if re.search(r"(^|[_\.])R2([_\.]|$)", p.name)]
        if len(r1s) == 1 and len(r2s) == 1:
            return FastqPair(sample_id=sample_id, r1=r1s[0], r2=r2s[0])
        # Try common Illumina pattern *_R1_001.fastq.gz etc
        r1s = [p for p in fqs if "_R1_" in p.name or p.name.endswith("_R1.fastq.gz")]
        r2s = [p for p in fqs if "_R2_" in p.name or p.name.endswith("_R2.fastq.gz")]
        if len(r1s) == 1 and len(r2s) == 1:
            return FastqPair(sample_id=sample_id, r1=r1s[0], r2=r2s[0])
        return None

    if subdirs:
        for d in subdirs:
            sid = d.name
            pair = find_pair_in_dir(d, sid)
            if pair is not None:
                pairs[sid] = pair

    # If no pairs from subdirs, attempt flat layout
    if not pairs:
        fqs = [p for p in root.iterdir() if p.is_file() and p.name.endswith(FASTQ_EXTS)]
        # group by prefix up to _R1/_R2
        grp: Dict[str, Dict[str, Path]] = {}
        for p in fqs:
            m = re.search(r"(.+?)(?:[_\.])R([12])(?:[_\.].*|$)", p.name)
            if not m:
                continue
            prefix = m.group(1)
            r = m.group(2)
            grp.setdefault(prefix, {})[r] = p
        for prefix, rr in grp.items():
            if "1" in rr and "2" in rr:
                sid = prefix
                pairs[sid] = FastqPair(sample_id=sid, r1=rr["1"], r2=rr["2"])

    if not pairs:
        raise SystemExit(
            f"[FATAL] No FASTQ pairs discovered under: {root}\n"
            f"        Expected either root/<sample>/*R1* & *R2* or root/*R1* & *R2*."
        )
    return pairs


def read_ff_from_analysis(analysis_dir: Path, sample_id: str) -> FfInfo:
    ff_file = analysis_dir / sample_id / "Output_FF" / f"{sample_id}.fetal_fraction.txt"
    gender_file = analysis_dir / sample_id / "Output_FF" / f"{sample_id}.gender.txt"

    if not ff_file.exists():
        raise FileNotFoundError(str(ff_file))
    if not gender_file.exists():
        raise FileNotFoundError(str(gender_file))

    # fetal_fraction.txt: 2-col TSV with rows like "YFF_2 <tab> 8.12"
    yff2 = 0.0
    mseqff = 0.0
    with ff_file.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 2:
                parts = re.split(r"\s+", line)
            if len(parts) < 2:
                continue
            k = parts[0].strip().strip('"')
            try:
                v = float(parts[1])
            except ValueError:
                continue
            if k == "YFF_2":
                yff2 = v
            elif k == "M-SeqFF":
                mseqff = v

    # gender.txt: pipeline JSON builder expects row[0]=="gd_2" and row[2] is XX/XY
    gender = "Unknown"
    with gender_file.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 3:
                parts = re.split(r"\s+", line)
            if len(parts) >= 3 and parts[0] == "gd_2":
                gender = "Female" if parts[2] == "XX" else "Male"
                break

    return FfInfo(gender=gender, yff2=yff2, mseqff=mseqff)


def read_ff_map(ff_map_path: Path) -> Dict[str, FfInfo]:
    """
    Read FF mapping file.
    Supports TSV/CSV with header. Required columns:
      - sample_id
      - gender (Male/Female/Unknown)  [optional but recommended]
      - yff2 (percent)               [optional]
      - mseqff (percent)             [optional]
    If gender is missing, it falls back to Unknown.
    If yff2/mseqff missing, value defaults to 0.0.
    """
    with ff_map_path.open("r", encoding="utf-8", errors="replace") as f:
        sniff = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sniff, delimiters=",\t")
        reader = csv.DictReader(f, dialect=dialect)
        cols = {c.strip().lower(): c for c in (reader.fieldnames or [])}
        need = {"sample_id"}
        missing = [c for c in need if c not in cols]
        if missing:
            raise SystemExit(
                f"[FATAL] --ff-map missing required columns: {missing}\n"
                f"        Found: {reader.fieldnames}"
            )
        out: Dict[str, FfInfo] = {}
        for row in reader:
            sid = (row.get(cols["sample_id"]) or "").strip()
            if not sid:
                continue
            gender = (row.get(cols.get("gender", ""), "") or "Unknown").strip() or "Unknown"
            yff2 = float((row.get(cols.get("yff2", ""), "") or 0.0) or 0.0)
            mseqff = float((row.get(cols.get("mseqff", ""), "") or 0.0) or 0.0)
            out[sid] = FfInfo(gender=gender, yff2=yff2, mseqff=mseqff)
        if not out:
            raise SystemExit("[FATAL] --ff-map parsed 0 rows.")
        return out


def seqtk_sample_to_gz(
    in_fq: Path,
    out_fq_gz: Path,
    seed: int,
    sample_arg: str,
    dry_run: bool,
) -> None:
    safe_mkdir(out_fq_gz.parent)
    cmd = (
        f"seqtk sample -s{seed} {sh_quote(str(in_fq))} {sample_arg} | "
        f"gzip -c > {sh_quote(str(out_fq_gz))}"
    )
    run_cmd(cmd, dry_run=dry_run)


def concat_gz(inputs: List[Path], out_gz: Path, dry_run: bool) -> None:
    safe_mkdir(out_gz.parent)
    in_list = " ".join(sh_quote(str(p)) for p in inputs)
    cmd = f"cat {in_list} > {sh_quote(str(out_gz))}"
    run_cmd(cmd, dry_run=dry_run)


def make_synthetic_one(
    pos: FastqPair,
    bg_pairs: List[FastqPair],
    ff0_method: str,
    ff0_percent: float,
    target_ff_percent: float,
    out_fastq_dir: Path,
    out_sample_id: str,
    total_pairs: int,
    seed: int,
    bg_oversample_factor: float,
    verify_counts: bool,
    dry_run: bool,
) -> Dict[str, object]:
    if ff0_percent <= 0:
        raise SystemExit(
            f"[FATAL] FF0 is <=0 for positive sample {pos.sample_id} "
            f"(method={ff0_method}, ff0={ff0_percent})."
        )
    if target_ff_percent <= 0:
        raise SystemExit("[FATAL] target FF must be > 0")

    pos_fraction = target_ff_percent / ff0_percent
    if pos_fraction > 1.0:
        return {
            "status": "SKIP",
            "reason": "target_ff_gt_ff0",
            "pos_sample": pos.sample_id,
            "target_ff": target_ff_percent,
            "ff0_method": ff0_method,
            "ff0": ff0_percent,
        }

    pos_reads = int(round(total_pairs * pos_fraction))
    pos_reads = max(0, min(total_pairs, pos_reads))
    bg_reads = total_pairs - pos_reads

    # temp working directory (inside output for easier debugging)
    work_dir = out_fastq_dir / "_tmp_artificial_lowff"
    safe_mkdir(work_dir)
    tmpd = Path(tempfile.mkdtemp(prefix=f"{out_sample_id}.", dir=str(work_dir)))

    try:
        # 1) sample exact reads from positive
        pos_r1_sub = tmpd / "pos.R1.fastq.gz"
        pos_r2_sub = tmpd / "pos.R2.fastq.gz"
        seqtk_sample_to_gz(pos.r1, pos_r1_sub, seed=seed, sample_arg=str(pos_reads), dry_run=dry_run)
        seqtk_sample_to_gz(pos.r2, pos_r2_sub, seed=seed, sample_arg=str(pos_reads), dry_run=dry_run)

        # 2) build background pool (oversample), then sample exact bg_reads
        bg_pool_r1_parts: List[Path] = []
        bg_pool_r2_parts: List[Path] = []

        if bg_reads > 0:
            # Determine a pool fraction from total background size
            # (counting reads is expensive; only do it once per run if verify_counts is requested)
            # We avoid exact per-sample quotas by oversampling a small fraction from each, then final exact sample.
            # If verify_counts is False, we still need total_bg_reads to compute a reasonable fraction.
            total_bg_reads = 0
            for b in bg_pairs:
                total_bg_reads += count_fastq_reads(b.r1)

            want_pool_reads = int(math.ceil(bg_reads * bg_oversample_factor))
            pool_frac = min(1.0, want_pool_reads / max(1, total_bg_reads))
            # Ensure pool_frac is not too tiny (seqtk can handle small, but keep meaningful)
            pool_frac = max(pool_frac, min(0.02, 1.0)) if total_bg_reads > 0 else 1.0

            for i, b in enumerate(bg_pairs):
                part_r1 = tmpd / f"bg{i+1}.R1.fastq.gz"
                part_r2 = tmpd / f"bg{i+1}.R2.fastq.gz"
                seqtk_sample_to_gz(b.r1, part_r1, seed=seed + 10 + i, sample_arg=f"{pool_frac:.8f}", dry_run=dry_run)
                seqtk_sample_to_gz(b.r2, part_r2, seed=seed + 10 + i, sample_arg=f"{pool_frac:.8f}", dry_run=dry_run)
                bg_pool_r1_parts.append(part_r1)
                bg_pool_r2_parts.append(part_r2)

            bg_pool_r1 = tmpd / "bg.pool.R1.fastq.gz"
            bg_pool_r2 = tmpd / "bg.pool.R2.fastq.gz"
            concat_gz(bg_pool_r1_parts, bg_pool_r1, dry_run=dry_run)
            concat_gz(bg_pool_r2_parts, bg_pool_r2, dry_run=dry_run)

            bg_r1_sub = tmpd / "bg.final.R1.fastq.gz"
            bg_r2_sub = tmpd / "bg.final.R2.fastq.gz"
            seqtk_sample_to_gz(bg_pool_r1, bg_r1_sub, seed=seed + 99, sample_arg=str(bg_reads), dry_run=dry_run)
            seqtk_sample_to_gz(bg_pool_r2, bg_r2_sub, seed=seed + 99, sample_arg=str(bg_reads), dry_run=dry_run)
        else:
            bg_r1_sub = None
            bg_r2_sub = None

        # 3) final concat to output layout
        out_dir = out_fastq_dir / out_sample_id
        safe_mkdir(out_dir)
        out_r1 = out_dir / f"{out_sample_id}_R1.fastq.gz"
        out_r2 = out_dir / f"{out_sample_id}_R2.fastq.gz"

        if bg_reads > 0 and bg_r1_sub and bg_r2_sub:
            concat_gz([pos_r1_sub, bg_r1_sub], out_r1, dry_run=dry_run)
            concat_gz([pos_r2_sub, bg_r2_sub], out_r2, dry_run=dry_run)
        else:
            # Pure positive (target == ff0, pos_fraction==1)
            concat_gz([pos_r1_sub], out_r1, dry_run=dry_run)
            concat_gz([pos_r2_sub], out_r2, dry_run=dry_run)

        r1_count = None
        r2_count = None
        if verify_counts and not dry_run:
            r1_count = count_fastq_reads(out_r1)
            r2_count = count_fastq_reads(out_r2)
            if r1_count != r2_count:
                raise SystemExit(
                    f"[FATAL] Pair count mismatch for {out_sample_id}: R1={r1_count}, R2={r2_count}"
                )
            if r1_count != total_pairs:
                eprint(
                    f"[WARN] Output read pairs != requested for {out_sample_id}: "
                    f"got {r1_count}, expected {total_pairs}"
                )

        return {
            "status": "OK",
            "synthetic_sample": out_sample_id,
            "pos_sample": pos.sample_id,
            "pos_r1": str(pos.r1),
            "pos_r2": str(pos.r2),
            "target_ff": target_ff_percent,
            "ff0_method": ff0_method,
            "ff0": ff0_percent,
            "pos_fraction": round(pos_fraction, 6),
            "total_pairs": total_pairs,
            "pos_pairs": pos_reads,
            "bg_pairs": bg_reads,
            "bg_pool_n": len(bg_pairs),
            "out_r1": str(out_r1),
            "out_r2": str(out_r2),
            "out_pairs_r1": r1_count,
            "out_pairs_r2": r2_count,
        }
    finally:
        # In dry-run mode, keep temp dir for inspection.
        if not dry_run:
            shutil.rmtree(tmpd, ignore_errors=True)


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Make artificial low-FF synthetic FASTQs by dilution with non-pregnant cfDNA pool."
    )
    p.add_argument("--positive-dir", required=True, type=Path, help="Directory with Positive sample FASTQ pairs")
    p.add_argument("--nonpreg-dir", required=True, type=Path, help="Directory with non-pregnant cfDNA FASTQ pairs (6 samples)")
    p.add_argument("--out-fastq-dir", required=True, type=Path, help="Output FASTQ directory (creates per-sample subdirs)")
    p.add_argument("--targets", default="3,3.5,4,5", help="Target FF percents, e.g. '3,3.5,4,5'")
    p.add_argument("--pairs", type=int, default=7_500_000, help="Total read pairs to output per synthetic sample (default: 7,500,000)")
    p.add_argument("--analysis-dir", type=Path, default=None, help="Pipeline analysis root to read FF0 (Output_FF) for positive samples")
    p.add_argument("--ff-map", type=Path, default=None, help="TSV/CSV mapping file (sample_id,gender,yff2,mseqff)")
    p.add_argument("--seed", type=int, default=100, help="Random seed for seqtk sampling")
    p.add_argument("--bg-oversample-factor", type=float, default=1.25, help="Oversample factor for background pool before final exact sampling")
    p.add_argument("--verify-counts", action="store_true", help="Verify final output read counts (slow)")
    p.add_argument("--dry-run", action="store_true", help="Print commands only, do not create outputs")

    args = p.parse_args(argv)

    check_dep("seqtk")
    check_dep("gzip")

    targets = parse_targets(args.targets)
    if args.pairs <= 0:
        raise SystemExit("[FATAL] --pairs must be > 0")
    if args.bg_oversample_factor < 1.0:
        raise SystemExit("[FATAL] --bg-oversample-factor must be >= 1.0")

    pos_pairs = discover_fastq_pairs(args.positive_dir)
    bg_pairs_map = discover_fastq_pairs(args.nonpreg_dir)
    bg_pairs = list(bg_pairs_map.values())
    if len(bg_pairs) < 1:
        raise SystemExit("[FATAL] No nonpreg FASTQ pairs discovered.")

    # FF0 info map
    ff_map: Dict[str, FfInfo] = {}
    if args.ff_map:
        ff_map = read_ff_map(args.ff_map)

    if args.analysis_dir is None and not ff_map:
        raise SystemExit(
            "[FATAL] Need FF0 source for positives.\n"
            "        Provide --analysis-dir (to read Output_FF) or --ff-map."
        )

    safe_mkdir(args.out_fastq_dir)
    manifest_path = args.out_fastq_dir / "synthetic_lowff_manifest.tsv"
    manifest_rows: List[Dict[str, object]] = []

    # Determine FF0 for each positive sample and generate
    for pos_id, pos in sorted(pos_pairs.items()):
        try:
            if pos_id in ff_map:
                ffinfo = ff_map[pos_id]
            elif args.analysis_dir is not None:
                ffinfo = read_ff_from_analysis(args.analysis_dir, pos_id)
            else:
                raise KeyError(pos_id)
        except Exception as e:
            eprint(f"[SKIP] {pos_id}: cannot read FF0 ({e})")
            manifest_rows.append({"status": "SKIP", "pos_sample": pos_id, "reason": f"ff0_missing:{e}"})
            continue

        ff0_method, ff0 = ffinfo.ff0_for_dilution()
        if ff0 <= 0:
            eprint(f"[SKIP] {pos_id}: FF0 <= 0 (method={ff0_method}, ff0={ff0})")
            manifest_rows.append({"status": "SKIP", "pos_sample": pos_id, "reason": "ff0_le_0", "ff0_method": ff0_method, "ff0": ff0})
            continue

        for t in targets:
            out_id = f"{pos_id}__FF{ff_tag(t)}"
            row = make_synthetic_one(
                pos=pos,
                bg_pairs=bg_pairs,
                ff0_method=ff0_method,
                ff0_percent=ff0,
                target_ff_percent=t,
                out_fastq_dir=args.out_fastq_dir,
                out_sample_id=out_id,
                total_pairs=args.pairs,
                seed=args.seed,
                bg_oversample_factor=args.bg_oversample_factor,
                verify_counts=args.verify_counts,
                dry_run=args.dry_run,
            )
            # enrich with gender / raw metrics
            row = dict(row)
            row["gender"] = ffinfo.gender
            row["yff2"] = ffinfo.yff2
            row["mseqff"] = ffinfo.mseqff
            manifest_rows.append(row)

            status = row.get("status")
            if status == "OK":
                print(f"[OK] {out_id} (pos={pos_id}, FF0={ff0_method}:{ff0} -> target={t}%)")
            else:
                print(f"[{status}] {pos_id} target={t}%: {row.get('reason', '')}")

    # Write manifest
    if manifest_rows:
        # stable column order: union of keys
        keys: List[str] = []
        seen = set()
        for r in manifest_rows:
            for k in r.keys():
                if k not in seen:
                    seen.add(k)
                    keys.append(k)
        if not args.dry_run:
            with manifest_path.open("w", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
                w.writeheader()
                for r in manifest_rows:
                    w.writerow({k: r.get(k, "") for k in keys})
            print(f"[INFO] Manifest written: {manifest_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

