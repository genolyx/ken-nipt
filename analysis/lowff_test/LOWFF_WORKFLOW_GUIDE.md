## Low-FF Artificial Sample Workflow (NIPT v3)

This document summarizes the **major scripts** and **how to run** the Low-FF artificial dilution experiment so it can be repeated later.

Goal:
- Take already-analyzed pregnant samples (High Risk / Normal)
- Create artificial `proper_paired.bam` files at target FF levels (e.g. 3.0/3.5/4.0/4.5/5.0) by diluting with **non‑pregnant cfDNA** (FF=0 background)
- Run the existing NIPT pipeline **starting from `proper_paired.bam`**
- Compare artificial vs original results (FN/TP/HR_MISMATCH/NoCall) and summarize patterns

---

## Directory conventions

Low-FF workspace:
- `analysis/lowff_test/manifest/` : manifest TSVs (what to generate/run)
- `src/util/lowff_script/`   : helper scripts for this workflow
- `analysis/lowff_test/config/`   : **test-only** config dir override (safe; does not touch production config)

Artificial BAMs (generated):
- `analysis/<work>/artificial/<sample_name>/<sample_name>.proper_paired.bam`

Pipeline outputs (artificial samples):
- `output/<work>/<sample_name>/<sample_name>.json`
- plus `analysis/<work>/<sample_name>/Output_*` folders

Where `<work>` is usually `lowff_test`.

---

## Safety / “production 영향” 주의

- **Do NOT edit** production config under `config/<lab>/pipeline_config.json` for testing.
- Use `--config-dir analysis/lowff_test/config` to mount a **separate** config directory into Docker.
  - Example (Cordlife test config): `analysis/lowff_test/config/cordlife/pipeline_config.json`
  - In this test config, `QC.YFF` and `QC.seqFF` were lowered (e.g. 2.0) for low-FF experiments.

---

## Prerequisites

- Host has `samtools` available (used during artificial BAM generation and read counting)
- Docker image exists: `nipt_docker_v1.2`
- Non‑pregnant donor BAMs exist (5 donors recommended), each has:
  - `<ID>.proper_paired.bam` (and ideally `.bai`)

Donor BAM dir can be either:
- `<bg-bam-dir>/<ID>.proper_paired.bam`
- `<bg-bam-dir>/<ID>/<ID>.proper_paired.bam`

---

## Key pipeline integration (important)

Main pipeline: `bin/scripts/nipt_pipeline.py`

Added/used options for this workflow:
- `--from_proper_paired` : start analysis from an existing `proper_paired.bam`
- `--proper_paired_bam <path>` : provide the input BAM path (container-visible path)

Docker runner wrapper (recommended): `src/util/lowff_script/run_nipt_from_bam.sh`
- Handles staging/mounting the input BAM so the container can read it
- Supports `--config-dir <dir>` for safe test config override

---

## Step-by-step workflow (recommended)

### A) Select pregnant source samples
Script: `src/util/lowff_script/lowff_select_samples.py`

Input:
- `data/refs/cordlife/reference_make/reference_sample_list_Cordlife_all.tsv` (aggregated table)

Examples:

```bash
# High risk, male fetus only (XY), exclude MD high risk/no-call (default)
python3 src/util/lowff_script/lowff_select_samples.py \
  --sample-list data/refs/cordlife/reference_make/reference_sample_list_Cordlife_all.tsv \
  --mode high_risk --gender male \
  --min-ff0 5.0 --min-mapping-rate 98.0 \
  --require-bam \
  --out analysis/lowff_test/manifest/selected_highrisk_male.tsv
```

```bash
# High risk, female fetus only (XX)
python3 src/util/lowff_script/lowff_select_samples.py \
  --sample-list data/refs/cordlife/reference_make/reference_sample_list_Cordlife_all.tsv \
  --mode high_risk --gender female \
  --min-ff0 5.0 --min-mapping-rate 98.0 \
  --require-bam \
  --out analysis/lowff_test/manifest/selected_highrisk_female.tsv
```

Notes:
- FF0 definition used for dilution planning:
  - XY: FF0 = `YFF_2`
  - XX: FF0 = `M-SeqFF`

---

### B) Create manifest TSV (what artificial samples to generate)
Script: `src/util/lowff_script/lowff_make_manifest.py`

Example:

```bash
python3 src/util/lowff_script/lowff_make_manifest.py \
  --preg-list analysis/lowff_test/manifest/selected_highrisk_male.tsv \
  --targets 3.0,3.5,4.0,4.5,5.0 \
  --pairs 7500000 \
  --seed 42 \
  --replicates 1 \
  --bg-mode pool_plus_per_donor \
  --bg-donors GNCI25100169,GNCI25100170,GNCI25100171,GNCI25100173,GNCI25100174 \
  --write-metadata \
  --out analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv
```

Output columns include:
- `sample_name`, `preg_bam`, `ff0`, `ff_target`, `pairs`, `seed`, `bg_mode`, `bg_donors`, `out_bam`, ...

Naming convention (example):
- `LF_M_<PREGID>_FF4p5_BGPOOL5_S42_P7p5M`
- `LF_F_<PREGID>_FF3_BG<DonorID>_S42_P7p5M`

---

### C) Generate artificial BAMs (+ optionally run pipeline) in parallel
Script: `src/util/lowff_script/lowff_run_parallel.py`

Common knobs:
- `--make-bams` : generate BAMs
- `--run-pipeline` : run full pipeline from BAM
- `--max-workers 10` : parallelism
- `--config-dir analysis/lowff_test/config` : use test config safely

Examples:

```bash
# Make BAMs only (recommended first)
python3 src/util/lowff_script/lowff_run_parallel.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --bg-bam-dir /home/ken/ken-nipt/analysis/2510 \
  --make-bams \
  --labcode cordlife --age 30 \
  --root /home/ken/ken-nipt --work lowff_test \
  --max-workers 10
```

```bash
# Run pipeline only (BAMs already exist)
python3 src/util/lowff_script/lowff_run_parallel.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --bg-bam-dir /home/ken/ken-nipt/analysis/2510 \
  --run-pipeline \
  --labcode cordlife --age 30 \
  --root /home/ken/ken-nipt --work lowff_test \
  --config-dir /home/ken/ken-nipt/analysis/lowff_test/config \
  --max-workers 10
```

---

### D) (Recommended) Check artificial BAM read counts before running pipeline
Script: `src/util/lowff_script/lowff_check_bam_reads.py`

```bash
python3 src/util/lowff_script/lowff_check_bam_reads.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --root /home/ken/ken-nipt --work lowff_test \
  --min-reads 10000000 --min-mapped-reads 9500000 \
  --max-workers 10 \
  --out /home/ken/ken-nipt/analysis/lowff_test/bam_read_check_male.tsv
```

---

### E) Rerun only samples missing QC outputs (qc.filter.txt)
Script: `src/util/lowff_script/lowff_rerun_missing_qc.py`

Use when many samples exist but some are NoCall due to missing `qc.filter.txt`.

```bash
python3 src/util/lowff_script/lowff_rerun_missing_qc.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --root /home/ken/ken-nipt --work lowff_test \
  --labcode cordlife --age 30 \
  --config-dir /home/ken/ken-nipt/analysis/lowff_test/config \
  --max-workers 10
```

---

### F) Recreate only “too short” BAMs (below QC reads)
Script: `src/util/lowff_script/lowff_regen_bams_below_qc.py`

```bash
python3 src/util/lowff_script/lowff_regen_bams_below_qc.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --root /home/ken/ken-nipt --work lowff_test \
  --min-reads 10000000 \
  --max-workers 10
```

---

### G) Aggregate results and generate final report
Script: `src/util/lowff_script/make_result.py`

```bash
python3 src/util/lowff_script/make_result.py \
  --manifest analysis/lowff_test/manifest/manifest_highrisk_male_sex.tsv \
  --root /home/ken/ken-nipt --work lowff_test \
  --out-table /home/ken/ken-nipt/analysis/lowff_test/result_table_male.tsv \
  --final-report /home/ken/ken-nipt/analysis/lowff_test/final_report_male.txt
```

Important definitions used in reporting:
- **FN**: original High Risk → artificial Low Risk
- **HR_MISMATCH**: original High Risk → artificial High Risk BUT **does not retain all original targets**
  - For “original disease retention” purpose, treat **FN + HR_MISMATCH** as failures.

---

### H) Summarize “pattern table” per original sample (preg_id)
Script: `src/util/lowff_script/lowff_summarize_patterns.py`

```bash
python3 src/util/lowff_script/lowff_summarize_patterns.py \
  --result-table /home/ken/ken-nipt/analysis/lowff_test/result_table_male.tsv \
  --out /home/ken/ken-nipt/analysis/lowff_test/pattern_summary_male.tsv
```

This produces one row per `preg_id`, with FF별로:
- pool outcome
- donor outcomes counts (TP/FN/HR_MISMATCH/NoCall)
- convenience columns like `min_ff_all_bg_tp`

---

## Script reference (major)

In `src/util/lowff_script/`:
- `lowff_select_samples.py`: select eligible pregnant sources from reference sample list TSV
- `lowff_make_manifest.py`: generate manifest TSV (+ optional metadata JSONs)
- `lowff_make_artificial_bam.py`: generate one artificial `proper_paired.bam` by dilution (samtools)
- `lowff_run_parallel.py`: parallel orchestrator (make BAMs / run pipeline / both)
- `run_nipt_from_bam.sh`: docker runner for pipeline `--from_proper_paired` (+ `--config-dir`)
- `lowff_check_bam_reads.py`: verify read counts against QC thresholds
- `lowff_rerun_missing_qc.py`: rerun only samples missing QC outputs
- `lowff_regen_bams_below_qc.py`: regenerate only BAMs below QC min reads
- `make_result.py`: compare artificial vs original JSONs and write `result_table.tsv` + `final_report.txt`
- `lowff_summarize_patterns.py`: produce per-`preg_id` “pattern summary” pivot table

Legacy / optional:
- `lowff_run_batch.py`: older sequential runner (kept for reference; prefer `lowff_run_parallel.py`)
- `lowff_run_pipeline_when_male_ready.py`: wait-for-male-then-run helper (not required if you run manually)

---

## Troubleshooting notes (high-signal)

- **All samples become NoCall**:
  - Check `analysis/<work>/<sample>/Output_QC/<sample>.qc.filter.txt`
  - If missing, use `lowff_rerun_missing_qc.py` (and ensure pipeline-from-BAM runs Qualimap + QC)

- **Artificial BAM is too small (e.g. ~1.5M reads)**:
  - Use `lowff_check_bam_reads.py` to confirm read counts
  - Use `lowff_regen_bams_below_qc.py` to regenerate only failing BAMs
  - Ensure `lowff_make_artificial_bam.py` has `--min_reads` aligned with QC (default 10M)

- **Docker can’t see input BAM**:
  - Always run pipeline via `src/util/lowff_script/run_nipt_from_bam.sh`
  - It stages/copies BAM into mounted `analysis/` when needed

