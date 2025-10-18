#!/bin/bash
set -euo pipefail

# MD-Only Pipeline Runner
# Runs MD analysis on existing BAM files without FASTQ processing

usage() {
    cat <<'USAGE'
Usage: run_md_pipeline.sh \
  -s <sample_id> \
  -l <labcode> \
  -root <root_directory> \
  -work <work_directory> \
  [--fetal_gender <M|F>] \
  [--no-log] [--detached] [-f|--force] [-h]

Examples:
  # Run MD analysis on sample
  ./run_md_pipeline.sh -s GNCI25080163 -l cordlife -root /home/ken/ken-nipt -work 2508

  # Force rerun
  ./run_md_pipeline.sh -s GNCI25080163 -l cordlife -root /home/ken/ken-nipt -work 2508 -f

  # Artificial sample with known fetal gender
  ./run_md_pipeline.sh -s 0001_1p36_FF05_10M -l cordlife -root /home/ken/ken-nipt -work md_test --fetal_gender M

Notes:
- Sample BAM file must exist: analysis/<work_dir>/<sample_id>/<sample_id>.proper_paired.bam (or .sorted.bam)
- Results will be saved in: analysis/<work_dir>/<sample_id>/Output_WC, Output_WCX, plots/
- --fetal_gender is optional. If not provided, gender will be detected from Output_FF/gender.txt or BAM
USAGE
    exit 1
}

# 기본값 설정
ENABLE_LOGGING=true
DETACHED_MODE=false
FORCE_EXECUTION=false
SAMPLE_ID=""
LABCODE=""
ROOT_DIR=""
WORK_DIR=""
FETAL_GENDER=""

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s) SAMPLE_ID="$2"; shift 2 ;;
        -l) LABCODE="$2"; shift 2 ;;
        -root) ROOT_DIR="$2"; shift 2 ;;
        -work) WORK_DIR="$2"; shift 2 ;;
        --fetal_gender) FETAL_GENDER="$2"; shift 2 ;;
        --no-log) ENABLE_LOGGING=false; shift ;;
        --detached) DETACHED_MODE=true; shift ;;
        -f|--force) FORCE_EXECUTION=true; shift ;;
        -h|--help) usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
done

# 필수 인자 확인
if [[ -z "$SAMPLE_ID" || -z "$LABCODE" || -z "$ROOT_DIR" || -z "$WORK_DIR" ]]; then
    echo "[ERROR] Missing required arguments"
    usage
fi

# 디렉토리 설정
HOST_ANALYSIS_DIR="$ROOT_DIR/analysis/$WORK_DIR"
HOST_OUTPUT_DIR="$ROOT_DIR/output/$WORK_DIR"
HOST_LOG_DIR="$ROOT_DIR/log"
HOST_DATA_DIR="$ROOT_DIR/data"
HOST_CONFIG_DIR="$ROOT_DIR/config"

# 디렉토리 생성
mkdir -p "$HOST_ANALYSIS_DIR"
mkdir -p "$HOST_OUTPUT_DIR"
mkdir -p "$HOST_LOG_DIR"
mkdir -p "$HOST_LOG_DIR/$WORK_DIR"

# Docker 설정
DOCKER_BIN="${DOCKER_BIN:-docker}"

# 로그 파일
RUN_LOG="$HOST_LOG_DIR/$WORK_DIR/${SAMPLE_ID}_md_run.log"

# 로깅 시작
if [ "$ENABLE_LOGGING" = true ]; then
    exec > >(tee -a "$RUN_LOG")
    exec 2>&1
fi

echo "========================================="
echo "MD Pipeline Runner"
echo "========================================="
echo "Sample ID: $SAMPLE_ID"
echo "Lab Code: $LABCODE"
echo "Work Dir: $WORK_DIR"
echo "Root Dir: $ROOT_DIR"
if [[ -n "$FETAL_GENDER" ]]; then
    echo "Fetal Gender: $FETAL_GENDER"
fi
echo "========================================="

# 마커 파일 확인
PIPELINE_MARKER="$HOST_ANALYSIS_DIR/${SAMPLE_ID}/${SAMPLE_ID}.md_pipeline_completed.marker"

if [[ "$FORCE_EXECUTION" = false && -f "$PIPELINE_MARKER" ]]; then
    echo "=== SKIPPING: MD analysis already completed. Use -f to force ==="
    exit 0
fi

# Force 모드: 마커 삭제
if [ "$FORCE_EXECUTION" = true ]; then
    echo "=== FORCE MODE: Removing previous marker and results ==="
    rm -f "$PIPELINE_MARKER"
fi

# BAM 파일 존재 확인
SAMPLE_DIR="$HOST_ANALYSIS_DIR/$SAMPLE_ID"
BAM_FILE=""

if [[ -f "$SAMPLE_DIR/${SAMPLE_ID}.proper_paired.bam" ]]; then
    BAM_FILE="$SAMPLE_DIR/${SAMPLE_ID}.proper_paired.bam"
elif [[ -f "$SAMPLE_DIR/${SAMPLE_ID}.sorted.bam" ]]; then
    BAM_FILE="$SAMPLE_DIR/${SAMPLE_ID}.sorted.bam"
else
    echo "[ERROR] No BAM file found for sample $SAMPLE_ID"
    echo "  Expected: $SAMPLE_DIR/${SAMPLE_ID}.proper_paired.bam or .sorted.bam"
    exit 1
fi

echo "BAM file found: ✓"

# Docker 실행
echo "=== Launching Docker container for MD analysis ==="
CONTAINER_NAME="md_${SAMPLE_ID}"

USER_UID="$(id -u)"
USER_GID="$(id -g)"
USER_NAME="${USER:-ken}"

# Docker arguments 구성
DOCKER_ARGS=(
    "--sample_id" "$SAMPLE_ID"
    "--work_dir" "$WORK_DIR"
    "--labcode" "$LABCODE"
    "--analysis_dir" "/Work/NIPT/analysis"
    "--output_dir" "/Work/NIPT/output"
    "--data_dir" "/Work/NIPT/data"
)

# Fetal gender 추가 (옵션)
if [[ -n "$FETAL_GENDER" ]]; then
    DOCKER_ARGS+=("--fetal_gender" "$FETAL_GENDER")
fi

CONTAINER_ID=$("$DOCKER_BIN" run --rm -d \
    --user "${USER_UID}:${USER_GID}" \
    --name "$CONTAINER_NAME" \
    -e TZ=Asia/Seoul \
    -e USER="$USER_NAME" \
    -e USERNAME="$USER_NAME" \
    -e HOME=/tmp \
    -e FONTCONFIG_PATH=/tmp \
    -v "$HOST_ANALYSIS_DIR:/Work/NIPT/analysis/$WORK_DIR" \
    -v "$HOST_OUTPUT_DIR:/Work/NIPT/output/$WORK_DIR" \
    -v "$HOST_LOG_DIR:/Work/NIPT/log" \
    -v "$HOST_DATA_DIR:/Work/NIPT/data" \
    -v "$HOST_CONFIG_DIR:/Work/NIPT/config" \
    -e PYTHON2="python2.7" \
    -e WC="/opt/wisecondor/wisecondor.py" \
    -e WCX="wisecondorx" \
    --entrypoint python3 \
    nipt_docker_v1.1 \
    /Work/NIPT/bin/scripts/md_pipeline.py \
    "${DOCKER_ARGS[@]}"
)

DOCKER_EXIT_CODE=$?

if [ $DOCKER_EXIT_CODE -ne 0 ]; then
    echo "[ERROR] Failed to launch Docker container"
    exit 1
fi

# Wait for container
echo "Waiting for MD analysis to complete..."
"$DOCKER_BIN" wait "$CONTAINER_NAME" 2>/dev/null || true

echo "========================================="
echo "MD Pipeline completed successfully!"
echo "========================================="
echo "Sample ID: $SAMPLE_ID"
echo "Results:"
echo "  - WC output: $HOST_ANALYSIS_DIR/$SAMPLE_ID/Output_WC/"
echo "  - WCX output: $HOST_ANALYSIS_DIR/$SAMPLE_ID/Output_WCX/"
echo "  - Plots: $HOST_ANALYSIS_DIR/$SAMPLE_ID/plots/"
echo "  - Log: $HOST_ANALYSIS_DIR/$SAMPLE_ID/${SAMPLE_ID}_md_analysis.log"
echo "========================================="
echo "Run log: $RUN_LOG"
echo "Analysis log: $HOST_ANALYSIS_DIR/$SAMPLE_ID/${SAMPLE_ID}_md_analysis.log"

