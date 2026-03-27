#!/bin/bash
set -euo pipefail

usage() {
  echo "Usage: $0 -s <sample_name> -l <labcode> -a <age> -root <root_directory> -work <work_directory> -b <proper_paired.bam> [--config-dir <dir>] [--detached] [-f]"
  exit 1
}

DETACHED_MODE=false
FORCE_EXECUTION=false
INPUT_BAM=""
CONFIG_DIR_OVERRIDE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -s) SAMPLE_NAME="$2"; shift 2 ;;
    -l) LABCODE="$2"; shift 2 ;;
    -a) AGE="$2"; shift 2 ;;
    -root) ROOT_DIR="$2"; shift 2 ;;
    -work) WORK_DIR="$2"; shift 2 ;;
    -b|--bam) INPUT_BAM="$2"; shift 2 ;;
    --config-dir) CONFIG_DIR_OVERRIDE="$2"; shift 2 ;;
    --detached) DETACHED_MODE=true; shift ;;
    -f|--force) FORCE_EXECUTION=true; shift ;;
    -h) usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "${SAMPLE_NAME-}" || -z "${LABCODE-}" || -z "${AGE-}" || -z "${ROOT_DIR-}" || -z "${WORK_DIR-}" || -z "${INPUT_BAM-}" ]]; then
  echo "[ERROR] Missing required arguments." >&2
  usage
fi

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "[FATAL] Input BAM not found: $INPUT_BAM" >&2
  exit 2
fi

# docker binary
DOCKER_BIN="${DOCKER_BIN:-}"
if [ -z "${DOCKER_BIN}" ]; then
  if [ -x /usr/bin/docker ]; then
    DOCKER_BIN="/usr/bin/docker"
  elif command -v docker >/dev/null 2>&1; then
    DOCKER_BIN="$(command -v docker)"
  else
    echo "[FATAL] docker not found. PATH=$PATH" >&2
    exit 127
  fi
fi

HOST_FASTQ_DIR="$ROOT_DIR/fastq/$WORK_DIR"
HOST_ANALYSIS_DIR="$ROOT_DIR/analysis/$WORK_DIR"
HOST_OUTPUT_DIR="$ROOT_DIR/output/$WORK_DIR"
HOST_LOG_DIR="$ROOT_DIR/log"
HOST_DATA_DIR="$ROOT_DIR/data"
HOST_CONFIG_DIR="${CONFIG_DIR_OVERRIDE:-$ROOT_DIR/config}"

if [[ -n "$CONFIG_DIR_OVERRIDE" && ! -d "$HOST_CONFIG_DIR" ]]; then
  echo "[FATAL] --config-dir does not exist: $HOST_CONFIG_DIR" >&2
  exit 2
fi

mkdir -p "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_OUTPUT_DIR" "$HOST_LOG_DIR"

# Decide the container-visible path for the input BAM.
# - If INPUT_BAM is already inside the mounted analysis dir, pass it directly (fast; no extra copy).
# - Otherwise, copy it into analysis/_proper_paired_inputs so the container can read it.
CONTAINER_PP_BAM=""
if [[ "$INPUT_BAM" == "$HOST_ANALYSIS_DIR/"* ]]; then
  REL_PATH="${INPUT_BAM#"$HOST_ANALYSIS_DIR/"}"
  CONTAINER_PP_BAM="/Work/NIPT/analysis/$REL_PATH"
  echo "[INFO] Using mounted input BAM directly: $INPUT_BAM -> $CONTAINER_PP_BAM"
else
  STAGED_DIR="$HOST_ANALYSIS_DIR/_proper_paired_inputs/$SAMPLE_NAME"
  mkdir -p "$STAGED_DIR"
  STAGED_BAM="$STAGED_DIR/${SAMPLE_NAME}.proper_paired.bam"
  echo "[INFO] Copying input BAM into mounted analysis dir: $STAGED_BAM"
  cp -f "$INPUT_BAM" "$STAGED_BAM"
  if [[ -f "${INPUT_BAM}.bai" ]]; then
    cp -f "${INPUT_BAM}.bai" "${STAGED_BAM}.bai" || true
  fi
  CONTAINER_PP_BAM="/Work/NIPT/analysis/_proper_paired_inputs/$SAMPLE_NAME/${SAMPLE_NAME}.proper_paired.bam"
fi

if [[ "$FORCE_EXECUTION" = true ]]; then
  rm -f "$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.pipeline_completed.marker" || true
  "$DOCKER_BIN" rm -f "$SAMPLE_NAME" 2>/dev/null || true
  rm -f "$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.json" || true
  rm -f "$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.html" || true
  rm -f "$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.output.tar" || true
fi

# Provide dummy FASTQ names (required by argparse, but not used in from_proper_paired mode)
FAKE_R1="${SAMPLE_NAME}_R1.fastq.gz"
FAKE_R2="${SAMPLE_NAME}_R2.fastq.gz"

USER_UID="$(id -u)"
USER_GID="$(id -g)"
USER_NAME="${USER:-ken}"

DOCKER_ARGS=(
  --sample_name "$SAMPLE_NAME"
  --fastq_r1 "$FAKE_R1"
  --fastq_r2 "$FAKE_R2"
  --labcode "$LABCODE"
  --age "$AGE"
  --from_proper_paired
  --proper_paired_bam "$CONTAINER_PP_BAM"
)

echo "=== Launching Docker container (from proper_paired.bam) ==="
CONTAINER_ID=$("$DOCKER_BIN" run -d \
  --user "${USER_UID}:${USER_GID}" \
  --name "$SAMPLE_NAME" \
  -e TZ=Asia/Seoul \
  -e USER="$USER_NAME" \
  -e USERNAME="$USER_NAME" \
  -e HOME=/tmp \
  -e FONTCONFIG_PATH=/tmp \
  -v "$HOST_FASTQ_DIR:/Work/NIPT/fastq" \
  -v "$HOST_ANALYSIS_DIR:/Work/NIPT/analysis" \
  -v "$HOST_LOG_DIR:/Work/NIPT/log" \
  -v "$HOST_DATA_DIR:/Work/NIPT/data" \
  -v "$HOST_CONFIG_DIR:/Work/NIPT/config" \
  -v "$HOST_OUTPUT_DIR:/Work/NIPT/output" \
  -v "$ROOT_DIR/bin/scripts/nipt_pipeline.py:/Work/NIPT/bin/nipt_pipeline.py" \
  -v "$ROOT_DIR/bin/scripts/modules:/Work/NIPT/bin/modules" \
  -e BWA2="bwa-mem2" \
  -e SAMTools="samtools" \
  -e PICARD="/Work/NIPT/bin/picard/picard.jar" \
  -e qualimap="qualimap" \
  -e QC.bwa_threads="8" \
  -e QC.samtools_threads="4" \
  -e QC.samtools_memory="6G" \
  -e QC.picard_memory="12G" \
  -e PYTHON2="python2.7" \
  -e WC="/opt/wisecondor/wisecondor.py" \
  -e WCX="wisecondorx" \
  -e WCFF="wisecondor-ff" \
  -e HMMcopy="/opt/conda/envs/nipt" \
  -e Rscript="Rscript" \
  -e DATA_DIR="/Work/NIPT/data" \
  -e ANALYSIS_DIR="/Work/NIPT/analysis" \
  -e OUTPUT_DIR="/Work/NIPT/output" \
  nipt_docker_v1.2 \
  "${DOCKER_ARGS[@]}"
)

if [[ -z "$CONTAINER_ID" ]]; then
  echo "[FATAL] Docker failed to start container" >&2
  exit 1
fi

if [ "$DETACHED_MODE" = true ]; then
  echo "Container is running in detached mode (ID: $CONTAINER_ID)"
  echo "Check progress with: docker logs -f $SAMPLE_NAME"
  exit 0
fi

echo "Waiting for container to complete..."
CONTAINER_EXIT_CODE="$("$DOCKER_BIN" wait "$SAMPLE_NAME")"

if [[ "$CONTAINER_EXIT_CODE" != "0" ]]; then
  echo "[ERROR] Container exited with code: $CONTAINER_EXIT_CODE" >&2
  echo "[INFO] Last 200 log lines:" >&2
  "$DOCKER_BIN" logs "$SAMPLE_NAME" 2>&1 | tail -200 >&2 || true
  echo "[INFO] Keeping container for debugging: $SAMPLE_NAME" >&2
  exit "$CONTAINER_EXIT_CODE"
fi

"$DOCKER_BIN" rm -f "$SAMPLE_NAME" >/dev/null 2>&1 || true
exit 0

