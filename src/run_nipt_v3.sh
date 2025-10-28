#!/bin/bash
set -euo pipefail

# 스크립트 사용법 출력 함수
usage() {
    echo "Usage: $0 -s <sample_name> -1 <fastq_r1> -2 <fastq_r2> -l <labcode> -a <age> -root <root_directory> -work <work_directory> [--no-log] [--detached] [-f] [-cf] [-ao] [-h]"
    exit 1
}

# 기본값 설정
ENABLE_LOGGING=true
DETACHED_MODE=false
FORCE_EXECUTION=false
CLEAN_FORCE=false
ALGORITHM_ONLY=false
FASTQ_R1=""
FASTQ_R2=""

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s) SAMPLE_NAME="$2"; shift 2 ;;
        -1) FASTQ_R1="$2"; shift 2 ;;
        -2) FASTQ_R2="$2"; shift 2 ;;
        -l) LABCODE="$2"; shift 2 ;;
        -a) AGE="$2"; shift 2 ;;
        -root) ROOT_DIR="$2"; shift 2 ;;
        -work) WORK_DIR="$2"; shift 2 ;;
        --no-log) ENABLE_LOGGING=false; shift ;;
        --detached) DETACHED_MODE=true; shift ;;
        -f|--force) FORCE_EXECUTION=true; shift ;;
        -cf|--clean_force) CLEAN_FORCE=true; shift ;;
        -ao|--algorithm_only) ALGORITHM_ONLY=true; shift ;;
        -h) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# docker 바이너리 탐색(절대경로 우선)
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

export PATH="/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:${PATH:-}"

# 필수 인자 기본 확인 (FASTQ는 AO 모드에선 제외)
if [[ -z "${SAMPLE_NAME-}" || -z "${LABCODE-}" || -z "${AGE-}" || -z "${ROOT_DIR-}" || -z "${WORK_DIR-}" ]]; then
    echo "[ERROR] Missing required arguments." >&2
    usage
fi

# full run 모드에서만 FASTQ 체크
if ! $ALGORITHM_ONLY; then
    if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" ]]; then
        echo "[ERROR] In full-run mode, -1 FASTQ_R1 and -2 FASTQ_R2 are required." >&2
        usage
    fi
else
    echo "[INFO] Algorithm-only mode (-ao): skipping FASTQ input check"
fi

# 디렉토리 설정
HOST_FASTQ_DIR="$ROOT_DIR/fastq/$WORK_DIR"
HOST_ANALYSIS_DIR="$ROOT_DIR/analysis/$WORK_DIR"
HOST_OUTPUT_DIR="$ROOT_DIR/output/$WORK_DIR"
HOST_LOG_DIR="$ROOT_DIR/log"
HOST_DATA_DIR="$ROOT_DIR/data"
HOST_CONFIG_DIR="$ROOT_DIR/config"

# 로그 설정
if [ "$ENABLE_LOGGING" = true ]; then
    LOG_DIR_FULL="$HOST_LOG_DIR/$WORK_DIR"
    RUN_LOG_FILE="$LOG_DIR_FULL/${SAMPLE_NAME}_run.log"
    mkdir -p "$LOG_DIR_FULL"
    exec > >(tee "$RUN_LOG_FILE") 2>&1
fi

# 마커 파일 확인
PIPELINE_MARKER="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.pipeline_completed.marker"
JSON_OUTPUT="$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.json"
HTML_OUTPUT="$HOST_OUTPUT_DIR//${SAMPLE_NAME}_report.html"
TAR_OUTPUT="$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.output.tar"

if [[ "$FORCE_EXECUTION" = false && "$CLEAN_FORCE" = false && -f "$PIPELINE_MARKER" ]]; then
    echo "=== SKIPPING: Already completed. Use -f to force ==="
    exit 0
fi

# clean_force 모드: 이전 결과만 지우고 마커도 삭제
if [ "$CLEAN_FORCE" = true ]; then
    echo "=== CLEAN_FORCE MODE: Removing previous results & marker ==="
    # 마커 삭제
    [ -f "$PIPELINE_MARKER" ] && rm -f "$PIPELINE_MARKER"

    # Docker 컨테이너 강제 삭제
    "$DOCKER_BIN" rm -f "$SAMPLE_NAME" 2>/dev/null || true

    # 분석 결과 디렉토리 내 clean 대상
    ANALYSIS_BASE="$HOST_ANALYSIS_DIR/$SAMPLE_NAME"
    rm -rf \
        "$ANALYSIS_BASE/Output_EZD"/* \
        "$ANALYSIS_BASE/Output_PRIZM"/* \
        "$ANALYSIS_BASE/Output_WC"/* \
        "$ANALYSIS_BASE/Output_WCX"/*

    OUTPUT_BASE="$HOST_OUTPUT_DIR/$SAMPLE_NAME"
    rm -rf \
        "$OUTPUT_BASE/Output_EZD"/* \
        "$OUTPUT_BASE/Output_PRIZM"/* \
        "$OUTPUT_BASE/Output_QC"/* \
        "$OUTPUT_BASE/Output_WC"/* \
        "$OUTPUT_BASE/Output_WCX"/*

    # 최종 JSON (output_dir) 삭제
    rm -f "$JSON_OUTPUT"
    rm -f "$HTML_OUTPUT"
    rm -f "$TAR_OUTPUT"
fi

if [ "$FORCE_EXECUTION" = true ]; then
    echo "=== FORCE MODE: Removing previous marker and container ==="
    [ -f "$PIPELINE_MARKER" ] && rm -f "$PIPELINE_MARKER"
    "$DOCKER_BIN" rm -f "$SAMPLE_NAME" 2>/dev/null || true
    rm -f "$JSON_OUTPUT"
    rm -f "$HTML_OUTPUT"
    rm -f "$TAR_OUTPUT"
fi

# 디렉토리 생성
#mkdir -p "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_DATA_DIR/bed" "$HOST_OUTPUT_DIR"
#chown -R ken:ken "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_DATA_DIR" "$HOST_OUTPUT_DIR"
#mkdir -p "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"
#chown -R ken:ken "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"
if [ ! -d "$HOST_ANALYSIS_DIR" ] || [ ! -d "$HOST_LOG_DIR" ] || [ ! -d "$HOST_OUTPUT_DIR" ]; then
  mkdir -p "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"
  chown -R ken:ken "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"
fi

if ! $ALGORITHM_ONLY; then
    if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" ]]; then
        echo "[ERROR] In full-run mode, -1 FASTQ_R1 and -2 FASTQ_R2 are required." >&2
        usage
    fi
else
    echo "[INFO] Algorithm-only mode: skipping FASTQ input check"
    # algorithm-only 모드면 fake FASTQ 이름을 sample_name 기반으로 설정
    FASTQ_R1="${SAMPLE_NAME}_R1.fastq.gz"
    FASTQ_R2="${SAMPLE_NAME}_R2.fastq.gz"
fi

# Docker 인자 설정
#OPTIONAL_ARGS=()
#[ "$ALGORITHM_ONLY" = true ] && OPTIONAL_ARGS+=("--algorithm_only")

DOCKER_ARGS=(--sample_name "$SAMPLE_NAME" --labcode "$LABCODE" --age "$AGE")
if ! $ALGORITHM_ONLY; then
    DOCKER_ARGS+=(--fastq_r1 "$FASTQ_R1" --fastq_r2 "$FASTQ_R2")
else
    # algorithm-only 모드에도 기본값으로 넘겨주기
    DOCKER_ARGS+=(--fastq_r1 "$FASTQ_R1" --fastq_r2 "$FASTQ_R2" --algorithm_only)
fi

USER_UID="$(id -u)"
USER_GID="$(id -g)"
USER_NAME="${USER:-ken}"

# Docker 실행
echo "=== Launching Docker container ==="
CONTAINER_ID=$("$DOCKER_BIN" run --rm -d \
    --user "${USER_UID}:${USER_GID}" \
    --name "$SAMPLE_NAME" \
    -e TZ=Asia/Seoul \
    -e USER="$USER_NAME" \
    -e USERNAME=="$USER_NAME" \
    -e HOME=/tmp \
    -e FONTCONFIG_PATH=/tmp \
    -v "$HOST_FASTQ_DIR:/Work/NIPT/fastq" \
    -v "$HOST_ANALYSIS_DIR:/Work/NIPT/analysis" \
    -v "$HOST_LOG_DIR:/Work/NIPT/log" \
    -v "$HOST_DATA_DIR:/Work/NIPT/data" \
    -v "$HOST_CONFIG_DIR:/Work/NIPT/config" \
    -v "$HOST_OUTPUT_DIR:/Work/NIPT/output" \
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
    nipt_docker_v1.1 \
    "${DOCKER_ARGS[@]}"
    #--sample_name "$SAMPLE_NAME" \
    #--fastq_r1 "$FASTQ_R1" \
    #--fastq_r2 "$FASTQ_R2" \
    #--labcode "$LABCODE" \
    #--age "$AGE" \
    #"${OPTIONAL_ARGS[@]}"
)

DOCKER_EXIT_CODE=$?

if [ $DOCKER_EXIT_CODE -ne 0 ]; then
    echo "Docker run failed. Exit code: $DOCKER_EXIT_CODE" >&2
    exit 1
fi

# wait 모드 처리
if [ "$DETACHED_MODE" = true ]; then
    echo "Container is running in detached mode (ID: $CONTAINER_ID)"
    echo "Check progress with: docker logs -f $SAMPLE_NAME"
else
    echo "Waiting for container to complete..."
    sleep 2
    "$DOCKER_BIN" wait "$SAMPLE_NAME"
    CONTAINER_EXIT_CODE=$?

    #COMPLETED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.completed"
    JSON_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.json"
    FAILED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.failed"
    PROGRESS_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}_progress.txt"
    ANALYSIS_LOG_FILE="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}_analysis.log"

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    if [ $CONTAINER_EXIT_CODE -eq 0 ]; then
        if [ -f "$JSON_FILE" ]; then
            echo "Pipeline completed successfully: $SAMPLE_NAME"
        else
            echo "Pipeline finished, but no completion marker found"
            echo "Possible error in internal pipeline. Marking as failed."
            echo "Failure timestamp: $(date)" > "$FAILED_FILE"
        fi
    else
        echo "Docker container exited with error code: $CONTAINER_EXIT_CODE"
        echo "Failure timestamp: $(date)" > "$FAILED_FILE"
    fi

    if [ "$ENABLE_LOGGING" = true ]; then
        echo "Logs: $RUN_LOG_FILE"
        [ -f "$ANALYSIS_LOG_FILE" ] && echo "Analysis log: $ANALYSIS_LOG_FILE"
    fi
fi

exit 0
