#!/bin/bash

# 스크립트 사용법 출력 함수
usage() {
    echo "Usage: $0 -s <sample_name> -1 <fastq_r1> -2 <fastq_r2> -l <labcode> -a <age> -root <root_directory> -work <work_directory> [--no-log] [--detached] [-f] [--algorithm_only] [-h]"
    exit 1
}

# 기본값 설정
ENABLE_LOGGING=true
DETACHED_MODE=false
FORCE_EXECUTION=false

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
        -ao|--algorithm_only) ALGORITHM_ONLY=true; shift ;;
        -h) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

# 필수 인자 확인
if [[ -z "$SAMPLE_NAME" || -z "$FASTQ_R1" || -z "$FASTQ_R2" || -z "$LABCODE" || -z "$AGE" || -z "$ROOT_DIR" || -z "$WORK_DIR" ]]; then
    echo "Error: Missing required arguments" >&2
    usage
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
if [ "$FORCE_EXECUTION" = false ] && [ -f "$PIPELINE_MARKER" ]; then
    echo "=== SKIPPING: Already completed. Use -f to force ==="
    exit 0
fi

if [ "$FORCE_EXECUTION" = true ]; then
    echo "=== FORCE MODE: Removing previous marker and container ==="
    [ -f "$PIPELINE_MARKER" ] && rm -f "$PIPELINE_MARKER"
    docker rm -f "$SAMPLE_NAME" 2>/dev/null || true
fi

# 디렉토리 생성
mkdir -p "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_DATA_DIR/bed" "$HOST_OUTPUT_DIR"
chown -R ken:ken "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_DATA_DIR" "$HOST_OUTPUT_DIR"

# FASTQ 확인
if [[ ! -f "$HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R1" || ! -f "$HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R2" ]]; then
    echo "FASTQ 파일이 존재하지 않습니다: $FASTQ_R1 또는 $FASTQ_R2" >&2
    exit 1
fi

# Docker 인자 설정
OPTIONAL_ARGS=()
[ "$ALGORITHM_ONLY" = true ] && OPTIONAL_ARGS+=("--algorithm_only")

# Docker 실행
echo "=== Launching Docker container ==="
CONTAINER_ID=$(docker run -d \
    --user "$(id -u):$(id -g)" \
    --name "$SAMPLE_NAME" \
    -e TZ=Asia/Seoul \
    -e USER=$(whoami) \
    -e USERNAME=$(whoami) \
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
    nipt_docker_v1.0 \
    --sample_name "$SAMPLE_NAME" \
    --fastq_r1 "$FASTQ_R1" \
    --fastq_r2 "$FASTQ_R2" \
    --labcode "$LABCODE" \
    --age "$AGE" \
    "${OPTIONAL_ARGS[@]}"
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
    docker wait "$SAMPLE_NAME"
    CONTAINER_EXIT_CODE=$?

    COMPLETED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.completed"
    FAILED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.failed"
    PROGRESS_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}_progress.txt"
    ANALYSIS_LOG_FILE="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}_analysis.log"

    if [ $CONTAINER_EXIT_CODE -eq 0 ]; then
        if [ -f "$COMPLETED_FILE" ]; then
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
