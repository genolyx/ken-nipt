#!/bin/bash
set -euo pipefail

# 스크립트 사용법 출력 함수
usage() {
    echo "Usage: $0 -s <sample_name> -1 <fastq_r1> -2 <fastq_r2> -l <labcode> -a <age> -root <root_directory> -work <work_directory> [--no-log] [--detached] [-f] [-rb] [-cf] [-ao] [-h]"
    exit 1
}

# 기본값 설정
ENABLE_LOGGING=true
DETACHED_MODE=false
FORCE_EXECUTION=false
CLEAN_FORCE=false
ALGORITHM_ONLY=false
REMOVE_BAMS=false
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
        -rb|--remove_bams) REMOVE_BAMS=true; shift ;;
        -h) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

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
    docker rm -f "$SAMPLE_NAME" 2>/dev/null || true

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
    docker rm -f "$SAMPLE_NAME" 2>/dev/null || true
    rm -f "$JSON_OUTPUT"
    rm -f "$HTML_OUTPUT"
    rm -f "$TAR_OUTPUT"
fi

# 디렉토리 생성
mkdir -p "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"
chown -R ken:genolyx "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR"

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
# fastq 파일 이름만 추출 (파이프라인이 파일 이름만 받음)
# create_symbolic_links 함수는 {FASTQ_DIR}/{sample_name}/{fastq_r1} 형태로 경로를 만듦
# 즉, /Work/NIPT/fastq/{sample_name}/파일이름 형태로 찾음

# Docker 추가 볼륨 배열 초기화
DOCKER_EXTRA_VOLUMES=()

if ! $ALGORITHM_ONLY; then
    # fastq 파일이 /data/fastq_backup에 있으면 실제 샘플 디렉토리를 직접 마운트
    if [[ "$FASTQ_R1" == /data/fastq_backup/* ]]; then
        # 실제 샘플 디렉토리 경로 추출 (예: /data/fastq_backup/2411/OPC241100002)
        SAMPLE_FASTQ_DIR=$(dirname "$FASTQ_R1")
        
        # 파일 이름만 추출
        FASTQ_R1_BASENAME=$(basename "$FASTQ_R1")
        FASTQ_R2_BASENAME=$(basename "$FASTQ_R2")
        
        echo "[INFO] Sample fastq directory: $SAMPLE_FASTQ_DIR"
        echo "[INFO] Fastq files: $FASTQ_R1_BASENAME, $FASTQ_R2_BASENAME"
        
        # Docker에 전달할 파일 이름 (basename만)
        CONTAINER_FASTQ_R1="$FASTQ_R1_BASENAME"
        CONTAINER_FASTQ_R2="$FASTQ_R2_BASENAME"
        
        # 실제 샘플 디렉토리를 Docker 볼륨으로 추가 마운트
        # /data/fastq_backup/2411/OPC241100002 -> /Work/NIPT/fastq/OPC241100002
        DOCKER_EXTRA_VOLUMES+=("-v" "$SAMPLE_FASTQ_DIR:/Work/NIPT/fastq/$SAMPLE_NAME")
    else
        # /data/fastq_backup가 아닌 경우는 그대로 사용
        CONTAINER_FASTQ_R1=$(basename "$FASTQ_R1")
        CONTAINER_FASTQ_R2=$(basename "$FASTQ_R2")
    fi
else
    # algorithm-only mode면 fake FASTQ 이름을 sample_name 기반으로 설정
    CONTAINER_FASTQ_R1="${SAMPLE_NAME}_R1.fastq.gz"
    CONTAINER_FASTQ_R2="${SAMPLE_NAME}_R2.fastq.gz"
fi

DOCKER_ARGS=(--sample_name "$SAMPLE_NAME" --labcode "$LABCODE" --age "$AGE")
if ! $ALGORITHM_ONLY; then
    # 파일 이름만 전달 (파이프라인이 {FASTQ_DIR}/{sample_name}/{filename} 형태로 경로를 만듦)
    DOCKER_ARGS+=(--fastq_r1 "$CONTAINER_FASTQ_R1" --fastq_r2 "$CONTAINER_FASTQ_R2")
else
    # algorithm-only 모드에도 기본값으로 넘겨주기
    DOCKER_ARGS+=(--fastq_r1 "$CONTAINER_FASTQ_R1" --fastq_r2 "$CONTAINER_FASTQ_R2" --algorithm_only)
fi

# Docker 실행
echo "=== Launching Docker container ==="
CONTAINER_ID=$(docker run -d \
    --security-opt seccomp=unconfined \
    --user "$(id -u):$(id -g)" \
    --name "$SAMPLE_NAME" \
    -e TZ=Asia/Seoul \
    -e USER=$(whoami) \
    -e USERNAME=$(whoami) \
    -e HOME=/tmp \
    -e FONTCONFIG_PATH=/tmp \
    -v "$HOST_FASTQ_DIR:/Work/NIPT/fastq" \
    -v "/data/fastq_backup:/Work/NIPT/fastq_backup" \
    "${DOCKER_EXTRA_VOLUMES[@]}" \
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
    nipt_docker_v1.2 \
    "${DOCKER_ARGS[@]}"
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
    echo "Note: Analysis directory cleanup will be handled by batch script after completion"
else
    echo "Waiting for container to complete..."
    sleep 2
    docker wait "$SAMPLE_NAME"
    CONTAINER_EXIT_CODE=$?

    COMPLETED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.completed"
    JSON_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.json"
    FAILED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.failed"
    PROGRESS_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}_progress.txt"
    ANALYSIS_LOG_FILE="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}_analysis.log"
    ANALYSIS_SAMPLE_DIR="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}"

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    if [ $CONTAINER_EXIT_CODE -eq 0 ]; then
        if [ -f "$JSON_FILE" ]; then
            echo "Pipeline completed successfully: $SAMPLE_NAME"
            
            # BAM 정리 모드가 켜져 있으면 cleanup_bam.sh 호출
            if [ "$REMOVE_BAMS" = true ]; then
                echo "=== REMOVE_BAMS: cleaning up BAM files ==="
                bash "$SCRIPT_DIR/cleanup_bam.sh" "$HOST_ANALYSIS_DIR" "$SAMPLE_NAME" "true"
            fi
            
            # v2 추가 기능: 분석 완료 후 analysis dir의 샘플 디렉토리 완전히 제거
            # output에 결과가 정상적으로 생성되었는지 확인 후 제거
            if [ -f "$JSON_FILE" ] && [ -d "$ANALYSIS_SAMPLE_DIR" ]; then
                echo "=== Removing analysis directory after successful completion ==="
                echo "Output JSON exists: $JSON_FILE"
                echo "Removing analysis directory: $ANALYSIS_SAMPLE_DIR"
                rm -rf "$ANALYSIS_SAMPLE_DIR"
                if [ $? -eq 0 ]; then
                    echo "Successfully removed analysis directory for $SAMPLE_NAME"
                else
                    echo "Warning: Failed to remove analysis directory for $SAMPLE_NAME"
                fi
            fi
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

