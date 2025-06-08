#!/bin/bash

# 스크립트 사용법 출력 함수
usage() {
    echo "Usage: $0 -s <sample_name> -1 <fastq_r1> -2 <fastq_r2> -l <labcode> -root <root_directory> -work <work_directory> [--no-log] [--detached]"
    echo "Options:"
    echo "  -s <sample_name>        Sample name"
    echo "  -1 <fastq_r1>           R1 FASTQ filename (not full path)"
    echo "  -2 <fastq_r2>           R2 FASTQ filename (not full path)"
    echo "  -l <labcode>            Labcode (cordlife, ucl, vn)"
    echo "  -root <root_directory>  Root directory (e.g. /home/ken/ken-nipt)"
    echo "  -work <work_directory>  Work directory name (e.g. 250430_01)"
    echo "  --no-log                Disable logging to file (output to stdout only)"
    echo "  --detached              Run in detached mode (for daemon use)"
    echo "  -f                      Force execution even if completed marker exists"
    echo "  -h                      Show this help message"
    exit 1
}

# 기본값 설정
ENABLE_LOGGING=true
DETACHED_MODE=false
FORCE_EXECUTION=false

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        -1)
            FASTQ_R1="$2"
            shift 2
            ;;
        -2)
            FASTQ_R2="$2"
            shift 2
            ;;
        -l)
            LABCODE="$2"
            shift 2
            ;;
        -root)
            ROOT_DIR="$2"
            shift 2
            ;;
        -work)
            WORK_DIR="$2"
            shift 2
            ;;
        --no-log)
            ENABLE_LOGGING=false
            shift
            ;;
        --detached)
            DETACHED_MODE=true
            shift
            ;;
        -f|--force)
            FORCE_EXECUTION=true
            shift
            ;;
        -h)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# 필수 인자 확인
if [ -z "$SAMPLE_NAME" ] || [ -z "$FASTQ_R1" ] || [ -z "$FASTQ_R2" ] || [ -z "$LABCODE" ] || [ -z "$ROOT_DIR" ] || [ -z "$WORK_DIR" ]; then
    echo "Error: Missing required arguments" 1>&2
    usage
fi

# 필수 디렉토리 설정
HOST_FASTQ_DIR="$ROOT_DIR/fastq/$WORK_DIR"
HOST_ANALYSIS_DIR="$ROOT_DIR/analysis/$WORK_DIR"
HOST_OUTPUT_DIR="$ROOT_DIR/output/$WORK_DIR"
HOST_LOG_DIR="$ROOT_DIR/log"
HOST_DATA_DIR="$ROOT_DIR/data"
HOST_CONFIG_DIR="$ROOT_DIR/config"

# 로그 설정 (인자 파싱 후, 디렉토리 생성 전에)
if [ "$ENABLE_LOGGING" = true ]; then
    # 로그 디렉토리와 파일 설정
    LOG_DIR_FULL="$HOST_LOG_DIR/$WORK_DIR"
    RUN_LOG_FILE="$LOG_DIR_FULL/${SAMPLE_NAME}_run.log"
    
    # 로그 디렉토리 생성
    mkdir -p "$LOG_DIR_FULL"
    
    # 모든 출력을 로그 파일과 화면에 동시 출력
    exec > >(tee "$RUN_LOG_FILE") 2>&1
    
    echo "=========================================="
    echo "=== NIPT Analysis Started ==="
    echo "Sample: $SAMPLE_NAME"
    echo "Work Dir: $WORK_DIR"
    echo "Run Log file: $RUN_LOG_FILE"
    echo "Timestamp: $(date)"
    echo "=========================================="
fi

# 디렉토리 생성 (없는 경우)
echo "Creating necessary directories..."
mkdir -p "$HOST_FASTQ_DIR"
mkdir -p "$HOST_ANALYSIS_DIR"
mkdir -p "$HOST_LOG_DIR"
mkdir -p "$HOST_DATA_DIR/bed"
mkdir -p "$HOST_OUTPUT_DIR"

chown -R ken:ken "$HOST_FASTQ_DIR"
chown -R ken:ken "$HOST_ANALYSIS_DIR"
chown -R ken:ken "$HOST_LOG_DIR"
chown -R ken:ken "$HOST_DATA_DIR"
chown -R ken:ken "$HOST_OUTPUT_DIR"


# FASTQ 파일 확인
echo "Checking FASTQ files..."
if [ ! -f "$HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R1" ]; then
    echo "Error: FASTQ file $HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R1 not found" 1>&2
    exit 1
fi

if [ ! -f "$HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R2" ]; then
    echo "Error: FASTQ file $HOST_FASTQ_DIR/$SAMPLE_NAME/$FASTQ_R2 not found" 1>&2
    exit 1
fi

echo "FASTQ files validated successfully."

# Docker 실행 정보 출력
echo ""
echo "=========================================="
echo "=== Docker Container Execution Info ==="
echo "Starting NIPT pipeline for sample: $SAMPLE_NAME"
echo "FASTQ R1: $FASTQ_R1"
echo "FASTQ R2: $FASTQ_R2"
echo "Labcode: $LABCODE"
echo "Root directory: $ROOT_DIR"
echo "Work directory: $WORK_DIR"
echo "FASTQ directory: $HOST_FASTQ_DIR"
echo "Config directory: $HOST_CONFIG_DIR"
echo "Analysis directory: $HOST_ANALYSIS_DIR"
echo "Output directory: $HOST_OUTPUT_DIR"
echo "=========================================="
echo ""

# Docker 실행 시작 시간 기록
echo "Docker execution started at: $(date)"

docker run \
    --user $(id -u):$(id -g) \
    --name "$SAMPLE_NAME" \
    -d \
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
    -v "$ROOT_DIR/config:/Work/NIPT/config" \
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
    nipt_docker_dev \
    --sample_name "$SAMPLE_NAME" \
    --fastq_r1 "$FASTQ_R1" \
    --fastq_r2 "$FASTQ_R2" \
    --labcode "$LABCODE" \
    --age 30 

# Docker 실행 결과 저장
DOCKER_EXIT_CODE=$?

# 실행 완료 시간 기록
echo ""
echo "Docker execution completed at: $(date)"

# Progress 파일 경로 설정
PROGRESS_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}_progress.txt"
COMPLETED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.completed"
FAILED_FILE="$HOST_OUTPUT_DIR/${SAMPLE_NAME}.failed"
ANALYSIS_LOG_FILE="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}_analysis.log"

# 여기에 추가: Docker 실행 전 완료 마커 체크
PIPELINE_MARKER="$HOST_ANALYSIS_DIR/${SAMPLE_NAME}/${SAMPLE_NAME}.pipeline_completed.marker"

# Force 옵션 확인 및 처리
if [ "$FORCE_EXECUTION" = false ]; then
    if [ -f "$PIPELINE_MARKER" ]; then
        echo "=========================================="
        echo "=== ALREADY COMPLETED ==="
        echo "Sample $SAMPLE_NAME has already been completed successfully."
        echo "Pipeline marker: $PIPELINE_MARKER"
        echo "Use -f option to force re-execution."
        echo "=========================================="
        exit 0
    fi
else
    echo "=========================================="
    echo "=== FORCE MODE ==="
    echo "Removing pipeline completion marker for $SAMPLE_NAME..."

    if [ -f "$PIPELINE_MARKER" ]; then
        rm -f "$PIPELINE_MARKER"
        echo "Removed: $PIPELINE_MARKER"
    else
        echo "No pipeline marker found: $PIPELINE_MARKER"
    fi

    # 기존 컨테이너 제거 (이름 충돌 방지)
    if docker ps -a --format "{{.Names}}" | grep -q "^${SAMPLE_NAME}$"; then
        docker rm -f "$SAMPLE_NAME" 2>/dev/null || true
        echo "Removed existing container: $SAMPLE_NAME"
    fi

    echo "Ready for fresh execution..."
    echo "=========================================="
fi

# 실행 상태 확인
if [ $DOCKER_EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "=== SUCCESS ==="
    
    if [ "$DETACHED_MODE" = true ]; then
        echo "Docker container started successfully for sample: $SAMPLE_NAME"
        echo "Container is running in background"
        echo "Use 'docker logs $SAMPLE_NAME' to monitor progress"
        
        if [ "$ENABLE_LOGGING" = true ]; then
            echo "Logs will be saved to: $RUN_LOG_FILE"
        fi
    else
        echo "NIPT pipeline started for sample: $SAMPLE_NAME"
        echo "Waiting for completion..."
        
        docker wait "$SAMPLE_NAME"
        CONTAINER_EXIT_CODE=$?

        if [ $CONTAINER_EXIT_CODE -eq 0 ]; then
            echo "NIPT pipeline completed successfully for sample: $SAMPLE_NAME" 
            # 완료 파일이 생성되지 않은 경우 (비정상 종료) 실패로 처리
            if [ ! -f "$COMPLETED_FILE" ]; then
                echo "Warning: Pipeline completed but no completion marker found"
                echo "Creating failure marker..."
                TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
                echo "Pipeline failed for $SAMPLE_NAME at $TIMESTAMP" > "$FAILED_FILE"
                echo "Reason: Docker exited successfully but pipeline did not complete normally" >> "$FAILED_FILE"
                
                # Progress 파일에도 기록
                if [ -f "$PROGRESS_FILE" ]; then
                    echo "Pipeline failed at $TIMESTAMP - Docker exited successfully but no completion marker" >> "$PROGRESS_FILE"
                fi
                
                echo "Check progress file: $PROGRESS_FILE"
                echo "Check log file: $RUN_LOG_FILE"
                exit 1
            fi
            
            if [ "$ENABLE_LOGGING" = true ]; then
                echo "Run log file saved: $RUN_LOG_FILE"
                if [ -f "$ANALYSIS_LOG_FILE" ]; then
                    echo "Analysis log file: $ANALYSIS_LOG_FILE"
                fi
            fi
        else
            echo "NIPT pipeline failed for sample: $SAMPLE_NAME"
            echo "Container exit code: $CONTAINER_EXIT_CODE"
            exit 1
        fi
    fi 
    echo "=========================================="
else
    echo "=========================================="
    echo "=== ERROR ==="
    echo "NIPT pipeline failed for sample: $SAMPLE_NAME"
    echo "Exit code: $DOCKER_EXIT_CODE"
    
    # 실패 파일 생성 (파이프라인에서 생성하지 못한 경우)
    if [ ! -f "$FAILED_FILE" ]; then
        TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
        echo "Pipeline failed for $SAMPLE_NAME at $TIMESTAMP" > "$FAILED_FILE"
        echo "Reason: Docker container exited with code $DOCKER_EXIT_CODE" >> "$FAILED_FILE"
        
        # Progress 파일에도 기록 (존재하는 경우)
        if [ -f "$PROGRESS_FILE" ]; then
            echo "Pipeline failed at $TIMESTAMP - Docker exit code: $DOCKER_EXIT_CODE" >> "$PROGRESS_FILE"
        else
            # Progress 파일도 없는 경우 생성
            echo "Pipeline started for $SAMPLE_NAME at unknown time" > "$PROGRESS_FILE"
            echo "Pipeline failed at $TIMESTAMP - Docker exit code: $DOCKER_EXIT_CODE" >> "$PROGRESS_FILE"
        fi
    fi
    
    if [ "$ENABLE_LOGGING" = true ]; then
        echo "Check run log file for details: $RUN_LOG_FILE"
        if [ -f "$ANALYSIS_LOG_FILE" ]; then
            echo "Check analysis log file: $ANALYSIS_LOG_FILE"
        fi
    fi
    echo "Check progress file: $PROGRESS_FILE"
    echo "=========================================="
    exit 1
fi

exit 0
