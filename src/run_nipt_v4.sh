#!/bin/bash
set -euo pipefail

########################################
# run_nipt_v4.sh  —  robust launcher (Ken 요청 반영)
# - 컨테이너 인자 단순화( --work/--force/--clean_force/--remove_bams 미전달 )
# - CLEAN_FORCE / FORCE 사전 정리 (호스트 측)
# - ANALYSIS/OUTPUT 경로는 $WORK_DIR/$SAMPLE_NAME 만 고려
# - detached 기본, non-detached 디버깅용
########################################

usage() {
  cat <<'USAGE'
Usage: run_nipt_v4.sh \
  -s <sample_name> \
  -1 <fastq_r1> \
  -2 <fastq_r2> \
  -l <labcode> \
  -a <age> \
  -root <root_directory> \
  -work <work_directory> \
  [--no-log] [--detached] [-f|--force] [-cf|--clean_force] [-ao|--algorithm_only] [-h]

Notes:
- In algorithm-only mode (-ao), FASTQ args are not required (names are still passed).
- <root_directory> example: /home/ken/ken-nipt
- <work_directory> example: 2507
USAGE
  exit 1
}

# ===== 기본값 =====
ENABLE_LOGGING=true
DETACHED_MODE=true
FORCE_EXECUTION=false
CLEAN_FORCE=false
ALGORITHM_ONLY=false

FASTQ_R1=""
FASTQ_R2=""
SAMPLE_NAME=""
LABCODE=""
AGE=""
ROOT_DIR=""
WORK_DIR=""

# 리소스/로그/ulimit (환경변수로 조정 가능)
CPUS="${CPUS:-8}"
MEMORY="${MEMORY:-24g}"
MEMORY_SWAP="${MEMORY_SWAP:-$MEMORY}"
NOFILE_ULIMIT="${NOFILE_ULIMIT:-65536}"
LOG_MAX_SIZE="${LOG_MAX_SIZE:-10m}"
LOG_MAX_FILE="${LOG_MAX_FILE:-3}"

# ===== 인자 파싱 =====
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s)    SAMPLE_NAME="$2"; shift 2 ;;
    -1)    FASTQ_R1="$2"; shift 2 ;;
    -2)    FASTQ_R2="$2"; shift 2 ;;
    -l)    LABCODE="$2"; shift 2 ;;
    -a)    AGE="$2"; shift 2 ;;
    -root) ROOT_DIR="$2"; shift 2 ;;
    -work) WORK_DIR="$2"; shift 2 ;;
    --no-log)     ENABLE_LOGGING=false; shift ;;
    --detached)   DETACHED_MODE=true; shift ;;          # 기본 true
    -f|--force)   FORCE_EXECUTION=true; shift ;;
    -cf|--clean_force) CLEAN_FORCE=true; shift ;;
    -ao|--algorithm_only) ALGORITHM_ONLY=true; shift ;;
    -h|--help) usage ;;
    *) echo "[ERROR] Unknown option: $1" >&2; usage ;;
  esac
done

# ===== 로깅 헬퍼 =====
ts() { date '+%Y-%m-%d %H:%M:%S'; }
log_info()  { $ENABLE_LOGGING && echo "[INFO  $(ts)] $*";  true; }
log_warn()  { $ENABLE_LOGGING && echo "[WARN  $(ts)] $*" >&2; true; }
log_error() { echo "[ERROR $(ts)] $*" >&2; }

# ===== 필수값 확인 =====
if [[ -z "${SAMPLE_NAME:-}" || -z "${LABCODE:-}" || -z "${AGE:-}" || -z "${ROOT_DIR:-}" || -z "${WORK_DIR:-}" ]]; then
  log_error "Missing required arguments (-s/-l/-a/-root/-work)"; usage
fi
if [[ "$ALGORITHM_ONLY" != true ]]; then
  if [[ -z "${FASTQ_R1:-}" || -z "${FASTQ_R2:-}" ]]; then
    log_error "FASTQ (-1/-2) are required unless -ao/--algorithm_only is set"; usage
  fi
fi

# ===== 호스트 경로 =====
HOST_ROOT="$ROOT_DIR"
HOST_FASTQ_DIR="$HOST_ROOT/fastq"
HOST_ANALYSIS_DIR="$HOST_ROOT/analysis"
HOST_LOG_DIR="$HOST_ROOT/log"
HOST_OUTPUT_DIR="$HOST_ROOT/output"

# 필요한 디렉토리 준비
mkdir -p "$HOST_FASTQ_DIR" "$HOST_ANALYSIS_DIR" "$HOST_LOG_DIR" "$HOST_OUTPUT_DIR" || true

# 샘플별 작업 경로 (A 타입만 사용)
ANALYSIS_BASE="$HOST_ANALYSIS_DIR/$WORK_DIR/$SAMPLE_NAME"
OUTPUT_BASE="$HOST_OUTPUT_DIR/$WORK_DIR/$SAMPLE_NAME"
mkdir -p "$ANALYSIS_BASE" "$OUTPUT_BASE" || true

# 마커/출력 파일 (경로는 환경에 맞게 조정 가능)
PIPELINE_MARKER="${PIPELINE_MARKER:-$ANALYSIS_BASE/.marker.done}"
JSON_OUTPUT="${JSON_OUTPUT:-$OUTPUT_BASE/${SAMPLE_NAME}.json}"
HTML_OUTPUT="${HTML_OUTPUT:-$OUTPUT_BASE/${SAMPLE_NAME}.html}"
TAR_OUTPUT="${TAR_OUTPUT:-$OUTPUT_BASE/${SAMPLE_NAME}.tar.gz}"

# ===== FASTQ 검사 (AO가 아니면) =====
if [[ "$ALGORITHM_ONLY" != true ]]; then
  for fq in "$FASTQ_R1" "$FASTQ_R2"; do
    if [[ -z "$fq" ]]; then
      log_error "FASTQ argument missing"; exit 1
    fi
    fq_host="$HOST_FASTQ_DIR/$WORK_DIR/$SAMPLE_NAME/$(basename "$fq")"
    if [[ ! -r "$fq_host" ]]; then
      log_error "FASTQ not found: $fq_host"
      exit 1
    fi
  done
else
  # AO에서도 이름은 컨테이너에 전달됨
  FASTQ_R1="${FASTQ_R1:-${SAMPLE_NAME}_R1.fastq.gz}"
  FASTQ_R2="${FASTQ_R2:-${SAMPLE_NAME}_R2.fastq.gz}"
fi

log_info "Launching pipeline"
log_info "sample=$SAMPLE_NAME lab=$LABCODE age=$AGE work=$WORK_DIR root=$HOST_ROOT detached=$DETACHED_MODE ao=$ALGORITHM_ONLY"
log_info "R1=$(basename "$FASTQ_R1") R2=$(basename "$FASTQ_R2")"

# ===== CLEAN_FORCE / FORCE 사전 정리 =====
if [[ "$CLEAN_FORCE" == true ]]; then
  echo "=== CLEAN_FORCE MODE: Removing previous results & marker ==="
  [ -f "$PIPELINE_MARKER" ] && rm -f -- "$PIPELINE_MARKER" || true
  docker rm -f "$SAMPLE_NAME" 2>/dev/null || true

  if [[ -d "$ANALYSIS_BASE" ]]; then
    rm -rf -- \
      "$ANALYSIS_BASE/Output_EZD"/* \
      "$ANALYSIS_BASE/Output_PRIZM"/* \
      "$ANALYSIS_BASE/Output_WC"/* \
      "$ANALYSIS_BASE/Output_WCX"/* 2>/dev/null || true
  fi

  if [[ -d "$OUTPUT_BASE" ]]; then
    rm -rf -- \
      "$OUTPUT_BASE/Output_EZD"/* \
      "$OUTPUT_BASE/Output_PRIZM"/* \
      "$OUTPUT_BASE/Output_QC"/* \
      "$OUTPUT_BASE/Output_WC"/* \
      "$OUTPUT_BASE/Output_WCX"/* 2>/dev/null || true
  fi

  rm -f -- "$JSON_OUTPUT" "$HTML_OUTPUT" "$TAR_OUTPUT" 2>/dev/null || true
fi

if [[ "$FORCE_EXECUTION" == true ]]; then
  echo "=== FORCE MODE: Removing previous marker and container ==="
  [ -f "$PIPELINE_MARKER" ] && rm -f -- "$PIPELINE_MARKER" || true
  docker rm -f "$SAMPLE_NAME" 2>/dev/null || true
  rm -f -- "$JSON_OUTPUT" "$HTML_OUTPUT" "$TAR_OUTPUT" 2>/dev/null || true
fi

# ===== 이미 실행 중이면 No-Op (idempotent) =====
if docker ps --format '{{.Names}}' | grep -qx "$SAMPLE_NAME"; then
  log_info "container '$SAMPLE_NAME' already running; noop"
  exit 0
fi
# 종료된 동명 컨테이너는 정리 (강제 삭제 실패해도 무시)
if docker ps -a --format '{{.Names}}' | grep -qx "$SAMPLE_NAME"; then
  docker rm -f "$SAMPLE_NAME" >/dev/null 2>&1 || true
fi

# ===== Docker 공통 플래그 =====
USER_UID="$(id -u)"
USER_GID="$(id -g)"
USER_NAME="${USER:-nipt}"

DOCKER_FLAGS=()
if [[ "$DETACHED_MODE" == true ]]; then
  #DOCKER_FLAGS+=( -d --rm )
  DOCKER_FLAGS+=( -d )
fi

DOCKER_FLAGS+=(
  --cpus="$CPUS"
  --memory="$MEMORY"
  --memory-swap="$MEMORY_SWAP"
  --ulimit "nofile=${NOFILE_ULIMIT}:${NOFILE_ULIMIT}"
  --log-driver=local
  --log-opt "max-size=${LOG_MAX_SIZE}"
  --log-opt "max-file=${LOG_MAX_FILE}"
  --user "${USER_UID}:${USER_GID}"
  --name "$SAMPLE_NAME"
  --label "nipt.sample=$SAMPLE_NAME"
  --label "nipt.work=$WORK_DIR"
  -e TZ=Asia/Taipei
  -e HOME=/tmp
  -e USER="$USER_NAME"
  -e USERNAME="$USER_NAME"
  -v "$HOST_FASTQ_DIR:/Work/NIPT/fastq:rw"
  -v "$HOST_ANALYSIS_DIR:/Work/NIPT/analysis:rw"
  -v "$HOST_LOG_DIR:/Work/NIPT/log:rw"
  -v "$HOST_OUTPUT_DIR:/Work/NIPT/output:rw"
)

DOCKER_IMAGE="${DOCKER_IMAGE:-nipt_docker_v1.0}"

# ===== 컨테이너에 넘길 인자(최소) =====
DOCKER_ARGS=( --sample_name "$SAMPLE_NAME" --labcode "$LABCODE" --age "$AGE" )
if [[ "$ALGORITHM_ONLY" == true ]]; then
  DOCKER_ARGS+=( --fastq_r1 "$FASTQ_R1" --fastq_r2 "$FASTQ_R2" --algorithm_only )
else
  DOCKER_ARGS+=( --fastq_r1 "$FASTQ_R1" --fastq_r2 "$FASTQ_R2" )
fi
# NOTE: --work / --force / --clean_force / --remove_bams 전달 안 함

echo "=== Launching Docker container ==="
set +e
CONTAINER_ID="$(docker run "${DOCKER_FLAGS[@]}" "$DOCKER_IMAGE" "${DOCKER_ARGS[@]}")"
rc=$?
set -e

if [[ "$DETACHED_MODE" == true ]]; then
  if [[ $rc -ne 0 || -z "${CONTAINER_ID:-}" ]]; then
    log_error "docker run (detached) failed (rc=$rc)."
    docker logs --tail=100 "$SAMPLE_NAME" 2>/dev/null || true
    exit 1
  fi

  # 즉시 종료 감지용 간단 체크
  sleep 1
  if ! docker ps --format '{{.Names}}' | grep -qx "$SAMPLE_NAME"; then
    log_error "Container exited immediately; tailing last 100 lines:"
    docker logs --tail=100 "$SAMPLE_NAME" 2>/dev/null || true
    EXIT_RC="$(docker inspect -f '{{.State.ExitCode}}' "$SAMPLE_NAME" 2>/dev/null || echo 1)"
    exit "${EXIT_RC:-1}"
  fi

  log_info "Detached start OK. name=$SAMPLE_NAME id=$CONTAINER_ID"
  echo "Hint: docker logs -f $SAMPLE_NAME"
  exit 0
else
  # 포그라운드 모드: 종료까지 wait
  if [[ $rc -ne 0 ]]; then
    log_error "docker run (foreground) failed (rc=$rc)."
    exit $rc
  fi

  if docker ps --format '{{.Names}}' | grep -qx "$SAMPLE_NAME"; then
    log_info "Container running; waiting for exit..."
    set +e
    docker wait "$SAMPLE_NAME" >/dev/null
    EXIT_RC="$(docker inspect -f '{{.State.ExitCode}}' "$SAMPLE_NAME" 2>/dev/null || echo 0)"
    set -e
    log_info "Container exited with rc=$EXIT_RC"
    # 디버깅 편의: 자동 정리
    docker rm -f "$SAMPLE_NAME" >/dev/null 2>&1 || true
    exit "${EXIT_RC:-0}"
  else
    log_warn "Container not running; printing last logs"
    docker logs --tail=200 "$SAMPLE_NAME" 2>/dev/null || true
    EXIT_RC="$(docker inspect -f '{{.State.ExitCode}}' "$SAMPLE_NAME" 2>/dev/null || echo 1)"
    exit "${EXIT_RC:-1}"
  fi
fi
