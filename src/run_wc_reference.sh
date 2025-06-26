#!/bin/bash
set -euo pipefail

# 사용법 출력
usage() {
    echo "Usage: $0 -t <target> -m <method> -r <root_dir> [-d docker_image]"
    echo "  -t <target>       : 예시) cordlife"
    echo "  -m <method>       : WC / WCX / WCFF"
    echo "  -r <root_dir>     : HOST base root directory (예: /home/ken/ken-nipt)"
    echo "  -d <docker_image> : (옵션) Docker image name (기본값: nipt_docker_v1.0)"
    exit 1
}

# 기본값
DOCKER_IMAGE="nipt_docker_v1.0"

# 인자 파싱
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t) TARGET="$2"; shift 2 ;;
        -m) METHOD="$2"; shift 2 ;;
        -r) ROOT_DIR="$2"; shift 2 ;;
        -d) DOCKER_IMAGE="$2"; shift 2 ;;
        -h) usage ;;
        *) echo "[ERROR] Unknown option: $1"; usage ;;
    esac
done

echo "[DEBUG] DOCKER_IMAGE=$DOCKER_IMAGE"

# 필수값 확인
if [[ -z "${TARGET-}" || -z "${METHOD-}" || -z "${ROOT_DIR-}" ]]; then
    echo "[ERROR] Missing required arguments."
    usage
fi

echo "[DEBUG] TARGET =$TARGET"
echo "[DEBUG] METHOD =$METHOD"

# 볼륨 경로 설정
HOST_DATA_DIR="${ROOT_DIR}/data"
HOST_BIN_DIR="${ROOT_DIR}/bin"  # generate_wcfamily_reference.sh

# Docker 실행
docker run --rm -it \
  --entrypoint /Work/NIPT/bin/generate_wcfamily_reference.sh \
  -v "${HOST_DATA_DIR}:/Work/NIPT/data" \
  -e WC="${WC:-/opt/wisecondor/wisecondor.py}" \
  -e WCX="${WCX:-wisecondorx}" \
  -e WCFF="${WCFF:-wisecondor-ff}" \
  -e PYTHON2="${PYTHON2:-python2.7}" \
  -e PYTHON3="${PYTHON3:-python3}" \
  "$DOCKER_IMAGE" \
  "$TARGET" "$METHOD"
