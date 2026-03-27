#!/bin/bash

# 기본값
IMAGE_NAME="nipt_docker_v1.3"
TAG="latest"
DOCKERFILE_PATH="./docker/Dockerfile"

# 사용법 안내
usage() {
  echo "Usage: $0 [-c] [-n IMAGE_NAME] [-t TAG]"
  echo "  -c              docker build --no-cache (캐시 없이 전체 재빌드)"
  echo "  -n IMAGE_NAME   Docker image name (default: nipt_docker_v1.3)"
  echo "  -t TAG          Docker tag version (default: latest)"
  exit 1
}

# 옵션 파싱
NO_CACHE=""
while getopts ":cn:t:" opt; do
  case ${opt} in
    c )
      NO_CACHE=1
      ;;
    n )
      IMAGE_NAME=$OPTARG
      ;;
    t )
      TAG=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# image:tag: use "${NAME}:$TAG" not "${NAME}:${TAG}" (bash treats ${NAME: as substring)
IMAGE_REF="${IMAGE_NAME}:$TAG"

echo "🛠️  Building Docker image → ${IMAGE_REF}"
echo "📁  Dockerfile path → ${DOCKERFILE_PATH}"

# Docker 빌드 실행
if [ -n "$NO_CACHE" ]; then
  echo "📦  --no-cache (no layer cache)"
  docker build --no-cache -t "${IMAGE_REF}" -f "${DOCKERFILE_PATH}" .
else
  docker build -t "${IMAGE_REF}" -f "${DOCKERFILE_PATH}" .
fi

# 결과 확인
if [ $? -eq 0 ]; then
  echo "✅ Docker image built successfully: ${IMAGE_REF}"
else
  echo "❌ Docker build failed"
  exit 1
fi
