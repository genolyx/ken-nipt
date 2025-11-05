#!/bin/bash

# 기본값
IMAGE_NAME="nipt_docker_v1.2"
TAG="latest"
DOCKERFILE_PATH="./docker/Dockerfile"

# 사용법 안내
usage() {
  echo "Usage: $0 [-n IMAGE_NAME] [-t TAG]"
  echo "  -n IMAGE_NAME   Docker image name (default: nipt-docker)"
  echo "  -t TAG          Docker tag version (default: latest)"
  exit 1
}

# 옵션 파싱
while getopts ":n:t:" opt; do
  case ${opt} in
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

echo "🛠️  Building Docker image → ${IMAGE_NAME}:${TAG}"
echo "📁  Dockerfile path → ${DOCKERFILE_PATH}"

# Docker 빌드 실행
#docker build --no-cache -t ${IMAGE_NAME}:${TAG} -f ${DOCKERFILE_PATH} .
docker build -t ${IMAGE_NAME}:${TAG} -f ${DOCKERFILE_PATH} .

# 결과 확인
if [ $? -eq 0 ]; then
  echo "✅ Docker image built successfully: ${IMAGE_NAME}:${TAG}"
else
  echo "❌ Docker build failed"
  exit 1
fi
