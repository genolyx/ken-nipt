#!/usr/bin/env bash
set -euo pipefail

# Docker 이미지 저장 및 배포 패키지 생성 스크립트
# Usage: ./create_deploy_package.sh [output_dir]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:-${SCRIPT_DIR}/deploy_package}"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
PACKAGE_NAME="Deploy_${TIMESTAMP}"

# 이미지 목록
IMAGES=(
    "roche_kapa_analysis:v2"
    "roche_web_ui:latest"
)

echo "=========================================="
echo "Docker 이미지 저장 및 배포 패키지 생성"
echo "=========================================="
echo "출력 디렉토리: ${OUTPUT_DIR}"
echo "패키지 이름: ${PACKAGE_NAME}"
echo ""

# 출력 디렉토리 생성
mkdir -p "${OUTPUT_DIR}/${PACKAGE_NAME}"
IMAGES_DIR="${OUTPUT_DIR}/${PACKAGE_NAME}/docker_images"
mkdir -p "${IMAGES_DIR}"

# Docker 이미지 저장
echo "1. Docker 이미지 저장 중..."
for image in "${IMAGES[@]}"; do
    image_name=$(echo "${image}" | cut -d':' -f1)
    image_tag=$(echo "${image}" | cut -d':' -f2)
    output_file="${IMAGES_DIR}/${image_name}_${image_tag}.tar"
    
    echo "  - ${image} → ${output_file}"
    
    # 이미지 존재 확인
    if ! docker image inspect "${image}" >/dev/null 2>&1; then
        echo "    ⚠️  경고: 이미지 ${image}를 찾을 수 없습니다. 건너뜁니다."
        continue
    fi
    
    # 이미지 저장
    docker save "${image}" -o "${output_file}"
    
    # 압축 (선택사항 - 용량 절약)
    echo "    압축 중..."
    gzip -f "${output_file}"
    
    echo "    ✅ 완료: ${output_file}.gz"
done

# 배포 스크립트 생성
echo ""
echo "2. 배포 스크립트 생성 중..."
DEPLOY_SCRIPT="${OUTPUT_DIR}/${PACKAGE_NAME}/deploy_images.sh"
cat > "${DEPLOY_SCRIPT}" << 'DEPLOY_EOF'
#!/usr/bin/env bash
set -euo pipefail

# Docker 이미지 로드 스크립트
# Usage: ./deploy_images.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGES_DIR="${SCRIPT_DIR}/docker_images"

echo "=========================================="
echo "Docker 이미지 로드"
echo "=========================================="
echo "이미지 디렉토리: ${IMAGES_DIR}"
echo ""

# Docker 설치 확인
if ! command -v docker &> /dev/null; then
    echo "❌ 오류: Docker가 설치되어 있지 않습니다."
    exit 1
fi

echo "Docker 버전: $(docker --version)"
echo ""

# 이미지 파일 찾기
image_files=($(find "${IMAGES_DIR}" -name "*.tar.gz" -o -name "*.tar" | sort))

if [ ${#image_files[@]} -eq 0 ]; then
    echo "❌ 오류: 이미지 파일을 찾을 수 없습니다."
    exit 1
fi

# 각 이미지 로드
for image_file in "${image_files[@]}"; do
    echo "로드 중: $(basename "${image_file}")"
    
    # .gz 파일인 경우 압축 해제
    if [[ "${image_file}" == *.gz ]]; then
        echo "  압축 해제 중..."
        gunzip -c "${image_file}" | docker load
    else
        docker load -i "${image_file}"
    fi
    
    if [ $? -eq 0 ]; then
        echo "  ✅ 완료"
    else
        echo "  ❌ 실패"
        exit 1
    fi
    echo ""
done

echo "=========================================="
echo "모든 이미지 로드 완료"
echo "=========================================="
echo ""
echo "로드된 이미지 목록:"
docker images | grep -E "roche_kapa_analysis|roche_web_ui" || true
DEPLOY_EOF

chmod +x "${DEPLOY_SCRIPT}"

# README 생성
echo "3. README 생성 중..."
README_FILE="${OUTPUT_DIR}/${PACKAGE_NAME}/README.md"
cat > "${README_FILE}" << README_EOF
# Docker 이미지 배포 패키지

생성일: $(date)

## 포함된 이미지

$(for image in "${IMAGES[@]}"; do echo "- ${image}"; done)

## 설치 방법

### 1. 패키지 압축 해제
\`\`\`bash
tar -xzf ${PACKAGE_NAME}.tar
cd ${PACKAGE_NAME}
\`\`\`

### 2. Docker 이미지 로드
\`\`\`bash
./deploy_images.sh
\`\`\`

### 3. 이미지 확인
\`\`\`bash
docker images | grep -E "roche_kapa_analysis|roche_web_ui"
\`\`\`

## 파일 구조

\`\`\`
${PACKAGE_NAME}/
├── docker_images/          # Docker 이미지 tar 파일들
│   ├── roche_kapa_analysis_v2.tar.gz
│   └── roche_web_ui_latest.tar.gz
├── deploy_images.sh        # 이미지 로드 스크립트
└── README.md              # 이 파일
\`\`\`

## 주의사항

- Docker가 설치되어 있어야 합니다.
- 이미지 로드에는 충분한 디스크 공간이 필요합니다.
- 이미지 크기: 약 3.2GB (압축 전)
README_EOF

# 패키지 크기 확인
echo ""
echo "4. 패키지 정보:"
TOTAL_SIZE=$(du -sh "${OUTPUT_DIR}/${PACKAGE_NAME}" | cut -f1)
echo "  총 크기: ${TOTAL_SIZE}"

# 최종 tar 파일 생성
echo ""
echo "5. 배포 패키지 압축 중..."
cd "${OUTPUT_DIR}"
tar -czf "${PACKAGE_NAME}.tar.gz" "${PACKAGE_NAME}"

FINAL_SIZE=$(du -sh "${PACKAGE_NAME}.tar.gz" | cut -f1)
echo "  ✅ 완료: ${PACKAGE_NAME}.tar.gz (${FINAL_SIZE})"

echo ""
echo "=========================================="
echo "배포 패키지 생성 완료!"
echo "=========================================="
echo "패키지 위치: ${OUTPUT_DIR}/${PACKAGE_NAME}.tar.gz"
echo ""
echo "타겟 시스템에서 사용 방법:"
echo "  1. tar -xzf ${PACKAGE_NAME}.tar.gz"
echo "  2. cd ${PACKAGE_NAME}"
echo "  3. ./deploy_images.sh"
echo ""


