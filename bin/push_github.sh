#!/bin/bash

# 간단한 설명 출력
echo "📥 Pulling latest changes from origin/main..."
git pull origin main || { echo "❌ git pull failed. Check for conflicts."; exit 1; }

# 수정 내용 확인
echo "📂 Checking modified files..."
git status

# 사용자 입력 받기
read -p "📝 Commit message: " msg

# 스테이징, 커밋, 푸시
git add .
git commit -m "$msg"
git push origin main

echo "✅ Done: changes pushed to GitHub!"
