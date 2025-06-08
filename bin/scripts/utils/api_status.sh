#!/bin/bash
# Usage: api_status.sh <sample_id> <stage> <status> <message> <api_url> <api_enabled>

sample_id=$1
stage=$2
status=$3
message=$4
api_url=$5
api_enabled=$6

# API 상태 업데이트
if [[ "${api_enabled}" == "true" ]]; then
    echo "[API] Updating status: ${sample_id} - ${stage} - ${status}" >&2
    
    curl -s -X POST "${api_url}/sample/status" \
        -H "Content-Type: application/json" \
        -d "{\"sample_id\": \"${sample_id}\", \"stage\": \"${stage}\", \"status\": \"${status}\", \"message\": \"${message}\"}" >/dev/null
    
    if [ $? -eq 0 ]; then
        echo "[API] Status update successful" >&2
    else
        echo "[API] Failed to update status" >&2
    fi
else
    echo "[API] Status updates disabled (would update ${sample_id} - ${stage} - ${status})" >&2
fi
