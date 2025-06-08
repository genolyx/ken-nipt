#!/bin/bash

INTERVAL=${1:-5}   # 기본: 5초 간격
OUTDIR=${2:-"./logs"}
mkdir -p "$OUTDIR"
LOGFILE="$OUTDIR/resource_log_$(date +%Y%m%d_%H%M%S).tsv"

# TSV 헤더
echo -e "timestamp\tcpu_usage_percent\tmem_used_mb\tmem_total_mb\tmemory_usage_percent" > "$LOGFILE"

# CPU 사용률 계산 함수 (1초 샘플링)
get_cpu_usage() {
    PREV=($(head -n1 /proc/stat)); sleep 1
    CUR=($(head -n1 /proc/stat))

    PREV_IDLE=${PREV[4]}
    CUR_IDLE=${CUR[4]}
    PREV_TOTAL=0; CUR_TOTAL=0
    for val in "${PREV[@]:1}"; do ((PREV_TOTAL += val)); done
    for val in "${CUR[@]:1}"; do ((CUR_TOTAL += val)); done

    DIFF_TOTAL=$((CUR_TOTAL - PREV_TOTAL))
    DIFF_IDLE=$((CUR_IDLE - PREV_IDLE))
    CPU_USAGE=$((100 * (DIFF_TOTAL - DIFF_IDLE) / DIFF_TOTAL))
    echo "$CPU_USAGE"
}

while true; do
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
    CPU=$(get_cpu_usage)

    read MEM_USED MEM_TOTAL <<< $(free -m | awk '/^Mem:/ {print $3, $2}')
    MEM_PCT=$(awk "BEGIN {printf \"%.1f\", 100 * $MEM_USED / $MEM_TOTAL}")

    echo -e "${TIMESTAMP}\t${CPU}\t${MEM_USED}\t${MEM_TOTAL}\t${MEM_PCT}" >> "$LOGFILE"
    sleep "$INTERVAL"
done
