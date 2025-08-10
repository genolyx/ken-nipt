#!/bin/bash

INTERVAL=${1:-5}
OUTDIR=${2:-"./logs"}
mkdir -p "$OUTDIR"
LOGFILE="$OUTDIR/resource_log_$(date +%Y%m%d_%H%M%S).tsv"

# TSV 헤더
echo -e "timestamp\tcpu_usage_percent\tmem_used_mb\tmem_total_mb\tmemory_usage_percent" > "$LOGFILE"

# CPU 사용률 계산 함수
get_cpu_usage() {
    # CPU stat 읽기
    read -r cpu prev_user prev_nice prev_system prev_idle prev_iowait prev_irq prev_softirq prev_steal _ < /proc/stat
    prev_idle_total=$((prev_idle + prev_iowait))
    prev_non_idle=$((prev_user + prev_nice + prev_system + prev_irq + prev_softirq + prev_steal))
    prev_total=$((prev_idle_total + prev_non_idle))

    sleep 1

    read -r cpu cur_user cur_nice cur_system cur_idle cur_iowait cur_irq cur_softirq cur_steal _ < /proc/stat
    cur_idle_total=$((cur_idle + cur_iowait))
    cur_non_idle=$((cur_user + cur_nice + cur_system + cur_irq + cur_softirq + cur_steal))
    cur_total=$((cur_idle_total + cur_non_idle))

    total_diff=$((cur_total - prev_total))
    idle_diff=$((cur_idle_total - prev_idle_total))

    if [ "$total_diff" -eq 0 ]; then
        echo 0
    else
        echo $(( (100 * (total_diff - idle_diff)) / total_diff ))
    fi
}

while true; do
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
    CPU=$(get_cpu_usage)

    read -r MEM_USED MEM_TOTAL <<< $(free -m | awk '/^Mem:/ {print $3, $2}')
    MEM_PCT=$(awk "BEGIN {printf \"%.1f\", 100 * $MEM_USED / $MEM_TOTAL}")

    echo -e "${TIMESTAMP}\t${CPU}\t${MEM_USED}\t${MEM_TOTAL}\t${MEM_PCT}" >> "$LOGFILE"
    sleep "$INTERVAL"
done
