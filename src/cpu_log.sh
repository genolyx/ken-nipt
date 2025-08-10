#!/bin/bash

while true; do
    TIMESTAMP=$(date "+%Y-%m-%d %H:%M:%S")
    read -r MEM_USED MEM_TOTAL <<< $(free -m | awk '/^Mem:/ {print $3, $2}')
    MEM_PCT=$(awk "BEGIN {printf \"%.1f\", 100 * $MEM_USED / $MEM_TOTAL}")
    echo -e "${TIMESTAMP}\t${MEM_USED}\t${MEM_TOTAL}\t${MEM_PCT}"
    sleep 5
done
