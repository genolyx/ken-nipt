#!/bin/bash
# Full batch run script with nohup

cd /home/ken/ken-nipt/src/batch_run

# Run batch processing with nohup
nohup python3 run_batch_manual.py \
    -l cordlife \
    -s sample_sheet.tsv \
    -m 2 \
    -i 5 \
    > run_batch_manual.log 2>&1 &

echo "Batch processing started in background"
echo "PID: $!"
echo "Log file: run_batch_manual.log"
echo "Monitor with: tail -f run_batch_manual.log"
