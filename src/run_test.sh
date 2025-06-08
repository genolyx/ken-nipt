#!/bin/bash

INPUT_TSV="test/250605_01_list_part.tsv"

# Skip header and process each line
tail -n +2 "$INPUT_TSV" | while IFS=$'\t' read -r SAMPLE_NAME WORK_DIR FQ1 FQ2 AGE LAB; do
  # Trim all fields (remove leading/trailing whitespace and carriage return)
  SAMPLE_NAME=$(echo "$SAMPLE_NAME" | tr -d '\r' | xargs)
  WORK_DIR=$(echo "$WORK_DIR" | tr -d '\r' | xargs)
  FQ1=$(echo "$FQ1" | tr -d '\r' | xargs)
  FQ2=$(echo "$FQ2" | tr -d '\r' | xargs)
  AGE=$(echo "$AGE" | tr -d '\r' | xargs)
  LAB=$(echo "$LAB" | tr -d '\r' | xargs)

  echo "Running for sample: $SAMPLE_NAME"
  bash run_nipt_dev.sh -s "$SAMPLE_NAME" -1 "$FQ1" -2 "$FQ2" -l "$LAB" -a "$AGE" -root /home/ken/ken-nipt -work "$WORK_DIR" --detached
done
