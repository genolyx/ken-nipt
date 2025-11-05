#!/usr/bin/env bash

# Merge mom_list.tsv and non_preg_mom_list.tsv into a single file
# Removes duplicates based on sample ID

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MOM_LIST1="mom_list.tsv"
MOM_LIST2="non_preg_mom_list.tsv"
OUTPUT="mom_list_combined.tsv"

if [[ ! -f "$MOM_LIST1" ]]; then
    echo "Error: $MOM_LIST1 not found"
    exit 1
fi

if [[ ! -f "$MOM_LIST2" ]]; then
    echo "Error: $MOM_LIST2 not found"
    exit 1
fi

echo "Merging mom lists..."
echo "  Input 1: $MOM_LIST1"
echo "  Input 2: $MOM_LIST2"
echo "  Output: $OUTPUT"
echo ""

# Get header from first file
head -1 "$MOM_LIST1" > "$OUTPUT"

# Combine and remove duplicates based on sample ID (column 2)
{
    tail -n +2 "$MOM_LIST1"
    tail -n +2 "$MOM_LIST2"
} | awk -F'\t' '
BEGIN { OFS="\t" }
{
    sample_id = $2
    if (!seen[sample_id]) {
        print
        seen[sample_id] = 1
    }
}' >> "$OUTPUT"

COUNT1=$(tail -n +2 "$MOM_LIST1" | wc -l)
COUNT2=$(tail -n +2 "$MOM_LIST2" | wc -l)
OUTPUT_COUNT=$(tail -n +2 "$OUTPUT" | wc -l)

echo "Summary:"
echo "  $MOM_LIST1: $COUNT1 samples"
echo "  $MOM_LIST2: $COUNT2 samples"
echo "  $OUTPUT: $OUTPUT_COUNT samples (after removing duplicates)"
echo ""
echo "✓ Merged mom list created: $OUTPUT"

