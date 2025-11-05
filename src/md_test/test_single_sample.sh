#!/usr/bin/env bash

# Test script: Create single artificial sample
# FF: 5%, Deletion: 5Mb

set -euo pipefail
export LC_ALL=C

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Performance settings
THREADS="${THREADS:-8}"
COMPRESS_LEVEL="${COMPRESS_LEVEL:-9}"  # Higher compression for smaller files

# Test parameters
MOM_BAM="/home/ken/ken-nipt/analysis/2510/GNCI25100169/GNCI25100169.proper_paired.bam"
FETUS_BAM="/home/ken/ken-nipt/analysis/2508/GNCI25080181/GNCI25080181.proper_paired.bam"
MD_BED="test_1p36_5Mb.bed"
FF_TARGET=5
READS=15000000  # 15M reads (eff_reads with MAPQ>=30)
MOM_IDX=1  # Mom 샘플 인덱스
FETUS_IDX=1  # Fetus 샘플 인덱스
GENDER=""  # 자동 감지 (또는 "M" 또는 "F"로 명시)
BASE_OUTDIR="/home/ken/ken-nipt/analysis/md_test"

# Generate sample_id to determine output directory
# Extract disease name from BED file
DISEASE_NAME=$(awk 'NR==1 {print $4}' "$MD_BED" 2>/dev/null || echo "")
if [[ -z "$DISEASE_NAME" ]]; then
    DISEASE_NAME=$(basename "$MD_BED" .bed)
fi
# Clean disease name (remove spaces, special chars)
DISEASE_NAME=$(echo "$DISEASE_NAME" | sed 's/[^a-zA-Z0-9_]//g' | tr '[:upper:]' '[:lower:]')

# Calculate deletion size
DEL_START=$(awk 'NR==1 {print $2}' "$MD_BED")
DEL_END=$(awk 'NR==1 {print $3}' "$MD_BED")
DEL_SIZE_MB=$(( (DEL_END - DEL_START) / 1000000 ))

# Format: {mom_idx}_{fetus_idx}_{disease}_{ff}FF_{reads}M_{del}Mb_{gender}
READS_M=$((READS / 1000000))
if [[ -z "$GENDER" ]]; then
    # Try to detect gender from fetus BAM (simple check)
    Y_READS=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$FETUS_BAM" chrY 2>/dev/null || echo "0")
    A_READS=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$FETUS_BAM" chr1 2>/dev/null || echo "0")
    if [[ "$A_READS" -gt 0 ]]; then
        # Use awk for comparison (no bc needed)
        RATIO=$(awk -v y="$Y_READS" -v a="$A_READS" 'BEGIN{printf "%.6f", y/a}')
        # Compare: if ratio > 0.0001 then Male
        if awk -v r="$RATIO" 'BEGIN{exit (r > 0.0001 ? 0 : 1)}'; then
            GENDER="M"
        else
            GENDER="F"
        fi
    else
        GENDER="F"  # Default to Female if can't detect
    fi
fi
GENDER=$(echo "$GENDER" | tr '[:lower:]' '[:upper:]')
GENDER="${GENDER:0:1}"  # Take first character (M/F)

SAMPLE_ID="${MOM_IDX}_${FETUS_IDX}_${DISEASE_NAME}_FF${FF_TARGET}_${READS_M}M_${DEL_SIZE_MB}Mb_${GENDER}"
OUTDIR="$BASE_OUTDIR/$SAMPLE_ID"

# Create output directory
mkdir -p "$OUTDIR"

# Create FF map
FF_MAP="$OUTDIR/ff_map.tsv"
cat > "$FF_MAP" <<EOF
GNCI25100169	0.00
GNCI25080181	27.8
EOF

echo "Sample ID: $SAMPLE_ID"
echo "Output Directory: $OUTDIR"
echo "Mom: $(basename $(dirname "$MOM_BAM"))"
echo "Fetus: $(basename $(dirname "$FETUS_BAM"))"
echo "Mom Index: ${MOM_IDX}"
echo "Fetus Index: ${FETUS_IDX}"
echo "FF Target: ${FF_TARGET}%"
echo "Reads: ${READS} (eff_reads with MAPQ>=30)"
echo "Deletion: 5Mb"
echo "BED: $MD_BED"
echo "Gender: ${GENDER} (will be used for sample_id)"
echo "Threads: ${THREADS}"
echo "Compression: Level ${COMPRESS_LEVEL}"
echo "=========================================="
echo ""

# Check if files exist
for f in "$MOM_BAM" "$FETUS_BAM" "$MD_BED"; do
    if [[ ! -f "$f" ]]; then
        echo "Error: File not found: $f"
        exit 1
    fi
done

# Run make_artificial.sh
echo "Running make_artificial.sh..."
echo ""

# Build make_artificial.sh command
MAKE_ARGS=(
    --mom_bam "$MOM_BAM"
    --fetus_bam "$FETUS_BAM"
    --ff_map "$FF_MAP"
    --md_bed "$MD_BED"
    --ff_target "$FF_TARGET"
    --reads "$READS"
    --mom_idx "$MOM_IDX"
    --fetus_idx "$FETUS_IDX"
    --outdir "$OUTDIR"
)

# Add gender (already detected or specified)
MAKE_ARGS+=(--gender "$GENDER")
# Add sample_id to ensure consistent naming
MAKE_ARGS+=(--sample_id "$SAMPLE_ID")

# Add performance settings as environment variables
THREADS="$THREADS" COMPRESS_LEVEL="$COMPRESS_LEVEL" \
./make_artificial.sh "${MAKE_ARGS[@]}"

echo ""
echo "=========================================="
echo "Test completed!"
echo "=========================================="

# Get sample_id from JSON file
JSON_FILE=$(find "$OUTDIR" -name "*.json" -type f | head -1)
if [[ -n "$JSON_FILE" ]]; then
    SAMPLE_ID=$(python3 -c "import json; print(json.load(open('$JSON_FILE'))['sample_id'])" 2>/dev/null || echo "")
    FINAL_BAM="$OUTDIR/${SAMPLE_ID}.proper_paired.bam"
    echo "Sample ID: $SAMPLE_ID"
    echo "JSON metadata: $JSON_FILE"
else
    echo "Warning: JSON metadata file not found"
    FINAL_BAM="$OUTDIR/output.bam"
fi

echo "Output BAM: $FINAL_BAM"
echo ""

# Check output
if [[ -f "$FINAL_BAM" ]]; then
    echo "✓ Output BAM created successfully"

    # Check file size
    SIZE=$(du -h "$FINAL_BAM" | cut -f1)
    echo "  Size: $SIZE"

    # Count reads in output (eff_reads with MAPQ>=30)
    OUTPUT_READS=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$FINAL_BAM" 2>/dev/null || echo "0")
    OUTPUT_PAIRS=$((OUTPUT_READS / 2))
    echo "  Reads: $OUTPUT_READS (target: $READS)"
    echo "  Pairs: $OUTPUT_PAIRS"

    echo ""
    echo "To verify deletion region, check the deletion ratios in the output above."
    echo "Expected: Deletion region ratio should be lower than upstream/downstream ratios."
else
    echo "✗ Output BAM not found!"
    exit 1
fi

