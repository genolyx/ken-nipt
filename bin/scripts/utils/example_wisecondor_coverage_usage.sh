#!/bin/bash
#
# Example usage of wisecondor_coverage_check.py for NIPT
#
# This script demonstrates how to use the Wisecondor coverage check tool
# for NIPT shallow depth WGS (0.2-0.3x average depth)

# Configuration
SAMPLE_ID="OPC241100001"
WORK_DIR="2411"
BAM_FILE="/Work/NIPT/analysis/${WORK_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
BED_FILE="/Work/NIPT/data/bed/common/MD_Target_8.bed"
OUTPUT_DIR="/Work/NIPT/analysis/${WORK_DIR}/${SAMPLE_ID}/coverage_analysis"
TOOL_PATH="/Work/NIPT/bin/scripts/utils/wisecondor_coverage_check.py"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "================================================================"
echo "Wisecondor Coverage Analysis for NIPT"
echo "================================================================"
echo "Sample ID: ${SAMPLE_ID}"
echo "BAM file: ${BAM_FILE}"
echo "BED file: ${BED_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "================================================================"
echo ""

# Check sequencing depth first
echo "Step 0: Checking sequencing depth..."
TOTAL_READS=$(samtools idxstats "${BAM_FILE}" | awk '{sum+=$3} END {print sum}')
DEPTH=$(echo "scale=3; ${TOTAL_READS} * 150 / 3000000000" | bc)
echo "Total mapped reads: ${TOTAL_READS}"
echo "Estimated average depth: ${DEPTH}x"
echo ""

if (( $(echo "${DEPTH} < 0.1" | bc -l) )); then
    echo "⚠ WARNING: Very low depth (<0.1x), analysis may be unreliable"
elif (( $(echo "${DEPTH} < 0.15" | bc -l) )); then
    echo "⚠ WARNING: Low depth (<0.15x), consider using relaxed parameters"
elif (( $(echo "${DEPTH} > 1.0" | bc -l) )); then
    echo "⚠ NOTE: High depth (>1.0x), this tool is designed for shallow depth"
else
    echo "✓ Depth is appropriate for NIPT analysis"
fi
echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 1: Basic usage (default parameters for NIPT)
echo "Example 1: Basic coverage check with default NIPT parameters"
echo "  - Bin size: 200kb (Wisecondor default)"
echo "  - Min reads per bin: 10 (sufficient for 0.2x depth)"
echo "  - Min bin fraction: 80% (good coverage required)"
echo "  - Min total reads: 1000"
echo ""

python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_default.bed" \
    --report

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 2: For low depth samples (0.1-0.15x)
if (( $(echo "${DEPTH} < 0.15" | bc -l) )); then
    echo "Example 2: Relaxed parameters for low depth samples"
    echo "  - Min reads per bin: 5 (relaxed)"
    echo "  - Min bin fraction: 70% (relaxed)"
    echo "  - Min total reads: 500 (relaxed)"
    echo ""

    python3 "${TOOL_PATH}" \
        --bam "${BAM_FILE}" \
        --bed "${BED_FILE}" \
        --output "${OUTPUT_DIR}/coverable_regions_relaxed.bed" \
        --min-reads-per-bin 5 \
        --min-bin-fraction 0.7 \
        --min-total-reads 500 \
        --report

    echo ""
    echo "----------------------------------------------------------------"
    echo ""
fi

# Example 3: Strict parameters for high-quality analysis
if (( $(echo "${DEPTH} >= 0.25" | bc -l) )); then
    echo "Example 3: Strict parameters for high-quality samples (≥0.25x)"
    echo "  - Min reads per bin: 15 (strict)"
    echo "  - Min bin fraction: 90% (strict)"
    echo "  - Min total reads: 1500 (strict)"
    echo ""

    python3 "${TOOL_PATH}" \
        --bam "${BAM_FILE}" \
        --bed "${BED_FILE}" \
        --output "${OUTPUT_DIR}/coverable_regions_strict.bed" \
        --min-reads-per-bin 15 \
        --min-bin-fraction 0.9 \
        --min-total-reads 1500 \
        --report

    echo ""
    echo "----------------------------------------------------------------"
    echo ""
fi

# Example 4: Check multiple MD targets
echo "Example 4: Check all MD targets"
echo ""

MD_TARGETS=(
    "MD_Target_1"
    "MD_Target_4"
    "MD_Target_5"
    "MD_Target_8"
    "MD_Target_15"
    "MD_Target_22"
)

MD_TARGET_DIR="${OUTPUT_DIR}/md_targets"
mkdir -p "${MD_TARGET_DIR}"

for TARGET in "${MD_TARGETS[@]}"; do
    TARGET_BED="/Work/NIPT/data/bed/common/${TARGET}.bed"
    
    if [ -f "${TARGET_BED}" ]; then
        echo "  Checking ${TARGET}..."
        python3 "${TOOL_PATH}" \
            --bam "${BAM_FILE}" \
            --bed "${TARGET_BED}" \
            --output "${MD_TARGET_DIR}/coverable_${TARGET}.bed" \
            --report 2>&1 | grep -E "(analyzable|Progress|Complete)" || true
    else
        echo "  Skipping ${TARGET} (BED file not found)"
    fi
done

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 5: Custom bin size (for different analysis)
echo "Example 5: Custom bin size (500kb for coarse analysis)"
echo ""

python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_500k.bed" \
    --bin-size 500000 \
    --report

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 6: Only output fully covered regions (no partial)
echo "Example 6: Only fully covered regions (no partial)"
echo ""

python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_full_only.bed" \
    --no-partial \
    --report

echo ""
echo "================================================================"
echo "Analysis complete!"
echo "================================================================"
echo ""
echo "Results saved to: ${OUTPUT_DIR}"
echo ""
echo "Output files:"
ls -lh "${OUTPUT_DIR}"/*.bed 2>/dev/null || echo "  No BED files generated"
echo ""
ls -lh "${OUTPUT_DIR}"/*_report.txt 2>/dev/null || echo "  No report files generated"
echo ""
echo "MD target results:"
ls -lh "${MD_TARGET_DIR}"/*.bed 2>/dev/null || echo "  No MD target files generated"
echo ""
echo "================================================================"
echo "Next steps:"
echo "================================================================"
echo ""
echo "1. Review the coverage reports:"
echo "   less ${OUTPUT_DIR}/coverable_regions_default_report.txt"
echo ""
echo "2. Check which regions are analyzable:"
echo "   grep 'PASS\\|FAIL' ${OUTPUT_DIR}/coverable_regions_default_report.txt"
echo ""
echo "3. Use coverable regions BED file for Wisecondor/WisecondorX analysis"
echo ""
echo "4. For low coverage samples, try relaxed parameters:"
echo "   --min-reads-per-bin 5 --min-bin-fraction 0.7"
echo ""
echo "5. Compare different parameter sets:"
echo "   diff ${OUTPUT_DIR}/coverable_regions_default.bed \\"
echo "        ${OUTPUT_DIR}/coverable_regions_relaxed.bed"
echo ""
echo "================================================================"
echo ""

