#!/bin/bash
#
# Example usage of md_coverage_check.py
#
# This script demonstrates how to use the MD coverage check tool
# to analyze BAM coverage for Wisecondor/WisecondorX analysis

# Configuration
SAMPLE_ID="your_sample_id"
WORK_DIR="2511"
BAM_FILE="/Work/NIPT/analysis/${WORK_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.proper_paired.bam"
BED_FILE="/Work/NIPT/data/bed/common/MD_Target_8.bed"
OUTPUT_DIR="/Work/NIPT/analysis/${WORK_DIR}/${SAMPLE_ID}/coverage_analysis"
TOOL_PATH="/Work/NIPT/bin/scripts/utils/md_coverage_check.py"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "================================================================"
echo "MD Coverage Analysis"
echo "================================================================"
echo "Sample ID: ${SAMPLE_ID}"
echo "BAM file: ${BAM_FILE}"
echo "BED file: ${BED_FILE}"
echo "Output directory: ${OUTPUT_DIR}"
echo "================================================================"
echo ""

# Example 1: Basic usage (default parameters)
echo "Example 1: Basic coverage check with default parameters"
python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_default.bed"

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 2: Strict parameters for high-quality analysis
echo "Example 2: Strict parameters (min 5x coverage, 95% threshold)"
python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_strict.bed" \
    --min-coverage 5 \
    --coverage-threshold 0.95 \
    --min-region-size 5000 \
    --report

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 3: For WisecondorX analysis with custom bin size
echo "Example 3: Custom bin size (500kb for different analysis)"
python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_500k.bed" \
    --bin-size 500000 \
    --report

echo ""
echo "----------------------------------------------------------------"
echo ""

# Example 4: Only output fully covered regions (no partial regions)
echo "Example 4: Only fully covered regions (no partial)"
python3 "${TOOL_PATH}" \
    --bam "${BAM_FILE}" \
    --bed "${BED_FILE}" \
    --output "${OUTPUT_DIR}/coverable_regions_full_only.bed" \
    --no-partial \
    --report

echo ""
echo "================================================================"
echo "Analysis complete!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "================================================================"
echo ""
echo "Next steps:"
echo "1. Review the coverage report: ${OUTPUT_DIR}/*_report.txt"
echo "2. Use coverable regions BED file for MD analysis"
echo "3. Run Wisecondor/WisecondorX with the filtered regions"
echo ""

