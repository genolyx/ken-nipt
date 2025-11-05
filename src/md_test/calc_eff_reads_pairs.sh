#!/usr/bin/env bash

# Calculate eff_reads and eff_pairs from BAM files
# 
# Usage:
#   # From file list:
#   ./calc_eff_reads_pairs.sh bam_list.txt
#
#   # From stdin:
#   cat bam_list.txt | ./calc_eff_reads_pairs.sh
#
#   # Direct arguments:
#   ./calc_eff_reads_pairs.sh file1.bam file2.bam file3.bam
#
# Output: TSV format (sample, eff_reads, eff_pairs)
#   GNCI25100169	14874986	7437493
#   GNCI25100170	15090124	7545062
#
# Notes:
# - eff_reads: properly paired + primary + non-supplementary + MAPQ >= 30
# - eff_pairs: eff_reads / 2

set -euo pipefail
export LC_ALL=C

usage() {
    cat <<'USAGE'
Usage: calc_eff_reads_pairs.sh [OPTIONS] [BAM_FILE...]

Options:
  -i, --input FILE    Read BAM file paths from FILE (one per line)
  -o, --output FILE   Output TSV file (default: stdout)
  -h, --help         Show this help message

Examples:
  # From file list:
  ./calc_eff_reads_pairs.sh -i bam_list.txt -o results.tsv

  # From stdin:
  cat bam_list.txt | ./calc_eff_reads_pairs.sh

  # Direct arguments:
  ./calc_eff_reads_pairs.sh file1.bam file2.bam file3.bam
USAGE
    exit 1
}

# Parse arguments
INPUT_FILE=""
OUTPUT_FILE=""
BAM_FILES=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage
            ;;
        *)
            BAM_FILES+=("$1")
            shift
            ;;
    esac
done

# Collect BAM files from input
if [[ -n "$INPUT_FILE" ]]; then
    if [[ ! -f "$INPUT_FILE" ]]; then
        echo "Error: Input file not found: $INPUT_FILE" >&2
        exit 1
    fi
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        BAM_FILES+=("$line")
    done < "$INPUT_FILE"
elif [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    # Read from stdin
    while IFS= read -r line || [[ -n "$line" ]]; do
        [[ -z "$line" || "$line" =~ ^# ]] && continue
        BAM_FILES+=("$line")
    done
fi

if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
    echo "Error: No BAM files provided" >&2
    usage
fi

# Open output file or use stdout
if [[ -n "$OUTPUT_FILE" ]]; then
    exec 3> "$OUTPUT_FILE"
else
    exec 3>&1
fi

# Print header
echo -e "sample\teff_reads\teff_pairs" >&3

# Process each BAM file
TOTAL=${#BAM_FILES[@]}
CURRENT=0

for bam in "${BAM_FILES[@]}"; do
    ((CURRENT++)) || true
    
    # Get sample name from directory
    if [[ "$bam" == *"/"* ]]; then
        sample=$(basename "$(dirname "$bam")")
    else
        sample=$(basename "$bam" .bam)
        sample="${sample%.proper_paired}"
        sample="${sample%.sorted}"
    fi
    
    # Check if file exists
    if [[ ! -f "$bam" ]]; then
        echo "Warning: File not found: $bam (skipping)" >&2
        echo -e "${sample}\tFILE_NOT_FOUND\tN/A" >&3
        continue
    fi
    
    # Calculate eff_reads
    eff_reads=$(samtools view -c -f 2 -F 256 -F 2048 -q 30 "$bam" 2>/dev/null || echo "ERROR")
    
    if [[ "$eff_reads" == "ERROR" ]]; then
        echo "Warning: Error counting reads in $bam" >&2
        echo -e "${sample}\tERROR\tN/A" >&3
        continue
    fi
    
    # Calculate eff_pairs
    eff_pairs=$((eff_reads / 2))
    
    # Output result
    echo -e "${sample}\t${eff_reads}\t${eff_pairs}" >&3
    
    # Progress (stderr)
    printf "[%d/%d] %s: eff_reads=%s, eff_pairs=%s\n" \
        "$CURRENT" "$TOTAL" "$sample" "$eff_reads" "$eff_pairs" >&2
done

# Close output file if opened
if [[ -n "$OUTPUT_FILE" ]]; then
    exec 3>&-
    echo "Results saved to: $OUTPUT_FILE" >&2
fi

echo "Done. Processed ${TOTAL} BAM files." >&2

