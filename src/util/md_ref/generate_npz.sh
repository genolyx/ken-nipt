#!/bin/bash

# Check if argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 Target"
    echo "Example: $0 ucl/cordlife"
    exit 1
fi

# Static paths
WCX_BIN="/usr/bin/miniconda3/bin/WisecondorX"
WC_PATH="/Work/NIPT/bin/wisecondor/wisecondor.py"
PYTHON2="/usr/bin/miniconda2/bin/python2"
PYTHON3="/usr/bin/miniconda3/bin/python3"

# Set target from argument
TARGET=$1

# Define base directories
BASE_DIR="/Work/NIPT/refs/${TARGET}"
DATA_DIR="${BASE_DIR}/data"
NPZ_WCX_DIR="${BASE_DIR}/npz/WCX"
NPZ_WC_DIR="${BASE_DIR}/npz/WC"

# Array of subdirectories
SAMPLE_TYPES=("orig" "fetus" "mom")
GENDERS=("F" "M")

# Create output directories if they don't exist
for type in "${SAMPLE_TYPES[@]}"; do
    for gender in "${GENDERS[@]}"; do
        mkdir -p "${NPZ_WCX_DIR}/${type}/${gender}"
        mkdir -p "${NPZ_WC_DIR}/${type}/${gender}"
    done
done

# Process each subdirectory
for type in "${SAMPLE_TYPES[@]}"; do
    for gender in "${GENDERS[@]}"; do
        # Define current working directory
        CURRENT_DIR="${DATA_DIR}/${type}/${gender}"
        
        # Check if directory exists
        if [ ! -d "$CURRENT_DIR" ]; then
            echo "Directory not found: $CURRENT_DIR"
            continue
        fi

        # Process each BAM file
        for bam_file in "${CURRENT_DIR}"/*.bam; do
            # Check if any bam files exist
            if [ ! -e "$bam_file" ]; then
                echo "No BAM files found in: $CURRENT_DIR"
                break
            fi

            # Get filename without path and extension
            filename=$(basename "$bam_file" .bam)
            
            echo "Processing: $bam_file"

            # WisecondorX conversion
            echo "Converting with WisecondorX..."
            ${WCX_BIN} convert --binsize 100000 "$bam_file" "${NPZ_WCX_DIR}/${type}/${gender}/${filename}.npz"
            
            # Original Wisecondor conversion
            echo "Converting with Wisecondor..."
            ${PYTHON2} ${WC_PATH} convert "$bam_file" "${NPZ_WC_DIR}/${type}/${gender}/${filename}.npz" -binsize 100000

            echo "Completed processing: $filename"
        done
    done
done

echo "Processing completed for target: ${TARGET}"
