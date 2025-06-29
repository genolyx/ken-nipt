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
NPZ_WCX_DIR="${BASE_DIR}/npz/WCX"
NPZ_WC_DIR="${BASE_DIR}/npz/WC"
WCX_REF_OUTPUT="${BASE_DIR}/MD/WCX_ref"
WC_REF_OUTPUT="${BASE_DIR}/MD/WC_ref"

# Set binsize
BINSIZE=200000

# Array of sample types
SAMPLE_TYPES=("orig" "fetus" "mom")
#SAMPLE_TYPES=("orig" "fetus")
#SAMPLE_TYPES=("orig")

# Process each sample type
for type in "${SAMPLE_TYPES[@]}"; do
    echo "Creating reference for ${type}"
    
    # Create output directories for reference files if they don't exist
    mkdir -p "${NPZ_WCX_DIR}/${type}"
    mkdir -p "${NPZ_WC_DIR}/${type}"
    
    # WisecondorX reference (combining M and F)
    echo "Creating WisecondorX reference for ${type}..."
    if [ "$type" == "mom" ]; then
        # For mom samples, specify female gender
        ${WCX_BIN} newref \
            "${NPZ_WCX_DIR}/${type}/F"/*.npz \
            "${NPZ_WCX_DIR}/${type}/M"/*.npz \
            "${WCX_REF_OUTPUT}/${type}_200k.npz" \
            --binsize ${BINSIZE} --yfrac 0 
    else
        # For other samples, use default settings
        ${WCX_BIN} newref \
            "${NPZ_WCX_DIR}/${type}/F"/*.npz \
            "${NPZ_WCX_DIR}/${type}/M"/*.npz \
            "${WCX_REF_OUTPUT}/${type}_200k.npz" \
            --binsize ${BINSIZE} --nipt
    fi
    
    # Original Wisecondor reference (combining M and F)
    echo "Creating Wisecondor reference for ${type}..."
    ${PYTHON2} ${WC_PATH} newref \
        "${NPZ_WC_DIR}/${type}/F"/*.npz \
        "${NPZ_WC_DIR}/${type}/M"/*.npz \
        "${WC_REF_OUTPUT}/${type}_200k.npz" \
        -binsize ${BINSIZE}
    
    echo "Completed reference creation for ${type}"

done

echo "Reference creation completed for target: ${TARGET}"
