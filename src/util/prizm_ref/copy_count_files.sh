#!/bin/bash

# Check if arguments are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <list_file> <gender>"
    echo "Example: $0 male_list.txt M"
    echo "         $0 female_list.txt F"
    exit 1
fi

# Get arguments
LIST_FILE=$1
GENDER=$2

# Validate gender argument
if [[ "$GENDER" != "M" && "$GENDER" != "F" ]]; then
    echo "Error: Gender must be 'M' or 'F'"
    exit 1
fi

# Check if list file exists
if [ ! -f "$LIST_FILE" ]; then
    echo "Error: List file '$LIST_FILE' not found"
    exit 1
fi

# Base directories
ORIG_DIR="/Work/NIPT/refs/cordlife/count/orig"
FETUS_DIR="/Work/NIPT/refs/cordlife/count/fetus"
FETUS1_DIR="/Work/NIPT/refs/cordlife/count/fetus1"
MOM_DIR="/Work/NIPT/refs/cordlife/count/mom"

# Create directories if they don't exist
mkdir -p "${ORIG_DIR}/${GENDER}"
mkdir -p "${FETUS_DIR}/${GENDER}"
mkdir -p "${FETUS1_DIR}/${GENDER}"
mkdir -p "${MOM_DIR}/${GENDER}"

# Process list file
while IFS= read -r sample_path; do
    # Skip empty lines
    [ -z "$sample_path" ] && continue

    echo "Processing: $sample_path"

    # Copy files to respective directories
    cp "${sample_path}"/Output_PRIZM/[!z]*10mb.txt "${ORIG_DIR}/${GENDER}/"
    cp "${sample_path}"/Output_PRIZM/[!z]*filter.10mb.txt "${FETUS_DIR}/${GENDER}/"
    cp "${sample_path}"/Output_PRIZM/[!z]*filter1.10mb.txt "${FETUS1_DIR}/${GENDER}/"
    cp "${sample_path}"/Output_PRIZM/[!z]*filter_out.10mb.txt "${MOM_DIR}/${GENDER}/"

    #cp "${sample_path}"/Output_PRIZM/*10mb.txt "${ORIG_DIR}/${GENDER}/"
    #cp "${sample_path}"/Output_PRIZM/*filter.10mb.txt "${FETUS_DIR}/${GENDER}/"
    #cp "${sample_path}"/Output_PRIZM/*filter1.10mb.txt "${FETUS1_DIR}/${GENDER}/"
    #cp "${sample_path}"/Output_PRIZM/*filter_out.10mb.txt "${MOM_DIR}/${GENDER}/"
done < "$LIST_FILE"

echo "Done!"
