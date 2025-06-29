#!/bin/bash

# Usage check: ensure two arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 combined_outliers.txt target_directory"
    exit 1
fi

OUTLIERS_FILE="$1"
TARGET_DIR="$2"

# Verify that the outliers file exists
if [ ! -f "$OUTLIERS_FILE" ]; then
    echo "Error: $OUTLIERS_FILE does not exist."
    exit 1
fi

# Verify that the target directory exists
if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: $TARGET_DIR does not exist or is not a directory."
    exit 1
fi

# Process each line in the outliers file
while IFS= read -r file; do
    # Skip empty lines
    if [ -n "$file" ]; then
        target_file="$TARGET_DIR/$file"
        if [ -f "$target_file" ]; then
            echo "Removing $target_file"
            rm "$target_file"
        else
            echo "File $target_file does not exist."
        fi
    fi
done < "$OUTLIERS_FILE"

echo "Deletion process complete."

