#!/bin/bash

# ===== Usage Check =====
if [ $# -ne 2 ]; then
    echo "Usage: $0 <target> <method>"
    echo "Method must be one of: WC, WCX, WCFF"
    echo "Example: $0 ucl/cordlife WCFF"
    exit 1
fi

TARGET=$1
METHOD=$2

# ===== Method Check =====
if [[ "$METHOD" != "WC" && "$METHOD" != "WCX" && "$METHOD" != "WCFF" ]]; then
    echo "[ERROR] Invalid method: $METHOD"
    exit 1
fi

# ===== Binaries =====
WCX_BIN="${WCX:-wisecondorx}"
WC_BIN="${WC:-/opt/wisecondor/wisecondor.py}"
WCFF_BIN="${WCFF:-wisecondor-ff}"

PYTHON2="${PYTHON2:-python2.7}"
PYTHON3="${PYTHON3:-python3}"

# ===== Constants =====
BASE_DIR="/Work/NIPT/data/refs/${TARGET}"
BINSIZE=200000

# ===== Paths =====
NPZ_DIR="${BASE_DIR}/gen_ref/${METHOD}"
REF_DIR="${BASE_DIR}/${METHOD}"
mkdir -p "${REF_DIR}"

# ===== Groups =====
#SAMPLE_TYPES=("orig" "fetus" "mom")
SAMPLE_TYPES=("fetus" "mom")

for type in "${SAMPLE_TYPES[@]}"; do
    echo "[INFO] Generating ${METHOD} reference for ${type}"

    # Ensure group dir exists
    if [ ! -d "${NPZ_DIR}/${type}/F" ]; then
        echo "[WARNING] Missing directory: ${NPZ_DIR}/${type}/F"
        continue
    fi

    if [ "$METHOD" == "WCX" ]; then
        if [ "$type" == "mom" ]; then
            ${WCX_BIN} newref \
                "${NPZ_DIR}/${type}/F"/*.npz \
                "${NPZ_DIR}/${type}/M"/*.npz \
                "${REF_DIR}/${type}_200k.npz" \
                --binsize ${BINSIZE} --yfrac 0.0
        else
            ${WCX_BIN} newref \
                "${NPZ_DIR}/${type}/F"/*.npz \
                "${REF_DIR}/${type}_F_200k.npz" \
                --binsize ${BINSIZE} --nipt --yfrac 0.0

            ${WCX_BIN} newref \
                "${NPZ_DIR}/${type}/M"/*.npz \
                "${REF_DIR}/${type}_M_200k.npz" \
                --binsize ${BINSIZE} --nipt --yfrac 0.0
        fi

    elif [ "$METHOD" == "WC" ]; then
        ${PYTHON2} ${WC_BIN} newref \
            "${NPZ_DIR}/${type}/F"/*.npz \
            "${NPZ_DIR}/${type}/M"/*.npz \
            "${REF_DIR}/${type}_200k.npz" \
            -binsize ${BINSIZE}

    elif [ "$METHOD" == "WCFF" ]; then
        ${WCFF_BIN} reference \
            -i "${NPZ_DIR}/${type}/F"/*.npz "${NPZ_DIR}/${type}/M"/*.npz \
            --binsize ${BINSIZE} \
            --refsize 100 \
            -o "${REF_DIR}/${type}_200k.npz"
    fi

    echo "[INFO] Completed ${METHOD} reference for ${type}"
done

echo "[DONE] Reference creation completed for target: ${TARGET}, method: ${METHOD}"
