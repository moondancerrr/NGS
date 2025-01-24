#!/usr/bin/env bash

DATA_DIR="./Data"

if [ ! -d "$DATA_DIR" ]; then
    echo "=============================================================="
    echo " Copying datasets from https://zenodo.org/record/6952995      "
    echo " -> Destination: $DATA_DIR                                   "
    echo "=============================================================="
    
    mkdir -p "$DATA_DIR"

    # Download Clover dataset (with -L)
    curl -L "https://zenodo.org/record/6952995/files/clover.tar.gz?download=1" \
        -o "$DATA_DIR/Clover_Data.tar.gz"
    tar -zxvf "$DATA_DIR/Clover_Data.tar.gz" -C "$DATA_DIR"

    # Download single-cell dataset (with -L)
    curl -L "https://zenodo.org/record/6952995/files/singlecell.tar.gz?download=1" \
        -o "$DATA_DIR/scrna_Data.tar.gz"
    tar -zxvf "$DATA_DIR/scrna_Data.tar.gz" -C "$DATA_DIR"

    rm -f "$DATA_DIR"/*.tar.gz
else
    echo "=========================================================="
    echo " Datasets folder already exists, no need to download it.  "
    echo "=========================================================="
fi
