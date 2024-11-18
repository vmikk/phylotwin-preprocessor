#!/bin/bash


## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -r H3_RESOLUTION -s SPECIES_KEY [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR]"
    echo "  -i INPUT_FILE      : Input Parquet file path"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -w OUTLIER_SCORES  : Tab-delimited file with H3 cell IDs and outlier scores"
    echo "  -s SPECIESKEY      : Species key for filtering"
    echo "  -b BASIS_OF_RECORD : Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    exit 1
}

## Initialize variables
INPUT_FILE=""
OUTPUT_FILE=""
OUTLIER_SCORES=""
SPECIES_KEY=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""
