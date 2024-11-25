#!/bin/bash

## Usage:
# remove_outliers_from_parquet.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -w 'ELKI_outlier_scores.txt.gz' \
#   -r 6 \
#   -s 2495667 \
#   -b "PRESERVED_SPECIMEN,MATERIAL_CITATION,MACHINE_OBSERVATION"



## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -w OUTLIER_SCORES -s SPECIES_KEY [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-e EXT_DIR] [-z COMPRESSION]"
    echo "  -i INPUT_FILE      : Input Parquet file path"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -w OUTLIER_SCORES  : Tab-delimited file with H3 cell IDs and outlier scores"
    echo "  -r H3_RESOLUTION   : H3 resolution (0-15)"
    echo "  -s SPECIESKEY      : Species key for filtering"
    echo "  -b BASIS_OF_RECORD : Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    echo "  -e EXT_DIR         : DuckDB extensions directory path (optional)"
    echo "  -z COMPRESSION     : ZSTD compression level (0-22) (optional, default: 10)"
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

## Validate input parameters
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$OUTLIER_SCORES" || -z "$SPECIES_KEY" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

if ! [[ "$SPECIES_KEY" =~ ^[0-9]+$ ]]; then
    echo -e "Error: Species key must be a positive integer!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

## `BASIS_OF_RECORD` should be a comma-separated list of valid values
VALID_BOR=(\
  "OBSERVATION" "OCCURRENCE" "MACHINE_OBSERVATION" "MATERIAL_SAMPLE" \
  "HUMAN_OBSERVATION" "MATERIAL_CITATION" "PRESERVED_SPECIMEN" \
  "FOSSIL_SPECIMEN" "LIVING_SPECIMEN")

if [[ -n "${BASIS_OF_RECORD}" ]]; then
    
    ## Split values into array
    ## NB! `read -a` works in bash, but not in zsh (in zsh, use `read -A`)
    IFS=',' read -r -a BOR <<< "${BASIS_OF_RECORD}"
    
    # Debug output
    # echo "Received values: ${BOR[@]}"
        
    ## Check each value
    for value in "${BOR[@]}"; do
        if [[ ! " ${VALID_BOR[@]} " =~ " ${value} " ]]; then
            echo "Error: Invalid basis of record value '${value}'"
            echo "Supported values are: ${VALID_BOR[*]}"
            exit 1
        fi
    done
fi

