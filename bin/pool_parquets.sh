#!/bin/bash

## Usage:
# pool_parquets.sh \
#   -i '/path/to/input/directory' \
#   -o '/path/to/output.parquet'

## TO DO:
# - Add hive partitioning? (by year?)
# - Extract all unique grid cells and get coordinates?

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-e EXT_DIR]"
    echo "  -i INPUT      : Input directory with Parquet files"
    echo "  -o OUTPUT     : Output Parquet file path"
    echo "  -t THREADS    : Number of CPU threads to use (optional)"
    echo "  -m MEMORY     : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR   : Temporary directory path (optional)"
    exit 1
}

## Initialize variables
INPUT=""
OUTPUT=""
THREADS=""
MEMORY=""
TEMP_DIR=""

## Parse command-line options
while getopts "i:o:t:m:x:" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done


## Validate input parameters
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input directory: $INPUT"
echo "Output file: $OUTPUT"
if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi
if [[ -n "$MEMORY" ]]; then
    echo "Memory: $MEMORY"
fi
if [[ -n "$TEMP_DIR" ]]; then
    echo "Temp directory: $TEMP_DIR"
fi


SQL_COMMAND=""

## Add configuration settings (if provided)
if [[ -n "$THREADS" ]]; then
    SQL_COMMAND+="
SET threads TO ${THREADS};
"
fi

if [[ -n "$MEMORY" ]]; then
    SQL_COMMAND+="
SET memory_limit = '${MEMORY}';
"
fi

if [[ -n "$TEMP_DIR" ]]; then
    SQL_COMMAND+="
PRAGMA temp_directory='${TEMP_DIR}';
"
fi


SQL_COMMAND+="
COPY (
    SELECT * FROM read_parquet('${INPUT}/*.parquet')
) 
TO '${OUTPUT}' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 14, ROW_GROUP_SIZE 200000, ROW_GROUPS_PER_FILE 1);
"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"


