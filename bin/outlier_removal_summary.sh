#!/bin/bash

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_DIR [-q THRESHOLD] -o OUTPUT_FILE [-t THREADS] [-m MEMORY] [-x TEMP_DIR]"
    echo "  -i INPUT_DIR      : Input Parquet file path"
    echo "  -q THRESHOLD       : Outlier score threshold (default: 0.5) (optional)"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    exit 1
}

## Initialize variables
INPUT_DIR=""
THRESHOLD="0.5"
OUTPUT_FILE=""
THREADS=""
MEMORY=""
TEMP_DIR=""

## Parse command-line options
while getopts "i:q:o:t:m:x:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        q) THRESHOLD="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done


## Start the SQL command
echo -e "\nPreparing SQL command"

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
CREATE TEMPORARY TABLE tbl AS SELECT * FROM read_csv('/dev/stdin',
      auto_detect = false,
      header = false,
      delim = '\t',
      columns = {
        'SpeciesKey':   'BIGINT',
        'H3Index':      'VARCHAR',
        'OutlierScore': 'DOUBLE'
      });
     COPY tbl TO '${OUTPUT_FILE}' (HEADER, DELIMITER '\t', COMPRESSION 'gzip');
"

## Function to add species key as the first column
add_specieskey() {
  zcat "$1" | awk -v id="$2" 'BEGIN { OFS="\t" } { print id, $0 }'
}
export -f add_specieskey

## Process all files in the input directory
find "${INPUT_DIR}" -name "*.txt.gz" \
  | parallel \
    -j 1 \
    --rpl '{s} s:.*/::; s/_.*//' \
    "add_specieskey {} {s}" \
  | duckdb :memory: "${SQL_COMMAND}"

echo -e "..Done\n"


