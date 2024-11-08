#!/bin/bash

## Script to prepare data for spatial outlier detection with ELKI

## Usage:
# ./h3_preprocessor.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -r 6 \
#   -s 2495667

## Main workflow:
# - filter records using the `specieskey` column
# - convert latitude/longitude coordinates to H3 cells
# - keep unique H3 grids
# - estimate the centroid of H3 grids

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -r H3_RESOLUTION -s SPECIES_KEY [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-d] [-c]"
    echo "  -i INPUT_FILE    : Input Parquet file path"
    echo "  -o OUTPUT_FILE   : Output Parquet file path"
    echo "  -r H3_RESOLUTION : H3 resolution (0-15)"
    echo "  -s SPECIES_KEY   : Species key for filtering"
    echo "  -t THREADS       : Number of CPU threads to use (optional)"
    echo "  -m MEMORY        : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR      : Temporary directory path (optional)"
    echo "  -d               : Save SQL script for debugging (optional, default: true)"
    echo "  -c               : Convert Parquet output to CSV (optional, default: true)"
    exit 1
}

## Initialize variables
INPUT_FILE=""
OUTPUT_FILE=""
H3_RESOLUTION=""
SPECIES_KEY=""
THREADS=""
MEMORY=""
TEMP_DIR=""
SAVE_SQL_SCRIPT=true
CONVERT_TO_CSV=true

## Parse command-line options
while getopts "i:o:r:s:t:m:x:dc" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        r) H3_RESOLUTION="$OPTARG" ;;
        s) SPECIES_KEY="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        d) SAVE_SQL_SCRIPT=true ;;
        c) CONVERT_TO_CSV=true ;;
        *) usage ;;
    esac
done

## Validate input parameters
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$H3_RESOLUTION" || -z "$SPECIES_KEY" ]]; then
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

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "H3 resolution: $H3_RESOLUTION"
echo "Species key: $SPECIES_KEY"
echo "Save SQL script: $SAVE_SQL_SCRIPT"
echo "Convert output to CSV: $CONVERT_TO_CSV"

if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi
if [[ -n "$MEMORY" ]]; then
    echo "Memory: $MEMORY"
fi
if [[ -n "$TEMP_DIR" ]]; then
    echo "Temp directory: $TEMP_DIR"
fi

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

-- Install and load the H3 extension
INSTALL h3 FROM community;
LOAD h3;

COPY (
    WITH
        inp AS (
            SELECT 
                h3_latlng_to_cell(decimallatitude, decimallongitude, ${H3_RESOLUTION})::UBIGINT AS h3_index,
                specieskey
            FROM 
                read_parquet('${INPUT_FILE}')
            WHERE specieskey = ${SPECIES_KEY}
        ),
        unique_grids AS (
            SELECT DISTINCT h3_index
            FROM inp
            WHERE h3_index IS NOT NULL
        ),
        grid_centroids AS (
            SELECT 
                h3_cell_to_lat(h3_index) AS LAT,
                h3_cell_to_lng(h3_index) AS LON,
                h3_h3_to_string(h3_index) as H3,
            FROM unique_grids
        )

    SELECT * FROM grid_centroids
) TO '${OUTPUT_FILE}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 14);
"

## Save the SQL command to a file (can be used for debugging)
if [ "$SAVE_SQL_SCRIPT" = true ]; then
    echo -e "Saving SQL command to a file"
    SQL_SCRIPT_FILE="${OUTPUT_FILE/.parquet/}.sql"
    echo "${SQL_COMMAND}" > "$SQL_SCRIPT_FILE"
fi

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"

## Convert Parquet output to CSV (for ELKI)
if [ "$CONVERT_TO_CSV" = true ]; then
    echo -e "\nConverting Parquet output to CSV"
    
    OUTPUT_FILE_CSV="${OUTPUT_FILE/.parquet/.csv.gz}"
    duckdb -c "COPY (
        SELECT * FROM read_parquet('${OUTPUT_FILE}')
    ) TO '${OUTPUT_FILE_CSV}' (FORMAT CSV, HEADER false, COMPRESSION 'gzip');"
fi

echo -e "\nDone"
