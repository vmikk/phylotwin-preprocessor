#!/bin/bash

## Script to prepare data for spatial outlier detection with ELKI

## Output format:
# LAT, LON, H3, Species

## Usage:
# ./h3_preprocessor.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -r 6 \
#   -s 2495667

## Main workflow:
# - filter records using the `species` column
# - optionally filter records using the `basis_of_record` column
# - convert latitude/longitude coordinates to H3 cells
# - keep unique H3 grids
# - estimate the centroid of H3 grids

## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -r H3_RESOLUTION -s SPECIES [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-d] [-c] [-e EXT_DIR] [-z COMPRESSION]"
    echo "  -i INPUT_FILE      : Input Parquet file path"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -r H3_RESOLUTION   : H3 resolution (0-15)"
    echo "  -s SPECIES         : Species for filtering"
    echo "  -b BASIS_OF_RECORD : Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    echo "  -d                 : Save SQL script for debugging (optional, default: true)"
    echo "  -c                 : Convert Parquet output to CSV (optional, default: true)"
    echo "  -e EXT_DIR         : DuckDB extensions directory path (optional)"
    echo "  -z COMPRESSION     : ZSTD compression level (0-22) (optional, default: 10)"
    exit 1
}

## Initialize variables
INPUT_FILE=""
OUTPUT_FILE=""
H3_RESOLUTION=""
SPECIES=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""
SAVE_SQL_SCRIPT=true
CONVERT_TO_CSV=true
EXT_DIR=""
COMPRESSION_LEVEL="10"

## Parse command-line options
while getopts "i:o:r:s:b:t:m:x:e:dcz:" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        r) H3_RESOLUTION="$OPTARG" ;;
        s) SPECIES="$OPTARG" ;;
        b) BASIS_OF_RECORD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        d) SAVE_SQL_SCRIPT=true ;;
        c) CONVERT_TO_CSV=true ;;
        e) EXT_DIR="$OPTARG" ;;
        z) COMPRESSION_LEVEL="$OPTARG" ;;
        *) usage ;;
    esac
done

## Validate input parameters
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$H3_RESOLUTION" || -z "$SPECIES" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

## Species key validation -- deprecated
# if ! [[ "$SPECIES_KEY" =~ ^[0-9]+$ ]]; then
#     echo -e "Error: Species key must be a positive integer!\n"
#     usage
# fi

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

## Validate compression level
if ! [[ "$COMPRESSION_LEVEL" =~ ^[0-9]+$ ]] || [ "$COMPRESSION_LEVEL" -lt 0 ] || [ "$COMPRESSION_LEVEL" -gt 22 ]; then
    echo -e "Error: Compression level must be an integer between 0 and 22!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "H3 resolution: $H3_RESOLUTION"
echo "Species: ${SPECIES}"
if [[ -n "$BASIS_OF_RECORD" ]]; then
    echo "Basis of record filter: $BASIS_OF_RECORD"
fi
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
if [[ -n "$EXT_DIR" ]]; then
    echo "DuckDB extensions directory: $EXT_DIR"
fi
echo "Parquet compression level (ZSTD): $COMPRESSION_LEVEL"

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

if [[ -n "$EXT_DIR" ]]; then
    SQL_COMMAND+="
SET extension_directory='${EXT_DIR}';
"
fi

SQL_COMMAND+="

-- Install and load the H3 extension
-- INSTALL h3 FROM community;
LOAD h3;

-- Main query
COPY (
    WITH
        inp AS (
            SELECT 
                h3_latlng_to_cell(decimallatitude, decimallongitude, ${H3_RESOLUTION})::UBIGINT AS h3_index,
                species
            FROM 
                read_parquet('${INPUT_FILE}')
            WHERE species = '${SPECIES}'"

# Add basis of record filter if specified
if [[ -n "$BASIS_OF_RECORD" ]]; then
    BASIS_LIST=$(echo "${BASIS_OF_RECORD}" | sed "s/,/','/g")
    BASIS_LIST="'$BASIS_LIST'"
    SQL_COMMAND+=" AND basisofrecord IN (${BASIS_LIST})"
fi

SQL_COMMAND+="),
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
                '${SPECIES}' as Species,
            FROM unique_grids
        )

    SELECT * FROM grid_centroids
) TO '${OUTPUT_FILE}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL ${COMPRESSION_LEVEL});
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
