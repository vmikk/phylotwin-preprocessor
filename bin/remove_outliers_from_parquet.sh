#!/bin/bash

## Usage:
# remove_outliers_from_parquet.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -w 'ELKI_outlier_scores.txt.gz' \
#   -r 6 \
#   -s "Anchoa_mitchilli" \
#   -b "PRESERVED_SPECIMEN,MATERIAL_CITATION,MACHINE_OBSERVATION"

## NB!
# - currently, only binary outlier scores are implemented
#   (threshold is hard-coded to 0.5)

## Main workflow:
# - filter records using the `species` column
# - optionally filter records using the `basis_of_record` column
# - exclude occurrences from H3 cell identified as outliers


## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -w OUTLIER_SCORES -s SPECIES [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-e EXT_DIR] [-z COMPRESSION]"
    echo "  -i INPUT_FILE      : Input Parquet file path (or directory with multiple files)"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -w OUTLIER_SCORES  : Tab-delimited file with H3 cell IDs and outlier scores"
    echo "  -r H3_RESOLUTION   : H3 resolution (0-15)"
    echo "  -s SPECIES         : Species for filtering"
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
H3_RESOLUTION=""
SPECIES=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""
EXT_DIR=""
COMPRESSION_LEVEL="10"

## Parse command-line options
while getopts "i:o:w:r:s:b:t:m:x:e:z:" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        w) OUTLIER_SCORES="$OPTARG" ;;
        r) H3_RESOLUTION="$OPTARG" ;;
        s) SPECIES="$OPTARG" ;;
        b) BASIS_OF_RECORD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        e) EXT_DIR="$OPTARG" ;;
        z) COMPRESSION_LEVEL="$OPTARG" ;;
        *) usage ;;
    esac
done

## Validate required input parameters
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$OUTLIER_SCORES" || -z "$SPECIES" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## Check that H3 resolution is a valid integer between 0 and 15
if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

## Check if outlier scores file exists
if [[ ! -e "$OUTLIER_SCORES" ]]; then
    echo -e "Error: Outlier scores file not found!\n"
    usage
fi

## Validate species key -- deprecated
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
echo "File with outlier scores: $OUTLIER_SCORES"
echo "Species: $SPECIES"
if [[ -n "$BASIS_OF_RECORD" ]]; then
    echo "Basis of record filter: $BASIS_OF_RECORD"
fi
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
    echo "Extensions directory: $EXT_DIR"
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

-- Load H3 grid cell IDs with outlier scores of interest
CREATE TEMP TABLE outlier_scores AS 
SELECT 
  CAST(column0 AS VARCHAR) AS h3_index, 
  CAST(column1 AS VARCHAR) AS species,
  CAST(column2 AS FLOAT) AS outlier_score 
FROM read_csv('${OUTLIER_SCORES}', header=false, delim = '\t');

-- Filter records and exclude outlier cells
COPY (
    WITH filtered_records AS (
        SELECT 
            *,
            h3_latlng_to_cell_string(decimallatitude, decimallongitude, ${H3_RESOLUTION})::VARCHAR AS h3_index
        FROM read_parquet('${INPUT_FILE}')
        WHERE species = '${SPECIES}'"

# Add basis of record filter if specified
if [[ -n "$BASIS_OF_RECORD" ]]; then
    BASIS_LIST=$(echo "${BASIS_OF_RECORD}" | sed "s/,/','/g")
    BASIS_LIST="'$BASIS_LIST'"
    SQL_COMMAND+=" AND basisofrecord IN (${BASIS_LIST})"
fi

SQL_COMMAND+="
    )
    SELECT r.* 
    FROM filtered_records r
    LEFT JOIN outlier_scores o ON r.h3_index = o.h3_index
    WHERE o.outlier_score < 0.5 OR o.outlier_score IS NULL
) TO '${OUTPUT_FILE}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL ${COMPRESSION_LEVEL});
"

## Save the SQL command to a file
echo -e "\nSaving SQL command to a file"
SQL_SCRIPT_FILE="${OUTPUT_FILE/.parquet/}.sql"
echo "${SQL_COMMAND}" > "$SQL_SCRIPT_FILE"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"
cat "${SQL_SCRIPT_FILE}" | duckdb

echo -e "\nDone"
