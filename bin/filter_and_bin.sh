#!/bin/bash

## Filter and bin GBIF records to H3 grids

## Usage:
# ./filter_and_bin.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -r 4 \
#   -s spnames.txt

## Main workflow:
# - filter records using the `species` column
# - optionally filter records using the `basis_of_record` column
# - bin and count occurrences per H3 cell, species, and collection year
# - aggregate data sources


## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT_FILE -r H3_RESOLUTION -s SPECIES [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-e EXT_DIR] [-z COMPRESSION]"
    echo "  -i INPUT           : Input Parquet file or directory path"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -r H3_RESOLUTION   : H3 resolution (0-15)"
    echo "  -s spnames.txt     : Text file with species names (one per line)"
    echo "  -b BASIS_OF_RECORD : Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    echo "  -e EXT_DIR         : DuckDB extensions directory path (optional)"
    echo "  -z COMPRESSION     : ZSTD compression level (0-22) (optional, default: 10)"
    exit 1
}

## Initialize variables
INPUT=""
OUTPUT_FILE=""
H3_RESOLUTION=""
SPECIES=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""
EXT_DIR=""
COMPRESSION_LEVEL="10"

## Parse command-line options
while getopts "i:o:r:s:b:t:m:x:e:z:" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
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

## Validate input parameters
if [[ -z "$INPUT" || -z "$OUTPUT_FILE" || -z "$H3_RESOLUTION" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

## If file with species names is not provided, create it from the input parquet files
if [ ! -r "$SPECIES" ]; then
    echo -e "\n..No file with species names provided, will create it from the input parquet files"
    
    SPECIES="${OUTPUT_FILE/.parquet/}__species_names.txt"
    echo "COPY (
        SELECT DISTINCT species
        FROM read_parquet('${INPUT}')
        WHERE species IS NOT NULL
        ORDER BY species
    ) TO '${SPECIES}' (HEADER FALSE);" | duckdb

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

## Check if the input is a directory
if [[ -d "$INPUT" ]]; then
    echo -e "\n..Specified input is a directory, will process all parquet files in the directory"
    INPUT="${INPUT}/*"
fi

## Validate compression level
if ! [[ "$COMPRESSION_LEVEL" =~ ^[0-9]+$ ]] || [ "$COMPRESSION_LEVEL" -lt 0 ] || [ "$COMPRESSION_LEVEL" -gt 22 ]; then
    echo -e "Error: Compression level must be an integer between 0 and 22!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input: $INPUT"
echo "Output file: $OUTPUT_FILE"
echo "H3 resolution: $H3_RESOLUTION"
echo "File with species names: $SPECIES"
NN=$(wc -l < "${SPECIES}")
echo "..Number of species names detected: $NN"
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

-- Create a table containing species keys of interest
CREATE TEMP TABLE species_names AS 
SELECT CAST(column0 AS VARCHAR) AS species 
FROM read_csv('${SPECIES}', header=false);

-- Bin and count occurrences per H3 cell
COPY (
    WITH
        inp AS (
            SELECT 
                h3_latlng_to_cell(decimallatitude, decimallongitude, ${H3_RESOLUTION})::UBIGINT AS h3_index,
                datasetkey, 
                species,
                year
            FROM 
                read_parquet('${INPUT}')
            WHERE species IN (SELECT species FROM species_names) "

# Add basis of record filter if specified
if [[ -n "$BASIS_OF_RECORD" ]]; then
    BASIS_LIST=$(echo "${BASIS_OF_RECORD}" | sed "s/,/','/g")
    BASIS_LIST="'$BASIS_LIST'"
    SQL_COMMAND+=" AND basisofrecord IN (${BASIS_LIST})"
fi

SQL_COMMAND+="),
        unique_grids AS (
            SELECT 
                h3_index,
                species,
                year,
                COUNT(*) as record_count,
                LIST(DISTINCT datasetkey) as dataset_keys
            FROM inp
            WHERE h3_index IS NOT NULL
            GROUP BY h3_index, species, year
        ),
        aggregated_data AS (
            SELECT 
                h3_h3_to_string(h3_index) as H3,
                species,
                year,
                record_count,
                dataset_keys
            FROM unique_grids
        )

    SELECT * FROM aggregated_data
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
