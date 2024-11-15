#!/bin/bash

## Filter and bin GBIF records to H3 grids

## Usage:
# ./filter_and_bin.sh \
#   -i '/path/to/GBIF/dump/*' \
#   -o '/path/to/output.parquet' \
#   -r 4 \
#   -s spkeys.txt

## Main workflow:
# - filter records using the `specieskey` column
# - optionally filter records using the `basis_of_record` column
# - bin and count occurrences per H3 cell, species, and collection year
# - aggregate data sources


## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_FILE -r H3_RESOLUTION -s SPECIES_KEY [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR]"
    echo "  -i INPUT_FILE      : Input Parquet file path"
    echo "  -o OUTPUT_FILE     : Output Parquet file path"
    echo "  -r H3_RESOLUTION   : H3 resolution (0-15)"
    echo "  -s spkeys.txt      : Text file with species keys (one per line)"
    echo "  -b BASIS_OF_RECORD : Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS         : Number of CPU threads to use (optional)"
    echo "  -m MEMORY          : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR        : Temporary directory path (optional)"
    exit 1
}

## Initialize variables
INPUT_FILE=""
OUTPUT_FILE=""
H3_RESOLUTION=""
SPECIES_KEYS=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""

## Parse command-line options
while getopts "i:o:r:s:b:t:m:x" opt; do
    case $opt in
        i) INPUT_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        r) H3_RESOLUTION="$OPTARG" ;;
        s) SPECIES_KEYS="$OPTARG" ;;
        b) BASIS_OF_RECORD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        *) usage ;;
    esac
done

## Validate input parameters
if [[ -z "$INPUT_FILE" || -z "$OUTPUT_FILE" || -z "$H3_RESOLUTION" || -z "$SPECIES_KEYS" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

if ! [[ "$H3_RESOLUTION" =~ ^[0-9]+$ ]] || [ "$H3_RESOLUTION" -lt 0 ] || [ "$H3_RESOLUTION" -gt 15 ]; then
    echo -e "Error: H3 resolution must be an integer between 0 and 15!\n"
    usage
fi

if [ ! -r "$SPECIES_KEYS" ]; then
    echo -e "Error: Species keys file '$SPECIES_KEYS' does not exist or is not readable!\n"
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

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input file: $INPUT_FILE"
echo "Output file: $OUTPUT_FILE"
echo "H3 resolution: $H3_RESOLUTION"
echo "File with species keys: $SPECIES_KEYS"
NN=$(wc -l < "${SPECIES_KEYS}")
echo "..Number of species keys detected: $NN"
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

-- Create a table containing species keys of interest
CREATE TEMP TABLE species_keys AS 
SELECT CAST(column0 AS BIGINT) AS specieskey 
FROM read_csv('${SPECIES_KEYS}', header=false);

-- Bin and count occurrences per H3 cell
COPY (
    WITH
        inp AS (
            SELECT 
                h3_latlng_to_cell(decimallatitude, decimallongitude, ${H3_RESOLUTION})::UBIGINT AS h3_index,
                datasetkey, 
                specieskey,
                year
            FROM 
                read_parquet('${INPUT_FILE}')
            WHERE specieskey IN (SELECT specieskey FROM species_keys) "

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
                specieskey,
                year,
                COUNT(*) as record_count,
                STRING_AGG(DISTINCT CAST(datasetkey AS VARCHAR), ',') as dataset_keys
            FROM inp
            WHERE h3_index IS NOT NULL
            GROUP BY h3_index, specieskey, year
        ),
        aggregated_data AS (
            SELECT 
                h3_h3_to_string(h3_index) as H3,
                specieskey,
                year,
                record_count,
                dataset_keys
            FROM unique_grids
        )

    SELECT * FROM aggregated_data
) TO '${OUTPUT_FILE}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 14);
"

## Save the SQL command to a file
echo -e "\nSaving SQL command to a file"
SQL_SCRIPT_FILE="${OUTPUT_FILE/.parquet/}.sql"
echo "${SQL_COMMAND}" > "$SQL_SCRIPT_FILE"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"
cat "${SQL_SCRIPT_FILE}" | duckdb

echo -e "\nDone"
