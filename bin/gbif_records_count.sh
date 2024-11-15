#!/bin/bash

## Count number of records per species in GBIF data

## Usage example:
# ./gbif_records_count.sh \
#  -i /path/to/GBIF/dump \
#  -o GBIF_record_counts.parquet


## Function to display usage information
usage() {
    echo "Usage: $0 -i INPUT -o OUTPUT [-b BASIS_OF_RECORD] [-t THREADS] [-m MEMORY] [-x TEMP_DIR] [-d] [-c]"
    echo "  -i INPUT          : Input directory containing Parquet files"
    echo "  -o OUTPUT         : Output Parquet file path"
    echo "  -b BASIS_OF_RECORD: Comma-separated list of basis of record values to include (optional)"
    echo "  -t THREADS        : Number of CPU threads to use (optional)"
    echo "  -m MEMORY         : Memory limit (e.g., '100GB') (optional)"
    echo "  -x TEMP_DIR       : Temporary directory path (optional)"
    echo "  -d                : Save SQL script for debugging (optional, default: true)"
    echo "  -c                : Convert Parquet output to CSV (optional, default: true)"
    exit 1
}

## Initialize variables
INPUT=""
OUTPUT=""
BASIS_OF_RECORD=""
THREADS=""
MEMORY=""
TEMP_DIR=""
SAVE_SQL_SCRIPT=true
CONVERT_TO_CSV=true

## Parse command-line options
while getopts "i:o:b:t:m:x:dc" opt; do
    case $opt in
        i) INPUT="$OPTARG" ;;
        o) OUTPUT="$OPTARG" ;;
        b) BASIS_OF_RECORD="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        m) MEMORY="$OPTARG" ;;
        x) TEMP_DIR="$OPTARG" ;;
        d) SAVE_SQL_SCRIPT=true ;;
        c) CONVERT_TO_CSV=true ;;
        *) usage ;;
    esac
done

## Validate input parameters
if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## View user-supplied parameters
echo -e "\nInput parameters:"
echo "Input directory: $INPUT"
echo "Output file: $OUTPUT"
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
-- Load the data from Parquet files, exclude records with missing coordinates
CREATE VIEW species_data AS
SELECT specieskey, decimallatitude, decimallongitude, kingdom, phylum, class, \"order\", family, genus, species
FROM read_parquet('${INPUT}/*')
WHERE decimallatitude IS NOT NULL AND decimallongitude IS NOT NULL"

# Add basis of record filter if specified
if [[ -n "$BASIS_OF_RECORD" ]]; then
    # Convert comma-separated string to SQL IN clause format
    BASIS_LIST=$(echo "${BASIS_OF_RECORD}" | sed "s/,/','/g")
    BASIS_LIST="'$BASIS_LIST'"
    SQL_COMMAND+=" AND basisofrecord IN (${BASIS_LIST})"
fi

SQL_COMMAND+=";

-- Count the number of records per unique 'specieskey' (raw counts)
CREATE TABLE raw_counts AS
SELECT 
    specieskey,
    -- Store one copy of taxonomy info per species
    MAX(kingdom) as kingdom,
    MAX(phylum) as phylum,
    MAX(class) as class,
    MAX(\"order\") as \"order\",
    MAX(family) as family,
    MAX(genus) as genus,
    MAX(species) as species,
    COUNT(*) AS NumRecords
FROM species_data
GROUP BY specieskey;

--  Count unique coordinates rounded to 2 decimal places
CREATE TABLE unique_coordinate_counts_2dp AS
SELECT specieskey, COUNT(DISTINCT CONCAT(ROUND(decimallatitude, 2), ROUND(decimallongitude, 2))) AS NumUniqCoordsRounded2dp
FROM species_data
GROUP BY specieskey;

-- Create a joined table that combines pure counts and unique coordinate counts
CREATE TABLE species_counts_combined AS
SELECT 
    r.specieskey,
    r.kingdom,
    r.phylum,
    r.class,
    r.\"order\",
    r.family,
    r.genus,
    r.species,
    r.NumRecords, 
    u2.NumUniqCoordsRounded2dp
FROM raw_counts r
JOIN unique_coordinate_counts_2dp u2 ON r.specieskey = u2.specieskey;

-- Export the joined result to a Parquet file
COPY species_counts_combined TO '${OUTPUT}' (FORMAT PARQUET, COMPRESSION 'ZSTD', COMPRESSION_LEVEL 14);

-- Clean up intermediate tables
DROP TABLE raw_counts;
DROP TABLE unique_coordinate_counts_2dp;
DROP VIEW species_data;
"

## Save the SQL command to a file (can be used for debugging)
if [ "$SAVE_SQL_SCRIPT" = true ]; then
    echo -e "Saving SQL command to a file"
    SQL_SCRIPT_FILE="${OUTPUT/.parquet/}.sql"
    echo "${SQL_COMMAND}" > "$SQL_SCRIPT_FILE"
fi

## Execute the SQL command
echo -e "\nExecuting DuckDB command"

duckdb -c "${SQL_COMMAND}"

## Convert Parquet output to CSV
if [ "$CONVERT_TO_CSV" = true ]; then
    echo -e "\nConverting Parquet output to CSV"
    
    OUTPUT_CSV="${OUTPUT/.parquet/.csv.gz}"
    duckdb -c "COPY (
        SELECT * FROM read_parquet('${OUTPUT}')
        ORDER BY kingdom, phylum, class, \"order\", family, genus, species, specieskey
    ) TO '${OUTPUT_CSV}' (FORMAT CSV, HEADER true, DELIMITER '\t', COMPRESSION 'gzip');"
fi

echo -e "\nDone"
