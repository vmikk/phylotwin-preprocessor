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

-- Function to process files in a directory
CREATE OR REPLACE FUNCTION process_tsv_files(directory_path STRING) AS $$
BEGIN
    -- Create a temporary table to store all results
    CREATE TEMP TABLE combined_data AS
    SELECT 
        H3Index,
        OutlierScore,
        regexp_replace(filename, '.*/(.*)\.[^.]*$', '\1') as source_file,
        CAST(regexp_extract(regexp_replace(filename, '.*/(.*)\.[^.]*$', '\1'), '^([0-9]+)_', 1) AS BIGINT) as specieskey
    FROM read_csv_auto(
        (SELECT list_concat(list_transform(
            glob(directory_path || '/*.txt.gz'),
            x -> { 'file': x, 'format': 'tsv', 'header': false, 
                   'columns': {'H3Index': 'VARCHAR', 'OutlierScore': 'DOUBLE'} }
        ))),
        filename=true
    );

    -- Write to parquet file
    COPY (
        SELECT 
            specieskey,
            H3Index,
            OutlierScore
        FROM combined_data
        ORDER BY specieskey, H3Index
    ) TO '${OUTPUT_FILE}' (FORMAT CSV, HEADER true, DELIMITER '\t', COMPRESSION 'gzip');

    -- Generate and print summary statistics
    SELECT 
        COUNT(*) as total_records,
        COUNT(*) FILTER (WHERE OutlierScore >= ${THRESHOLD}) as outliers,
        COUNT(*) FILTER (WHERE OutlierScore < ${THRESHOLD}) as non_outliers,
        MIN(OutlierScore) as min_score,
        MAX(OutlierScore) as max_score,
        AVG(OutlierScore) as avg_score
    FROM combined_data;

    -- Clean up
    DROP TABLE combined_data;
END;
$$;

-- Execute the function
CALL process_tsv_files('${INPUT_DIR}');
"

## Save the SQL command to a file
echo -e "\nSaving SQL command to a file"
echo "${SQL_COMMAND}" > "count_outliers.sql"

## Execute the SQL command
echo -e "\nExecuting DuckDB command"
cat count_outliers.sql | duckdb

echo -e "\nDone"

