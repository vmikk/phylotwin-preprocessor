#!/bin/bash

## Outlier removal summary
# outlier_removal_summary.sh \
#   -i scores/ \
#   -o outlier_scores.txt.gz \
#   -q 0.5 -t 2


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

## Validate input parameters
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_FILE" ]]; then
    echo -e "Error: Missing required parameters!\n"
    usage
fi

## Threads should be a positive integer
if [[ -n "$THREADS" && "$THREADS" -le 0 ]]; then
    echo -e "Error: Threads must be a positive integer!\n"
    usage
fi

echo -e "\nInput parameters:"
echo "Input directory: $INPUT_DIR"
echo "Outlier threshold: $THRESHOLD"
echo "Output file: $OUTPUT_FILE"
if [[ -n "$THREADS" ]]; then
    echo "Threads: $THREADS"
fi
if [[ -n "$MEMORY" ]]; then
    echo "Memory: $MEMORY"
fi
if [[ -n "$TEMP_DIR" ]]; then
    echo "Temp directory: $TEMP_DIR"
fi


####################################### Merge outlier scores

echo -e "\nMerging outlier scores\n"

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


####################################### Outlier summary

echo -e "\nOutlier summary\n"

SUMMARY_COMMAND=""
## Add configuration settings (if provided)
if [[ -n "$THREADS" ]]; then
    SUMMARY_COMMAND+="
SET threads TO ${THREADS};
"
fi

if [[ -n "$MEMORY" ]]; then
    SUMMARY_COMMAND+="
SET memory_limit = '${MEMORY}';
"
fi

if [[ -n "$TEMP_DIR" ]]; then
    SUMMARY_COMMAND+="
PRAGMA temp_directory='${TEMP_DIR}';
"
fi

SUMMARY_COMMAND+="

-- Reading the combined outlier scores
CREATE TABLE inp AS SELECT * FROM read_csv('${OUTPUT_FILE}', 
      auto_detect = false,
      delim = '\t', 
      header = true,
      columns = {
        'SpeciesKey':   'BIGINT',
        'H3Index':      'VARCHAR',
        'OutlierScore': 'DOUBLE'
      });

-- Create a table with summary statistics in long format
CREATE TEMP TABLE summary AS (

    -- Overall summary
    SELECT 'Total_records' as Variable, CAST(COUNT(*) AS VARCHAR) as Value
    FROM inp
    UNION ALL
    SELECT 'Unique_H3_cells', CAST(COUNT(DISTINCT H3Index) AS VARCHAR)
    FROM inp
    UNION ALL
    SELECT 'Unique_species', CAST(COUNT(DISTINCT SpeciesKey) AS VARCHAR)
    FROM inp
    UNION ALL

    -- Total Outliers
    SELECT 'Total_outliers', CAST(COUNT(*) AS VARCHAR)
    FROM inp
    WHERE OutlierScore > ${THRESHOLD}
    UNION ALL

    -- Average number of outliers per species
    SELECT 'Average_number_of_outliers_per_species', CAST(ROUND(AVG(column_count), 3) AS VARCHAR)
    FROM (
        SELECT SpeciesKey, COUNT(*) AS column_count
        FROM inp
        WHERE OutlierScore > ${THRESHOLD}
        GROUP BY SpeciesKey
    ) t
    UNION ALL

    -- Score statistics
    SELECT 'Minimum_outlier_score', CAST(MIN(OutlierScore) AS VARCHAR)
    FROM inp
    UNION ALL
    SELECT 'Q1_outlier_score', CAST(PERCENTILE_CONT(0.25) WITHIN GROUP (ORDER BY OutlierScore) AS VARCHAR)
    FROM inp
    UNION ALL
    SELECT 'Median_outlier_score', CAST(PERCENTILE_CONT(0.5) WITHIN GROUP (ORDER BY OutlierScore) AS VARCHAR)
    FROM inp
    UNION ALL
    SELECT 'Q3_outlier_score', CAST(PERCENTILE_CONT(0.75) WITHIN GROUP (ORDER BY OutlierScore) AS VARCHAR)
    FROM inp
    UNION ALL
    SELECT 'Maximum_outlier_score', CAST(MAX(OutlierScore) AS VARCHAR)
    FROM inp
);

-- Output both to file and console
COPY (
    SELECT * FROM summary
) TO '${OUTPUT_FILE%.txt.gz}_summary.txt' (HEADER, DELIMITER '\t');

-- Also display to console
SELECT * FROM summary;
"

duckdb -c "${SUMMARY_COMMAND}"

echo -e "..Done\n"
