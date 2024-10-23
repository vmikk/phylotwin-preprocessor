#!/bin/bash

# Main workflow:
# - filter records using the `specieskey` column
# - convert latitude/longitude coordinates to H3 cells
# - keep unique H3 grids
# - estimate the centroid of H3 grids

VARINPP="$1"
VAROUTP="$2"
VARH3RES="$3"
VARSPEC="$4"
SAVE_SQL_SCRIPT="${5:-true}" # Save commands sent to DuckDB (e.g., for debugging purpose)

## Start the SQL command
SQL_COMMAND="
-- Install and load the H3 extension
INSTALL h3 FROM community;
LOAD h3;

COPY (
    WITH
        inp AS (
            SELECT 
                h3_latlng_to_cell(decimallatitude, decimallongitude, ${VARH3RES})::UBIGINT AS h3_index,
                specieskey
            FROM 
                read_parquet('${VARINPP}')
            WHERE specieskey = ${VARSPEC}
        ),
        unique_grids AS (
            SELECT DISTINCT h3_index
            FROM inp
            WHERE h3_index IS NOT NULL
        ),
        grid_centroids AS (
            SELECT 
                h3_index,
                h3_cell_to_lat(h3_index) AS LAT,
                h3_cell_to_lng(h3_index) AS LON
            FROM unique_grids
        )

    SELECT * FROM grid_centroids
) TO '${VAROUTP}' (FORMAT 'parquet', COMPRESSION 'ZSTD', COMPRESSION_LEVEL 14);
"

## Save the SQL command to a file if requested
if [ "${SAVE_SQL_SCRIPT}" = true ]; then
    echo "${SQL_COMMAND}" > "${VAROUTP/.parquet/}.sql"
fi

## Execute the SQL command
duckdb -c "${SQL_COMMAND}"
