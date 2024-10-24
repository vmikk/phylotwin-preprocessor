#!/bin/bash

# Isolation Forest: Isolation-Based Anomaly Detection
# HySortOD: Hypercube-Based Outlier Detection

## Path to ELKI
ELKI="${HOME}/bin/elki-bundle-0.8.0.jar"

## Parameters
INPUT=""
OUTPUT=""
METHOD=""
INDEXTYPE=""

## Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --output)
            OUTPUT="$2"
            shift 2
            ;;
        --method)
            METHOD="$2"
            shift 2
            ;;
        --indextype)
            INDEXTYPE="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

## Validate required arguments
if [[ -z $INPUT || -z $OUTPUT || -z $METHOD || -z $INDEXTYPE ]]; then
    echo "Usage: $0 --input <input_file> --output <output_file> --method <METHOD> --indextype <MTree|RStarTree>"
    exit 1
fi

## Validate method
VALID_METHODS=("IsolationForest" "HySortOD")
if [[ ! " ${VALID_METHODS[@]} " =~ " ${METHOD} " ]]; then
    echo "Error: Invalid method. Supported methods are: ${VALID_METHODS[*]}"
    exit 1
fi

## Validate nearest neighbor index type
if [[ $INDEXTYPE != "MTree" && $INDEXTYPE != "RStarTree" ]]; then
    echo "Error: Index type must be either 'MTree' or 'RStarTree'"
    exit 1
fi

## Substitute full class names for index types
if [[ $INDEXTYPE == "MTree" ]]; then
    INDEXTYPEF="tree.metrical.mtreevariants.mtree.MTreeFactory"
    ADDITIONAL_PARAMS="-mtree.distancefunction geo.LatLngDistance"
fi
if [[ $INDEXTYPE == "RStarTree" ]]; then
    INDEXTYPEF="tree.spatial.rstarvariants.rstar.RStarTreeFactory"
    ADDITIONAL_PARAMS="-spatial.bulkstrategy SpatialSortBulkSplit -rtree.bulk.spatial-sort HilbertSpatialSorter"
fi


echo -e "\nRunning ${METHOD} outlier detection with the following parameters:"
echo "Input:       ${INPUT}"
echo "Output:      ${OUTPUT}"
echo "Index type:  $INDEXTYPE"

ALGORITHM_PARAMS=""
case $METHOD in
    IsolationForest)
        # Isolation Forest: Isolation-Based Anomaly Detection
        ALGORITHM_PARAMS="-algorithm outlier.density.IsolationForest -iforest.numtrees 100 -iforest.subsample 256 -iforest.seed 42"
        ;;
    HySortOD)
        # HySortOD: Hypercube-Based Outlier Detection
        ALGORITHM_PARAMS="-algorithm outlier.density.HySortOD -hysortod.b 5 -hysortod.minsplit 100"
        ;;
esac

echo "Algorithm parameters: ${ALGORITHM_PARAMS}"
echo "Index parameters: ${ADDITIONAL_PARAMS}"
echo -e "\n"

NUMRECORDS=$(zcat "${INPUT}" | wc -l)
echo "Number of records in input file: ${NUMRECORDS}"
echo -e "\n"

## Run outlier detection
echo -e "\nRunning outlier detection...\n"

java -jar "${ELKI}" \
      KDDCLIApplication \
      $ALGORITHM_PARAMS \
      -pagefile.pagesize 1024 \
      -dbc.in    "${INPUT}" \
      -db.index  "$INDEXTYPEF" \
      $ADDITIONAL_PARAMS \
      -evaluator NoAutomaticEvaluation \
      -resulthandler tutorial.outlier.SimpleScoreDumper \
      | gzip > "${OUTPUT}"

echo -e "\nOutlier detection completed\n"

