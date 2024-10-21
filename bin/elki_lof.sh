#!/bin/bash


## Path to ELKI
ELKI="${HOME}/bin/elki-bundle-0.8.0.jar"

## Parameters
INPUT=""
OUTPUT=""
METHOD=""
GEOMODEL=""
INDEXTYPE=""
K=""

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
        --geomodel)
            GEOMODEL="$2"
            shift 2
            ;;
        --indextype)
            INDEXTYPE="$2"
            shift 2
            ;;
        --k)
            K="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

## Validate required arguments
if [[ -z $INPUT || -z $OUTPUT || -z $METHOD || -z $GEOMODEL || -z $INDEXTYPE ]]; then
    echo "Usage: $0 --input <input_file> --output <output_file> --method <METHOD> --geomodel <geomodel> --indextype <MTree|RStarTree> [--epsilon <epsilon> --minpts <minpts> --k <k>]"
    exit 1
fi

## Validate method
VALID_METHODS=("LOF" "ALOCI" "COF" "FlexibleLOF" "INFLO" "KDEOS" "LDF" "LDOF" "LOCI" "LoOP" "SimplifiedLOF" "ParallelSimplifiedLOF" "SimpleKernelDensityLOF")
if [[ ! " ${VALID_METHODS[@]} " =~ " ${METHOD} " ]]; then
    echo "Error: Invalid method. Supported methods are: ${VALID_METHODS[*]}"
    exit 1
fi

## Validate nearest neighbor index type
if [[ $INDEXTYPE != "MTree" && $INDEXTYPE != "RStarTree" ]]; then
    echo "Error: Index type must be either 'MTree' or 'RStarTree'"
    exit 1
fi

## Validate geometry model
if [[ $GEOMODEL != "WGS84SpheroidEarthModel" && $GEOMODEL != "SphericalVincentyEarthModel" && $GEOMODEL != "SphericalHaversineEarthModel" ]]; then
    echo "Error: Geomodel must be either 'WGS84SpheroidEarthModel' or 'SphericalVincentyEarthModel' or 'SphericalHaversineEarthModel'"
    exit 1
fi

## Substitute full class names for index types
if [[ $INDEXTYPE == "MTree" ]]; then
    INDEXTYPEF="tree.metrical.mtreevariants.mtree.MTreeFactory"
    ADDITIONAL_PARAMS="-mtree.distancefunction geo.LatLngDistance"
fi
if [[ $INDEXTYPE == "RStarTree" ]]; then
    INDEXTYPEF="tree.spatial.rstarvariants.rstar.RStarTreeFactory"
    ADDITIONAL_PARAMS="-rtree.bulk.spatial-sort HilbertSpatialSorter"
fi


echo -e "\nRunning ${METHOD} outlier detection with the following parameters:"
echo "Input:       ${INPUT}"
echo "Output:      ${OUTPUT}"
echo "Geomodel:    $GEOMODEL"
echo "Index type:  $INDEXTYPE"
echo "K:           $K"

ALGORITHM_PARAMS=""
case $METHOD in
    LOF)
        # LOF: Local Outlier Factor
        ALGORITHM_PARAMS="-algorithm outlier.lof.LOF -lof.k ${K}"
        ;;
    ALOCI)
        # Approximate LOCI: Fast Outlier Detection Using the "approximate Local Correlation Integral"
        ALGORITHM_PARAMS="-algorithm outlier.lof.ALOCI"
        ;;
    COF)
        # COF: Connectivity-based Outlier Factor
        ALGORITHM_PARAMS="-algorithm outlier.lof.COF -cof.k ${K}"
        ;;
    INFLO)
        # INFLO: Influenced Outlierness Factor
        ALGORITHM_PARAMS="-algorithm outlier.lof.INFLO -inflo.k ${K}"
        ;;
    LoOP)
        # LoOP: Local Outlier Probabilities
        ALGORITHM_PARAMS="-algorithm outlier.lof.LoOP -loop.kcomp ${K}"
        ;;
    SimplifiedLOF)
        # SimplifiedLOF: Simplified Local Outlier Factor (does not use the reachability distance)
        ALGORITHM_PARAMS="-algorithm outlier.lof.SimplifiedLOF -lof.k ${K}"
        ;;
    SimpleKernelDensityLOF)
        # SimpleKernelDensityLOF: uses a simple kernel density estimation instead of the local reachability density
        ALGORITHM_PARAMS="-algorithm outlier.lof.SimpleKernelDensityLOF -lof.k ${K}"
        ;;
    VarianceOfVolume)
        # Variance of Volume for outlier detection (The volume is estimated by the distance to the k-nearest neighbor, then the variance of volume is computed)
        ALGORITHM_PARAMS="-algorithm outlier.lof.VarianceOfVolume -vov.k ${K}"
        ;;
esac


echo "Algorithm parameters: ${ALGORITHM_PARAMS}"
echo -e "\n"

NUMRECORDS=$(zcat "${INPUT}" | wc -l)
echo "Number of records in input file: ${NUMRECORDS}"
echo -e "\n"

## Run outlier detection
echo -e "\nRunning outlier detection...\n"

java -jar "${ELKI}" \
      KDDCLIApplication \
      $ALGORITHM_PARAMS \
      -algorithm.distancefunction geo.LatLngDistance \
      -pagefile.pagesize 1024 \
      -spatial.bulkstrategy SpatialSortBulkSplit \
      -geo.model "$GEOMODEL" \
      -dbc.in    "${INPUT}" \
      -db.index  "$INDEXTYPEF" \
      $ADDITIONAL_PARAMS \
      -evaluator NoAutomaticEvaluation \
      -resulthandler tutorial.outlier.SimpleScoreDumper \
      | gzip > "${OUTPUT}"

echo -e "\nOutlier detection completed\n"


## Default result handler
## NB! in the output rows are doubled! (unorted and sorted?)
# -resulthandler ResultWriter -out.gzip \
# -out "${OUTPUT}"
# 
# parse_lofs (){
#   zcat "$1" | awk '/^#/ { next } {
#     printf "%s\t%s\t%s\t%s\t", $1, $2, $3, $4
#     match($0, /=-?[0-9]+(\.[0-9]+)?$/)
#     if (RSTART > 0) {
#         print substr($0, RSTART+1, RLENGTH-1)
#     } else {
#         print ""
#     }}'
# }
# 
# parse_lofs "${OUTPUT}" > "${OUTPUT}".tsv
