#!/bin/bash

## Usage:
# ./elki_clustering.sh \
#   --input input.csv \
#   --output output.csv.gz \
#   --method OPTICS \
#   --geomodel WGS84SpheroidEarthModel \
#   --epsilon 100000 \
#   --minpts 5 \
#   --indextype RStarTree

# Epsilon is in meters

## Path to ELKI
ELKI="${HOME}/bin/elki-bundle-0.8.0.jar"

## Parameters
INPUT=""
OUTPUT=""
METHOD=""
GEOMODEL=""
EPSILON=""
MINPTS=""
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
        --geomodel)
            GEOMODEL="$2"
            shift 2
            ;;
        --epsilon)
            EPSILON="$2"
            shift 2
            ;;
        --minpts)
            MINPTS="$2"
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
if [[ -z $INPUT || -z $OUTPUT || -z $METHOD || -z $GEOMODEL || -z $EPSILON || -z $MINPTS || -z $INDEXTYPE ]]; then
    echo "Usage: $0 --input <input_file> --output <output_file> --method <OPTICS|DBSCAN> --geomodel <geomodel> --epsilon <epsilon> --minpts <minpts> --indextype <MTree|RStarTree>"
    exit 1
fi

## Validate method
if [[ $METHOD != "OPTICS" && $METHOD != "DBSCAN" ]]; then
    echo "Error: Method must be either OPTICS or DBSCAN"
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
    # ADDITIONAL_PARAMS=""
fi


echo -e "\nRunning ${METHOD} clustering with the following parameters:"
echo "Input:       ${INPUT}"
echo "Output:      ${OUTPUT}"
echo "Geomodel:    $GEOMODEL"
echo "Epsilon:     $EPSILON"
echo "MinPts:      $MINPTS"
echo "Index type:  $INDEXTYPE"
echo -e "\n"

## Run OPTICS clustering
if [[ $METHOD == "OPTICS" ]]; then

    java -jar "${ELKI}" \
      KDDCLIApplication \
      -algorithm clustering.optics.OPTICSXi \
      -opticsxi.algorithm OPTICSHeap \
      -algorithm.distancefunction geo.LatLngDistance \
      -pagefile.pagesize    1024 \
      -spatial.bulkstrategy SortTileRecursiveBulkSplit \
      -geo.model "$GEOMODEL" \
      -dbc.in    "${INPUT}" \
      -db.index  "$INDEXTYPEF" \
      "$ADDITIONAL_PARAMS" \
      -optics.epsilon "$EPSILON" \
      -optics.minpts  "$MINPTS" \
      -opticsxi.xi     0.006 \
      -resulthandler ResultWriter -out.gzip \
      -out "${OUTPUT}"

fi

echo -e "\nClustering completed\n"

