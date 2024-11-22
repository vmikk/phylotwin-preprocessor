#!/bin/bash

## Different methods for outlier detection:
# - LOF-like or KNN-based methods
# - Clustering-based methods
# - Distance-based method (DBOutlierScore)
# - Spatial methods (SLOM, SOF, SLZ)

# Distance (`-d`) is in meters

## Usage:
# ./elki_outlier.sh \
#   --input input.csv \
#   --output output.csv.gz \
#   --method LOF \
#   --geomodel WGS84SpheroidEarthModel \
#   --indextype RStarTree \
#   --k 5

## NB!
##  Every JVM creates a temporary performance instrumentation file in `/tmp/hsperfdata_<user>/<pid>`
##  With Singularity, there could be PID collisions when running multiple instances of the same container
##  --> use `-XX:-UsePerfData` or `-XX:+PerfDisableSharedMem` with Java

## Path to ELKI
# - check in PATH, /usr/local/bin, and ~/bin
ELKI_JAR="elki-bundle-0.8.0.jar"
if which "${ELKI_JAR}" >/dev/null 2>&1; then
    ELKI=$(which "${ELKI_JAR}")
elif [ -f "/usr/local/bin/${ELKI_JAR}" ]; then
    ELKI="/usr/local/bin/${ELKI_JAR}"
elif [ -f "${HOME}/bin/${ELKI_JAR}" ]; then
    ELKI="${HOME}/bin/${ELKI_JAR}"
else
    echo "Error: Could not find ${ELKI_JAR} in PATH, /usr/local/bin, or ${HOME}/bin"
    echo "Please ensure ${ELKI_JAR} is installed in one of these locations"
    exit 1
fi

## Parameters
INPUT=""
OUTPUT=""
METHOD=""
GEOMODEL=""
INDEXTYPE=""
K=""   # number of nearest neighbors
D=""   # distance threshold (for DBOutlierScore and DBSCANOutlierDetection)

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
        --d)
            D="$2"
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
    echo "Usage: $0 --input <input_file> --output <output_file> --method <METHOD> --geomodel <geomodel> --indextype <MTree|RStarTree> [--k <k> --d <d>]"
    exit 1
fi

## Validate method
VALID_METHODS=(\
   "LOF" "ALOCI" "COF" "FlexibleLOF" "INFLO" "KDEOS" "LDF" \
   "LDOF" "LOCI" "LoOP" "SimplifiedLOF" "SimpleKernelDensityLOF" \
   "VarianceOfVolume" "KNNOutlier" "ODIN" "SOS" "KNNDD" "KNNSOS" \
   "KNNWeightOutlier" "LocalIsolationCoefficient" "DBOutlierScore" \
   "DBSCANOutlierDetection" "EMOutlier" "KMeansOutlierDetection" "KMeansMinusMinusOutlierDetection" "GLOSH" \
   "SLOM" "SOF" "CTLuScatterplotOutlier" \
   "DWOF" "COP" "SimpleCOP" "OPTICSOF" \
   "TrivialAverageCoordinateOutlier")

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

## Validate k and d (both should be non-negative)
## TODO.. but not for all methods
# if [[ $K -lt 0 || $D -lt 0 ]]; then
#     echo "Error: K and D must be non-negative"
#     exit 1
# fi

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
echo "Geomodel:    $GEOMODEL"
echo "Index type:  $INDEXTYPE"
echo "K:           $K"
echo "D:           $D"

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
    LDOF)
        # LDOF:Local Distance-Based Outlier Factor
        ALGORITHM_PARAMS="-algorithm outlier.lof.LDOF -ldof.k ${K}"
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
    KNNOutlier)
        # Outlier Detection based on the distance of an object to its k nearest neighbor
        ALGORITHM_PARAMS="-algorithm outlier.distance.KNNOutlier -knno.k ${K}"
        ;;
    ODIN)
        # ODIN: Outlier detection using k-nearest neighbour graph
        ALGORITHM_PARAMS="-algorithm outlier.distance.ODIN -odin.k ${K}"
        ;;
    SOS)
        # SOS: Stochastic Outlier Selection
        ALGORITHM_PARAMS="-algorithm outlier.distance.SOS -sos.perplexity 4.5"
        ;;
    KNNDD)
        # KNNDD: Nearest Neighbor Data Description (original k = 1)
        ALGORITHM_PARAMS="-algorithm outlier.distance.KNNDD -knndd.k ${K}"
        ;;
    KNNSOS)
        # KNNSOS: k-Nearest-Neighbor Stochastic Outlier Selection (default k = 16; perplexity = k/3)
        ALGORITHM_PARAMS="-algorithm outlier.distance.KNNSOS -sos.k ${K}"
        ;;
    KNNWeightOutlier)
        # KNNWeightOutlier: Outlier detection based on the sum of distances of an object to its k nearest neighbors
        ALGORITHM_PARAMS="-algorithm outlier.distance.KNNWeightOutlier -knnwod.k ${K}"
        ;;
    LocalIsolationCoefficient)
        # Local Isolation Coefficient
        ALGORITHM_PARAMS="-algorithm outlier.distance.LocalIsolationCoefficient -lic.k ${K}"
        ;;
    DBOutlierScore)
        # Distance Based Outlier Score (Compute percentage of neighbors in the given neighborhood with size d)
        ALGORITHM_PARAMS="-algorithm outlier.distance.DBOutlierScore -dbod.d ${D}"
        ;;
    DBSCANOutlierDetection)
        # DBSCAN Outlier Detection: Outlier Detection based on the Generalized DBSCAN clustering
        ALGORITHM_PARAMS="-algorithm outlier.clustering.DBSCANOutlierDetection -gdbscan.neighborhood EpsilonNeighborPredicate -gdbscan.core-model -dbscan.epsilon ${D} -dbscan.minpts ${K}"
        ;;
    EMOutlier)
        # EM Outlier: Outlier Detection based on the generic EM clustering
        ALGORITHM_PARAMS="-algorithm outlier.clustering.EMOutlier -em.k ${K}"
        ;;
    # KMeansOutlierDetection)
    #     # Outlier detection by using k-means clustering
    #     # NB! should be used with squared Euclidean distance only!
    #     ALGORITHM_PARAMS="-algorithm outlier.clustering.KMeansOutlierDetection -kmeans.k ${K}"
    #     ;;
    KMeansMinusMinusOutlierDetection)
        # KMeansMinusMinusOutlierDetection: (k-means--) A Unified Approach to Clustering and Outlier Detection
        ALGORITHM_PARAMS="-algorithm outlier.clustering.KMeansMinusMinusOutlierDetection -kmeans.k ${K}"
        ;;
    GLOSH)
        # GLOSH: Global-Local Outlier Scores from Hierarchies
        ALGORITHM_PARAMS="-algorithm outlier.clustering.GLOSH -hdbscan.minPts ${K} -hdbscan.minclsize 1"
        ;;
    SLOM)
        # SLOM: Spatial local outlier measure 
        ALGORITHM_PARAMS=" -algorithm outlier.spatial.SLOM -neighborhood PrecomputedKNearestNeighborNeighborhood -neighborhood.distancefunction geo.LatLngDistance -neighborhood.k ${K}"
        ;;
    SOF)
        # SOF: Spatial outlier factor
        ALGORITHM_PARAMS="-algorithm outlier.spatial.SOF -neighborhood PrecomputedKNearestNeighborNeighborhood -neighborhood.distancefunction geo.LatLngDistance -neighborhood.k ${K}"
        ;;
    CTLuScatterplotOutlier)
        # CTLuScatterplotOutlier (aka SLZ): Spatial Outlier Detection Algorithm using linear regression of attributes and the mean of their neighbors
        ALGORITHM_PARAMS="-algorithm outlier.spatial.CTLuScatterplotOutlier -neighborhood PrecomputedKNearestNeighborNeighborhood -neighborhood.distancefunction geo.LatLngDistance -neighborhood.k ${K}"
        ;;
    DWOF)
        # DWOF: Dynamic Window Outlier Factor
        ALGORITHM_PARAMS="-algorithm outlier.DWOF -dwof.k ${K} -dwof.delta 1.1"
        ;;
    COP)
        # COP: Correlation Outlier Probability
        ALGORITHM_PARAMS="-algorithm outlier.COP -cop.k ${K}"
        ;;
    SimpleCOP)
        # SimpleCOP: Correlation Outlier Probability
        ALGORITHM_PARAMS="-algorithm outlier.SimpleCOP -cop.k 5"
        ;;
    OPTICSOF)
        # OPTICS-OF: Identifying Local Outliers
        ALGORITHM_PARAMS="-algorithm outlier.OPTICSOF -optics.minpts ${K}"
        ;;
    TrivialAverageCoordinateOutlier)
        # TrivialAverageCoordinateOutlier: Trivial method that takes the average of all dimensions as outlier score
        ALGORITHM_PARAMS="-algorithm outlier.trivial.TrivialAverageCoordinateOutlier"
        ;;
esac

    ##### To be added

    # ParallelSimplifiedLOF)
    #     # ParallelSimplifiedLOF: Parallel implementation of Simplified-LOF Outlier detection using processors
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.ParallelSimplifiedLOF"
    #     ;;
    # ParallelLOF)
    #     # ParallelLOF: Parallel implementation of the LOF Algorithm
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.ParallelLOF"
    #     ;;

    # -algorithm outlier.clustering.NoiseAsOutliers ++
    #     -algorithm optics.OPTICSXi / -algorithm dbscan.LSDBC / .. etc

    # CBLOF)
    #     # CBLOF: Cluster-based local outlier factor
    #     ALGORITHM_PARAMS="-algorithm outlier.clustering.CBLOF -cblof.alpha ? -cblof.beta ? -kmeans.k ${K}"
    #     # alpha = The ratio of the size that separates the large clusters from the small clusters.
    #     # beta = The minimal ratio between two consecutive clusters (when ordered descending by size) at which the boundary between the large and small clusters is set.
    #     ;;


    ##### Algorithms requiring additional parameters:

    # LOCI)
    #     # LOCI:  Local Correlation Integral (outlier = LOCI score > 3)
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.LOCI -loci.rmax ${D} -loci.nmin 20 -loci.alpha 0.5"
    #     ;;

    # FlexibleLOF)
    #     # FlexibleLOF: Local Outlier Factor with additional options
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.FlexibleLOF -lof.krefer ${K} -lof.kreach ${K2}"
    #     ;;

    # KDEOS)
    #     # KDEOS: Kernel Density Estimator Outlier Score
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.KDEOS -kdeos.k.min ${K} -kdeos.k.max ${K2}"  # -kdeos.kernel.minbw
    #     ;;

    # LDF)
    #     # LDF: Outlier Detection with Kernel Density Functions
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.LDF -ldf.k ${K} -ldf.h ${H}"
    #     ;;

    # OnlineLOF)
    #     # OnlineLOF: Incremental version of the LOF Algorithm, supports insertions and removals
    #     ALGORITHM_PARAMS="-algorithm outlier.lof.OnlineLOF ..."
    #     ;;


echo "Algorithm parameters: ${ALGORITHM_PARAMS}"
echo "Index parameters: ${ADDITIONAL_PARAMS}"
echo -e "\n"

NUMRECORDS=$(zcat "${INPUT}" | wc -l)
echo "Number of records in input file: ${NUMRECORDS}"
echo -e "\n"

## Run outlier detection
echo -e "\nRunning outlier detection...\n"

java \
  -XX:-UsePerfData \
  -jar "${ELKI}" \
      KDDCLIApplication \
      $ALGORITHM_PARAMS \
      -algorithm.distancefunction geo.LatLngDistance \
      -pagefile.pagesize 1024 \
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
