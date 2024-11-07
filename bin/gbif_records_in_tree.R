#!/usr/bin/env Rscript

## Count the number of GBIF records represented and not represented in phylogenetic trees

## Usage:
# Rscript bin/gbif_records_in_tree.R \
#   --gbif   data/gbif_records_summary.tsv \
#   --tree   phylogenies/tree.newick \
#   --output gbif_records_in_tree.tsv

## Format of tip labels in the tree:
# SpeciesKey___SpeciesName (e.g., "1925221___Lysandra_coridon")


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("data.table")
load_pckg("ape")
load_pckg("optparse")

# load_pckg("plyr")
# load_pckg("ggplot2")
# load_pckg("patchwork")

# theme_set(theme_bw())


## Define the option parser
option_list <- list(
    make_option(c("-g", "--gbif"),
        type = "character", default = NULL,
        help = "Input - GBIF record summary per species (TSV format)"),
    make_option(c("-t", "--tree"),
        type = "character", default = NULL,
        help = "Input - Tree (Newick format)"),
    make_option(c("-s", "--sourcetree"),
        type = "character", default = NULL,
        help = "Input - Initial tree, before species name matching (Newick format), optional"),
    make_option(c("-o", "--outputstats"),
        type = "character", default = NULL,
        help = "Output file with summary statistics (TSV format)"),
    make_option(c("-x", "--occurrence_threshold"),
        type = "numeric", default = 100,
        help = "Threshold for the number of occurrences"),
    make_option(c("-L", "--output_large_occ"),
        type = "character", default = NULL,
        help = "Output file with a list of species keys with large number of occurrences (TSV format)"),
    make_option(c("-l", "--output_low_occ"),
        type = "character", default = NULL,
        help = "Output file with a list of species keys with low number of occurrences (TSV format)"),
    make_option(c("-v", "--verbose"),
        action = "store_true", default = TRUE,
        help = "Print verbose output")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Input parameters
GBIF         <- opt$gbif
TREE         <- opt$tree
SOURCETREE   <- opt$sourcetree
OUTPUTSTATS  <- opt$outputstats
OCC_THRESHOLD    <- opt$occurrence_threshold
OUTPUT_LARGE_OCC <- opt$output_large_occ
OUTPUT_LOW_OCC   <- opt$output_low_occ
VERBOSE          <- opt$verbose

cat("\nInput parameters:\n")
cat("..GBIF record summary:", GBIF, "\n")
cat("..Phylogenetic tree:", TREE, "\n")
cat("..Initial tree:", SOURCETREE, "\n")
cat("..Output with summary statistics:", OUTPUTSTATS, "\n")
cat("..Occurrence threshold:", OCC_THRESHOLD, "\n")
cat("..Output with a list of species keys with large number of occurrences:", OUTPUT_LARGE_OCC, "\n")
cat("..Output with a list of species keys with low number of occurrences:", OUTPUT_LOW_OCC, "\n")
cat("..Verbose:", VERBOSE, "\n")

## Load GBIF records
cat("\nLoading GBIF record summary\n")
CNT <- fread(file = GBIF, sep = "\t")     # count of the number of records per species

## Load phylogenetic tree
cat("Loading phylogenetic tree\n")
TRE <- read.tree(file = TREE)

## Load initial tree
if(! is.null(SOURCETREE)){
  cat("Loading initial tree (prior to species name matching)\n")
  TRE0 <- read.tree(file = SOURCETREE)
}

## Parse species names and GBIF specieskeys
SPP <- data.table(TipLabel = TRE$tip.label)
SPP[ , c("SpeciesKey", "Species") := tstrsplit(x = TipLabel, split = "___", keep = 1:2) ]

## Subset the GBIF records to the species in the tree
CNTS <- CNT[ specieskey %in% SPP$SpeciesKey ]

## Summary
RES <- data.table(
  Tree                       = basename(TREE),
  SpeciesInTree              = nrow(SPP),                                   # species that passed name matching and are in the tree
  SpeciesInGBIF              = nrow(CNTS),                                  # species that passed name matching and have occurrences in GBIF
  SpeciesNotInGBIF           = sum(! SPP$SpeciesKey %in% CNT$specieskey),   # species that passed name matching but are not in GBIF occurrences
  NumUniqueOccurrences       = sum(CNTS$NumRecords),
  NumNonRedundantOccurrences = sum(CNTS$NumUniqCoordsRounded2dp) )

## Estimate percentages
RES[ , Pct_SpeciesInGBIF    := round(SpeciesInGBIF / SpeciesInTree * 100, digits = 2) ]
RES[ , Pct_SpeciesNotInGBIF := round(SpeciesNotInGBIF / SpeciesInTree * 100, digits = 2) ]

setcolorder(x = RES, neworder = c(
    "Tree", "SpeciesInTree", "SpeciesInGBIF", "SpeciesNotInGBIF",
    "Pct_SpeciesInGBIF", "Pct_SpeciesNotInGBIF"))

## Add initial tree info
if(! is.null(SOURCETREE)){
  
  RES[ , SpeciesInSourceTree      := length(TRE0$tip.label) ]
  RES[ , OverallTreeTipsNotInGBIF := SpeciesInSourceTree - SpeciesInGBIF ]  # includes matched and unmatched species

  # NB. `OverallTreeTipsNotInGBIF` might be overestimated because there could be multiple occurrences of the same species in the initial tree

  ## Estimate percentages
  RES[ , Pct_NameMatchedSpecies       := round(SpeciesInTree / SpeciesInSourceTree * 100, digits = 2) ]
  RES[ , Pct_OverallTreeTipsNotInGBIF := round(OverallTreeTipsNotInGBIF / SpeciesInSourceTree * 100, digits = 2) ]

  setcolorder(
    x = RES,
    neworder = c(
      "Tree", "SpeciesInSourceTree", "SpeciesInTree",
      "SpeciesInGBIF", "SpeciesNotInGBIF",
      "OverallTreeTipsNotInGBIF",
      "Pct_NameMatchedSpecies", "Pct_SpeciesInGBIF", "Pct_SpeciesNotInGBIF",
      "Pct_OverallTreeTipsNotInGBIF",
      "NumUniqueOccurrences", "NumNonRedundantOccurrences"
    )
  )

}

if(VERBOSE){
  cat("\nSummary:\n")
  print(RES$Tree)

  if(! is.null(SOURCETREE)){
    print(RES[ , .(SpeciesInTree, SpeciesInGBIF, SpeciesNotInGBIF) ], row.names = FALSE)
    print(RES[ , .(Pct_SpeciesInGBIF, Pct_SpeciesNotInGBIF)], row.names = FALSE)
    print(RES[ , .(SpeciesInSourceTree, OverallTreeTipsNotInGBIF, Pct_NameMatchedSpecies, Pct_OverallTreeTipsNotInGBIF)], row.names = FALSE)
    print(RES[ , .(NumUniqueOccurrences, NumNonRedundantOccurrences)], row.names = FALSE)
  } else {
    print(RES[ , .(SpeciesInTree, SpeciesInGBIF, SpeciesNotInGBIF) ], row.names = FALSE)
    print(RES[ , .(Pct_SpeciesInGBIF, Pct_SpeciesNotInGBIF)], row.names = FALSE)
    print(RES[ , .(NumUniqueOccurrences, NumNonRedundantOccurrences)], row.names = FALSE)
  }

}

## Save summary statistics
if(! is.null(OUTPUTSTATS)){
  fwrite(x = RES, file = OUTPUTSTATS, sep = "\t")
}


## Check occurrence thresholds
if(! is.null(OCC_THRESHOLD)){

  ## Split species into those with large and small number of occurrences
  CNTS[ , OccGroup := fifelse(NumUniqCoordsRounded2dp >= OCC_THRESHOLD, "Large", "Small") ]

  cat("..Number of species with large number of occurrences:", sum(CNTS$OccGroup %in% "Large"), "\n")
  cat("..Number of species with small number of occurrences:", sum(CNTS$OccGroup %in% "Small"), "\n")

  ## Save the results 
  if(! is.null(OUTPUT_LARGE_OCC)){
    fwrite(x = CNTS[ OccGroup == "Large", .(specieskey, OccGroup, NumUniqCoordsRounded2dp) ], file = OUTPUT_LARGE_OCC, sep = "\t")
  }

  if(! is.null(OUTPUT_LOW_OCC)){
    fwrite(x = CNTS[ OccGroup == "Small", .(specieskey, OccGroup, NumUniqCoordsRounded2dp) ], file = OUTPUT_LOW_OCC, sep = "\t")
  }

}

cat("\nDone\n")
