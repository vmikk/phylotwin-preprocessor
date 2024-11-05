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
    make_option(c("-o", "--output"),
        type = "character", default = NULL,
        help = "Output file (TSV format)"),
    make_option(c("-v", "--verbose"),
        action = "store_true", default = TRUE,
        help = "Print verbose output")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Input parameters
GBIF    <- opt$gbif
TREE    <- opt$tree
OUTPUT  <- opt$output
VERBOSE <- opt$verbose

cat("\nInput parameters:\n")
cat("..GBIF record summary:", GBIF, "\n")
cat("..Phylogenetic tree:", TREE, "\n")
cat("..Output file:", OUTPUT, "\n")
cat("..Verbose:", VERBOSE, "\n")

## Load GBIF records
cat("\nLoading GBIF record summary\n")
CNT <- fread(file = GBIF, sep = "\t")     # count of the number of records per species

## Load phylogenetic tree
cat("Loading phylogenetic tree\n")
TRE <- read.tree(file = TREE)

## Parse species names and GBIF specieskeys
SPP <- data.table(TipLabel = TRE$tip.label)
SPP[ , c("SpeciesKey", "Species") := tstrsplit(x = TipLabel, split = "___", keep = 1:2) ]

## Subset the GBIF records to the species in the tree
CNTS <- CNT[ specieskey %in% SPP$SpeciesKey ]

## Summary
RES <- data.table(
  Tree                       = basename(TREE),
  SpeciesInTree              = nrow(SPP),
  SpeciesInGBIF              = nrow(CNTS),
  SpeciesNotInGBIF           = sum(! SPP$SpeciesKey %in% CNT$specieskey),
  NumUniqueOccurrences       = sum(CNTS$NumRecords),
  NumNonRedundantOccurrences = sum(CNTS$NumUniqCoordsRounded2dp) )

if(VERBOSE){
  print(RES)
}

## Save the results
fwrite(x = RES, file = OUTPUT, sep = "\t")

