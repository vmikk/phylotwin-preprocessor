#!/usr/bin/env Rscript

## Pool species occurrence counts (from different phylogenetic trees)

## Usage:
# Rscript bin/pool_occurrence_counts.R \
#   --inpdir        occ_counts/ \
#   --output_large  pooled_occ_counts_large.tsv \
#   --output_small  pooled_occ_counts_small.tsv


## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("data.table")
load_pckg("optparse")
load_pckg("plyr")


## Define the option parser
option_list <- list(
    make_option(c("-i", "--inpdir"),
        type = "character", default = NULL,
        help = "Input - Directory with species occurrence counts (TSV format)"),
    make_option(c("-o", "--output_large"),
        type = "character", default = NULL,
        help = "Output file with a list of species with large number of occurrences"),
    make_option(c("-s", "--output_small"),
        type = "character", default = NULL,
        help = "Output file with a list of species with small number of occurrences")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Input parameters
INPDIR <- opt$inpdir
OUTPUT_LARGE <- opt$output_large
OUTPUT_SMALL <- opt$output_small

cat("\nInput parameters:\n")
cat("..Input directory:", INPDIR, "\n")
cat("..Output [large]:", OUTPUT_LARGE, "\n")
cat("..Output [small]:", OUTPUT_SMALL, "\n")

## List files in the input directory
fls_large <- list.files(path = INPDIR, pattern = "*_large.txt", full.names = TRUE)
fls_small <- list.files(path = INPDIR, pattern = "*_small.txt", full.names = TRUE)

cat("\nNumber of files [large]:", length(fls_large), "\n")
cat(paste(fls_large, collapse = "\n"))
cat("\n")

cat("\nNumber of files [small]:", length(fls_small), "\n")
cat(paste(fls_small, collapse = "\n"))
cat("\n")

## Load data
cat("\nLoading data:\n")

load_fls <- function(x){
  # x = file name
  res <- fread(x)
  res[ , FileID := sub(pattern = ".txt$", replacement = "", x = basename(x)) ]
  return(res)
}

occ_large <- alply(.data = fls_large, .margins = 1, .fun = load_fls)
occ_small <- alply(.data = fls_small, .margins = 1, .fun = load_fls)

occ_large <- rbindlist(occ_large)
occ_small <- rbindlist(occ_small)

cat("\nKeeping only unique species:\n")

occ_large <- unique(occ_large, by = "specieskey")
occ_small <- unique(occ_small, by = "specieskey")

cat("Number of unique specieskeys [large]:", nrow(occ_large), "\n")
cat("Number of unique specieskeys [small]:", nrow(occ_small), "\n")

cat("\nExporting results:\n")

fwrite(x = occ_large[, .(specieskey) ], file = OUTPUT_LARGE, sep = "\t", row.names = FALSE, col.names = FALSE)
fwrite(x = occ_small[, .(specieskey) ], file = OUTPUT_SMALL, sep = "\t", row.names = FALSE, col.names = FALSE)
