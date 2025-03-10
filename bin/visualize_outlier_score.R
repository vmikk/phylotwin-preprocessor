#!/usr/bin/env Rscript

## Visualize outlier scores

## Function to load packages
load_pckg <- function(pkg = "data.table"){
    suppressPackageStartupMessages( library(package = pkg, character.only = TRUE) )
    cat(paste(pkg, packageVersion(pkg), "\n"))
}

cat("Loading packages:\n")
load_pckg("data.table")
load_pckg("plyr")
load_pckg("sf")
load_pckg("rnaturalearth")
load_pckg("ggplot2")
load_pckg("patchwork")
load_pckg("optparse")

theme_set(theme_bw())


## Define the option parser
option_list <- list(
    make_option(c("-i", "--input"),
        type = "character", default = NULL,
        help = "Input coordinates (CSV format)"),
    make_option(c("-o", "--outlier"),
        type = "character", default = NULL,
        help = "Outlier scores (TSV format)"),
    make_option(c("-t", "--threshold"),
        type = "numeric", default = NULL,
        help = "Fixed threshold for outlier scores"),
    make_option(c("-p", "--plot"),
        type = "character", default = NULL,
        help = "Output plots (PNG format)"),
    make_option(c("-c", "--cropmap"),
        type = "logical", default = FALSE,
        help = "Crop map to the bounding box of the coordinates")
)

## Parse the command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

## Input parameters
INPUT     <- opt$input
OUTLIER   <- opt$outlier
THRESHOLD <- opt$threshold
PLOT      <- opt$plot
CROP      <- opt$cropmap

cat("\nInput parameters:\n")
cat("..Input coordinates:", INPUT, "\n")
cat("..Outlier scores:", OUTLIER, "\n")
cat("..Fixed threshold:", THRESHOLD, "\n")
cat("..Output plots:", PLOT, "\n")
cat("..Crop map:", CROP, "\n")

## Load world map
cat("\n..Loading world map\n")
world <- ne_countries(scale = "medium", returnclass = "sf")

## Load coordinates
cat("..Loading coordinates\n")
CRD <- fread(file = INPUT, header = F, col.names = c("Latitude", "Longitude", "ID"))

## Load outlier scores
cat("..Loading outlier scores\n")
SCR <- fread(file = OUTLIER, header = F, col.names = c("ID", "OutlierScore"))

## Combine coordinates and outlier scores
cat("..Combining coordinates and outlier scores\n\n")
CRD$OutlierScore <- SCR$OutlierScore

## Score summary
cat("..Score summary\n")
print(summary(CRD$OutlierScore))

## Find bounding box of the coordinates
if(CROP == TRUE){
  cat("..Finding bounding box of the coordinates\n")
  bbox <- list(
    x = c(min(CRD$Longitude), max(CRD$Longitude)),
    y = c(min(CRD$Latitude),  max(CRD$Latitude))
  )
}

## Plotting function
## with optional outlier bubbles (circle around point proportional to the score)
show_occ <- function(x, outlier_mode = NULL, threshold = NULL, threshold_multiplier = 5, verbose = TRUE){
  # x <- copy(CRD)
  # outlier_mode <- c("low", "medium", "high")  # quantile-based (ArcGIS-style) 
  # outlier_mode <- "threshold"                 # fixed outlier threshold
  # threshold_multiplier <- 5  # since scores can be too hight to display, restrict max values

  x <- copy(x)  # to avoid modifications to the original data (outside function)

  ## Quantile-based thresholds for outlier score (as in ArcGIS)
  if(!is.null(outlier_mode)){

    if(!outlier_mode %in% "threshold"){

      q3  <- quantile(x = x$OutlierScore, probs = 0.75, na.rm = TRUE)
      iqr <- IQR(x$OutlierScore, na.rm = TRUE)
      if(outlier_mode %in% "low")   { iqr <- iqr * 2   }
      if(outlier_mode %in% "medium"){ iqr <- iqr * 1.5 }
      thrsh <- q3 + iqr
    
    } else { # fixed user-provided threshold
      thrsh <- threshold
    }

    if(any(is.infinite(x$OutlierScore))){
      mx <- as.numeric(thrsh * 2)
        x[ is.infinite(OutlierScore), OutlierScore := mx ]
    }
    x[ , Outlier := fifelse(OutlierScore < thrsh, 0, 1, na=NA) ]
    
    ## Restrict max values of thresholds to show on a map
    x[ , OutlierScore := OutlierScore ]
    x[ OutlierScore > thrsh * threshold_multiplier, OutlierScore := thrsh * threshold_multiplier ]

  } else {     # no outlier mode
    x[ , Outlier := 0 ]
  }

  x$Outlier <- factor(x$Outlier, levels = c(0, 1))

  ## Stats
  if(verbose == TRUE & !is.null(outlier_mode)){
    nr <- nrow(x)
    no <- nrow(x[ Outlier == 1 ])
    cat("\nOutlier mode:", outlier_mode, "\n")
    cat("..Number of records:", nr, "\n")
    cat("..Number of outliers:", no, "\n")
    cat("..Percentage of outliers:", round(no / nr * 100, 2), "%\n")
  }

  ## Convert to sf object
  xx <- sf::st_as_sf(x,
    coords = c("Longitude", "Latitude"),
    crs = 4326)

  pp <- ggplot() + 
    geom_sf(data = world)

  ## Show outlier bubbles
  if( !is.null(outlier_mode) ){

    pp <- pp +
      geom_sf(data = xx, aes(size = OutlierScore, color = Outlier), shape = 21) + 
      scale_color_manual(values = c("blue", "red")) +
      geom_sf(data = xx, aes(color = Outlier), size = 1) +
      labs(size = "Outlier score")
  
  } else {

    pp <- pp + 
      geom_sf(data = xx, color = "blue")

  }

  return(pp)
}

## Plot
cat("\n..Plotting\n")
if(is.null(THRESHOLD)){

  p1 <- show_occ(CRD, outlier_mode = "low")
  p2 <- show_occ(CRD, outlier_mode = "medium")
  p3 <- show_occ(CRD, outlier_mode = "high")

  ## Arrange plots in a row
  cat("\n..Arranging plots\n")
  pp <- p1 + p2 + p3 + plot_layout(ncol = 1)

  if(CROP == TRUE){
    pp <- pp + coord_sf(xlim = bbox$x, ylim = bbox$y, expand = TRUE)
  }

  ## Export plots
  cat("..Exporting plots\n")
  ggsave(
    filename = PLOT,
    plot = pp,
    width = 15, height = 22,
    units = "in", dpi = 300)

} else {
  ## Plot with fixed outlier threshold

  p <- show_occ(CRD, outlier_mode = "threshold", threshold = THRESHOLD)

  if(CROP == TRUE){
    p <- p + coord_sf(xlim = bbox$x, ylim = bbox$y, expand = TRUE)
  }

  cat("..Exporting plot\n")
  ggsave(
    filename = PLOT,
    plot = p,
    width = 15, height = 7.5,
    units = "in", dpi = 300)

}

cat("\nDone.\n")

