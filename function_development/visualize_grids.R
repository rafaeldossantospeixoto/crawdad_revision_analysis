
# Create function ---------------------------------------------------------

#' Visualize grids and clusters
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' and clusters.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param gridSize numeric; size of the grid to plot
#' @param square if false, make hexagonal grid (default TRUE)
#' @param ofInterest character vector; a vector of specific clusters to visualize
#' @param pointSize numeric; size of points
#' @param alpha numeric; transparency of points
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
#' vizAllClusters(cells)
#' }
#' 
#' @export
vizGrids <- function(cells, scale, shuffledList = NULL,
                     totalPermutations = NULL, permutation = 1, 
                     square = square,
                     ofInterest = NULL, pointSize = 1, alpha = 0.5){

  ## the total number of permutations is used to calculate the offset, so if the
  ## user asks for a different permutation, but does not provide it, the code
  ## will not work
  if ((permutation != 1) & (is.null(totalPermutations))){
    stop('provide the total number of permutations')
  }
  
  ## create grids
  if (permutation == 1) {
    ## for no permutations 
    ## create grid with no offset
    grid <- sf::st_make_grid(cells, cellsize = scale, square = square)
  } else {
    ## for permutations 
    ## calculate offset
    offset <- -seq(from = 0, to = scale, by = scale/totalPermutations)[permutation]
    ## get bounding box with min and max coordinates
    bbox <- sf::st_bbox(cells$geometry)
    ## create grid with offset
    grid <- sf::st_make_grid(cells, cellsize = scale, 
                             offset = c(bbox[['xmin']] + offset, 
                                        bbox[['ymin']] + offset),
                             square = square)
  }
  
  
  ## get the coordinates of the centers of the grid tiles to add the tile IDs
  grid_coords_centroids <- as.data.frame(sf::st_coordinates(sf::st_centroid(grid)))
  grid_coords_centroids$name <- as.character(rownames(grid_coords_centroids))
  ## create plot
  crawdad::vizClusters(cells = cells, ofInterest, pointSize, alpha) + 
    ## add in the grid information on top of the plot
    ggplot2::geom_sf(data = grid, fill = NA) +
    ggplot2::geom_text(data = grid_coords_centroids, 
                       ggplot2::aes(X, Y, label = name)) +
    ggplot2::guides(color = "none") +
    ggplot2::theme(axis.line = ggplot2::element_blank(),       # Remove axis lines
          axis.text = ggplot2::element_blank(),       # Remove axis text
          axis.ticks = ggplot2::element_blank(),      # Remove axis ticks
          axis.title = ggplot2::element_blank(),      # Remove axis titles
          panel.background = ggplot2::element_blank(), # Remove panel background,
          axis.text.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.title.y = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(color = "black", linewidth = 1, fill = NA))
}

library(crawdad)
## load the spleen data of the pkhl sample 
data('pkhl')
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
## plot
vizGrids(cells, scale = 500, pointSize = .1)
vizClusters(cells, pointSize = .1)



# Apply -------------------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

data(sim)
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)
dat <- readRDS('running_code/processed_data/dat_sim_50.RDS')

vizGrids(sim, 300)
