
# Create function ---------------------------------------------------------

#' Visualize grids and clusters
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' used to create the null background.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param scale numeric; size of the scale to plot
#' @param permutation numeric; the number of the permutation of interest.
#' @param totalPermutations numeric; the total number of permutations used to 
#' shuffle the data.
#' @param square boolean; if true, create a squared grid, if false, make 
#' hexagonal grid (default TRUE)
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
vizGrids <- function(cells, scale,
                     permutation = 1, totalPermutations = NULL, 
                     square = TRUE,
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
                       ggplot2::aes(X, Y, label = name))
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
# library(tidyverse)
ncores <- 7

data(sim)
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)
dat <- readRDS('running_code/processed_data/dat_sim_50.RDS')

vizGrids(sim, scale =  300)
vizGrids(sim, scale =  300, permutation = 3, totalPermutations = 10)
