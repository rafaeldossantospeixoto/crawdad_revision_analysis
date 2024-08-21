
#' Visualize grids and clusters after shuffling
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' and the shuffled null background.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param shuffledList list; cell type labels shuffled at different scales 
#' (output from makeShuffledCells())
#' @param scale numeric; size of the scale to plot
#' @param permutation numeric; the number of the permutation of interest.
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
vizShuffledGrids <- function(cells, shuffledList, scale, 
                             permutation = 1, 
                             square = TRUE,
                             ofInterest = NULL, pointSize = 1, alpha = 0.5){
  
  ## define colors to be consistent with the vizCluster function
  cluster_cols <- rainbow(n = length(unique(cells$celltypes)))
  names(cluster_cols) <- unique(cells$celltypes)
  cluster_cols['other'] <- '#E6E6E6'
  
  ## assing the shuffled labels to the celltypes column
  shuffled_labels <- shuffledList[[as.character(scale)]][[as.character(permutation)]]
  cells$celltypes <- as.factor(shuffled_labels)
  ## determine the total number of permutations to perform the offsetting
  totalPermutations <- length(shuffledList[[1]])
  
  ## use vizGrids to plot
  vizGrids(cells = cells, scale = scale, 
           permutation = permutation, totalPermutations = totalPermutations,
           square = square,
           ofInterest = ofInterest, pointSize = pointSize, alpha = alpha) +
    ## NA to gray
    ggplot2::scale_color_manual(values = cluster_cols, na.value = "#E6E6E6")
  
}



# Test --------------------------------------------------------------------

library(crawdad)

data(sim)
## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = sim[,c("x", "y")],
                        celltypes = sim$celltypes)
## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(250, 500, by=50),
                                           perms = 3,
                                           ncores = 5,
                                           seed = 1,
                                           verbose = TRUE)

vizGrids(cells, scale =  300, permutation = 2, totalPermutations = 3)
vizShuffledGrids(cells = cells, shuffledList = shuffle.list, scale = 300, 
                  permutation = 2)
