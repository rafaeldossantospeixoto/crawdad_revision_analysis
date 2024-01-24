
#' Visualize grids and clusters
#' 
#' @description Uses the cells sf object and size of grid to visualize the grids 
#' and clusters.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param shuffledList list; cell type labels shuffled at different scales 
#' (output from makeShuffledCells())
#' @param gridSize numeric; size of the grid to plot
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
vizShuffledLabels <- function(cells, shuffledList, gridSize, permutation = 1, 
                              ofInterest = NULL,
                              pointSize = 1, alpha = 0.5){
  
  shuffled_labels <- shuffledList[[as.character(gridSize)]][[as.character(permutation)]]
  cells$celltypes <- as.factor(shuffled_labels)
  
  vizGrids(cells, gridSize, ofInterest, pointSize, alpha)
  
}



library(crawdad)

data(sim)
## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = sim[,c("x", "y")],
                        celltypes = sim$celltypes)
## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 3,
                                           ncores = 5,
                                           seed = 1,
                                           verbose = TRUE)

vizShuffledLabels(cells, shuffle.list, 300, 1)
