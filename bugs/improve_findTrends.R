
# CRAWDAD -----------------------------------------------------------------

library(crawdad)
library(tidyverse)
## load the spleen data of the pkhl sample 
data('pkhl')
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## shuffle cells to create null background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = 7,
                                            seed = 1,
                                            verbose = TRUE)


# Testing findTrends --------------------------------------------------------------

new_dist <- 100 # * (377.4038462 / 1000)
ncores <- 7


## findTrends --------------------------------------------------------------


### ngbh 100 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)
    
d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {
      
# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 24.11239 secs



### ngbh 37 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = cells
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 24.76162 secs
## Time difference of 14.42391 mins

      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      # if(removeDups){
      #   # message("number of neighbor cells before: ", nrow(neigh.cells))
      #   neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
      #   # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      # }


### call evaluateSignificance -----------------------------------------------
      
## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
## I think a bottle neck originally was waiting for certain cell types to finish
## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
## I could also split up the permutations, but then each cell type for each scale is done one by one
## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
results <- evaluateSignificance(cells = cells,
                                randomcellslist = shuffle.list,
                                trueNeighCells = neigh.cells,
                                cellBuffer = ref.buffer_union,
                                ncores = ncores,
                                removeDups = removeDups,
                                returnMeans = returnMeans)
# return(results)
# names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
  if(!is.null(subset.list)){
    
    ## load in the subset file if it exists, or make it and probably be a good
    ## idea to save it, too
    if(!is.list(subset.list)){
      stop(paste0("`subset.list` is not a list. You can build this using: `binomialTestMatrix()` then `selectSubsets()`"))
    }
    
    combo_ids <- names(subset.list)
    d <- dist
    
    if(verbose){
      message("Calculating trends for each subset in `subset.list` with respect to the cell types in `cells$celltypes`")
    }
    
    ## initialize list
    results.all <- list()
    
    ## for each subset of cells, evaluate significance against each reference cell type across scales
    for(i in combo_ids){
      
      if(verbose){
        message(i)
      }
      
      ## get area around the subset cells to identify neighbors
      ref.buffer <- sf::st_buffer(cells[subset.list[[i]], ], d)
      ## union polygons to avoid memory overflow
      ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer_union)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      # if(removeDups){
      #   # message("number of neighbor cells before: ", nrow(neigh.cells))
      #   neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
      #   # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      # }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer_union,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans)
      results.all[[i]] <- results
      
      rm(ref.buffer)
      rm(neigh.cells)
      rm(results)
      gc(verbose = FALSE, reset = TRUE)
    }
    
  }
  
  ## return results
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  return(results.all)





## evaluateSignificance ----------------------------------------------------

# evaluateSignificance <- function(cells,
#                                  randomcellslist,
#                                  trueNeighCells,
#                                  cellBuffer,
#                                  ncores = 1,
#                                  removeDups = TRUE,
#                                  returnMeans = TRUE){

cells = cells
randomcellslist = shuffle.list
trueNeighCells = neigh.cells
cellBuffer = ref.buffer_union
ncores = ncores
removeDups = removeDups
returnMeans = returnMeans
  
allcells <- cells
trueNeighCells <- trueNeighCells
cellBuffer <- cellBuffer
  
## for each scale:
# results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
cellsAtRes <- randomcellslist[[1]]
      
## Iterate through each permutation of a given scale and 
## produce the scores for each neighbor cell type. 
## Scores for a given permutation added as a row to a dataframe.
## Take the column mean of the scores for each neighbor cell type across permutations.
## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
# scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
randomcellslabels <- cellsAtRes[[1]]

randomcells <- allcells
randomcells$celltypes <- as.factor(randomcellslabels)
sf::st_agr(randomcells) <- "constant"

it <- Sys.time()
bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer) #$geometry)
Sys.time() - it
## Time difference of 14.40688 mins
        
## remove duplicate neighbor cells to prevent them from being counted multiple times
## and inflate the Z scores
# if(removeDups){
#   # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
#   bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
#   # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
# }
        
## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
it <- Sys.time()
y1 <- table(trueNeighCells$celltypes)
y2 <- table(bufferrandomcells$celltypes)
n1 <- length(trueNeighCells$celltypes)
n2 <- length(bufferrandomcells$celltypes)
p1 <- y1/n1
p2 <- y2/n2
p <- (y1+y2)/(n1+n2)
Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
Sys.time() - it
        
# return(Z)
      
# ## returning the data.frame of Z scores for each permutation (row)
# return(scores)
    
# ## will be list of data.frame in this case
# ## each data.frame is a scale
# names(results) <- names(randomcellslist)
  
# return(results)


# Testing findTrends original -----------------------------------------------------

## ngbh 37 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

## calculate the zscore for the cell-type pairs at different scales
results <- findTrends(cells,
                               dist = new_dist,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 15.52 mins
## with $geometry
results_geometry <- findTrends(cells,
                               dist = new_dist,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 9.14 mins

## ngbh 1000 ----------------------------------------------------------------

## test if there will be memory overflow when there are many cells inside the 
## neighborhood
new_dist <- 3000
ncores <- 7

## calculate the zscore for the cell-type pairs at different scales
results <- findTrends(cells,
                      dist = new_dist,
                      shuffle.list = shuffle.list,
                      ncores = 7,
                      verbose = TRUE,
                      returnMeans = FALSE)
## Time was 68.31 mins

# st_join -----------------------------------------------------------------

## findTrends --------------------------------------------------------------

### ngbh 100 ----------------------------------------------------------------

new_dist <- 100 #* (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
# it <- Sys.time()
# ref.buffer_union <- sf::st_union(ref.buffer)
# Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer)
Sys.time() - it
## Time difference of 14.75686 secs



### ngbh 37 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = cells
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
# it <- Sys.time()
# ref.buffer_union <- sf::st_union(ref.buffer)
# Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer)
Sys.time() - it
## Time difference of 2.124815 secs






# Simplify ----------------------------------------------------------------

## it is not great, but it is better, need to compare the results
## st_simplify(original_polygons, dTolerance = 0.01)

## findTrends --------------------------------------------------------------


### ngbh 37 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
plot(ref.buffer_union)
ref.buffer_union <- st_simplify(ref.buffer_union, dTolerance = 10)
plot(ref.buffer_union)
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 34.18229 secs






# Convert sfc to polygon --------------------------------------------------
## nope

new_dist <- 100 # * (377.4038462 / 1000)
ncores <- 7


## findTrends --------------------------------------------------------------


### ngbh 100 ----------------------------------------------------------------

new_dist <- 100 # * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
class(ref.buffer_union)
typeof(ref.buffer_union)
object.size(ref.buffer_union)/1024/1024 
## 0.1 bytes
## cast to polygon
ref.buffer_union <- st_cast(ref.buffer_union, "POLYGON")
class(ref.buffer_union)
typeof(ref.buffer_union)
object.size(ref.buffer_union)/1024/1024 
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 24.11239 secs




# Precision findTrends --------------------------------------------------------------
## nope

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7


## findTrends --------------------------------------------------------------

### ngbh 100 ----------------------------------------------------------------

new_dist <- 100 # * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
sf::st_precision(cells) <- .001
st_precision(ref.buffer) <- .001
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_precision_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 0.6888142 secs

## compare intersections
dim(neigh.cells_precision_100)
## [1] 154446      2
dim(neigh.cells_100)
## [1] 148331      2


### ngbh 37 ----------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]
# results.all <- lapply(levels(celltypes), function(ct) {

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
sf::st_precision(cells) <- .001
st_precision(ref.buffer) <- .001
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
# get the different types of neighbor cells that are within "d" of the ref cells
it <- Sys.time()
neigh.cells_precision_100 <- sf::st_intersection(cells, ref.buffer_union)
Sys.time() - it
## Time difference of 0.6888142 secs

## compare intersections
dim(neigh.cells_precision_100)
## [1] 154446      2
dim(neigh.cells_100)
## [1] 148331      2


## Testing precision -------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
ncores <- 7

cells = crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dist = new_dist
shuffle.list = shuffle.list
ncores = ncores
verbose = TRUE
returnMeans = FALSE

## define the atributes of the geometry
## https://r-spatial.github.io/sf/reference/sf.html
sf::st_agr(cells) <- "constant"

celltypes <- factor(cells$celltypes)

d <- dist

ct <- levels(celltypes)[1]

# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
## set precision
sf::st_precision(cells) <- 10
st_precision(ref.buffer) <- 10
## union polygons to avoid memory overflow
it <- Sys.time()
ref.buffer_union <- sf::st_union(ref.buffer)
Sys.time() - it
plot(ref.buffer_union)



