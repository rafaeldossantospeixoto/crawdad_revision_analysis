library(crawdad)
library(tidyverse)
# source('bugs/testing_parallelization.R')

# Find the error ----------------------------------------------------------


# results <- crawdad::findTrends(cells_sim,
#                                dist = 100,
#                                shuffle.list = shuffle.list_sim,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)

## run findTrends
cells = cells_sim
dist = 100
ncores = 1
shuffle.list = shuffle.list_sim
subset.list = NULL
verbose = TRUE
removeDups = TRUE
returnMeans = FALSE
  
sf::st_agr(cells) <- "constant"
  
# if(is.null(subset.list)){
    
celltypes <- factor(cells$celltypes)
    
d <- dist
    
# results.all <- lapply(levels(celltypes), function(ct) {
ct <- levels(celltypes)[1]
ct
      
# get polygon geometries of reference cells of "celltype" up to defined distance "dist"
# use this to assess neighbors within "d" um of each cell
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
# > object.size(ref.buffer)/1024/1024
# 199.6 bytes

## union does fix it, is it better than intersection or join
df_pts <- data.frame(x = c(0, 0, 150, 150, 450), y = c(0, 150, 0, 150, 75))
pts <- crawdad::toSF(pos = df_pts, celltypes = rep('A', nrow(df_pts)))
circles <- sf::st_buffer(pts, 100)$geometry
plot(circles)
# removes intersections and creates multipolygon
plot(sf::st_union(circles))
# keeps intersections and creates mutipolygon
plot(sf::st_combine(circles))
## apply union and intersect
tmp <- sf::st_union(ref.buffer$geometry)
it <- Sys.time()
neigh.cells <- sf::st_intersection(cells, tmp)
Sys.time() - it ## Time difference of 2.117847 secs
## test intersection and then join
plot(rbind(cells[1:10, ], ref.buffer[1:10, ])$geometry)
tmp_int <- sf::st_intersection(cells[1:10, ], ref.buffer[1:10, ]$geometry)
tmp_int
plot(rbind(tmp_int, ref.buffer[1:10, ])$geometry)
## intersection shows the points inside all the polygons, but does not tell 
## which polygon it insersects to
# sf::st_join(cells[1:10, ], ref.buffer[1:10, ]$geometry)
## Intersections that involve points using st_intersection() is the same as 
## using st_join() because points are length-less and area-less (nothing to cut).



# Create new findTrends ---------------------------------------------------

findTrends <- function(cells,
                       dist = 100,
                       ncores = 1,
                       shuffle.list,
                       subset.list = NULL,
                       verbose = TRUE,
                       removeDups = TRUE,
                       returnMeans = TRUE){
  
  if(!is.list(shuffle.list)){
    stop("`shuffle.list` is not a list. You can make this using `makeShuffledCells()`")
  }
  
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  sf::st_agr(cells) <- "constant"
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  if(is.null(subset.list)){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    celltypes <- factor(cells$celltypes)
    
    d <- dist
    results.all <- lapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      # use this to assess neighbors within "d" um of each cell
      ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
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
      return(results)
    }) 
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
  }
  
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
  
}





evaluateSignificance <- function(cells,
                                 randomcellslist,
                                 trueNeighCells,
                                 cellBuffer,
                                 ncores = 1,
                                 removeDups = TRUE,
                                 returnMeans = TRUE){
  
  allcells <- cells
  trueNeighCells <- trueNeighCells
  cellBuffer <- cellBuffer
  
  ## if true, will return a data.frame
  if(returnMeans){
    
    ## for each scale:
    results <- do.call(rbind, BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer) #$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        # if(removeDups){
        #   # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
        #   bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
        #   # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        # }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the mean Z score across permutations for the given scale
      return(colMeans(scores))
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
    
    ## otherwise, returns a list
  } else {
    
    ## for each scale:
    results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer) #$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        # if(removeDups){
        #   # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
        #   bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
        #   # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        # }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the data.frame of Z scores for each permutation (row)
      return(scores)
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
    
    ## will be list of data.frame in this case
    ## each data.frame is a scale
    names(results) <- names(randomcellslist)
    
  }
  
  rm(allcells)
  rm(trueNeighCells)
  rm(cellBuffer)
  gc(verbose = FALSE, reset = TRUE)
  
  return(results)
}


# Test new functions ------------------------------------------------------


## seq ---------------------------------------------------------------------

data('seq')
ggplot(seq, aes(x, -y)) + geom_point(size = .1)

ncores <- 5
cells_seq <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_seq <- crawdad:::makeShuffledCells(cells_seq,
                                                scales = scales,
                                                perms = 3,
                                                ncores = ncores,
                                                seed = 1,
                                                verbose = TRUE)
## Time was 1.52 mins
## large list 42 MB

## old function
results <- crawdad::findTrends(cells_seq,
                               dist = 100,
                               shuffle.list = shuffle.list_seq,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 1.23 mins

## run new functions and try again
results_new <- findTrends(cells_seq,
                               dist = 100,
                               shuffle.list = shuffle.list_seq,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 13.95 mins
## 11.34146 times longer

## check if results match
all.equal(results, results_new)
## TRUE


## slide -------------------------------------------------------------------

data('slide')
ggplot(slide, aes(x, -y)) + geom_point(slide = .1)

ncores <- 5
cells_slide <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_slide <- crawdad:::makeShuffledCells(cells_slide,
                                                scales = scales,
                                                perms = 3,
                                                ncores = ncores,
                                                seed = 1,
                                                verbose = TRUE)
## Time was 0.83 mins
## large list 21.9 MB

## old function
results <- crawdad::findTrends(cells_slide,
                               dist = 100,
                               shuffle.list = shuffle.list_slide,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 0.92 mins

## run new functions and try again
results_new <- findTrends(cells_slide,
                          dist = 100,
                          shuffle.list = shuffle.list_slide,
                          ncores = ncores,
                          verbose = TRUE,
                          returnMeans = FALSE)
## Time was 3.88 mins
## 4.217391 longer

## check if results match
all.equal(results, results_new)
## TRUE


## pkhl --------------------------------------------------------------------

data('pkhl')
ggplot(pkhl, aes(x, -y)) + geom_point(size = .1)

ncores <- 5
cells_pkhl <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_pkhl <- crawdad:::makeShuffledCells(cells_pkhl,
                                                  scales = scales,
                                                  perms = 3,
                                                  ncores = ncores,
                                                  seed = 1,
                                                  verbose = TRUE)
## Time was 1.1 mins
## large list 21.9 MB

## old function
# results <- crawdad::findTrends(cells_pkhl,
#                                dist = 100,
#                                shuffle.list = shuffle.list_pkhl,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
## Not working

## run new functions and try again
results_new <- findTrends(cells_pkhl,
                          dist = 100,
                          shuffle.list = shuffle.list_pkhl,
                          ncores = ncores,
                          verbose = TRUE,
                          returnMeans = FALSE)
## Time was 5.05 mins

## fewer cores
results_new <- findTrends(cells_pkhl,
                          dist = 100,
                          shuffle.list = shuffle.list_pkhl,
                          ncores = 2,
                          verbose = TRUE,
                          returnMeans = FALSE)
## Time was 6.06 mins

