## This file is for running CRAWDAD on MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(here)

# Run method --------------------------------------------------------------

## based on the tutorial https://jef.works/CRAWDAD/

slices <- seq(1,3)
replicates <- seq(1,3)

## define the scales to analyze the data
scales <- seq(100, 1000, by=100)

## define the number of cores for parallelization
ncores <- parallel::detectCores() - 2

for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Running CRAWDAD on ", dataset_name))
    
    ## load dataset
    spe <- readRDS(file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))
    
    ## convert dataframe to sf data.frame
    cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)
    
    ## visualize
    # crawdad::vizEachCluster(cells = cells,
    #                         coms = as.factor(cells$celltypes),
    #                         s = 2)
    
    ## shuffle cells to create null background
    shuffle.list <- crawdad:::makeShuffledCells(cells,
                                                scales = scales,
                                                perms = 3,
                                                ncores = ncores,
                                                seed = 1,
                                                verbose = TRUE)
    saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
    ## calculate the zscore for the cell-type pairs at different scales
    ## error: Error in FUN(X[[i]], ...) : object 'neigh.cells' not found
    results <- crawdad::findTrends(cells = cells,
                                   dist = 100,
                                   shuffle.list = shuffle.list,
                                   ncores = ncores,
                                   verbose = TRUE,
                                   returnMeans = FALSE)
    saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))
  }
}
