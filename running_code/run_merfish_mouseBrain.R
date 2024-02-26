## This file is for running CRAWDAD on MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(here)
library(ggplot2)

# Run method --------------------------------------------------------------

## based on the tutorial https://jef.works/CRAWDAD/

slices <- seq(1,3)
replicates <- seq(1,3)

## define the scales to analyze the data
scales <- seq(100, 1000, by=100)

## define the number of cores for parallelization
ncores <- parallel::detectCores() - 2

## whether to use original or cleaned celltypes ("original" or "cleaned")
ct_type <- "cleaned"
ct_type

for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Running CRAWDAD on ", dataset_name))
    
    ## load dataset
    spe <- readRDS(file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))
    
    ## convert dataframe to sf data.frame
    if (ct_type == "original") {
      # use original cell types
      cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)
    } else if (ct_type == "cleaned") {
      # use cleaned cell types
      cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype_cleaned)
    }
    
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
    if (ct_type == "original") {
      saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
    } else if (ct_type == "cleaned") {
      saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells_ct_cleaned.RDS")))
    }
    ## calculate the zscore for the cell-type pairs at different scales
    ## error: Error in FUN(X[[i]], ...) : object 'neigh.cells' not found
    results <- crawdad::findTrends(cells = cells,
                                   dist = 100,
                                   shuffle.list = shuffle.list,
                                   ncores = ncores,
                                   verbose = TRUE,
                                   returnMeans = FALSE)
    if (ct_type == "original") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))
    } else if (ct_type == "cleaned") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    }
  }
}


# Plot --------------------------------------------------------------------

slices <- seq(1,3)
replicates <- seq(1,3)

for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Plotting CRAWDAD results for ", dataset_name))
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))
    
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    ## calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    ## summary visualization
    vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, dot.sizes = c(2, 8)) +
      theme(axis.text.x = element_text(angle = 35, h = 0)) +
      ggtitle(dataset_name)
    ggsave(filename = here("running_code", "plots", paste0(dataset_name, "_summary.png")), width = 12, height = 8, dpi = 300)
  }
}
