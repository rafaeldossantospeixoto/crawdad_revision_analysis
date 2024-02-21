## This file is for running CRAWDAD on 10X Genomics Xenium human breast cancer dataset

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(here)

dataset_name <- "xenium_humanBreastCancer"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "running_code/processed_data/xenium_humanBreastCancer_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()

ct_labels <- colData(spe)$celltype

# Run method --------------------------------------------------------------

## based on the tutorial https://jef.works/CRAWDAD/

## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)

## visualize
# crawdad::vizEachCluster(cells = cells,
#                         coms = as.factor(cells$celltypes),
#                         s = 2)

## define the scales to analyze the data
scales <- seq(100, 1000, by=100)
## shuffle cells to create null background
ncores <- parallel::detectCores() - 2
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


# Plot --------------------------------------------------------------------

findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))

dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, dot.sizes = c(3, 12)) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
ggsave(filename = here("running_code", "plots", paste0(dataset_name, "_summary.pdf")), dpi = 300)
