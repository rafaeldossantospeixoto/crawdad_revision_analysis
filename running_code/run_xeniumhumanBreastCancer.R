## This file is for running CRAWDAD on 10X Genomics Xenium human breast cancer dataset


# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)

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
## define the scales to analyze the data
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
## shuffle cells to create null background
ncores <- parallel::detectCores() - 2
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
## error: Error in FUN(X[[i]], ...) : object 'neigh.cells' not found
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = 7,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))


# Bug fix -----------------------------------------------------------------

cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)
ct <- "B_Cells"
d <- 100
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
