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


