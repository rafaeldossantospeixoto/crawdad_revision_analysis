
# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

dfs <- readRDS('simulating_data/null_sim/cells_nullsim.RDS')

vizClusters(dfs[[1]])
