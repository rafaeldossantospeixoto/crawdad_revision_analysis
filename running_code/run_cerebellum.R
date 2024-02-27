# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

data(slide)
## it should be aprox 2904 and 2425 microns, not 4628.55 and 3800.465
max(slide$x) - min(slide$x)
max(slide$y) - min(slide$y)

## convert to sf
slide <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)

scales <- seq(100, 1000, by=100)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(slide,
                                            scales = scales,
                                            perms = 5,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(slide,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_slide_50.RDS')
