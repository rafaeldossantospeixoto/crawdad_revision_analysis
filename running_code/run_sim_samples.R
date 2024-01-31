
# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 5

data(sim)

## convert to sf
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)

scales <- seq(100, 1000, by=100)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(sim,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 100
results_100 <- crawdad::findTrends(sim,
                                   dist = 100,
                                   shuffle.list = shuffle.list,
                                   ncores = ncores,
                                   verbose = TRUE,
                                   returnMeans = FALSE)
dat_100 <- crawdad::meltResultsList(results_100, withPerms = T)

## find trends, dist 50
results_50 <- crawdad::findTrends(sim,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)

## multiple-test correction
ntests <- length(unique(dat_100$reference)) * length(unique(dat_100$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

## viz 100
vizColocDotplot(dat_100, reorder = TRUE, zsig.thresh = zsig, zscore.limit = zsig*2) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
vizColocDotplot(dat_50, reorder = TRUE, zsig.thresh = zsig, zscore.limit = zsig*2) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
