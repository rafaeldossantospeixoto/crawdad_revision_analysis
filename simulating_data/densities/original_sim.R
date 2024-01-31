library(crawdad)
library(tidyverse)

ncores <- 5
data(sim)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = sim[,c("x", "y")],
                        celltypes = sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 3,
                                           ncores = ncores,
                                           seed = 1,
                                           verbose = TRUE)

## find trends, passing background as parameter
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE, 
                               returnMeans = FALSE)
## convert results to data.frame
dat <- crawdad::meltResultsList(results, withPerms = T)

## multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

vizColocDotplot(dat, reorder = F, zsig.thresh = zsig, zscore.limit = zsig*2) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        range = c(1, 41))


# Up sampling ----------------------------------------------------------

vizClusters(cells)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup')

up_sim <- rbind(sim, bc_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = up_sim[,c("x", "y")],
                        celltypes = up_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 3,
                                           ncores = ncores,
                                           seed = 1,
                                           verbose = TRUE)

## find trends, passing background as parameter
results <- crawdad::findTrends(cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE, 
                               returnMeans = FALSE)
## convert results to data.frame
dat <- crawdad::meltResultsList(results, withPerms = T)

## multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

vizColocDotplot(dat, reorder = F, zsig.thresh = zsig, zscore.limit = zsig*2) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        range = c(1, 41))