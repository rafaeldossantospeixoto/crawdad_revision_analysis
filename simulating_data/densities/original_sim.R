
library(crawdad)
library(tidyverse)

ncores <- 5


# Original ----------------------------------------------------------------


data(sim)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = sim[,c("x", "y")],
                        celltypes = sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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

data(sim)

## duplicate (2x)
vizClusters(cells)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup2')

up_sim <- rbind(sim, bc_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = up_sim[,c("x", "y")],
                        celltypes = up_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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


# Up 3 sampling ----------------------------------------------------------

data(sim)

## duplicate (2x)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup2')
table(bc_df$celltypes)

## add 1 dup
bc_df2 <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df2) <- paste0(rownames(bc_df2), 'dup3')
table(bc_df2$celltypes)

up_sim <- rbind(sim, bc_df, bc_df2)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = up_sim[,c("x", "y")],
                        celltypes = up_sim$celltypes)
vizClusters(cells)
hist(sim$x, breaks = seq(0, 2000, by = 25))
hist(up_sim$x, breaks = seq(0, 2000, by = 25))

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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


# Up 4 sampling ----------------------------------------------------------

vizClusters(cells)
data(sim)
## duplicate (2x)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup2')
up_sim <- rbind(sim, bc_df)

## duplicate again (2x2x)
bc_df <- up_sim[which(up_sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup4')
up_sim <- rbind(up_sim, bc_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = up_sim[,c("x", "y")],
                        celltypes = up_sim$celltypes)
vizClusters(cells)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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


# Up 8 sampling ----------------------------------------------------------

vizClusters(cells)

## duplicate (2x)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup2')
up_sim <- rbind(sim, bc_df)

## duplicate again (2x2x)
bc_df <- up_sim[which(up_sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup4')
up_sim <- rbind(up_sim, bc_df)

## duplicate again (2x2x2x)
bc_df <- up_sim[which(up_sim$celltypes %in% c('B', 'C')), ]
rownames(bc_df) <- paste0(rownames(bc_df), 'dup8')
up_sim <- rbind(up_sim, bc_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = up_sim[,c("x", "y")],
                        celltypes = up_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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


# Down sampling others ---------------------------------------------------------

vizClusters(cells)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
a_df <- sim[which(sim$celltypes %in% c('A')), ]
d_df <- sim[which(sim$celltypes %in% c('D')), ]

a_df <- a_df %>% 
  slice_sample(prop = .5)
d_df <- d_df %>% 
  slice_sample(prop = .5)

down_sim <- rbind(bc_df, a_df, d_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = down_sim[,c("x", "y")],
                        celltypes = down_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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

# Down 1/4 sampling others ---------------------------------------------------------

vizClusters(cells)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
a_df <- sim[which(sim$celltypes %in% c('A')), ]
d_df <- sim[which(sim$celltypes %in% c('D')), ]

a_df <- a_df %>% 
  slice_sample(prop = .25)
d_df <- d_df %>% 
  slice_sample(prop = .25)

down_sim <- rbind(bc_df, a_df, d_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = down_sim[,c("x", "y")],
                        celltypes = down_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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


# Down 1/8 sampling others ---------------------------------------------------------

vizClusters(cells)
bc_df <- sim[which(sim$celltypes %in% c('B', 'C')), ]
a_df <- sim[which(sim$celltypes %in% c('A')), ]
d_df <- sim[which(sim$celltypes %in% c('D')), ]

a_df <- a_df %>% 
  slice_sample(prop = .125)
d_df <- d_df %>% 
  slice_sample(prop = .125)

down_sim <- rbind(bc_df, a_df, d_df)

## convert to sp::SpatialPointsDataFrame
cells <- crawdad:::toSF(pos = down_sim[,c("x", "y")],
                        celltypes = down_sim$celltypes)

## generate background
shuffle.list <- crawdad::makeShuffledCells(cells,
                                           scales = seq(100, 1000, by=50),
                                           perms = 5,
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