library(crawdad)
library(tidyverse)
ncores <- 28


# pkhl --------------------------------------------------------------------

data('pkhl')
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")],
                       celltypes = pkhl$celltypes)

ggplot(pkhl, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(pkhl$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 25.74 mins in my pc
## Time was 70.47 mins in Easley with 28 cores
saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_pkhl_50.RDS')


## find trends, passing background as parameter
## changed distance to 50
results <- crawdad::findTrends(cells,
                               dist = 50,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 107.72 mins in my pc
## Time was 104.3 mins in Easley with 28 cores
## The number of cores does not seem to help much

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
saveRDS(dat, 'running_code/processed_data/spleen/dat_pkhl_50.RDS')




## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.68mins
# head(binomMat)
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_pkhl_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_pkhl_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_pkhl_50.RDS')
## Time to compute was 0.17mins
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_pkhl_50.RDS')

results.subsets <- crawdad::findTrends(cells,
                                       dist = 50,
                                       shuffle.list = shuffle.list,
                                       subset.list = subset.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE)
## 8.0865 hours to run
## subsets
dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
saveRDS(dats, file = 'running_code/processed_data/spleen/dats_pkhl_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_pkhl_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zsig.thresh = zsig, zscore.limit = zsig*2, 
                     reorder = TRUE, dot.sizes = c(2, 14)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  # scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p



# xxcd --------------------------------------------------------------------

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
head(xxcd)

cells <- crawdad::toSF(pos = xxcd[,c("x", "y")],
                       celltypes = xxcd$celltypes)

ggplot(xxcd, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(xxcd$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 24.92 mins

## changed distance to 50
results <- crawdad::findTrends(cells,
                               dist = 50,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 607.45 mins

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
saveRDS(dat, 'running_code/processed_data/dat_xxcd_50.RDS')