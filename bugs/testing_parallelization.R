library(crawdad)
library(tidyverse)


# seqFISH ----------------------------------------------------------------

data('seq')
## [1] 19416     3
max(seq$x) - min(seq$x) # 5069.662
max(seq$y) - min(seq$y) # 6984.168
## it should be aprox 1000 microns, not 5000 or 6000
ggplot(seq, aes(x, -y)) + geom_point(size = .1)
# > object.size(seq)/1024/1024
# 0.5 bytes

ncores <- 5
cells_seq <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_seq <- crawdad:::makeShuffledCells(cells_seq,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 1.5 mins
## large list 42 MB
results <- crawdad::findTrends(cells_seq,
                               dist = 100,
                               shuffle.list = shuffle.list_seq,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 1.25 mins
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
zsig <- correctZBonferroni(dat)
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))


# slideSeq ----------------------------------------------------------------

data('slide')
max(slide$x) - min(slide$x) # 4628.55
max(slide$y) - min(slide$y) # 3800.465
## it should be aprox 2904 and 2425 microns, not 4628.55 or 3800.465
ggplot(slide, aes(x, -y)) + geom_point(size = .1)
# > object.size(slide)/1024/1024
# 0.9 bytes


# spleen ----------------------------------------------------------------

data('pkhl')
## [1] 154446      3
max(pkhl$x) - min(pkhl$x) # 249.2555
max(pkhl$y) - min(pkhl$y) # 240.3526
# > object.size(pkhl)/1024/1024
# 3.5 bytes

ncores <- 5
cells_pkhl <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_pkhl <- crawdad:::makeShuffledCells(cells_pkhl,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 1.11 mins
## large list 333.6 MB
# results <- crawdad::findTrends(cells,
#                                dist = 100,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
## Stop in B cells and runs out of memory, even on Easley
## B cells is not the largest in number of cells, 3322, it should not take
## this long


# New simulation ----------------------------------------------------------

## create a new simulation with as many cells as the pkhl dataset to see if the
## memory overflow is due to the number of cells/larger shuffle.list



bkg <- simulate_background(size = 154446, cts = c("A"), prob = c(1), seed = 1, 
                            scale = 250)
df <- simulate_circles(bkg, 
                 locs = list(250*c(.25,.25), 250*c(.25,.75), 250*c(.75,.25)), 
                 radii = list(list(inner = 10, outer = 50),
                              list(inner = 15, outer = 50),
                              list(inner = 25, outer = 70)), 
                 cts = list(list(inner = 'B', outer = 'C'),
                            list(inner = 'B', outer = 'C'),
                            list(inner = 'C', outer = 'B')), 
                 probs = list(list(inner = 0.5, outer = 0.5),
                              list(inner = 0.5, outer = 0.5),
                              list(inner = 0.5, outer = 0.5)))
ggplot(df, aes(x, y, color = type)) + geom_point(size = .1)

ncores <- 5
cells_sim <- crawdad::toSF(pos = df[,c("x", "y")], celltypes = df$type)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_sim <- crawdad:::makeShuffledCells(cells_sim,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 1.17 mins
## large list 333.6 MB same as pkhl
results <- crawdad::findTrends(cells_sim,
                               dist = 100,
                               shuffle.list = shuffle.list_sim,
                               ncores = ncores, ## number of cores doesn't matter
                               verbose = TRUE,
                               returnMeans = FALSE)

