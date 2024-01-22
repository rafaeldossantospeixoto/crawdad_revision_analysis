library(crawdad)
library()

data('pkhl')
## [1] 154446      3
max(pkhl$x) - min(pkhl$x) # 249.2555
max(pkhl$y) - min(pkhl$y) # 240.3526
# > object.size(pkhl)/1024/1024
# 3.5 bytes

ncores <- 5
conv_microns <- 377.4038462 / 10 ## wrong
cells_pkhl <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
scales <- c(25, 50, 75, 100, 150, 200, 250) 
shuffle.list_pkhl <- crawdad:::makeShuffledCells(cells_pkhl,
                                                 scales = scales,
                                                 perms = 3,
                                                 ncores = ncores,
                                                 seed = 1,
                                                 verbose = TRUE)
## Time was 0.98 mins
## large list 233.5 MB

results <- crawdad::findTrends(cells_pkhl,
                               dist = 100/conv_microns,
                               shuffle.list = shuffle.list_pkhl,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Stop in B cells and runs out of memory, even on Easley
## B cells is not the largest in number of cells, 3322, it should not take
## this long