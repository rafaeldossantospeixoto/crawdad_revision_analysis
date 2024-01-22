library(crawdad)
library(tidyverse)


# Get filenames -----------------------------------------------------------

lfs <- list.files('simulating_data/symmetries/')
lfs <- unlist(str_extract_all(lfs, "\\b\\w+\\.csv\\b"))


# CRAWDAD -----------------------------------------------------------------

ncores <- 7
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

dats <- list()
for (file_name in lfs) {
  df <- read.csv(paste0('simulating_data/symmetries/', file_name), 
                 row.names = 1)
  cells <- crawdad::toSF(pos = df[,c("x", "y")], celltypes = df$celltype)
  shuffle.list <- crawdad:::makeShuffledCells(cells,
                                              scales = scales,
                                              perms = 3,
                                              ncores = ncores,
                                              seed = 1,
                                              verbose = TRUE)
  results <- crawdad::findTrends(cells,
                                 dist = 50,
                                 shuffle.list = shuffle.list,
                                 ncores = ncores,
                                 verbose = TRUE,
                                 returnMeans = FALSE)
  dat <- crawdad::meltResultsList(results, withPerms = TRUE)
  dats[[file_name]] <- dat
}

## Time was 38.81 mins
saveRDS(dats, 'simulating_data/symmetries/results/dats.RDS')
dats <- readRDS('simulating_data/symmetries/results/dats.RDS')


# Run CRAWDAD -------------------------------------------------------------

file_name <- lfs[3]
df <- read.csv(paste0('simulating_data/symmetries/', file_name), 
               row.names = 1)
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype)) +
  scale_color_manual(values = rainbow(3))
dat <- dats[[file_name]]
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))
