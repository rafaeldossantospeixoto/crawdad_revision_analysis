library(crawdad)
library(tidyverse)
library(gridExtra)
library(stringr)

# Get filenames -----------------------------------------------------------

lfs <- list.files('simulating_data/densities/')
lfs <- unlist(str_extract_all(lfs, "\\b\\w+\\.csv\\b"))
lfs

# CRAWDAD -----------------------------------------------------------------

ncores <- 15
scales <- c(25, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)

dats <- list()
for (file_name in lfs) {
  df <- read.csv(paste0('simulating_data/densities/', file_name), 
                 row.names = 1)
  cells <- crawdad::toSF(pos = df[,c("x", "y")], celltypes = df$celltype)
  shuffle.list <- crawdad:::makeShuffledCells(cells,
                                              scales = scales,
                                              perms = 3,
                                              ncores = ncores,
                                              seed = 1,
                                              verbose = TRUE)
  gc()
  results <- crawdad::findTrends(cells,
                                 dist = 50,
                                 shuffle.list = shuffle.list,
                                 ncores = ncores,
                                 verbose = TRUE,
                                 returnMeans = FALSE)
  gc()
  dat <- crawdad::meltResultsList(results, withPerms = TRUE)
  dats[[file_name]] <- dat
}

## Time was 38.81 mins
saveRDS(dats, 'simulating_data/densities/results/dats.RDS')
dats <- readRDS('simulating_data/densities/results/dats.RDS')

# Run Figures -------------------------------------------------------------

dev.off()
for (file_name in lfs) {
  print(file_name)
  df <- read.csv(paste0('simulating_data/densities/', file_name), 
                 row.names = 1)
  p_spat <- df %>% 
    ggplot() +
    geom_point(aes(x, y, color = celltype), size = .01) +
    scale_color_manual(values = rainbow(3)) +
    labs(title = file_name)
  pdf(paste0('simulating_data/densities/results/',
             'spat_',
             str_extract(file_name, '.*(?=\\.csv)'),
             '.pdf'))
  p_spat
  dev.off()
  dat <- dats[[file_name]]
  ## calculate the zscore for the multiple-test correction
  zsig <- correctZBonferroni(dat)
  ## summary visualization
  p <- vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  pdf(paste0('simulating_data/densities/results/',
             'dotplot_',
             str_extract(file_name, '.*(?=\\.csv)'),
             '.pdf'),
      height = 3.5, width = 4)
  p
  dev.off()
  ## trends
  p1 <- dat %>% 
    filter(reference == c('A','B')) %>% 
    filter(neighbor == 'C') %>% 
    vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
  p2 <- dat %>% 
    filter(reference == 'C') %>% 
    filter(neighbor == c('A','B')) %>% 
    vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
  pdf(paste0('simulating_data/densities/results/',
             'trends_',
             str_extract(file_name, '.*(?=\\.csv)'),
             '.pdf'),
      height = 7, width = 8)
  grid.arrange(p1,p2)
  dev.off()
}



