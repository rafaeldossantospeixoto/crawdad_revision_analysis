# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

data(seq)

## convert to sf
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                        celltypes = seq$celltypes)

scales <- seq(100, 1000, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(seq,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(seq,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_seq_50.RDS')


# Paper figures -----------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zsig.thresh = zsig, 
                zscore.limit = zsig*2, 
                dot.sizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = TRUE, zsig.thresh = zsig, 
                     zscore.limit = zsig*2,
                     dot.sizes = c(1, 8)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
pdf('running_code/paper_figures/embryo_dotplot_crawdad.pdf',
    height = 7.5, width = 9)
p
dev.off()
