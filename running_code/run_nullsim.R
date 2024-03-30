
# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

cells <- readRDS('simulating_data/null_sim/cells_nullsim_s0.RDS')

vizClusters(cells)

## reduced to 500
scales <- seq(100, 500, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_nullsim_50.RDS')

## multiple-test correction
zsig <- correctZBonferroni(dat_50)

## viz 50
vizColocDotplot(dat_50, zscore.limit = zsig*2)
vizColocDotplot(dat_50, zsig.thresh = zsig, zscore.limit = zsig*2)



# Paper figures -----------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_nullsim_50.RDS')

## error because no dot with significance
zsig <- correctZBonferroni(dat_50)
# vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
#                 zscoreLimit = zsig*2, 
#                 dotSizes = c(2, 14)) +
#   theme(legend.position='right',
#         axis.text.x = element_text(angle = 45, h = 0))



## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = F, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1))
p
pdf('running_code/paper_figures/nullsim_dotplot_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()



## Trends ------------------------------------------------------------------

## will plot using same design for squidpy and ripleys
ref_cts <- c(rep('A', 4), LETTERS[2:4])
neigh_cts <- c(LETTERS[1:4], rep('A', 3))
for (i in 1:length(ref_cts)) {
  p <- dat_50 %>% 
    filter(reference == ref_cts[i]) %>% 
    filter(neighbor == neigh_cts[i]) %>% 
    vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
  print(p)
  # pdf(sprintf('running_code/paper_figures/nullsim/trend_ref%s_neigh%s_crawdad.pdf',
  #             ref_cts[i], neigh_cts[i]),
  #     height = 4, width = 5)
  # print(p)
  # dev.off()
}