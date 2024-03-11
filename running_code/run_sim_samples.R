
# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

data(sim)

## convert to sf
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)

scales <- seq(100, 1000, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(sim,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

# ## find trends, dist 100
# results_100 <- crawdad::findTrends(sim,
#                                    dist = 100,
#                                    shuffle.list = shuffle.list,
#                                    ncores = ncores,
#                                    verbose = TRUE,
#                                    returnMeans = FALSE)
# dat_100 <- crawdad::meltResultsList(results_100, withPerms = T)

## find trends, dist 50
results_50 <- crawdad::findTrends(sim,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_sim_50.RDS')

## multiple-test correction
ntests <- length(unique(dat_50$reference)) * length(unique(dat_50$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

## viz 50
vizColocDotplot(dat_50, reorder = F, 
                zsig.thresh = zsig, zscore.limit = zsig*2,
                dot.sizes = c(10, 40)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



# Paper figures -----------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zsig.thresh = zsig, 
                zscore.limit = zsig*2, 
                dot.sizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = F, 
                     zsig.thresh = zsig, zscore.limit = zsig*2,
                     dot.sizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1))
p
pdf('running_code/paper_figures/sim_dotplot_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()



## Trends ------------------------------------------------------------------

p <- dat_50 %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'B') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/sim_trend_refC_neighB_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()

p <- dat_50 %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'B') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/sim_trend_refA_neighB_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()

p <- dat_50 %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/sim_trend_refB_neighC_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()

p <- dat_50 %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'A') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/sim_trend_refB_neighA_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()


# Compare -----------------------------------------------------------------

dat1 <- dat_100
dat2 <- dat_50

vizDiffZscores(dat1, dat2, scale.thresh = 200, reorder = F, dot.size = 50)
vizDiffZscores(dat2, dat1, scale.thresh = 200, reorder = F, dot.size = 50)

vizDiffScales(dat1, dat2, zscore.thresh = 1.96, reorder = F, dot.size = 50)
vizDiffScales(dat2, dat1, zscore.thresh = 1.96, reorder = F, dot.size = 50)
