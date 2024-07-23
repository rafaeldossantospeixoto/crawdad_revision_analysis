library(crawdad)
library(tidyverse)
ncores <- 7


# Extend sim --------------------------------------------------------------

data(sim)
simbr <- sim %>% 
  mutate(celltypes = 'E',
         x = x + 2000,
         y = y) 
simtr <- sim %>% 
  mutate(celltypes = 'E',
         x = x + 2000,
         y = y + 2000) 
simtl <- sim %>% 
  mutate(celltypes = 'E',
         x = x,
         y = y + 2000) 
ext_sim <- bind_rows(sim, simbr, simtr, simtl)

p <- ext_sim %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes), size = 0.1) +
  scale_color_manual(values = c(rainbow(4), '#ffc0cb')) + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/extsim/spat_viz.pdf')
p
dev.off()


# CRAWDAD -----------------------------------------------------------------

## convert to sf
cells <- crawdad:::toSF(pos = ext_sim[,c("x", "y")],
                      celltypes = ext_sim$celltypes)

scales <- seq(100, 1000, by=50)

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
saveRDS(dat_50, 'running_code/processed_data/dat_extsim_50.RDS')

## multiple-test correction
zsig <- correctZBonferroni(dat_50)

## viz 50
vizColocDotplot(dat_50, reorder = F, 
                zSigThresh = zsig, zScoreLimit = zsig*2,
                dotSizes = c(1, 10)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



# Paper figures -----------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_extsim_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zSigThresh = zsig, 
                zScoreLimit = zsig*2, 
                dotSizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = F, 
                     zSigThresh = zsig, zScoreLimit = zsig*2,
                     dotSizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1)) +
  coord_equal()
p
pdf('running_code/paper_figures/extsim/dotplot_crawdad.pdf', height = 5.5, width = 5.5)
p
dev.off()



# Export data -------------------------------------------------------------

write.csv(ext_sim, 'running_code/squidpy/exported_data/extsim.csv')
