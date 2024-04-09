# Run sim samples ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

data(slide)
## it should be aprox 2904 and 2425 microns, not 4628.55 and 3800.465
max(slide$x) - min(slide$x)
max(slide$y) - min(slide$y)

## convert to sf
slide <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)

scales <- slide(100, 1000, by=100)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(slide,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(slide,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_slide_50.RDS')

dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zsig.thresh = zsig, 
                zscore.limit = zsig*2, 
                dot.sizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



# Paper figures -----------------------------------------------------------

set.seed(42)
dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
                zscoreLimit = zsig*2, 
                dotSizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(1, 9)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
pdf('running_code/paper_figures/cerebellum_dotplot_crawdad.pdf',
    height = 7, width = 8)
p
dev.off()



## Spatial plot ------------------------------------------------------------

data(slide)
cells <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)

## tried to reproduce the same colors as in the preprint paper, but could not
## so I decided to manually select the color of the specific cell types in 
## illustrator and pass them as argument of the plot function
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts))
# ct_colors <- setNames(tail(SteppedSequential5Steps, length(interest_cts)),
#                       interest_cts)
# ct_colors <- setNames(sample(rainbow(length(interest_cts))), interest_cts)
# saveRDS(ct_colors, 'running_code/processed_data/colors_slide.RDS')
ct_colors <- readRDS('running_code/processed_data/colors_slide.RDS')


p <- vizClusters(cells, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/cerebellum/spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



### Selected cell types -----------------------------------------------------

## Purkinje
interest_cts <- c('Purkinje', 
                  'Bergmann',
                  'Oligodendrocytes')
## these colors were alpha = .5 in the paper, how to convert them back?
interest_ct_colors <- ct_colors[interest_cts]
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/cerebellum/spatial_plot_purkinje.pdf',
    height = 7, width = 12)
p
dev.off()