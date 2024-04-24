library(crawdad)
library(tidyverse)
ncores <- 7

# Run cerebellum ---------------------------------------------------------

## Spatial plot ------------------------------------------------------------

data(slide)
cells <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)

ct_colors <- readRDS('running_code/processed_data/colors_slide.RDS')

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

p <- vizClusters(cells, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p


### Selected cell types -----------------------------------------------------

## UBCs and Granule
interest_cts <- c('UBCs', 
                  'Granule')
interest_ct_colors <- ct_colors[interest_cts]
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = 1) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_ubcs_granule.pdf',
    height = 7, width = 12)
p
dev.off()


## Trend plot --------------------------------------------------------------

set.seed(42)
dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')
zsig <- correctZBonferroni(dat_50)

### Selected cell types -----------------------------------------------------

## UBCs and Granule
p <- dat_50 %>% 
  filter(reference == 'UBCs') %>% 
  filter(neighbor == 'Granule') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refUBCs_neighGranule.pdf',
    height = 4, width = 5)
p
dev.off()
p <- dat_50 %>% 
  filter(reference == 'Granule') %>% 
  filter(neighbor == 'UBCs') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refGranule_neighUBCs.pdf',
    height = 4, width = 5)
p
dev.off()





# Sim ---------------------------------------------------------------------


