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
  # coord_fixed() + 
  theme_minimal()
p


### Selected cell types -----------------------------------------------------

## UBCs and Granule
interest_cts <- c('UBCs', 
                  'Granule')
interest_ct_colors <- ct_colors[interest_cts]

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'UBCs', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_ubcs_granule_refubcs.pdf',
    height = 7, width = 12)
p
dev.off()

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Granule', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_ubcs_granule_refgranule.pdf',
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
pdf('running_code/paper_figures/asymmetry/trend_refubcs_neighgranule.pdf',
    height = 4, width = 5)
p
dev.off()
p <- dat_50 %>% 
  filter(reference == 'Granule') %>% 
  filter(neighbor == 'UBCs') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refgranule_neighubcs.pdf',
    height = 4, width = 5)
p
dev.off()





# Sim ---------------------------------------------------------------------

sim <- readRDS('running_code/processed_data/sim_asymmetrical_imbalance_B.RDS')


## Spatial plot ------------------------------------------------------------

cells <- crawdad:::toSF(pos = sim[,c("x", "y")],
                        celltypes = sim$celltypes)

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

p <- vizClusters(cells, alpha = 1, pointSize = 1) +
  coord_fixed() + 
  theme_minimal()
p


### Selected cell types -----------------------------------------------------

## UBCs and Granule
interest_cts <- c('B', 
                  'C')
interest_ct_colors <- c("#FF0000", "#80FF00", "#00FFFF", "#8000FF")
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = 1) + 
  scale_color_manual(values = interest_ct_colors[2:3], na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_b_c.pdf',
    height = 7, width = 12)
p
dev.off()


## Trend plot --------------------------------------------------------------

set.seed(42)
dat_50 <- readRDS('running_code/processed_data/dat_50_sim_asymmetrical_imbalance_B.RDS')
zsig <- correctZBonferroni(dat_50)

### Selected cell types -----------------------------------------------------

## UBCs and Granule
p <- dat_50 %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refB_neighC.pdf',
    height = 4, width = 5)
p
dev.off()
p <- dat_50 %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'B') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refC_neighB.pdf',
    height = 4, width = 5)
p
dev.off()




