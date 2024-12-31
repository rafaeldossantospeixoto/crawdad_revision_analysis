library(crawdad)
library(tidyverse)
ncores <- 7



# Run cerebellum ---------------------------------------------------------

## Spatial plot ------------------------------------------------------------

data(seq)
cells <- crawdad:::toSF(pos = seq[,c("x", "y")],
                        cellTypes = seq$celltypes)

ct_colors <- readRDS('running_code/processed_data/colors_seq.RDS')

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
interest_cts <- c('Presomitic mesoderm', 
                  'Spinal cord')
interest_ct_colors <- ct_colors[interest_cts]

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Presomitic mesoderm', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_presomitic_spinal_refpresomitic.pdf',
    height = 7, width = 12)
p
dev.off()

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Spinal cord', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p
pdf('running_code/paper_figures/asymmetry/spatial_plot_presomitic_spinal_refspinal.pdf',
    height = 7, width = 12)
p
dev.off()



## Trend plot --------------------------------------------------------------

set.seed(42)
dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')
zsig <- correctZBonferroni(dat_50)

### Selected cell types -----------------------------------------------------

## Presomitic mesoderm and Spinal cord
p <- dat_50 %>% 
  filter(reference == 'Presomitic mesoderm') %>% 
  filter(neighbor == 'Spinal cord') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refpresomitic_neighspinal.pdf',
    height = 4, width = 5)
p
dev.off()

p <- dat_50 %>% 
  filter(reference == 'Spinal cord') %>% 
  filter(neighbor == 'Presomitic mesoderm') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
p
pdf('running_code/paper_figures/asymmetry/trend_refspinal_neighpresomitic.pdf',
    height = 4, width = 5)
p
dev.off()

