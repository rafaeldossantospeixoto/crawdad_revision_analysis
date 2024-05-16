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

vizClusters(cells, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  theme_minimal()


### Selected cell types -----------------------------------------------------

## Granule
interest_cts <- c('Granule')
interest_ct_colors <- ct_colors[interest_cts]

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Granule', dist = 100) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p


for (d in c(10, 50, 100)){
  p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .001,
                   ref = 'Granule', dist = d) + 
    scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
    theme_minimal()
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'slide_granule_', d, '.pdf'),
      height = 2, width = 4)
  print(p)
  dev.off()
}



# Run embryo ---------------------------------------------------------

## Spatial plot ------------------------------------------------------------

data(seq)
cells <- crawdad:::toSF(pos = seq[,c("x", "y")],
                        celltypes = seq$celltypes)

ct_colors <- readRDS('running_code/processed_data/colors_seq.RDS')

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

vizClusters(cells, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  theme_minimal()



### Selected cell types -----------------------------------------------------

## Endothelium or Haematoendothelial progenitors
interest_cts <- c('Endothelium')
interest_ct_colors <- ct_colors[interest_cts]

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Endothelium', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p


for (d in c(10, 50, 100)){
  p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                   ref = 'Endothelium', dist = d) + 
    scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
    theme_minimal()
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'seq_endothelium_', d, '.pdf'),
      height = 7, width = 12)
  print(p)
  dev.off()
}


# Run spleen ---------------------------------------------------------

## Spatial plot ------------------------------------------------------------

data(pkhl)
cells <- crawdad:::toSF(pos = pkhl[,c("x", "y")],
                        celltypes = pkhl$celltypes)

ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

vizClusters(cells, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  theme_minimal()



### Selected cell types -----------------------------------------------------

## Neutrophils/Monocytes
interest_cts <- c('Neutrophils/Monocytes')
interest_ct_colors <- ct_colors[interest_cts]

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                 ref = 'Neutrophils/Monocytes', dist = 50) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
p


for (d in c(10, 50, 100)){
  p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .5,
                   ref = 'Neutrophils/Monocytes', dist = d) + 
    scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
    theme_minimal()
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'pkhl_neutrophilsmonocytes_', d, '.pdf'),
      height = 7, width = 12)
  print(p)
  dev.off()
}

