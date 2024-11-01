
# All cts -----------------------------------------------------------------

library(crawdad)
library(tidyverse)
library(SpatialExperiment)

## S2R2
spe <- readRDS('running_code/processed_data/merfish_mouseBrain_s2_r2.RDS')
cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                       celltypes = colData(spe)$celltype_merged)

# ct_colors <- rainbow(14)
# names(ct_colors) <- unique(cells$celltypes)
# ## switch colors to match figures with only cts of interest
# ## Inhib neurons were cyan before
# ct_colors['Inhibitory Interneurons'] <- ct_colors['GABAergic Estrogen-Receptive Neurons']
# ct_colors['GABAergic Estrogen-Receptive Neurons'] <- 'maroon' # B03060
# ct_colors['Excitatory Neurons'] <- 'cyan' ## 00FFFF
# saveRDS(ct_colors, 'running_code/processed_data/colors_merfish.RDS')
ct_colors <- readRDS('running_code/processed_data/colors_merfish.RDS')

p <- vizClusters(cells, alpha = 1, pointSize = 0.001) +
  scale_color_manual(values = ct_colors)
print(p)
pdf(paste0('running_code/paper_figures/brain/',
           'merfish_brains_legend', '.pdf'),
    height = 7, width = 7)
print(p)
dev.off()

## S1R*
for (r in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s1_r',
                        r, '.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors) +
    ggplot2::guides(colour = "none")
  print(p)
  pdf(paste0('running_code/paper_figures/brain/',
             'merfish_brains_s1r', r, '.pdf'),
      height = 10, width = 10)
  print(p)
  dev.off()
}

## S2R*
for (r in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s2_r',
                        r, '.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors) +
    ggplot2::guides(colour = "none")
  print(p)
  pdf(paste0('running_code/paper_figures/brain/',
             'merfish_brains_s2r', r, '.pdf'),
      height = 10, width = 10)
  print(p)
  dev.off()
}

## S3R*
for (r in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s3_r',
                        r, '.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors) +
    ggplot2::guides(colour = "none")
  print(p)
  pdf(paste0('running_code/paper_figures/brain/',
             'merfish_brains_s3r', r, '.pdf'),
      height = 10, width = 10)
  print(p)
  dev.off()
}



# Selected ----------------------------------------------------------------

library(crawdad)
library(tidyverse)
library(SpatialExperiment)

ct_colors <- readRDS('running_code/processed_data/colors_merfish.RDS')
ct_colors <- c('GABAergic Estrogen-Receptive Neurons' = 'maroon', 
               'Excitatory Neurons' = 'cyan', 
               'other' = '#E6E6E6')

ref <- 'GABAergic Estrogen-Receptive Neurons' 
nei <- 'Excitatory Neurons'

## S1R1
spe <- readRDS('running_code/processed_data/merfish_mouseBrain_s1_r1.RDS')
cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                       celltypes = colData(spe)$celltype_merged)
p <- vizClusters(cells, ofInterest = c(ref, nei), alpha = 1, pointSize = 0.001) +
  scale_color_manual(values = ct_colors) +
  ggplot2::guides(colour = "none")
print(p)
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_s1r1.pdf'),
    height = 10, width = 10)
print(p)
dev.off()

## S*R1
for (r in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s1_r',
                        r, '.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, ofInterest = c(ref, nei), alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors) +
    ggplot2::guides(colour = "none")
  print(p)
  pdf(paste0('function_development/comparing_samples/paper_figures/',
             'merfish_brains_s1r', r, '.pdf'),
      height = 10, width = 10)
  print(p)
  dev.off()
}

## S1R*
for (s in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s',
                        s, '_r1.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, ofInterest = c(ref, nei), alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors) +
    ggplot2::guides(colour = "none")
  print(p)
  pdf(paste0('function_development/comparing_samples/paper_figures/',
             'merfish_brains_s', s, 'r1.pdf'),
      height = 10, width = 10)
  print(p)
  dev.off()
}
