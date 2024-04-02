library(crawdad)
library(tidyverse)


# Sim ---------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/results_data/dat_sim_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')
dat_rk <- readRDS('running_code/processed_data/dat_sim_ripleys.RDS')

dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)

data(sim)
cells <- crawdad::toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
p <- vizClusters(cells)
p
pdf('function_development/comparing_methods/paper_figures/sim_viz.pdf')
p
dev.off()

## Visualize specific pairs ------------------------------------------------

# ## B and C
# dat_sp %>% 
#   filter(reference == 'B') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=distance, y=probability))
# dat_50 %>% 
#   group_by(reference, neighbor, scale) %>% 
#   summarize(Z = mean(Z)) %>% 
#   filter(reference == 'B') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z))
# 
# ## D and C
# dat_sp %>% 
#   filter(reference == 'D') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=distance, y=probability))
# dat_50 %>% 
#   group_by(reference, neighbor, scale) %>% 
#   summarize(Z = mean(Z)) %>% 
#   filter(reference == 'D') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z))
# 
# ## A and C
# dat_sp %>% 
#   filter(reference == 'A') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=distance, y=probability))
# dat_50 %>% 
#   group_by(reference, neighbor, scale) %>% 
#   summarize(Z = mean(Z)) %>% 
#   filter(reference == 'A') %>% 
#   filter(neighbor == 'C') %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z))



## Highlight specific pairs ------------------------------------------------

selected_cts <- c('D', 'C')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)


### Save plots --------------------------------------------------------------

reference_ct <- 'B'

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/sim_crawdad.pdf', 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/sim_squidpy.pdf', 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/sim_ripleys.pdf', 
    height = 4, width = 5)
p
dev.off()




# Slide -------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/results_data/dat_slide_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')
dat_rk <- readRDS('running_code/processed_data/dat_cerebellum_ripleys.RDS')

dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)

# data(slide)
# cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
# p <- vizClusters(cells)
# p
# pdf('function_development/comparing_methods/paper_figures/cerebellum_viz.pdf')
# p
# dev.off()
# 
# data(slide)
# cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
# p <- vizClusters(cells, ofInterest = c('Purkinje', 'Bergmann', 'Oligodendrocytes'))
# p
# pdf('function_development/comparing_methods/paper_figures/cerebellum_viz_selected.pdf')
# p
# dev.off()

## Visualize specific pairs ------------------------------------------------

# ## Purkinje and Bergmann
# dat_sp %>% 
#   filter(reference == 'Purkinje') %>% 
#   filter(neighbor == 'Bergmann') %>% 
#   ggplot() + 
#   geom_line(aes(x=distance, y=probability))
# 
# dat_rk %>% 
#   filter(reference == 'Purkinje') %>% 
#   filter(neighbor == 'Bergmann') %>% 
#   ggplot() + 
#   geom_line(aes(x=radius, y=score))
# 
# dat_50 %>% 
#   filter(reference == 'Purkinje') %>% 
#   filter(neighbor == 'Bergmann') %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z)) +
#   ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.6, linetype = "dotted")
# 
# ## Astrocytes and Microglia
# dat_sp %>% 
#   filter(reference == 'Astrocytes') %>% 
#   filter(neighbor == 'Microglia') %>% 
#   ggplot() + 
#   geom_line(aes(x=distance, y=probability))
# 
# dat_50 %>% 
#   filter(reference == 'Astrocytes') %>% 
#   filter(neighbor == 'Microglia') %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z))


## Highlight specific pairs ------------------------------------------------

selected_cts <- c('Oligodendrocytes', 'Bergmann')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)



### Different plots -------------------------------------------------------

# ## all cells in legend
# dat_50 %>% 
#   filter(reference == 'Purkinje') %>% 
#   mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z, color = neighbor)) +
#   scale_color_manual(values = ordered_colors) +
#   theme_bw()
# 
# ## selected cells in legend
# dat_50 %>% 
#   filter(reference == 'Purkinje') %>% 
#   mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
#   mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
#   ggplot() + 
#   geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_neighbor), 
#             size = .5) +
#   scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
#   scale_size_manual(name = 'Neighbor', values = c(.5, .5, 3)) +
#   theme_bw()
# 
# ## highlighting one of them
# dat_50 %>% 
#   filter(reference == 'Purkinje') %>% 
#   mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
#   ggplot(aes(x = scale, y = Z)) + 
#   geom_line(aes(colour = neighbor == "Bergmann", size = neighbor == "Bergmann", group = neighbor)) +
#   scale_color_manual(name = "neighbor", labels = c("Other", "Bergmann"), values = c("gray", "red")) +
#   scale_size_manual(name = "neighbor", labels = c("Other", "Bergmann"), values = c(0.5, 1)) +
#   theme_bw()


### Save plots --------------------------------------------------------------

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
pdf('function_development/comparing_methods/paper_figures/cerebellum_crawdad.pdf', 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/cerebellum_squidpy.pdf', 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = selected_neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/cerebellum_ripleys.pdf', 
    height = 4, width = 5)
p
dev.off()


# Seq -------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/results_data/dat_seq_squidpy.csv', 
                   row.names = 1)
dat_sp <- dat_sp %>% 
  filter(distance <= 750)

dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')
dat_50 <- dat_50 %>% 
  filter(scale <= 750)
dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)

# data(seq)
# cells <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
# p <- vizClusters(cells)
# p
# pdf('function_development/comparing_methods/paper_figures/embryo_viz.pdf')
# p
# dev.off()
# 
# data(seq)
# cells <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
# p <- vizClusters(cells, ofInterest = c('Endothelium', 'Spinal cord','Haematoendothelial progenitors'))
# p
# pdf('function_development/comparing_methods/paper_figures/embryo_viz_selected.pdf')
# p
# dev.off()


## Highlight Endothelium ------------------------------------------------

reference_ct <- 'Endothelium'
dat_rk <- readRDS(paste0("running_code/processed_data/dat_embryo_ripleys_",
                         sub(" ", "_", reference_ct),
                         ".RDS"))
dat_rk <- dat_rk %>% 
  filter(radius <= 750)

selected_cts <- c('Forebrain/Midbrain/Hindbrain','Haematoendothelial progenitors')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)

### Save plots --------------------------------------------------------------

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == reference_ct) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
p
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_crawdad_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == reference_ct) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_squidpy_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == reference_ct) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = selected_neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_ripleys_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()


## Highlight Intermediate mesoderm ------------------------------------------------

reference_ct <- 'Intermediate mesoderm'
dat_rk <- readRDS(paste0("running_code/processed_data/dat_embryo_ripleys_",
                         sub(" ", "_", reference_ct),
                         ".RDS"))

selected_cts <- c('Cardiomyocytes','Lateral plate mesoderm')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)

### Save plots --------------------------------------------------------------

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == reference_ct) %>% 
  filter(scale <= 750) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_crawdad_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == reference_ct) %>% 
  filter(distance <= 750) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_squidpy_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == reference_ct) %>% 
  filter(radius <= 750) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = selected_neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf(paste0('function_development/comparing_methods/paper_figures/embryo_ripleys_',
           sub(" ", "_", reference_ct),
           ".pdf"), 
    height = 4, width = 5)
p
dev.off()




# NullSim ------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/results_data/dat_nullsim_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_nullsim_50.RDS')
dat_rk <- readRDS('running_code/processed_data/dat_nullsim_ripleys.RDS')

dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)

# data(sim)
# cells <- crawdad::toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
# p <- vizClusters(cells)
# p
# pdf('function_development/comparing_methods/paper_figures/nullsim_viz.pdf')
# p
# dev.off()


## Highlight specific pairs ------------------------------------------------

selected_cts <- c('B', 'C')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)


### Save plots --------------------------------------------------------------

reference_ct <- 'A'

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/nullsim_crawdad.pdf', 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/nullsim_squidpy.pdf', 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/nullsim_ripleys.pdf', 
    height = 4, width = 5)
p
dev.off()




# Spleen ------------------------------------------------------------------

## pkhl -------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/results_data/spleen/dat_pkhl_squidpy.csv', row.names = 1)
# dat_50 <- readRDS('running_code/processed_data/dat_nullsim_50.RDS')
# dat_rk <- readRDS('running_code/processed_data/dat_nullsim_ripleys.RDS')

dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)

# data(sim)
# cells <- crawdad::toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
# p <- vizClusters(cells)
# p
# pdf('function_development/comparing_methods/paper_figures/nullsim_viz.pdf')
# p
# dev.off()


## Highlight specific pairs ------------------------------------------------

selected_cts <- c('Sinusoidal cells', 'Fol B cells')
all_cts <- unique(as.character(dat_sp$reference))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)


### Save plots --------------------------------------------------------------

reference_ct <- 'CD4 Memory T cells'

## CRAWDAD
p <- dat_50 %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_cts), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = rainbow(4)) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.3, linetype = "dashed") + 
  ggplot2::geom_hline(yintercept = -zsig, color = "black", size = 0.3, linetype = "dashed") + 
  labs(title = 'CRAWDAD') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/pkhl_crawdad.pdf', 
    height = 4, width = 5)
p
dev.off()

## Squidpy
p <- dat_sp %>% 
  filter(reference == reference_ct) %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = 'Squidpy Co-occurrence') + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/pkhl_squidpy.pdf', 
    height = 4, width = 5)
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == reference_ct) %>% 
  # mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  # mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=radius, y=score, group = neighbor, color = neighbor), 
            size = .5) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/pkhl_ripleys.pdf', 
    height = 4, width = 5)
p
dev.off()
