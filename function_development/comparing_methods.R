library(crawdad)
library(tidyverse)


# Sim ---------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/dat_sim_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')

zsig <- 3

# dat_sp %>% 
#   filter(reference == 'A') %>% 
#   filter(neighbor == 'B') %>% 
#   vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)


## CRAWDAD vs Squidpy ------------------------------------------------------

## B and C
dat_sp %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))



## D and C
dat_sp %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))



## A and C
dat_sp %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))


# Slide -------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/dat_slide_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')
dat_rk <- readRDS('running_code/processed_data/dat_cerebellum_ripleys.RDS')

dat_50 <- dat_50 %>% group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z))

zsig <- correctZBonferroni(dat_50)


## Visualize specific pairs ------------------------------------------------

## Purkinje and Bergmann
dat_sp %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Bergmann') %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability))

dat_rk %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Bergmann') %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=score))

dat_50 %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Bergmann') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z)) +
  ggplot2::geom_hline(yintercept = zsig, color = "black", size = 0.6, linetype = "dotted")

## Astrocytes and Microglia
dat_sp %>% 
  filter(reference == 'Astrocytes') %>% 
  filter(neighbor == 'Microglia') %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=probability))

dat_50 %>% 
  filter(reference == 'Astrocytes') %>% 
  filter(neighbor == 'Microglia') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))


## Highlight specific pairs ------------------------------------------------

selected_cts <- c('Oligodendrocytes', 'Bergmann')
all_cts <- unique(as.character(dat_50$neighbor))
ordered_cts <- c(all_cts[!all_cts %in% selected_cts], selected_cts)
selected_colors <- c('blue', 'red')
all_colors <- c(rep('lightgray', times = length(all_cts) - length(selected_cts)), selected_colors)
ordered_colors <- setNames(all_colors, ordered_cts)



### Different plots -------------------------------------------------------

## all cells in legend
dat_50 %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, color = neighbor)) +
  scale_color_manual(values = ordered_colors) +
  theme_bw()

## selected cells in legend
dat_50 %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  scale_size_manual(name = 'Neighbor', values = c(.5, .5, 3)) +
  theme_bw()

## highlighting one of them
dat_50 %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  ggplot(aes(x = scale, y = Z)) + 
  geom_line(aes(colour = neighbor == "Bergmann", size = neighbor == "Bergmann", group = neighbor)) +
  scale_color_manual(name = "neighbor", labels = c("Other", "Bergmann"), values = c("gray", "red")) +
  scale_size_manual(name = "neighbor", labels = c("Other", "Bergmann"), values = c(0.5, 1)) +
  theme_bw()


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
pdf('function_development/comparing_methods/paper_figures/cerebellum_crawdad.pdf')
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
pdf('function_development/comparing_methods/paper_figures/cerebellum_squidpy.pdf')
p
dev.off()

## Ripleys
p <- dat_rk %>% 
  filter(reference == 'Purkinje') %>% 
  mutate(neighbor = factor(neighbor, levels = ordered_cts)) %>% 
  mutate(selected_neighbor = fct_other(neighbor, keep = selected_cts)) %>% 
  ggplot() + 
  geom_line(aes(x=distance, y=score, group = neighbor, color = selected_neighbor), 
            size = .5) +
  scale_color_manual(name = 'Neighbor', values = c(selected_colors, 'lightgray')) +
  labs(title = "Ripley's K (isotropic-corrected minus theoretical)") + 
  theme_bw()
p
pdf('function_development/comparing_methods/paper_figures/cerebellum_ripleys.pdf')
p
dev.off()
