library(crawdad)
library(tidyverse)


# Paper figures -----------------------------------------------------------

data(seq)
dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')
zsig <- correctZBonferroni(dat_50)

ref_ct <- 'Intermediate mesoderm'
ngb_ct <- 'Lateral plate mesoderm'

ref_ct <- 'Endothelium'
ngb_ct <- 'Haematoendothelial progenitors'

filtered_dat <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  filter(reference == ref_ct)
ggplot() + 
  geom_line(data = filter(filtered_dat, neighbor == ngb_ct), 
            aes(x = scale, y = Z),
            linewidth = 2) +
  geom_line(data = filtered_dat,
            aes(x = scale, y = Z, color = neighbor)) + 
  scale_color_manual(values = rainbow(length(unique(seq$celltypes)))) + 
  geom_hline(yintercept = zsig, linetype = 2) + 
  geom_hline(yintercept = -zsig, linetype = 2) +
  labs(title = paste(ref_ct, 'vs', ngb_ct))
  


# Correlation of self-enrichment with number --------------------------------

ncells_df <- as.data.frame(table(seq$celltypes))
names(ncells_df) <- c('celltype', 'ncells')

df <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE)) %>% 
  filter(reference == neighbor) %>%
  left_join(ncells_df, by = c('reference' = 'celltype')) 
df %>% 
  ggplot() +
  geom_point(aes(x = ncells, y = scale), size = 5) + 
  geom_point(aes(x = ncells, y = scale, color = reference), size = 3) + 
  scale_color_manual(values = rainbow(length(unique(df$reference)))) + 
  guides(color = guide_legend(ncol = 1)) +
  labs(title = 'Self-enrichment vs number of cells of that celltype')



# Correlation self-enrichment with density ----------------------------------

d <- 50
data(seq)
cells <- toSF(pos = seq[, c('x', 'y')], celltypes = seq$celltypes)

areas <- sapply(unique(cells$celltypes), function(ct){
  ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
  ref.buffer_union <- sf::st_union(ref.buffer)
  plot(ref.buffer_union)
  area <- sf::st_area(ref.buffer_union)
  return(area)
})
df_areas <- data.frame(area = areas, 
                       celltype = unique(cells$celltypes))
df <- df %>%
  left_join(df_areas, by = c('reference' = 'celltype')) %>% 
  mutate(density = ncells/area)

df %>% 
  ggplot() +
  geom_point(aes(x = density, y = scale), size = 5) + 
  geom_point(aes(x = density, y = scale, color = reference), size = 3) + 
  scale_color_manual(values = rainbow(length(unique(df$reference)))) + 
  guides(color = guide_legend(ncol = 1)) +
  labs(title = 'Self-enrichment vs density of that celltype')



# Correlation of n of sig relationships -----------------------

df <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE)) %>% 
  group_by(reference) %>%
  summarise(nsig = n()) %>% 
  left_join(ncells_df, by = c('reference' = 'celltype')) %>%
  left_join(df_areas, by = c('reference' = 'celltype')) %>% 
  mutate(density = ncells/area)

df %>% ggplot() +
  geom_point(aes(x = ncells, y = nsig), size = 5) + 
  geom_point(aes(x = ncells, y = nsig, color = reference), size = 3) + 
  geom_smooth(aes(x = ncells, y = nsig), method = "lm", se = T) +
  ggpubr::stat_cor(aes(x = ncells, y = nsig), method = "pearson") + 
  scale_color_manual(values = rainbow(length(unique(df$reference)))) + 
  guides(color = guide_legend(ncol = 1)) +
  labs(title = 'Number of sig relationships vs number of cells')

df %>% ggplot() +
  geom_point(aes(x = density, y = nsig), size = 5) + 
  geom_point(aes(x = density, y = nsig, color = reference), size = 3) + 
  geom_smooth(aes(x = density, y = nsig), method = "lm", se = T) + 
  ggpubr::stat_cor(aes(x = density, y = nsig), method = "pearson") + 
  scale_color_manual(values = rainbow(length(unique(df$reference)))) + 
  guides(color = guide_legend(ncol = 1)) +
  labs(title = 'Number of sig relationships vs density')
