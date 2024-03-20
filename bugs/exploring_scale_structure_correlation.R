library(crawdad)
library(tidyverse)

set.seed(42)
ncores <- 7

# Original sim ------------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')
zsig <- correctZBonferroni(dat_50)

df_original <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE))
df_original <- df_original %>% 
  mutate(pair = paste0(reference, neighbor),
         mode = 'original')
df_original

data(sim)
filtered_sim <- sim %>% 
  filter(x < 1000, 
         y < 1000)
cells <- toSF(pos = filtered_sim[, c('x', 'y')], 
              celltypes = filtered_sim$celltypes)
vizClusters(cells, pointSize = 3)

filtered_sim %>% 
  filter(celltypes == 'A') %>% 
  pull(x) %>% 
  range()


# New sims ----------------------------------------------------------------

simulate_data <- function(ncells = 2000, size = 1000, 
                          inner_radius = .1, outter_radius = .3){
  x <- runif(ncells, min = 0, max = 1)
  y <- runif(ncells, min = 0, max = 1)
  
  df <- data.frame(x = x, y = y, celltypes = 'D')
  
  df <- df %>% 
    mutate(celltypes = case_when( (sqrt((x-.5)^2 + (y-.5)^2) < outter_radius) & 
                                    (sqrt((x-.5)^2 + (y-.5)^2) > inner_radius) ~ 'A',
                                  sqrt((x-.5)^2 + (y-.5)^2) < inner_radius ~ 'mixed',
                                  T ~ 'D'))
  mixed_idx <- which(df$celltypes == 'mixed')
  b_idx <- sample(mixed_idx, size = round(length(mixed_idx)/2))
  
  df[mixed_idx, 'celltypes'] <- 'C'
  df[b_idx, 'celltypes'] <- 'B'
  
  df <- df %>% 
    mutate(x = x * size,
           y = y * size)
  
  full_df <- rbind(df,
                   df %>% mutate(x = x + size),
                   df %>% mutate(y = y + size),
                   df %>% mutate(x = x + size, y = y + size))
  
  return(full_df)
}

ncell <- 2000
size <- 1000
inner_radius <- .1
outter_radius <- .3
df <- simulate_data()
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))



## Increase 100 ------------------------------------------------------------

ncell <- 2000
size <- 1000
inner_radius <- .2
outter_radius <- .4
df <- simulate_data(ncell, size, inner_radius, outter_radius)
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))

cells <- toSF(pos = df[, c('x', 'y'),], celltypes = df$celltypes)
scales <- seq(100, 1000, by=50)
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 5,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
results_50 <- crawdad::findTrends(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)

df_larger_r100 <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE))
df_larger_r100 <- df_larger_r100 %>% 
  mutate(pair = paste0(reference, neighbor),
         mode = 'larger_r100')
df_larger_r100

df %>% 
  filter(x < 1000, y < 1000) %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))


## Increase 200 ------------------------------------------------------------

ncell <- 2000
size <- 1000
inner_radius <- .3
outter_radius <- .5
df <- simulate_data(ncell, size, inner_radius, outter_radius)
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))

cells <- toSF(pos = df[, c('x', 'y'),], celltypes = df$celltypes)
scales <- seq(100, 1000, by=50)
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 5,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
results_50 <- crawdad::findTrends(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)

df_larger_r200 <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE))
df_larger_r200 <- df_larger_r200 %>% 
  mutate(pair = paste0(reference, neighbor),
         mode = 'larger_r200')
df_larger_r200

df %>% 
  filter(x < 1000, y < 1000) %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))




## Increase 100 inner -------------------------------------------------------

ncell <- 2000
size <- 1000
inner_radius <- .2
outter_radius <- .3
df <- simulate_data(ncell, size, inner_radius, outter_radius)
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))

cells <- toSF(pos = df[, c('x', 'y'),], celltypes = df$celltypes)
scales <- seq(100, 1000, by=50)
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 5,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
results_50 <- crawdad::findTrends(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)

df_larger_ic <- dat_50 %>% 
  group_by(scale, reference, neighbor) %>% 
  summarise(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE))
df_larger_ic <- df_larger_ic %>% 
  mutate(pair = paste0(reference, neighbor),
         mode = 'larger_ic')

df %>% 
  filter(x < 1000, y < 1000) %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(4))



# Compare all -------------------------------------------------------------

df_all <- rbind(df_original,
                df_larger_r100,
                df_larger_r200,
                df_larger_ic)
head(df_all)
saveRDS(df_all, 'bugs/processed_data/df_scale_structure_correlation.RDS')

df_changed <- df_all %>% 
  group_by(pair) %>% 
  summarise(all_same = n_distinct(scale) == T)
  
df_all %>% 
  left_join(df_changed) %>% 
  mutate(mode = factor(mode, levels = c('larger_ic', 'original',
                                        'larger_r100', 'larger_r200'))) %>% 
  filter(all_same == F) %>% 
  filter(reference == 'A') %>% 
  ggplot() + 
  geom_line(aes(mode, scale, group = pair, color = pair)) + 
  geom_point(aes(mode, scale, color = pair))

df_all %>% 
  left_join(df_changed) %>% 
  mutate(mode = factor(mode, levels = c('larger_ic', 'original',
                                        'larger_r100', 'larger_r200'))) %>% 
  filter(all_same == F) %>% 
  ggplot() + 
  geom_line(aes(mode, scale, group = pair, color = pair)) + 
  geom_point(aes(mode, scale, color = pair)) + 
  scale_color_manual(values = rainbow(10)) + 
  facet_wrap('reference')
