library(crawdad)
library(tidyverse)


# Read data ---------------------------------------------------------------

## read auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')


## Check nans --------------------------------------------------------------

## check for nans
auc[is.na(auc$auc),]
## most are Medium Spiny Neurons and s1r3

## viz some nans
readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  filter(reference == 'Ependymal Cells',
         neighbor == 'Medium Spiny Neurons') %>% 
  vizTrends(zSigThresh = 1.96)
readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  filter(reference == 'Medium Spiny Neurons',
         neighbor == 'Oligodendrocytes') %>% 
  vizTrends()
## there are nan Z-scores for some of the scales, so the auc is nan


## Filter shared pairs -----------------------------------------------------

create_id_column <- function(df){
  df_ids <- df %>%
    dplyr::mutate(id = paste0('s', slice, 'r', replicate))
  return(df_ids)
}
create_pair_column <- function(df){
  df_pairs <- df %>%
    dplyr::mutate(pair = paste0(paste0(reference, ' - ', neighbor)))
  return(df_pairs)
}
## filter data based on common values
filter_shared_pairs <- function(df_pairs){
  shared_pairs <- df_pairs %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarize(n_unique = dplyr::n_distinct(id)) %>% 
    dplyr::filter(n_unique == dplyr::n_distinct(df_pairs$id)) %>% 
    dplyr::pull(pair)
  
  filtered_df <- df_pairs %>% 
    filter(pair %in% shared_pairs)
  
  return(filtered_df)
}

auc_pairs <- auc %>% 
  create_id_column() %>% 
  create_pair_column()

## check the number of pairs for each sample
auc_pairs %>% 
  group_by(id) %>% 
  summarize(unique_pairs = n_distinct(pair))

filtered_auc <- filter_shared_pairs(auc_pairs)

  

# Reduced dimension -------------------------------------------------------


## Substitute NaNs for 0 ---------------------------------------------------

## I will substitute the NaNs for 0, just to try
## read auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')
auc <- auc %>% 
  mutate(sample = paste0('s', slice, 'r', replicate), 
         pair = paste0(paste0(reference, ' - ', neighbor))) %>% 
  select(c(pair, sample, auc)) %>% 
  replace_na(list(auc = 0))
auc_mtx <- auc %>% 
  pivot_wider(names_from = sample, values_from = auc) %>% 
  select(!pair) %>% 
  as.matrix() %>% 
  t()
## there are missing values in the matrix too because cell type pairs are missing
auc_mtx[is.na(auc_mtx)] <- 0

## calculate pca
pca <- prcomp(auc_mtx)
pcs <- pca$x[, 1:2] %>% 
  as.data.frame() %>% 
  mutate(sample = rownames(pca$x)) %>% 
  mutate(slice = sapply(rownames(pca$x), 
                        FUN = function(x) str_split(x, 'r')[[1]][1]))

## plot
pcs %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = slice)) + 
  scale_color_manual(values = c('#009440', '#ffcb00', '#302681'))
## it worked



## Drop pairs with NaNs ---------------------------------------------------

## read auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')
auc <- auc %>% 
  mutate(sample = paste0('s', slice, 'r', replicate), 
         pair = paste0(paste0(reference, ' - ', neighbor))) %>% 
  select(c(pair, sample, auc))
auc_mtx <- auc %>% 
  pivot_wider(names_from = sample, values_from = auc) %>% 
  select(!pair) %>% 
  drop_na() %>% ## drop nas
  as.matrix() %>% 
  t()

## calculate pca
pca <- prcomp(auc_mtx)
pcs <- pca$x[, 1:2] %>% 
  as.data.frame() %>% 
  mutate(sample = rownames(pca$x)) %>% 
  mutate(slice = sapply(rownames(pca$x), 
                        FUN = function(x) str_split(x, 'r')[[1]][1]))

p <- pcs %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = slice)) + 
  scale_color_manual(values = rainbow(3)) +
  theme_bw() +
  coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_pca.pdf'),
    height = 5, width = 7)
p
dev.off()
## still works



## Standardize and drop pairs with NaNs ------------------------------------

## The way to go!

## read and process auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')
processed_auc <- auc %>% 
  drop_na() %>% 
  create_id_column() %>% 
  create_pair_colum() %>% 
  filter_shared_pairs()

## check unique pairs
processed_auc %>% 
  group_by(pair) %>%
  summarize(groups_count = n_distinct(id))

## create matrix
processed_auc <- processed_auc %>% 
  select(c(pair, id, auc))
auc_mtx <- processed_auc %>% 
  pivot_wider(names_from = id, values_from = auc) %>% 
  select(!pair) %>% 
  drop_na() %>% ## drop nas
  as.matrix() %>% 
  t()

## normalize
auc_mtx <- scale(auc_mtx)
apply(auc_mtx, 2, var)

## calculate pca
pca <- prcomp(auc_mtx)
pcs <- pca$x[, 1:2] %>% 
  as.data.frame() %>% 
  mutate(id = rownames(pca$x)) %>% 
  mutate(slice = sapply(rownames(pca$x), 
                        FUN = function(x) str_split(x, 'r')[[1]][1]))

## plot
p <- pcs %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = slice)) + 
  scale_color_manual(values = c('#009440', '#ffcb00', '#302681')) +
  theme_bw() +
  coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_pca_standardized.pdf'),
    height = 5, width = 7)
p
dev.off()





# Variance ----------------------------------------------------------------

## read and process auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')
processed_auc <- auc %>% 
  drop_na() %>% 
  create_id_column() %>%
  create_pair_column() %>%
  filter_shared_pairs()

## all samples
processed_auc %>% 
  group_by(reference, neighbor) %>% 
  summarize(variance = var(auc)) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance)) +
  scale_radius(range = c(1, 10)) +
  theme_bw() +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 1))

## plotting parameters
all_cts <- unique(c(as.character(processed_auc$reference), 
                    as.character(processed_auc$neighbor)))

## slice
p <- processed_auc %>% 
  filter(slice == 1) %>% 
  group_by(reference, neighbor) %>% 
  summarize(variance = var(auc)) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance), color = '#006437') +
  scale_radius(range = c(1, 15),
               limits=c(1, 4e7),
               breaks=c(1e7, 2e7, 3e7, 4e7)) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(limits = all_cts, position = 'top') +
  ggplot2::scale_y_discrete(limits = all_cts) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=0)) +
  labs(size = 'mean scale') +
  ggplot2::coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_auc_variance_s1.pdf'),
    height = 5, width = 7)
p
dev.off()

## replicate
p <- processed_auc %>% 
  filter(replicate == 1) %>% 
  group_by(reference, neighbor) %>% 
  summarize(variance = var(auc)) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance), color = '#006437') +
  scale_radius(range = c(1, 15),
               limits=c(1, 4e7),
               breaks=c(1e7, 2e7, 3e7, 4e7)) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(limits = all_cts, position = 'top') +
  ggplot2::scale_y_discrete(limits = all_cts) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=0)) +
  labs(size = 'mean scale') +
  ggplot2::coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_auc_variance_r1.pdf'),
    height = 5, width = 7)
p
dev.off()




# Trends ------------------------------------------------------------------

ref <- 'GABAergic Estrogen-Receptive Neurons' 
nei <- 'Excitatory Neurons'

## S1R* --------------------------------------------------------------------

## S1R*
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r1')
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r2_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r2')
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r3')
dats <- bind_rows(dat1, dat2, dat3)

zsig <- correctZBonferroni(dats)

p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c(ggdark::invert_color('#009440'), 
                                ggdark::invert_color('#ffcb00'), 
                                ggdark::invert_color('#302681')))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_trends_s1.pdf'),
    height = 5, width = 7)
p
dev.off()



## S*R1 --------------------------------------------------------------------

## S*R1
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r1')
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s2_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's2r1')
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s3_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's3r1')
dats <- bind_rows(dat1, dat2, dat3)

zsig <- correctZBonferroni(dats)

p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('#009440', '#ffcb00', '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_trends_r1.pdf'),
    height = 5, width = 7)
p
dev.off()






# Spatial visualization ---------------------------------------------------

library(SpatialExperiment)
ref <- 'GABAergic Estrogen-Receptive Neurons' 
nei <- 'Excitatory Neurons'
ct_colors <- c('GABAergic Estrogen-Receptive Neurons' = 'maroon', 
               'Excitatory Neurons' = 'cyan', 
               'other' = '#E6E6E6')

## S1R1
spe <- readRDS('running_code/processed_data/merfish_mouseBrain_s1_r1.RDS')
cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                       celltypes = colData(spe)$celltype_merged)
p <- vizClusters(cells, ofInterest = c(ref, nei), alpha = 1, pointSize = 0.001) +
  scale_color_manual(values = ct_colors)
print(p)
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_s1r1.pdf'),
    height = 5, width = 7)
print(p)
dev.off()

## S*R1
for (r in 1:3) {
  spe <- readRDS(paste0('running_code/processed_data/merfish_mouseBrain_s1_r',
                        r, '.RDS'))
  cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), 
                         celltypes = colData(spe)$celltype_merged)
  p <- vizClusters(cells, ofInterest = c(ref, nei), alpha = 1, pointSize = 0.001) +
    scale_color_manual(values = ct_colors)
  print(p)
  pdf(paste0('function_development/comparing_samples/paper_figures/',
             'merfish_brains_s1r', r, '.pdf'),
      height = 5, width = 7)
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
    scale_color_manual(values = ct_colors)
  print(p)
  pdf(paste0('function_development/comparing_samples/paper_figures/',
             'merfish_brains_s', s, 'r1.pdf'),
      height = 5, width = 7)
  print(p)
  dev.off()
}
