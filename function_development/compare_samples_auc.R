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

create_colums <- function(auc){
  auc_pairs <- auc %>%
    mutate(sample = paste0('s', slice, 'r', replicate), 
           pair = paste0(paste0(reference, ' - ', neighbor)))
  
  return(auc_pairs)
}

auc_pairs <- create_colums(auc)

## check the number of pairs for each sample
auc_pairs %>% 
  group_by(sample) %>% 
  summarize(unique_pairs = n_distinct(pair))

## filter data based on common values
filter_shared_pairs <- function(auc_pairs){
  shared_pairs <- auc_pairs %>% 
    group_by(pair) %>%
    summarize(groups_count = n_distinct(sample)) %>% 
    filter(groups_count == n_distinct(auc_pairs$sample)) %>% 
    pull(pair)
  
  filtered_auc <- auc_pairs %>% 
    filter(pair %in% shared_pairs)
  
  return(filtered_auc)
}

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
  scale_color_manual(values = rainbow(3))
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
  create_colums() %>% 
  filter_shared_pairs()

## check unique pairs
processed_auc %>% 
  group_by(pair) %>%
  summarize(groups_count = n_distinct(sample))

## create matrix
processed_auc <- processed_auc %>% 
  select(c(pair, sample, auc))
auc_mtx <- processed_auc %>% 
  pivot_wider(names_from = sample, values_from = auc) %>% 
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
  mutate(sample = rownames(pca$x)) %>% 
  mutate(slice = sapply(rownames(pca$x), 
                        FUN = function(x) str_split(x, 'r')[[1]][1]))

## plot
p <- pcs %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = slice)) + 
  scale_color_manual(values = rainbow(3)) +
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
  create_colums() %>%
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
  theme_bw() +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 1)) + 
  coord_equal()
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
  theme_bw() +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 1)) +
  coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_auc_variance_r1.pdf'),
    height = 5, width = 7)
p
dev.off()
