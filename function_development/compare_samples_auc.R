library(crawdad)
library(tidyverse)


# Read data ---------------------------------------------------------------

## read auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')

## check for nans
auc[is.na(auc$auc),]
## most are Medium Spiny Neurons and s1r3

## viz some nans
readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  filter(reference == 'Ependymal Cells',
         neighbor == 'Medium Spiny Neurons') %>% 
  vizTrends()
readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  filter(reference == 'Medium Spiny Neurons',
         neighbor == 'Oligodendrocytes') %>% 
  vizTrends()
## there are nan Z-scores for some of the scales, so the auc is nan




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


# Drop NaNs ---------------------------------------------------------------

## read auc data
auc <- readRDS('running_code/outputs/merfish_mouseBrain_diff_auc_dist_50.RDS')
auc <- auc %>% 
  drop_na()

## all samples
auc %>% 
  group_by(reference, neighbor) %>% 
  summarize(variance = var(auc)) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance)) +
  scale_radius(range = c(1, 10)) +
  theme_bw() +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 1))

## slice and replicate
auc %>% 
  filter(replicate == 3) %>% 
  group_by(reference, neighbor) %>% 
  summarize(variance = var(auc)) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance)) +
  scale_radius(range = c(1, 10),
               limits=c(1, 4e7),
               breaks=c(1e7, 2e7, 3e7, 4e7)) +
  theme_bw() +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 1))
