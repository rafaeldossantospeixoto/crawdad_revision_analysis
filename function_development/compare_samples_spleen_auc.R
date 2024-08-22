
# Description -------------------------------------------------------------

## Draft compare samples using AUC and run comparison for spleen


library(crawdad)
library(tidyverse)


# Read data ---------------------------------------------------------------

samples <- c('pkhl', 'xxcd', 'pbvn', 'fsld', 'ngpl', 'ksfb')
patients <- c('vnkn', 'vnkn', 'zwnt', 'zwnt', 'kgnj', 'kgnj')
auc_list <- list()

for (sn in 1:length(samples)) {
  
  sample <- samples[sn]
  patient <- patients[sn]
  dat <- readRDS(paste0('running_code/processed_data/spleen/dat_',
                        sample, '_50.RDS'))
  
  ## define pairs
  pairs <- distinct(dat[c("neighbor", "reference")])
  
  ## compute AUC for each pair
  ## could use map_dfc instead of do.call rbind
  auc_sample <- do.call(rbind, lapply(seq(nrow(pairs)), function(i) {
    summarized_dat <- dat %>%
      filter(reference == pairs$reference[i]) %>%
      filter(neighbor == pairs$neighbor[i]) %>%
      group_by(neighbor, scale, reference) %>%
      summarise(mean = mean(Z),
                sd = sd(Z))
    return(data.frame(sample = sample,
                      patient = patient,
                      reference = pairs$reference[i], 
                      neighbor = pairs$neighbor[i], 
                      auc = pracma::trapz(summarized_dat$scale, summarized_dat$mean)))
  }))
  
  auc_list[[sn]] <- auc_sample
}

auc <- bind_rows(auc_list)


# Process data ------------------------------------------------------------

## create ct pair and id colunms
df_pairs <- auc %>% 
  dplyr::mutate(id = paste0('p_', patient, '_s_', sample)) %>% 
  dplyr::mutate(pair = paste0(paste0(reference, ' - ', neighbor))) 

## calculate the ct pairs that are shared across all samples
shared_pairs <- df_pairs %>% 
  dplyr::group_by(pair) %>%
  dplyr::summarize(n_unique = dplyr::n_distinct(id)) %>% 
  dplyr::filter(n_unique == dplyr::n_distinct(df_pairs$id)) %>% 
  dplyr::pull(pair)

## filter data based on shared pairs
processed_auc <- df_pairs %>% 
  dplyr::filter(pair %in% shared_pairs)

## create matrix
auc_mtx <- processed_auc %>% 
  select(c(pair, id, auc)) %>% 
  pivot_wider(names_from = id, values_from = auc) %>% 
  select(!pair) %>% 
  drop_na() %>% ## drop nas
  as.matrix() %>% 
  t()

## normalize
auc_mtx <- scale(auc_mtx)
apply(auc_mtx, 2, mean)

## calculate pca
pca <- prcomp(auc_mtx)
pcs <- pca$x[, 1:2] %>% 
  as.data.frame() %>% 
  mutate(id = rownames(pca$x)) %>% 
  mutate(patient = sapply(rownames(pca$x), 
                        FUN = function(x) str_split(x, '_')[[1]][2]))



# Visualize ---------------------------------------------------------------


## PCA ---------------------------------------------------------------------

## plot
p <- pcs %>% 
  ggplot() + 
  geom_point(aes(x = PC1, y = PC2, color = patient)) + 
  scale_color_manual(values = c('#009440', '#ffcb00', '#302681')) +
  theme_bw() +
  coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_pca_standardized.pdf'),
    height = 5, width = 7)
p
dev.off()



## Variance ----------------------------------------------------------------

## all samples
p <- processed_auc %>% 
  group_by(reference, neighbor) %>%
  filter(reference != 'indistinct',
         neighbor != 'indistinct') %>% 
  summarize(variance = (var(auc))) %>%
  ggplot() +
  geom_point(aes(x = reference, y = neighbor, size = variance), color = '#006437') +
  scale_radius(range = c(1, 10)) +
  ggplot2::theme_bw() +
  ggplot2::scale_x_discrete(position = 'top') +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=0)) +
  ggplot2::coord_equal()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_auc_variance_all_samples.pdf'),
    height = 5, width = 7)
p
dev.off()

processed_auc %>% 
  group_by(reference, neighbor) %>% 
  filter(reference != 'indistinct',
         neighbor != 'indistinct') %>% 
  summarize(variance = var(auc)) %>% 
  arrange(desc(variance)) %>% 
  print(n = 15)
# reference             neighbor             variance
# 1 Sinusoidal cells      Fol B cells      680558136.
# 2 Ki67 proliferating    Fol B cells      669102209.
# 3 CD8 Memory T cells    Fol B cells      506100950.
# 4 Blood endothelial     Fol B cells      367240977.
# 5 B cells, red pulp     Fol B cells      336635815.
# 6 Neutrophils/Monocytes Fol B cells      195936557.
# 11 Podoplanin            Sinusoidal cells      142024678.
# 12 Sinusoidal cells      CD4 Memory T cells    129769372.
# 13 Fol B cells           CD4 Memory T cells    121904619.

processed_auc %>% 
  group_by(reference, neighbor) %>% 
  filter(reference != 'indistinct',
         neighbor != 'indistinct') %>% 
  summarize(variance = var(auc)) %>% 
  arrange((variance)) %>% 
  print(n = 15)
# 1 Neutrophils/Monocytes Podoplanin          755749.
# 2 CD4 Memory T cells    Macrophages        1409475.
# 3 Neutrophils/Monocytes B cells, red pulp  1934633.
# 4 Macrophages           Podoplanin         2090387.
# 5 Neutrophils/Monocytes Ki67 proliferating 2202569.
# 6 Neutrophils/Monocytes Blood endothelial  2205917.
# 7 Macrophages           B cells, red pulp  2219647.
# 8 Macrophages           Myeloid cells      2237449.
# 9 CD4 Memory T cells    Ki67 proliferating 2359841.
# 10 Myeloid cells         Blood endothelial  2390182.


## Trends ------------------------------------------------------------------

samples <- c('pkhl', 'xxcd', 'pbvn', 'fsld', 'ngpl', 'ksfb')
patients <- c('vnkn', 'vnkn', 'zwnt', 'zwnt', 'kgnj', 'kgnj')
dat_list <- list()

## read data
for (sn in 1:length(samples)) {
  sample <- samples[sn]
  patient <- patients[sn]
  dat <- readRDS(paste0('running_code/processed_data/spleen/dat_',
                        sample, '_50.RDS')) %>% 
    mutate(sample = sample,
           patient = patient,
           id = paste0('p_', patient, '_s_', sample))
  dat_list[[sn]] <- dat
}
dats <- bind_rows(dat_list)

zsig <- correctZBonferroni(dats)



### High variance -----------------------------------------------------------

ref <- 'Sinusoidal cells' 
nei <- 'Fol B cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top1.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'Ki67 proliferating' 
nei <- 'Fol B cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top2.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'CD8 Memory T cells' 
nei <- 'Fol B cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top3.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'Blood endothelial' 
nei <- 'Fol B cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top4.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'B cells, red pulp' 
nei <- 'Fol B cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top5.pdf'),
    height = 5, width = 7)
p
dev.off()

# 11 Podoplanin            Sinusoidal cells      142024678.
# 12 Sinusoidal cells      CD4 Memory T cells    129769372.
# 13 Fol B cells           CD4 Memory T cells    121904619.
ref <- 'Podoplanin' 
nei <- 'Sinusoidal cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top11.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'Sinusoidal cells' 
nei <- 'CD4 Memory T cells'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_top12.pdf'),
    height = 5, width = 7)
p
dev.off()


### Low variance -----------------------------------------------------------

# 1 Neutrophils/Monocytes Podoplanin          755749.
# 2 CD4 Memory T cells    Macrophages        1409475.
# 3 Neutrophils/Monocytes B cells, red pulp  1934633.
# 4 Macrophages           Podoplanin         2090387.
# 5 Neutrophils/Monocytes Ki67 proliferating 2202569.

ref <- 'Neutrophils/Monocytes' 
nei <- 'Podoplanin'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_bottom1.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'CD4 Memory T cells' 
nei <- 'Macrophages'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_bottom2.pdf'),
    height = 5, width = 7)
p
dev.off()

ref <- 'Neutrophils/Monocytes' 
nei <- 'B cells, red pulp'
p <- dats %>% 
  filter(reference == ref,
         neighbor == nei) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) + 
  scale_color_manual(values = c('p_vnkn_s_pkhl' = '#009440', 
                                'p_vnkn_s_xxcd' = '#009440', 
                                'p_zwnt_s_pbvn' = '#ffcb00',
                                'p_zwnt_s_fsld' = '#ffcb00',
                                'p_kgnj_s_ngpl' = '#302681',
                                'p_kgnj_s_ksfb' = '#302681'))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'spleen_trends_var_bottom3.pdf'),
    height = 5, width = 7)
p
dev.off()



## Plot spatial ------------------------------------------------------------

## Top 2 and top 4
## fsld, xxcd, ngpl
colors_subset <- c('Fol B cells' = 'yellow', 'Blood endothelial' = 'red', 
                   'Ki67 proliferating' = 'blue', 'other' = '#E6E6E6')

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
cells <- crawdad::toSF(pos = xxcd[,c("x", "y")],
                       celltypes = xxcd$celltypes) %>% 
  mutate(celltypes = case_when(celltypes %in% names(colors_subset) ~ celltypes,
                               T ~ 'other'))
ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells, ofInterest = names(colors_subset), 
                 alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'xxcd_spatial.pdf'),
    height = 7, width = 12)
p
dev.off()



fsld <- read.csv2(file = '../CRAWDAD/data/spleen/FSLD.meta.csv.gz', row.names = 1)
cells <- crawdad::toSF(pos = fsld[,c("x", "y")],
                       celltypes = fsld$celltypes) %>% 
  mutate(celltypes = case_when(celltypes %in% names(colors_subset) ~ celltypes,
                               T ~ 'other'))
ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))
p <- vizClusters(cells, ofInterest = names(colors_subset), 
                 alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'fsld_spatial.pdf'),
    height = 7, width = 12)
p
dev.off()



ngpl <- read.csv2(file = '../CRAWDAD/data/spleen/NGPL.meta.csv.gz', row.names = 1)
cells <- crawdad::toSF(pos = ngpl[,c("x", "y")],
                       celltypes = ngpl$celltypes) %>% 
  mutate(celltypes = case_when(celltypes %in% names(colors_subset) ~ celltypes,
                               T ~ 'other'))
ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))
p <- vizClusters(cells, ofInterest = names(colors_subset), 
                 alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'ngpl_spatial.pdf'),
    height = 7, width = 12)
p
dev.off()