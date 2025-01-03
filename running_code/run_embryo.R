# Run embryo ---------------------------------------------------------

library(crawdad)
library(tidyverse)
library(gridExtra)
ncores <- 7

data(seq)
seq %>% ggplot() +
  geom_point(aes(x, y, color = celltypes), size = .5) + 
  facet_wrap('celltypes')

seq %>% 
  filter(celltypes %in% c('Presomitic mesoderm', 'Spinal cord')) %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltypes), size = .5)

## convert to sf
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                      cellTypes = seq$celltypes)
cells <- seq


## CRAWDAD -----------------------------------------------------------------

scales <- seq(100, 1000, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(seq,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(seq,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_seq_50.RDS')



# Exploring results -------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zSigThresh = zsig, 
                zScoreLimit = zsig*2, 
                dotSizes = c(2, 12), symmetrical = T) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



## Plot spat -------------------------------------------------------------


data(seq)
ref_ct <- unique(seq$celltypes)[1]
ngb_ct <- unique(seq$celltypes)[2]
pair <- c(ref_ct, ngb_ct)
for (ref_ct in unique(seq$celltypes)) {
  plots <- list()
  for (ngb_ct in unique(seq$celltypes)) {
    pair <- c(ref_ct, ngb_ct)
    plots[[ngb_ct]] <- seq %>% 
      mutate(selected = case_when(celltypes %in% pair ~ celltypes, 
                                  T ~ 'other')) %>% 
      ggplot() +
      geom_point(aes(x, y, color = selected), size = .025) +
      scale_color_manual(values = c('red', 'blue', 'lightgray') %>% 
                           `names<-`(c(ref_ct, ngb_ct, 'other')))
  }
  g <- gridExtra::arrangeGrob(grobs = plots)
  ggsave(paste0('running_code/exploring_embryo/spat_ref_',
                gsub('[ /]', '_', ref_ct), '.png'),
         g,
         height = 20, width = 30)
}


# Trends ------------------------------------------------------------------

## convert to sf
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                      celltypes = seq$celltypes)

## convert to sf
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                      celltypes = seq$celltypes)

all_cts <- unique(dat_50$reference)

for (ref_ct in all_cts) {
  ## spat
  p1 <- vizClusters(seq, ofInterest = c(ref_ct), alpha = 1) + 
    theme_void() +
    theme(legend.position = "none")
    
  
  ## trends
  sig_cts <- dat_50 %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(abs(Z) > zsig) %>% 
    pull(neighbor) %>% 
    unique() %>% 
    as.character()
  not_sig_cts <- all_cts[! all_cts %in% sig_cts]
  
  filtered_dat <- dat_50 %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(neighbor %in% sig_cts)
  # mutate(selected_neighbor = fct_other(neighbor, keep = sig_cts)) %>% 
  p2 <- filtered_dat %>% ggplot() +
    geom_line(aes(scale, Z, color = neighbor)) +
    scale_color_manual(values = rainbow(length(sig_cts))) +
    geom_hline(yintercept = zsig, linetype = 2) + 
    geom_hline(yintercept = -zsig, linetype = 2) +
    xlim(c(50, 1250)) + 
    ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                                 label = neighbor), 
                             data = filtered_dat %>%
                               filter(scale == max(scale)),
                             hjust = 'right', box.padding = 0.5,
                             nudge_x = .5, nudge_y = .5,
                             segment.alpha = .5,
                             xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
    labs(title = ref_ct)
  
  p <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1/3, 2/3))
  
  print(p)
  # png(paste0('running_code/exploring_embryo/trends_ref_', 
  #            gsub('[ /]', '_', ref_ct), '.png'),
  #     height = 700, width = 1400)
  # print(p)
  # dev.off()
}


## Dotplot -----------------------------------------------------------------


p <- vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(1, 8),
                     mutual = T) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
pdf('running_code/exploring_embryo/embryo_dotplot_crawdad.pdf',
    height = 7.5, width = 9)
p
dev.off()

  


# Paper figures -----------------------------------------------------------

set.seed(42)
dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')

zsig <- correctZBonferroni(dat_50)
p <- vizColocDotplot(dat_50, reorder = TRUE, zSigThresh = zsig, 
                zScoreLimit = zsig*2, dotSizes = c(1, 10)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
pg <- ggplot_build(p)
ct_order <- pg$layout$panel_params[[1]]$x$get_labels()

dat_50 <- dat_50 %>% 
  filter(scale <= 750)
vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
                zscoreLimit = zsig*2, mutual = T, dotSizes = c(1, 10)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  # scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

## Dotplot -----------------------------------------------------------------


p <- dat_50 %>% 
  filter(neighbor != 'Low quality',
         reference != 'Low quality') %>% 
  vizColocDotplot(reorder = TRUE, zsigThresh = zsig, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(1, 8)) +
  coord_fixed() + 
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
pdf('running_code/paper_figures/embryo/dotplot_crawdad.pdf',
    height = 7.5, width = 9)
p
dev.off()


## Spatial plot ------------------------------------------------------------

## flip to look like in the paper
data(seq)
seq$y <- -seq$y
cells <- crawdad:::toSF(pos = seq[,c("x", "y")],
                        celltypes = seq$celltypes)

## tried to reproduce the same colors as in the preprint paper, but could not
## so I decided to manually select the color of the specific cell types in 
## illustrator and pass them as argument of the plot function
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'Low quality']))
# ct_colors <- setNames(tail(SteppedSequential5Steps, length(interest_cts)), 
#                       interest_cts) 
# ct_colors <- setNames(sample(rainbow(length(interest_cts))), interest_cts)
# saveRDS(ct_colors, 'running_code/processed_data/colors_seq.RDS')
ct_colors <- readRDS('running_code/processed_data/colors_seq.RDS')

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/embryo/spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



### Selected cell types -----------------------------------------------------

## Endothelium
interest_cts <- c('Endothelium', 
                  'Haematoendothelial progenitors',
                  'Forebrain/Midbrain/Hindbrain')
## these colors were alpha = .5 in the paper, how to convert them back?
interest_ct_colors <- ct_colors[interest_cts]
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/embryo/spatial_plot_endothelium.pdf',
    height = 7, width = 12)
p
dev.off()
  
## Intermediate
interest_cts <- c('Intermediate mesoderm', 
                  'Lateral plate mesoderm',
                  'Cardiomyocytes')
## these colors were alpha = .5 in the paper, how to convert them back?
interest_ct_colors <- ct_colors[interest_cts]
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  coord_fixed() + 
  theme_minimal()
p
pdf('running_code/paper_figures/embryo/spatial_plot_intermediatemesoderm.pdf',
    height = 7, width = 12)
p
dev.off()



## Trend plot --------------------------------------------------------------

p <- dat_50 %>% 
  filter(reference == 'Lateral plate mesoderm') %>% 
  filter(neighbor == 'Intermediate mesoderm') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig) 
p
pdf('running_code/paper_figures/embryo/trend_refLateral_neighIntermediate_crawdad.pdf',
    height = 4, width = 6)
p
dev.off()

p <- dat_50 %>% 
  filter(reference == 'Cardiomyocytes') %>% 
  filter(neighbor == 'Intermediate mesoderm') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig)
p
pdf('running_code/paper_figures/embryo/trend_refCardiomyocytes_neighIntermediate_crawdad.pdf',
    height = 4, width = 6)
p
dev.off()


# Verify for proofing -------------------------------------------------------

dat_50 %>% 
  mutate(id = 'embryo') %>% 
  filter(reference == 'Intermediate mesoderm') %>% 
  filter(neighbor %in% c('Lateral plate mesoderm', 'Cardiomyocytes')) %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, zSigThresh = zsig)
