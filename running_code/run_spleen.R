library(crawdad)
library(tidyverse)
ncores <- 7


# pkhl --------------------------------------------------------------------

data('pkhl')
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")],
                       celltypes = pkhl$celltypes)

ggplot(pkhl, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(pkhl$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)

# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 25.74 mins in my pc
# ## Time was 70.47 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_pkhl_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_pkhl_50.RDS')

## find trends, passing background as parameter
## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 107.72 mins in my pc
# ## Time was 104.3 mins in Easley with 28 cores
# ## The number of cores does not seem to help much
# 
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_pkhl_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_pkhl_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.68mins
# head(binomMat)
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_pkhl_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_pkhl_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_pkhl_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_pkhl_50.RDS')

cells$celltypes <- as.factor(cells$celltypes)
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_pkhl_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_pkhl_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- c('Sinusoidal cells' = '#FF0080',
               'Myeloid cells' = '#0000FF',
               'Neutrophils/Monocytes' = '#8000FF',
               'Blood endothelial' = '#FF8000',
               'CD8 Memory T cells' = '#80FF00',
               'Macrophages' = '#0080FF',
               'Fol B cells' = '#00FF00',
               'B cells, red pulp' = '#FF0000',
               'Ki67 proliferating' = '#00FFFF',
               # 'indistinct' = '#00FF80',
               'CD4 Memory T cells' = '#FFFF00',
               'Podoplanin' = '#FF00FF')
saveRDS(ct_colors, 'running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/pkhl_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0,
                       breaks = c(-zsig*2, 0, zsig*2),
                       limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/pkhl_dotplot.pdf', height = 8, width = 6.7)
p
dev.off()

# pg <- ggplot_build(p)
# ct_order <- pg$layout$panel_params[[1]]$x$get_labels()
# saveRDS(ct_order, 'running_code/processed_data/ct_order_spleen.RDS')
ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')



## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
## got the colors from the umap image and corrected the 50% opacy
## details in the xxcd figure
## the other color is the Fol B
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/pkhl_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                  'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/pkhl_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()





# xxcd --------------------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
head(xxcd)
range(xxcd$x)
range(xxcd$y)

cells <- crawdad::toSF(pos = xxcd[,c("x", "y")],
                       celltypes = xxcd$celltypes)

ggplot(xxcd, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(xxcd$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 66.27 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_xxcd_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_xxcd_50.RDS')

## find trends, passing background as parameter
## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 101.13 mins in Easley with 28 cores
# ## The number of cores does not seem to help much
# 
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_xxcd_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_xxcd_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.64mins
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_xxcd_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_xxcd_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_xxcd_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_xxcd_50.RDS')

##
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_xxcd_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_xxcd_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/xxcd_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        breaks = c(-zsig*2, 0, zsig*2),
                        limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'bottom') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  coord_fixed() + 
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 1),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/xxcd_dotplot.pdf', height = 8, width = 6.7)
p
dev.off()





## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
## selected CD4 Color '#ff0080' and got 30% shades above and below
## https://maketintsandshades.com/#ff0080
## the other color is the Fol B
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/xxcd_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                   'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/xxcd_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()









# fsld --------------------------------------------------------------------

fsld <- read.csv2(file = '../CRAWDAD/data/spleen/FSLD.meta.csv.gz', row.names = 1)
head(fsld)
range(fsld$x)
range(fsld$y)

cells <- crawdad::toSF(pos = fsld[,c("x", "y")],
                       celltypes = fsld$celltypes)

ggplot(fsld, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(fsld$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
## Time was 66.27 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_fsld_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_fsld_50.RDS')

## find trends, passing background as parameter
## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
## Time was 101.13 mins in Easley with 28 cores
## The number of cores does not seem to help much

# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_fsld_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_fsld_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.64mins
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_fsld_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_fsld_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_fsld_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_fsld_50.RDS')

# ##
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_fsld_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_fsld_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/fsld_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        breaks = c(-zsig*2, 0, zsig*2),
                        limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  coord_fixed() + 
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/fsld_dotplot.pdf', height = 6.7, width = 8)
p
dev.off()





## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/fsld_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                   'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/fsld_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()




# pbvn --------------------------------------------------------------------

pbvn <- read.csv2(file = '../CRAWDAD/data/spleen/PBVN.meta.csv.gz', row.names = 1)
head(pbvn)
range(pbvn$x)
range(pbvn$y)

cells <- crawdad::toSF(pos = pbvn[,c("x", "y")],
                       celltypes = pbvn$celltypes)

ggplot(pbvn, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(pbvn$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 66.27 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_pbvn_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_pbvn_50.RDS')

## find trends, passing background as parameter
## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
## Time was 101.13 mins in Easley with 28 cores
## The number of cores does not seem to help much

# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_pbvn_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_pbvn_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.64mins
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_pbvn_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_pbvn_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_pbvn_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_pbvn_50.RDS')

##
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_pbvn_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_pbvn_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/pbvn_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        breaks = c(-zsig*2, 0, zsig*2),
                        limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  coord_fixed() + 
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/pbvn_dotplot.pdf', height = 6.7, width = 8)
p
dev.off()





## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/pbvn_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                   'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/pbvn_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()





# ksfb --------------------------------------------------------------------

ksfb <- read.csv2(file = '../CRAWDAD/data/spleen/KSFB.meta.csv.gz', row.names = 1)
head(ksfb)
range(ksfb$x)
range(ksfb$y)

cells <- crawdad::toSF(pos = ksfb[,c("x", "y")],
                       celltypes = ksfb$celltypes)

ggplot(ksfb, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(ksfb$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 66.27 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_ksfb_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_ksfb_50.RDS')

# ## find trends, passing background as parameter
# ## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 101.13 mins in Easley with 28 cores
# ## The number of cores does not seem to help much
# 
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_ksfb_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_ksfb_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 25.64mins
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_ksfb_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_ksfb_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_ksfb_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_ksfb_50.RDS')

# ##
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_ksfb_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_ksfb_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/ksfb_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        breaks = c(-zsig*2, 0, zsig*2),
                        limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  coord_fixed() + 
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/ksfb_dotplot.pdf', height = 6.7, width = 8)
p
dev.off()





## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/ksfb_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                   'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/ksfb_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()





# ngpl --------------------------------------------------------------------

ngpl <- read.csv2(file = '../CRAWDAD/data/spleen/NGPL.meta.csv.gz', row.names = 1)
head(ngpl)
range(ngpl$x)
range(ngpl$y)

cells <- crawdad::toSF(pos = ngpl[,c("x", "y")],
                       celltypes = ngpl$celltypes)

ggplot(ngpl, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(ngpl$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 66.27 mins in Easley with 28 cores
# saveRDS(shuffle.list, 'running_code/processed_data/spleen/shufflelist_ngpl_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/spleen/shufflelist_ngpl_50.RDS')

## find trends, passing background as parameter
## changed distance to 50
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 180 mins in Easley with 7 cores
# ## The number of cores does not seem to help much
# 
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/spleen/dat_ngpl_50.RDS')
dat <- readRDS('running_code/processed_data/spleen/dat_ngpl_50.RDS')





## Subset analysis ---------------------------------------------------------

## changed neighborhood to 50
# binomMat <- crawdad::binomialTestMatrix(cells,
#                                         neigh.dist = 50,
#                                         ncores = ncores,
#                                         verbose = TRUE)
# ## Time to compute was 104.53mins in E7
# saveRDS(binomMat, file = 'running_code/processed_data/spleen/binomMat_ngpl_50.RDS')
binomMat <- readRDS('running_code/processed_data/spleen/binomMat_ngpl_50.RDS')

# subset.list <- crawdad::selectSubsets(binomMat,
#                                       cells$celltypes,
#                                       sub.type = "near",
#                                       sub.thresh = 0.05,
#                                       ncores = ncores,
#                                       verbose = TRUE)
# ## Time to compute was 0.17mins
# saveRDS(subset.list, file = 'running_code/processed_data/spleen/subsetlist_ngpl_50.RDS')
subset.list <- readRDS('running_code/processed_data/spleen/subsetlist_ngpl_50.RDS')

# ##
# results.subsets <- crawdad::findTrends(cells,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        subset.list = subset.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE)
# ## Time was 715.76 mins
# ## subsets
# dats <- crawdad::meltResultsList(results.subsets, withPerms = TRUE)
# saveRDS(dats, file = 'running_code/processed_data/spleen/dats_ngpl_50.RDS')
dats <- readRDS('running_code/processed_data/spleen/dats_ngpl_50.RDS')

## Multiple-test correction
ntestss <- length(unique(dats$reference)) * length(unique(dats$neighbor))
psigs <- 0.05/ntestss
zsigs <- round(qnorm(psigs/2, lower.tail = F), 2)



## Paper figures -----------------------------------------------------------

ct_order <- readRDS('running_code/processed_data/ct_order_spleen.RDS')

## Spatial plot
all_cts <- unique(cells$celltypes)
interest_cts <- sort(as.character(all_cts[all_cts != 'indistinct']))
ct_colors <- readRDS('running_code/processed_data/colors_spleen.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/spleen/ngpl_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zSigThresh = zsig, zScoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0,
                        breaks = c(-zsig*2, 0, zsig*2),
                        limits=c(-zsig*2, zsig*2)) + 
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = rev(c(1750, 1000, 500, 100)),
                        limits = c(1750, 100),
                        range = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  coord_fixed() + 
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/spleen/ngpl_dotplot.pdf', height = 6.7, width = 8)
p
dev.off()





## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
colors_subset <- c('#b3005a', '#ff4da6', '#00FF80')

idx_ref_all <- as.character(which(cells$celltypes == ct_ref))
idx_ngb_all <- as.character(which(cells$celltypes == ct_ngb))
idx_ngb_near <- subset.list[[ct_subset]]
idx_ngb_notnear <- setdiff(idx_ngb_all, idx_ngb_near)


## bar plot
## percentage of subsets
df_pct <- data.frame(near = c(length(idx_ngb_near)/length(idx_ngb_all)),
                     notnear = c(length(idx_ngb_notnear)/length(idx_ngb_all)),
                     celltypes = c(ct_ngb)) %>% 
  pivot_longer(!celltypes, names_to = 'condition')

## bar
p <- df_pct %>% ggplot() + 
  geom_bar(aes(y=value, x=celltypes, fill=condition),
           position="fill", stat="identity") +
  scale_fill_manual(values = colors_subset[1:2]) +
  geom_text(aes(label = round(value, 2), 
                x = celltypes, 
                y = c(value[2] + value[1]/2,
                      value[2]/2)),
            color = 'white', size = 10) +
  theme_void()  + 
  theme(legend.position="none") +
  annotate("text", x=1, y=-.035, label = length(idx_ngb_all),
           size = 10)
p
pdf('running_code/paper_figures/spleen/ngpl_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ 'other'))

colors_subset <- c('Near' = '#b3005a', 'Not near' = '#ff4da6', 
                   'Reference' = '#00FF80', 'other' = '#E6E6E6')

ordered_cts <- names(sort(table(cells_subset$celltypes), decreasing = T))
cells_subset <- cells_subset %>% 
  arrange(match(celltypes, ordered_cts))


p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset) +
  theme_void()
p
pdf('running_code/paper_figures/spleen/ngpl_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()

