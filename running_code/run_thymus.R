library(crawdad)
library(tidyverse)
ncores <- 7


# vhck --------------------------------------------------------------------

vhck <- read.csv2(file = 'running_code/data/thymus/VHCK.meta.csv.gz', 
                  row.names = 1)
head(vhck)
paste(range(vhck$x), range(vhck$y))

ggplot(vhck, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(vhck$celltypes)))) +
  theme_void()

cells <- crawdad::toSF(pos = vhck[,c("x", "y")],
                       celltypes = vhck$celltypes)

scales <- seq(100, 1750, by=50)

# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 25.99 mins
# saveRDS(shuffle.list, 'running_code/processed_data/thymus/shufflelist_vhck_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/thymus/shufflelist_vhck_50.RDS')

## find trends, passing background as parameter
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 152.16 mins
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/thymus/dat_vhck_50.RDS')
dat <- readRDS('running_code/processed_data/thymus/dat_vhck_50.RDS')


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
saveRDS(ct_colors, 'running_code/processed_data/colors_thymus.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/thymus/vhck_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zsigThresh = zsig, zscoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  # scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/thymus/vhck_dotplot.pdf', height = 8, width = 6.7)
p
dev.off()

# pg <- ggplot_build(p)
# ct_order <- pg$layout$panel_params[[1]]$x$get_labels()
# saveRDS(ct_order, 'running_code/processed_data/ct_order_thymus.RDS')
ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')



## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
## selected CD4 Color '#FFFF00' and got 30% shades above and below
## https://maketintsandshades.com/#FFFF00
## the other color is the Fol B
colors_subset <- c('#b3b300', '#ffff4d', '#00FF00')

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
pdf('running_code/paper_figures/thymus/vhck_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ NA))
crawdad::vizClusters(cells = cells_subset)
p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset, na.value = '#E6E6E6') +
  theme_void()
p
pdf('running_code/paper_figures/thymus/vhck_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



# vhck --------------------------------------------------------------------

vhck <- read.csv2(file = 'running_code/data/thymus/VHCK.meta.csv.gz', 
                  row.names = 1)
head(vhck)
paste(range(vhck$x), range(vhck$y))

ggplot(vhck, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(vhck$celltypes)))) +
  theme_void()

cells <- crawdad::toSF(pos = vhck[,c("x", "y")],
                       celltypes = vhck$celltypes)

scales <- seq(100, 1750, by=50)

# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 25.99 mins
# saveRDS(shuffle.list, 'running_code/processed_data/thymus/shufflelist_vhck_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/thymus/shufflelist_vhck_50.RDS')

## find trends, passing background as parameter
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 152.16 mins
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/thymus/dat_vhck_50.RDS')
dat <- readRDS('running_code/processed_data/thymus/dat_vhck_50.RDS')


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
saveRDS(ct_colors, 'running_code/processed_data/colors_thymus.RDS')

p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/thymus/vhck_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

p <- vizColocDotplot(dat_filtered, zsigThresh = zsig, zscoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  # scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/thymus/vhck_dotplot.pdf', height = 8, width = 6.7)
p
dev.off()

# pg <- ggplot_build(p)
# ct_order <- pg$layout$panel_params[[1]]$x$get_labels()
# saveRDS(ct_order, 'running_code/processed_data/ct_order_thymus.RDS')
ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')



## Subsets
ct_ngb <- 'CD4 Memory T cells'
ct_ref <- 'Fol B cells'
ct_subset <- paste(ct_ngb, 'near', ct_ref, sep = '_')
## selected CD4 Color '#FFFF00' and got 30% shades above and below
## https://maketintsandshades.com/#FFFF00
## the other color is the Fol B
colors_subset <- c('#b3b300', '#ffff4d', '#00FF00')

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
pdf('running_code/paper_figures/thymus/vhck_subset_bar.pdf', height = 6, width = 2)
p
dev.off()


## visualize the subset only
cells_subset <- cells %>% 
  mutate(celltypes = case_when(row_number() %in% idx_ref_all ~ 'Reference', 
                               row_number() %in% idx_ngb_near ~ 'Near',
                               row_number() %in% idx_ngb_notnear ~ 'Not near',
                               T ~ NA))
crawdad::vizClusters(cells = cells_subset)
p <- vizClusters(cells_subset, alpha = 1, pointSize = .01) +
  scale_color_manual(values = colors_subset, na.value = '#E6E6E6') +
  theme_void()
p
pdf('running_code/paper_figures/thymus/vhck_subset_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



