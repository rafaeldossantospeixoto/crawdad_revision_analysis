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
ct_colors <- c(# 'indistinct' = '#00FF80',
               'Blood endothelial' = '#FF8000',
               'cortical and medullary thymocytes' = 'lightblue',
               'cortical thymocytes/epithelial cells' = 'salmon',
               'cortical thymocytes/DP T cells' = 'yellow',
               "Hassal's corpuscles" = 'darkred',
               "Medullar thymic epithelial cells/B cells" = 'blue',
               "thymic epithelial cells/B cells" = 'darkviolet',
               "B cells/CD4 T cells" = 'darkblue' 
               )
saveRDS(ct_colors, 'running_code/processed_data/colors_thymus.RDS')

ordered_cts <- names(sort(table(cells$celltypes), decreasing = T))
cells <- cells %>% 
  arrange(match(celltypes, ordered_cts))

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



## Checking symmetry -------------------------------------------------------

sort(table(cells$celltypes))

vizColocDotplot(dat_filtered, zsigThresh = zsig, zscoreLimit = zsig*2,
                reorder = TRUE, mutual = T, dotSizes = c(2, 14)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')

vhck %>% 
  filter(celltypes == c("cortical thymocytes/epithelial cells",
                        'thymic epithelial cells/B cells')) %>%
  ggplot(aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.5) +
  scale_color_manual(values=c('red', 'blue')) +
  theme_void()





# ktjk --------------------------------------------------------------------

ktjk <- read.csv2(file = 'running_code/data/thymus/KTJK.meta.csv.gz', 
                  row.names = 1)
head(ktjk)
paste(range(ktjk$x), range(ktjk$y))

ggplot(ktjk, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(ktjk$celltypes)))) +
  theme_void()

cells <- crawdad::toSF(pos = ktjk[,c("x", "y")],
                       celltypes = ktjk$celltypes)

scales <- seq(100, 1750, by=50)

# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# ## Time was 46.95 mins
# saveRDS(shuffle.list, 'running_code/processed_data/thymus/shufflelist_ktjk_50.RDS')
shuffle.list <- readRDS('running_code/processed_data/thymus/shufflelist_ktjk_50.RDS')

## find trends, passing background as parameter
# results <- crawdad::findTrends(cells,
#                                dist = 50,
#                                shuffle.list = shuffle.list,
#                                ncores = ncores,
#                                verbose = TRUE,
#                                returnMeans = FALSE)
# ## Time was 354.76 mins
# dat <- crawdad::meltResultsList(results, withPerms = TRUE)
# saveRDS(dat, 'running_code/processed_data/thymus/dat_ktjk_50.RDS')
dat <- readRDS('running_code/processed_data/thymus/dat_ktjk_50.RDS')


## Paper figures -----------------------------------------------------------

ct_colors <- saveRDS('running_code/processed_data/colors_thymus.RDS')
p <- vizClusters(cells, ofInterest = interest_cts, alpha = 1, pointSize = .01) +
  scale_color_manual(values = ct_colors, na.value = '#00FF80') +
  theme_void()
p
pdf('running_code/paper_figures/thymus/ktjk_spatial_plot.pdf',
    height = 7, width = 12)
p
dev.off()



## Dotplot
zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')
p <- vizColocDotplot(dat_filtered, zsigThresh = zsig, zscoreLimit = zsig*2, 
                     reorder = TRUE, dotSizes = c(2, 14)) +
  # scale_x_discrete(limits = ct_order, position = 'top') +
  # scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p
pdf('running_code/paper_figures/thymus/ktjk_dotplot.pdf', height = 8, width = 6.7)
p
dev.off()



