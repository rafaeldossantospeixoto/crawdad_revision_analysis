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

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 25.74 mins

## find trends, passing background as parameter
## changed distance to 50
results <- crawdad::findTrends(cells,
                               dist = 50,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 107.72 mins

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
saveRDS(dat, 'running_code/processed_data/dat_pkhl_50.RDS')

zsig <- correctZBonferroni(dat)

## filter indistinct cells
dat_filtered <- dat %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct')

## reorder to be the same as the paper in biorxiv?
ct_order <- c('Podoplanin', 'CD4 Memory T cells', 'Fol B cells',
              'Macrophages', 'CD8 Memory T cells', 'Ki67 proliferating',
              'Myeloid cells', 'B cells, red pulp', 'Blood endothelial',
              'Sinusoidal cells', 'Neutrophils/Monocytes')
p <- vizColocDotplot(dat_filtered, zsig.thresh = zsig, zscore.limit = zsig*2, 
                     reorder = TRUE, dot.sizes = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')
p



# xxcd --------------------------------------------------------------------

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
head(xxcd)

cells <- crawdad::toSF(pos = xxcd[,c("x", "y")],
                       celltypes = xxcd$celltypes)

ggplot(xxcd, aes(x=x, y=y, col=celltypes)) + 
  geom_point(size=0.2, alpha=0.5) +
  scale_color_manual(values=rainbow(length(unique(xxcd$celltypes)))) +
  theme_void()

scales <- seq(100, 1750, by=50)

shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 24.92 mins
saveRDS(shuffle.list, 'tmp/shufflelist.RDS')
