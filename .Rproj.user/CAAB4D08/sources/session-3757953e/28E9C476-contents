library(crawdad)
library(tidyverse)
## load the spleen data of the pkhl sample 
data('pkhl')
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)


# Compare neighborhoods ---------------------------------------------------


#' Compare neighborhoods
#' 
#' For a reference cell type and different neighborhood distances, plot the 
#' proportion of each neigbor cell type inside the neighborhood of the reference
#' cell type at those distances.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types anotated in the celltypes column and point positions in the 
#' geometry column
#' @param reference character; reference cell type to define neihborhood around
#' @param distances numeric vector; distances used to define the neighborhoods
#' 
compareNeighborhoods <- function(cells, reference, distances) {
  df_distances <- data.frame()
  ## parallelize this loop
  for (distance in distances) {
    ## create a circle around each reference cell 
    buffer <- sf::st_buffer(cells[cells$celltypes == reference,], distance) 
    ## merge the circles into a neighborhood (can take some time to compute)
    neighborhood <- sf::st_union(buffer)
    ## calculate cells inside the neighborhood
    neighbor_cells <- sf::st_intersection(cells, neighborhood)
    ## calculate proportions
    proportions <- as.vector(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
    names(proportions) <- names(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
    ## create a dataframe with this information
    df_distance <- data.frame(celltypes = names(proportions), 
                              proportions = proportions,
                              distances = distance)
    df_distances <- rbind(df_distances, df_distance)
  }
  df_distances %>% 
    dplyr::filter(celltypes != reference) %>% 
    ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = distances, y = proportions, 
                                    color = celltypes)) +
    ggplot2::scale_color_manual(values = grDevices::rainbow(length(unique(cells$celltypes)) - 1)) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 2))) + 
    ggplot2::theme_minimal()
}

## Podoplanin
reference <- 'Podoplanin'
distances <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
it <- Sys.time()
plt_distances <- compareNeighborhoods(cells, reference, distances)
Sys.time() - it
## Time difference of 7.18728 mins
plt_distances

## use geom_text_repel
end_points <- df_distances %>%
  group_by(celltypes) %>%
  summarize(distances = last(distances), proportions = last(proportions))
plt_distances + 
  ggrepel::geom_text_repel(data = end_points, 
                           ggplot2::aes(x = distances, y = proportions, label = celltypes), 
                           box.padding = 0.5)

## Fol B cells
reference <- 'Fol B cells'
distances <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
it <- Sys.time()
plt_distances <- compareNeighborhoods(cells, reference, distances)
Sys.time() - it
## Time difference of 27.20221 mins
plt_distances

## CD8 Memory T cells
reference <- 'CD8 Memory T cells'
distances <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
it <- Sys.time()
plt_distances <- compareNeighborhoods(cells, reference, distances)
Sys.time() - it
## Time difference of 26.67759 mins
plt_distances

## Sinusoidal cells
reference <- 'Sinusoidal cells'
distances <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
it <- Sys.time()
plt_distances <- compareNeighborhoods(cells, reference, distances)
Sys.time() - it
## Time difference of 21.73519 mins
plt_distances


# Run CRAWDAD -------------------------------------------------------------

ncores <- 7
cells_pkhl <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
scales <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
shuffle.list_pkhl <- crawdad:::makeShuffledCells(cells_pkhl,
                                                 scales = scales,
                                                 perms = 3,
                                                 ncores = ncores,
                                                 seed = 1,
                                                 verbose = TRUE)
## Time was 3.2 mins
## large list 333.6 MB


## 100 ---------------------------------------------------------------------

results_pkhl_100 <- crawdad::findTrends(cells_pkhl,
                               dist = 100,
                               shuffle.list = shuffle.list_pkhl,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
## Time was 38.81 mins
saveRDS(results_pkhl_100, 'function_development/spleen/results_pkhl_100.RDS')

dat <- crawdad::meltResultsList(results_pkhl_100, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))


## 100 pixels --------------------------------------------------------------

new_dist <- 100 * (377.4038462 / 1000)
results_pkhl_100 <- crawdad::findTrends(cells_pkhl,
                                        dist = new_dist,
                                        shuffle.list = shuffle.list_pkhl,
                                        ncores = ncores,
                                        verbose = TRUE,
                                        returnMeans = FALSE)
## Time was too long, I gave up
saveRDS(results_pkhl_100, 'function_development/spleen/results_pkhl_37.RDS')

