library(crawdad)
# library(tidyverse)

# Compare neighborhoods ---------------------------------------------------


#' Calculate proportions of cell types inside the neighborhood
#' 
#' For a reference cell type and neighborhood distance, calculate the 
#' proportion of each neighbor cell type inside the neighborhood of the reference
#' cell type at that distance.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types annotated in the celltypes column and point positions in the 
#' geometry column
#' @param ref character; reference cell type to define neighborhood around
#' @param dist numeric vector; distances used to define the neighborhoods
#' 
#' @return named vector; 
calculateProportions <- function(cells, ref, dist) {
  ## create a circle around each reference cell 
  neighborhood <- sf::st_buffer(cells[cells$celltypes == ref,], dist) 
  ## merge the circles into a neighborhood (can take some time to compute)
  # neighborhood <- sf::st_union(buffer)
  ## calculate cells inside the neighborhood
  neighbor_cells <- sf::st_intersection(cells, neighborhood)
  ## remove duplicates
  self_cells <- cells[cells$celltypes == ref, ]
  neighbor_cells <- neighbor_cells[setdiff(rownames(neighbor_cells), 
                                           rownames(self_cells)), ]
  ## hack to accommodate self cells that are neighbors of another self cell
  neighbor_cells <- neighbor_cells[intersect(rownames(neighbor_cells), 
                                             c(rownames(cells), 
                                               paste0(rownames(self_cells), '.1'))), ]
  ## calculate proportions
  proportions <- as.vector(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
  names(proportions) <- names(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
  
  return(proportions)
}



#' Compare neighborhoods
#' 
#' For a reference cell type and different neighborhood distances, plot the 
#' proportion of each neighbor cell type inside the neighborhood of the reference
#' cell type at those distances.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types annotated in the celltypes column and point positions in the 
#' geometry column
#' @param dist numeric; distance used to define the neighborhood
#' @param dotSize numeric; size of the dot
#' 
#' @return ggplot2 plot; the proportion of each neighbor cell type for the 
#' reference cell types, given a distance
#' 
#' @export
plotProportions <- function(cells, dist, dotSize = 5) {
  
  ## for each cell type
  celltypes <- unique(cells$celltypes)
  props <- lapply(celltypes, calculateProportions, cells = cells, dist = dist)
  df <- dplyr::bind_rows(props) %>% 
    dplyr::mutate(reference = celltypes) %>% 
    tidyr::pivot_longer(!reference, names_to = 'neighbor', values_to = 'proportion')
  
  df %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = reference, y = neighbor, 
                                     color = proportion), size = dotSize) +
    ggplot2::scale_colour_gradient(low = 'white', high = '#006437',
                                   limits = c(0, 100)) +
    ggplot2::scale_x_discrete(position = 'top') +
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position = 'right',
                   axis.text.x = ggplot2::element_text(angle = 45, h = 0))
}





## Test --------------------------------------------------------------------

data('slide')
cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
ref <- 'Bergmann'
dist <- 50

calculateProportions(cells, ref, dist)

plotProportions(cells, dist = 10)
plotProportions(cells, dist = 50)
plotProportions(cells, dist = 100)
plotProportions(cells, dist = 250)



# Paper figures -----------------------------------------------------------


## Cerebellum --------------------------------------------------------------

data('slide')
cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)

for (d in c(10, 50, 100, 250)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'slide_', d, '.pdf'),
      height = 7, width = 8)
  print(p)
  dev.off()
}



## Embryo ------------------------------------------------------------------

data('seq')
cells <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)

for (d in c(10, 50, 100, 250)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'seq_', d, '.pdf'),
      height = 7, width = 8)
  print(p)
  dev.off()
}



## Sim ------------------------------------------------------------------

data('sim')
cells <- crawdad::toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)

for (d in c(10, 50, 100, 250)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'sim_', d, '.pdf'),
      height = 7, width = 8)
  print(p)
  dev.off()
}



## Pkhl ------------------------------------------------------------------

data('pkhl')
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)

for (d in c(10, 50, 100)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'pkhl_', d, '.pdf'),
      height = 7, width = 8)
  print(p)
  dev.off()
}








# Deprecated --------------------------------------------------------------

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