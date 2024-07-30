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
  proportions <- as.vector(round(100*(table(neighbor_cells$celltypes)) / 
                                   (table(cells$celltypes)), 2))
  names(proportions) <- names(round(100*table(neighbor_cells$celltypes) / 
                                      table(cells$celltypes), 2))
  
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
plotProportions <- function(cells, dist) {
  
  ## for each cell type
  celltypes <- unique(cells$celltypes)
  props <- lapply(celltypes, calculateProportions, cells = cells, dist = dist)
  df <- data.frame(proportions = unlist(props))
  
  df %>% ggplot2::ggplot(ggplot2::aes(x = proportions)) + 
    ggplot2::geom_histogram(color='#006437', fill='white', bins = 100) +
    ggplot2::theme_bw()
  
}



## Test --------------------------------------------------------------------

data('slide')
ct_colors <- readRDS('running_code/processed_data/colors_slide.RDS')
cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)

ref <- 'Bergmann'
dist <- 50
calculateProportions(cells, ref, dist)
plotProportions(cells, dist)

plotProportions(cells, 10)
plotProportions(cells, 50)
plotProportions(cells, 75)
plotProportions(cells, 100)
plotProportions(cells, 250)



# Paper figures -----------------------------------------------------------



## Cerebellum --------------------------------------------------------------

data('slide')
cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)

for (d in c(10, 50, 100)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'histograms/',
             'slide_', d, '.pdf'),
      height = 5, width = 7)
  print(p)
  dev.off()
}



## Embryo ------------------------------------------------------------------

data('seq')
cells <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)

for (d in c(10, 50, 100)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'histograms/',
             'seq_', d, '.pdf'),
      height = 5, width = 7)
  print(p)
  dev.off()
}



## Sim ------------------------------------------------------------------

data('sim')
cells <- crawdad::toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)

for (d in c(10, 50, 100)){
  p <- plotProportions(cells, dist = d)
  print(p)
  pdf(paste0('function_development/compare_neighborhoods/paper_figures/',
             'histograms/',
             'sim_', d, '.pdf'),
      height = 5, width = 7)
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
             'histograms/',
             'pkhl_', d, '.pdf'),
      height = 5, width = 7)
  print(p)
  dev.off()
}



