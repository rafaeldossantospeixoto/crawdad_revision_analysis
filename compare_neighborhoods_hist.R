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
  # proportions <- as.vector(round(100*sum(table(neighbor_cells$celltypes)) / 
  #                                  sum(table(cells$celltypes)), 2))
  proportions <- as.vector(round(100*(table(neighbor_cells$celltypes)) / 
                                   (table(cells$celltypes)), 2))
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
    ggplot2::scale_x_discrete(position = 'bottom') +
    ggplot2::theme_bw() + 
    ggplot2::theme(legend.position = 'right',
                   axis.text.x = ggplot2::element_text(angle = 90, h = 1)) + 
    ggplot2::coord_equal()
}





## Test --------------------------------------------------------------------

data('slide')
ct_colors <- readRDS('running_code/processed_data/colors_slide.RDS')
cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)

ref <- 'Bergmann'
dist <- 50
calculateProportions(cells, ref, dist)

cts <- unique(slide$celltypes)
plotHistPct <- function(cells, references, dist) {
  proportions <- lapply(references, function(ct) {
    all_props <- calculateProportions(cells, ct, dist)
    props <- all_props[names(all_props) != ct]
    as.numeric(props)
  })
  print(unlist(proportions))
  hist(unlist(proportions), breaks = 100, main = dist)
}

plotHistPct(cells, cts, dist)


plotHistPct(cells, cts, 10)
plotHistPct(cells, cts, 50)
plotHistPct(cells, cts, 75)
plotHistPct(cells, cts, 100)
plotHistPct(cells, cts, 250)
