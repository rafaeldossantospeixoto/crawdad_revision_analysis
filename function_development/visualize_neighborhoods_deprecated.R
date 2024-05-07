library(crawdad)
library(tidyverse)
## load the spleen data of the pkhl sample 
data('pkhl')
max(pkhl$x) - min(pkhl$x) # 3550.238
max(pkhl$y) - min(pkhl$y) # 3423.43
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)


# Visualize neighborhood --------------------------------------------------

vizAllClusters(pkhl, coms = pkhl$celltypes)

#' Visualize Neighborhood
#' 
#' Given a reference cell type and a neighborhood distance, plot the 
#' neighborhood using the spatial coordinates and also show the percentage of 
#' each cell type within the neighborhood.
#' 
#' @param cells sf data.frame; as produced by crawdad::toSF function: cells with 
#' cell types anotated in the celltypes column and point positions in the 
#' geometry column
#' @param reference character; reference cell type to define neihborhood around
#' @param distance numeric; distance used to define the neighborhood
#' 
vizNeighborhood <- function(cells, reference, distance){
  message('This function might take a minute to calculate the proportions')
  
  ## create a circle around each reference cell 
  buffer <- sf::st_buffer(cells[cells$celltypes == reference,], distance) 
  ## merge the circles into a neighborhood (can take some time to compute)
  neighborhood <- sf::st_union(buffer)
  ## calculate cells inside the neighborhood
  neighbor_cells <- sf::st_intersection(cells, neighborhood)
  ## calculate proportions
  proportions <- as.vector(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
  names(proportions) <- names(round(100*table(neighbor_cells$celltypes)/table(cells$celltypes), 2))
  ## create labels with proportions
  proportion_labels <- paste0(format(proportions, nsmall = 2), '% ', names(proportions))
  names(proportion_labels)  <- names(proportions)
  
  ## plot neighborhood and cells
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = cells, ggplot2::aes(color = celltypes), size = .01) + 
    ggplot2::scale_color_manual(values = grDevices::rainbow(length(unique(cells$celltypes))),
                                labels = proportion_labels) +
    ggplot2::geom_sf(data = neighborhood, fill = NA, color = 'black', linewidth = 1) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_minimal()
}

reference <- 'Podoplanin'
distance <- 100
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 39.80471 secs

reference <- 'Podoplanin'
distance <- 50
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 59.10783 secs

reference <- 'Fol B cells'
distance <- 100
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 1.343152 mins

reference <- 'Fol B cells'
distance <- 50
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 4.432859 mins

reference <- 'CD8 Memory T cells'
distance <- 100
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 20.93948 secs

reference <- 'CD8 Memory T cells'
distance <- 50
it <- Sys.time()
vizNeighborhood(cells, reference, distance)
Sys.time() - it
## Time difference of 1.51143 mins
