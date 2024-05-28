#' Visualize clusters in the same plot
#' 
#' @description Uses the cells sf object to visualize the clusters together in 
#' the same plot.
#' 
#' @param cells sf object; spatial (x and y) coordinates and celltypes column
#' @param ofInterest character vector; a vector of specific clusters to visualize
#' @param pointSize numeric; size of points
#' @param alpha numeric; transparency of points
#' @param ref character; reference cell type to draw the neighborhood 
#' around. If NULL, it will not create the neighborhood (default: NULL) 
#' @param dist numeric; distance to define neighbor cells with respect to each 
#' reference cell. If NULL, it will not create the neighborhood (default: NULL)
#' @param lineWidth numeric; width of neighborhood line
#' 
#' @return plot
#' 
#' @examples
#' \dontrun{
#' data(slide)
#' cells <- crawdad::toSF(pos = slide[,c("x", "y")], celltypes = slide$celltypes)
#' vizClusters(cells)
#' }
#' 
#' @export
vizClusters <- function(cells, ofInterest = NULL,
                        pointSize = 1, alpha = 0.5,
                        ref = NULL, dist = NULL, lineWidth = 0.1){
  
  ## if cells are a data.frame with "x" and "y" cell coordinate columns
  if( class(cells)[1] %in% c("data.frame", "matrix") ){
    stop('Use an sf object created by the crawdad::toSF function.')
  }
  
  ## define colors
  cluster_cols <- rainbow(n = length(unique(cells$celltypes)))
  names(cluster_cols) <- unique(cells$celltypes)
  cluster_cols['other'] <- '#E6E6E6'
  
  ## separate cells of interest to plot on top of others
  if(!is.null(ofInterest)){
    cells <- cells %>% 
      dplyr::mutate(celltypes = 
                      dplyr::case_when((!cells$celltypes %in% ofInterest) ~ 'other',
                                       T ~ celltypes))
    # cells$celltypes <- droplevels(cells$celltypes)
  }
  
  ## order cell types based on abundance
  ordered_cts <- c(names(sort(table(cells$celltypes), decreasing = T)))
  ordered_cts <- c('other', ordered_cts[ordered_cts != 'other'])
  cells <- cells %>% 
    dplyr::arrange(match(celltypes, ordered_cts))
  
  ## plot
  plt <- ggplot2::ggplot() +
    ## plot other cells
    ggplot2::geom_sf(data = cells, 
                     ggplot2::aes(color = celltypes), 
                     size = pointSize, alpha = alpha) +
    ## NA to gray
    ggplot2::scale_color_manual(values = cluster_cols, na.value = "#E6E6E6")
  plt
  
  if( (!is.null(ref)) & (!is.null(dist)) ) {
    ## create a circle around each reference cell 
    buffer <- sf::st_buffer(cells[cells$celltypes == ref,], dist) 
    ## merge the circles into a neighborhood (can take some time to compute)
    neighborhood <- sf::st_union(buffer)
    ## add to plot
    plt <- plt +
      ggplot2::geom_sf(data = neighborhood, fill = NA, 
                       color = 'black', linewidth = lineWidth)
  }
  
  ## add labels
  plt <- plt + 
    ggplot2::labs(x = "x",
                  y = "y") +
    ggplot2::theme_minimal() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=2), 
                                                   ncol = 2))
    # ggplot2::coord_equal() ## geom_sf seems to be equal already
  
  return(plt)
  
}




# Testing -----------------------------------------------------------------

data(slide)
cells <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)
ct_colors <- readRDS('running_code/processed_data/colors_slide.RDS')

## UBCs and Granule
interest_cts <- c('UBCs', 
                  'Granule')
interest_ct_colors <- ct_colors[interest_cts]
vizClusters(cells, ofInterest = interest_cts, ref = 'UBCs', dist = 50,
            alpha = 1, pointSize = 1) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()
vizClusters(cells, ofInterest = interest_cts, ref = 'Granule', dist = 50,
            alpha = 1, pointSize = 1) + 
  scale_color_manual(values = interest_ct_colors, na.value = '#E6E6E6') +
  theme_minimal()

