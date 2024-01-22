library(crawdad)
library(tidyverse)

# Create dataset --------------------------------------------------------

## similar to a real dataset (spleen)
## spleen is 3.5 x 3.5 mm, with 150k cells and structures of .5 mm
tissue_size <- 4 * 1e3 ## 4 mm or 4000 microns
circle_diameter <- 1e3
## number of cells is estimated based on sleen density
## spleen has 45k per 1 mm in size
number_of_cells <- 45 * tissue_size

centroids <- list(c(tissue_size/4, tissue_size/4),
                  c(tissue_size*3/4, tissue_size/4),
                  c(tissue_size*3/4, tissue_size*3/4),
                  c(tissue_size/4, tissue_size*3/4),
                  c(tissue_size/2, tissue_size/2))

simulate_symmetries <- function(tissue_size, number_of_cells, circle_diameter, 
                               centroids) {
  
  dfs <- list()
  
  number_of_dfs <- length(centroids)
  ## creates a list of Boolean vectors that will decide whether a circle should
  ## have both cell types or not
  both_celltypes <- lapply(seq_len(number_of_dfs), function(x) {
    c(rep(TRUE, times = total_length - x + 1), rep(FALSE, times = x - 1))
  })
  
  for (ndf in 1:number_of_dfs) {
    ## define background
    df <- crawdad:::simulate_background(size = number_of_cells, scale = tissue_size)
    
    for (ncd in 1:length(centroids)) {
      ## simulate circle
      x_center <- centroids[[ncd]][1]
      y_center <- centroids[[ncd]][2]
      radius <- circle_diameter/2
      cells_circle <- rownames(df[((df$x-x_center)^2 + (df$y - y_center)^2 < radius^2),])
      df[cells_circle, 'type'] <- 'B'
      if (both_celltypes[[ndf]][ncd]) {
        cells_type_c <- sample(cells_circle, 
                               size = round(length(cells_circle) * .5), ## 50 
                               replace = FALSE)
        df[cells_type_c, 'type'] <- 'C'
      }
    }
    
    colnames(df) <- c('x', 'y', 'celltype')
    dfs[[ndf]] <- df
    
  }
  return(dfs)
}

dfs <- simulate_symmetries(tissue_size, number_of_cells, circle_diameter, 
                           centroids)

# Visualize ---------------------------------------------------------------

dfs[[1]] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .1) + 
  scale_color_manual(values = rainbow(3))

dfs[[2]] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .1) + 
  scale_color_manual(values = rainbow(3))

dfs[[5]] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .1) + 
  scale_color_manual(values = rainbow(3))

# Save --------------------------------------------------------------------

for (i in 1:length(dfs)) {
  print(table(dfs[[i]]$celltype))
  write.csv(dfs[[i]], paste0('simulating_data/symmetries/symmetry_', 
                             i-1, '.csv'))
}

