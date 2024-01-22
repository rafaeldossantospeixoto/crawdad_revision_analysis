library(crawdad)
library(tidyverse)

# Create dataset --------------------------------------------------------

## similar to a real dataset (spleen)
## spleen is 3.5 x 3.5 mm, with 150k cells and structures of .5 mm
tissue_size <- 2 * 1e3 ## 5 mm or 5000 microns
circle_diameter <- 1e3
## number of cells is estimated based on sleen density
## spleen has 45k per 1 mm in size
number_of_cells <- 45 * tissue_size

proportions <- c(.5, .4, .3, .2, .1, .05, .01)

simulate_densities <- function(tissue_size, number_of_cells, circle_diameter, 
                               proportions) {
  dfs <- list()
  for (i in 1:length(proportions)) {
    ## define background
    df <- crawdad:::simulate_background(size = number_of_cells, scale = tissue_size)
    
    ## simulate circle
    x_center <- tissue_size/2
    y_center <- tissue_size/2
    radius <- circle_diameter/2
    cells_circle <- rownames(df[((df$x-x_center)^2 + (df$y - y_center)^2 < radius^2),])
    df[cells_circle, 'type'] <- 'B'
    cells_type_c <- sample(cells_circle, 
                           size = round(length(cells_circle) * proportions[i]), 
                           replace = FALSE)
    df[cells_type_c, 'type'] <- 'C'
    colnames(df) <- c('x', 'y', 'celltype')
    dfs[[i]] <- df
  }
  return(dfs)
}

dfs <- simulate_densities(tissue_size, number_of_cells, circle_diameter, 
                          proportions)

# Visualize ---------------------------------------------------------------

dfs[[1]] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .1) + 
  scale_color_manual(values = rainbow(3))

dfs[[7]] %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype), size = .1) + 
  scale_color_manual(values = rainbow(3))



# Save --------------------------------------------------------------------

for (i in 1:length(dfs)) {
  print(table(dfs[[i]]$celltype))
  write.csv(dfs[[i]], paste0('simulating_data/densities/proportion_', 
                             str_replace(100*proportions[[i]], '\\.', ''),
                             '.csv'))
}

