library(crawdad)
library(tidyverse)


# Load data ---------------------------------------------------------------

df <- read.csv('simulating_data/symmetries/symmetry_0.csv', row.names = 1)
dim(df)
head(df)
df %>% 
  ggplot() +
  geom_point(aes(x, y, color = celltype)) +
  scale_color_manual(values = rainbow(3))
## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = df[,c("x", "y")], celltypes = df$celltype)


# Compare neighborhoods ---------------------------------------------------

## run function from compare_neighborhoods.R

## Podoplanin
reference <- 'C'
distances <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
it <- Sys.time()
plt_distances <- compareNeighborhoods(cells, reference, distances)
Sys.time() - it
## Time difference of 7.18728 mins
plt_distances

