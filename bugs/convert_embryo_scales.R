library(crawdad)
library(tidyverse)

setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/crawdad_revision_analysis")

# seqFISH ----------------------------------------------------------------

data('seq')
max(seq$x) - min(seq$x) # 5069.662
max(seq$y) - min(seq$y) # 6984.168
## it should be aprox 1000+ and 1500 microns, not 5000 or 6000

data <- seq
max_x <- data[which.max(data$x), ]
min_x <- data[which.min(data$x), ]
max_y <- data[which.max(data$y), ]
min_y <- data[which.min(data$y), ]

indices_to_mutate <- c(which.max(data$x), which.min(data$x),
                       which.max(data$y), which.min(data$y))
seq %>% 
  mutate(extremes = ifelse(row_number() %in% indices_to_mutate, T, F)) %>% 
  ggplot(aes(x, -y, color = extremes)) + 
  geom_point(size = .5) + 
  scale_color_manual(values = c('white', 'darkblue'))

## investigate data from the publication
meta_seq <- readRDS('bugs/seqfish_embryo/metadata.Rds')
head(meta_seq)
table(meta_seq$embryo)
meta_seq <- meta_seq %>% 
  filter(embryo == "embryo1")
names(meta_seq)
max(meta_seq$x_global) - min(meta_seq$x_global) # 157.1595
max(meta_seq$y_global) - min(meta_seq$y_global) # 216.5092
max(meta_seq$x_global_affine) - min(meta_seq$x_global_affine) # 5.069662
max(meta_seq$y_global_affine) - min(meta_seq$y_global_affine) # 6.984168


# Update scale ------------------------------------------------------------

data('seq')
max(seq$x) - min(seq$x) # 5069.662
max(seq$y) - min(seq$y) # 6984.168
## it should be aprox 1000+ and 1500 microns, not 5000 or 6000

x_res <- 1067.278 / (max(seq$x) - min(seq$x)) # 1067.27784044
y_res <- 1578.216 / (max(seq$y) - min(seq$y)) # 1578.21592795

mean(c(x_res, y_res))

seq$x <- seq$x * x_res
seq$y <- seq$y * y_res
paste(min(seq$x), max(seq$x), min(seq$y), max(seq$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

save(seq, file = "bugs/data/seq.rda")
