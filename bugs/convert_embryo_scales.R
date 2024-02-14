library(crawdad)
library(tidyverse)

setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/crawdad_revision_analysis")

# seqFISH ----------------------------------------------------------------

data('seq')
max(seq$x) - min(seq$x) # 5069.662
max(seq$y) - min(seq$y) # 6984.168
## it should be aprox 1000 microns, not 5000 or 6000
ggplot(seq, aes(x, -y)) + geom_point(size = .1)

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