library(crawdad)
library(tidyverse)

file1 <- readRDS('bugs/slideseq_cerebellum/cerebellum/F_GRCm38.81.P60Cerebellum_ALT.filtered.raw.dge.RDS')
dim(file1)
file1[1:7,1:7]

file1 <- readRDS('bugs/slideseq_cerebellum/cerebellum/F_GRCm38.81.P60Cerebellum_ALT.cell_cluster_outcomes.RDS')
dim(file1)
head(file1)

file1 <- read.csv('bugs/slideseq_cerebellum/Puck_180819_11/BeadMapping_10-17_0913/BeadLocationsForR.csv')
dim(file1)
head(file1)
max(file1$xcoord) - min(file1$xcoord) # 3886.176
max(file1$ycoord) - min(file1$ycoord) # 4366.586
file1 %>% 
  ggplot() +
  geom_point(aes(x=xcoord, y=ycoord), size = .1)

file1 <- read.csv('bugs/slideseq_cerebellum/Puck_180819_11/BeadMapping_10-17_0913/BeadLocationsForRUncropped.csv')
dim(file1)
head(file1)
max(file1$xcoord) - min(file1$xcoord) # 4203.508
max(file1$ycoord) - min(file1$ycoord) # 4366.586
file1 %>% 
  ggplot() +
  geom_point(aes(x=xcoord, y=ycoord), size = .1)

file1 <- read.csv('bugs/slideseq_cerebellum/Puck_180819_12/BeadLocationsForR.csv')
dim(file1)
head(file1)
max(file1$xcoord) - min(file1$xcoord) # 5798.119
max(file1$ycoord) - min(file1$ycoord) # 5616.116
file1 %>% 
  ggplot() +
  geom_point(aes(x=xcoord, y=ycoord), size = .1)

setwd("~/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/crawdad_revision_analysis")
setwd("/Users/rafaeldossantospeixoto/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/CRAWDAD")

# slideSeq ----------------------------------------------------------------

data('slide')
max(slide$x) - min(slide$x) # 4628.55
max(slide$y) - min(slide$y) # 3800.465
## it should be aprox 2904 and 2425 microns, not 4628.55 and 3800.465
ggplot(slide, aes(x, -y)) + geom_point(size = .1)
