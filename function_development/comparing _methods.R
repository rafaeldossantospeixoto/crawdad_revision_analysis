library(crawdad)
library(tidyverse)


# Sim ---------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/dat_sim_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')

zsig <- 3

# dat_sp$scale <- round(dat_sp$scale, digits = 0)
dat_sp %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'B') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)


## CRAWDAD vs Squidpy ------------------------------------------------------

## B and C
dat_sp %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))



## D and C
dat_sp %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))



## A and C
dat_sp %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'C') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))


# Slide -------------------------------------------------------------------

dat_sp <- read.csv('running_code/squidpy/dat_slide_squidpy.csv', row.names = 1)
dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')

zsig <- 3

# dat_sp$scale <- round(dat_sp$scale, digits = 0)
dat_sp %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'B') %>% 
  vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)



## CRAWDAD vs Squidpy ------------------------------------------------------

## Purkinje and Bergmann
dat_sp %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Bergmann') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Bergmann') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

## Purkinje and Oligodendrocytes
dat_sp %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Oligodendrocytes') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'Purkinje') %>% 
  filter(neighbor == 'Oligodendrocytes') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))


## Astrocytes and Microglia
dat_sp %>% 
  filter(reference == 'Astrocytes') %>% 
  filter(neighbor == 'Microglia') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))

dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == 'Astrocytes') %>% 
  filter(neighbor == 'Microglia') %>% 
  ggplot() + 
  geom_line(aes(x=scale, y=Z))