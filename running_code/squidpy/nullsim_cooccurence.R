library(tidyverse)
library(crawdad)

dfs <- read.csv('running_code/squidpy/results_data/dfs_nullsim_squidpy.csv', row.names = 1)
head(dfs)
unique(dfs$id)

## calculate auc
dfs_auc <- dfs %>% 
  group_by(id, reference, neighbor) %>% ## calculate auc
  summarise(auc = pracma::trapz(distance, probability)) 

## select which neighbor has the minimum auc for each reference
dfs_auc %>% 
  mutate(pair = paste(reference, neighbor)) %>% ## calculate min and max auc
  group_by(id, reference) %>% 
  filter(auc == min(auc)) %>% 
  ungroup() %>% 
  mutate(right = case_when((reference == 'A') & (neighbor == 'B') ~ T,
                           (reference == 'B') & (neighbor == 'A') ~ T,
                           (reference == 'C') & (neighbor == 'D') ~ T,
                           (reference == 'D') & (neighbor == 'C') ~ T,
                           T ~ F)) %>% 
  pull(right) %>% 
  sum()
## [1] 36

## select which neighbor has the minimum auc for each reference
dfs_auc %>% 
  mutate(pair = paste(reference, neighbor)) %>% ## calculate min and max auc
  group_by(id, reference) %>% 
  filter(auc == max(auc)) %>% 
  ungroup() %>% 
  mutate(right = case_when((reference == 'A') & (neighbor == 'A') ~ T,
                           (reference == 'B') & (neighbor == 'B') ~ T,
                           (reference == 'C') & (neighbor == 'C') ~ T,
                           (reference == 'D') & (neighbor == 'D') ~ T,
                           T ~ F)) %>% print(n=40)
  pull(right) %>% 
  sum()
## [1] 33

## verify one case
dfs %>% 
  filter(id == 3) %>% 
  filter(reference == 'B') %>% 
  ggplot() + 
  geom_line(aes(distance, probability, color = neighbor)) + 
  labs(title = 'Dataframe 3 - Reference B') +
  scale_color_manual(values = rainbow(4)) + 
  theme_minimal()


## verify one case
dfs %>% 
  filter(id == 4) %>% 
  filter(reference == 'A') %>% 
  ggplot() + 
  geom_line(aes(distance, probability, color = neighbor)) + 
  labs(title = 'Dataframe 4 - Reference A') +
  scale_color_manual(values = rainbow(4)) + 
  theme_minimal()

## spatial viz
dfs_cells <- readRDS('simulating_data/null_sim/cells_nullsim.RDS')
table(dfs_cells[[4]]$celltypes)
dfs_cells[[4]] %>% 
  vizClusters(ofInterest = c('A', 'D'), alpha = 1)