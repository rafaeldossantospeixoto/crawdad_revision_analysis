
# Dataset seed 1 ----------------------------------------------------------

library(crawdad)
library(tidyverse)

df1 <- readRDS("simulating_data/null_sim/seed_1/df1.RDS")
df2 <- readRDS("simulating_data/null_sim/seed_1/df2.RDS")

color_names <- rainbow(4)


## df1 ---------------------------------------------------------------------

p <- df1 %>% 
  ggplot() + 
  geom_point(aes(x, y, color = values)) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') +
                        # breaks = c(-3, 0, 3)) + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/nullsim/viz_df1_values.pdf')
p
dev.off()

p <- df1 %>% 
  ggplot() + 
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = color_names[1:2]) + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/nullsim/viz_df1_cts.pdf')
p
dev.off()

## df2 ---------------------------------------------------------------------

p <- df2 %>% 
  ggplot() + 
  geom_point(aes(x, y, color = values)) + 
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/nullsim/viz_df2_values.pdf')
p
dev.off()

p <- df2 %>% 
  ggplot() + 
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = color_names[3:4]) + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/nullsim/viz_df2_cts.pdf')
p
dev.off()


## both --------------------------------------------------------------------

p <- rbind(df1, df2) %>% 
  ggplot() + 
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = color_names) + 
  theme_bw() + 
  coord_equal()
p
pdf('running_code/paper_figures/nullsim/viz_dfs_cts.pdf')
p
dev.off()




## Run CRAWDAD ---------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

cells <- readRDS('simulating_data/null_sim/cells_nullsim_s1.RDS')

vizClusters(cells)

## reduced to 500
scales <- seq(100, 500, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)

## find trends, dist 50
results_50 <- crawdad::findTrends(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_nullsim_50.RDS')

## multiple-test correction
zsig <- correctZBonferroni(dat_50)

## viz 50
vizColocDotplot(dat_50, zScoreLimit = zsig*2)
vizColocDotplot(dat_50, zSigThresh = zsig, zScoreLimit = zsig*2)



# Paper figures -----------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_nullsim_50.RDS')

## error because no dot with significance
zsig <- correctZBonferroni(dat_50)
# vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
#                 zscoreLimit = zsig*2, 
#                 dotSizes = c(2, 14)) +
#   theme(legend.position='right',
#         axis.text.x = element_text(angle = 45, h = 0))



## Dotplot -----------------------------------------------------------------

p <- vizColocDotplot(dat_50, reorder = F, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1))
p
pdf('running_code/paper_figures/nullsim_dotplot_crawdad.pdf',
    height = 4, width = 5)
p
dev.off()



## Trends ------------------------------------------------------------------

## will plot using same design for squidpy and ripleys
ref_cts <- c(rep('A', 4), LETTERS[2:4])
neigh_cts <- c(LETTERS[1:4], rep('A', 3))
for (i in 1:length(ref_cts)) {
  p <- dat_50 %>% 
    filter(reference == ref_cts[i]) %>% 
    filter(neighbor == neigh_cts[i]) %>% 
    vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
  print(p)
  # pdf(sprintf('running_code/paper_figures/nullsim/trend_ref%s_neigh%s_crawdad.pdf',
  #             ref_cts[i], neigh_cts[i]),
  #     height = 4, width = 5)
  # print(p)
  # dev.off()
}



# Datasets 1 to 10 --------------------------------------------------------


## Run CRAWDAD -------------------------------------------------------------


library(crawdad)
library(tidyverse)

ncores <- 7
scales <- seq(100, 500, by=50)

dfs <- readRDS('simulating_data/null_sim/cells_nullsim.RDS')
dats <- list()

for (i in 1:length(dfs)) {
  
  print(paste('Analyzing dataset', i))
  
  cells <- dfs[[i]]
  
  ## generate background
  shuffle.list <- crawdad:::makeShuffledCells(cells,
                                              scales = scales,
                                              perms = 10,
                                              ncores = ncores,
                                              seed = i,
                                              verbose = TRUE)
  
  ## find trends, dist 50
  results_50 <- crawdad::findTrends(cells,
                                    dist = 50,
                                    shuffle.list = shuffle.list,
                                    ncores = ncores,
                                    verbose = TRUE,
                                    returnMeans = FALSE)
  dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
  dat_50$id <- i
  
  dats[[i]] <- dat_50
}

saveRDS(dats, 'running_code/processed_data/nullsim/dats.RDS')





## Evaluate performance ----------------------------------------------------

dat <- bind_rows(dats)

## Bonferroni
zsig <- correctZBonferroni(dat)
df <- define_relationship_type(dat, zSigThresh = zsig)

df %>% 
  group_by(reference, neighbor) %>% 
  summarize(enrichment = sum(enrichment), 
            depletion = sum(depletion))

## 1.96
zsig <- 1.96
df <- define_relationship_type(dat, zSigThresh = zsig)

df %>% 
  group_by(reference, neighbor) %>% 
  summarize(enrichment = sum(enrichment), 
            depletion = sum(depletion))




## AUC ---------------------------------------------------------------------

dat <- bind_rows(dats)

## calculate auc
dat_auc <- dat %>% 
  group_by(id, scale, reference, neighbor) %>% ## mean across perms
  summarise(Z = mean(Z)) %>% 
  group_by(id, reference, neighbor) %>% ## calculate auc
  summarise(auc = pracma::trapz(scale, Z)) 

## select which neighbor has the minimum auc for each reference
dat_auc %>% 
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
## [1] 38

## select which neighbor has the minimum auc for each reference
dat_auc %>% 
  mutate(pair = paste(reference, neighbor)) %>% ## calculate min and max auc
  group_by(id, reference) %>% 
  filter(auc == max(auc)) %>% 
  ungroup() %>% 
  mutate(right = case_when((reference == 'A') & (neighbor == 'A') ~ T,
                           (reference == 'B') & (neighbor == 'B') ~ T,
                           (reference == 'C') & (neighbor == 'C') ~ T,
                           (reference == 'D') & (neighbor == 'D') ~ T,
                           T ~ F)) %>%
  pull(right) %>% 
  sum()
## [1] 38

## verify one case
dats[[3]] %>% 
  group_by(scale, reference, neighbor) %>% ## mean across perms
  summarise(Z = mean(Z)) %>% 
  filter(reference == 'B') %>% 
  ggplot() + 
  geom_line(aes(scale, Z, color = neighbor)) + 
  labs(title = 'Dataframe 3 - Reference B') +
  scale_color_manual(values = rainbow(4)) + 
  theme_minimal()
dfs[[3]] %>% 
  vizClusters(ofInterest = c('A', 'B'), alpha = 1)
dfs[[3]] %>% 
  vizClusters(ofInterest = c('D', 'B'), alpha = 1)






# D = 30 --------------------------------------------------------


## Run CRAWDAD -------------------------------------------------------------


library(crawdad)
library(tidyverse)

ncores <- 7
scales <- seq(100, 500, by=50)

dfs <- readRDS('simulating_data/null_sim/cells_nullsim.RDS')
dats <- list()

for (i in 1:length(dfs)) {
  
  print(paste('Analyzing dataset', i))
  
  cells <- dfs[[i]]
  
  ## generate background
  shuffle.list <- crawdad:::makeShuffledCells(cells,
                                              scales = scales,
                                              perms = 5,
                                              ncores = ncores,
                                              seed = i,
                                              verbose = TRUE)
  
  ## find trends, dist 30
  results <- crawdad::findTrends(cells,
                                 dist = 30,
                                 shuffle.list = shuffle.list,
                                 ncores = ncores,
                                 verbose = TRUE,
                                 returnMeans = FALSE)
  dat <- crawdad::meltResultsList(results, withPerms = T)
  dat$id <- i
  
  dats[[i]] <- dat
}

saveRDS(dats, 'running_code/processed_data/nullsim/dats_30.RDS')
dats <- readRDS('running_code/processed_data/nullsim/dats_30.RDS')



## Evaluate performance ----------------------------------------------------

dat <- bind_rows(dats)

## Bonferroni
zsig <- correctZBonferroni(dat)
df <- define_relationship_type(dat, zSigThresh = zsig)

df %>% 
  group_by(reference, neighbor) %>% 
  summarize(enrichment = sum(enrichment), 
            depletion = sum(depletion))

## 1.96
zsig <- 1.96
df <- define_relationship_type(dat, zSigThresh = zsig)

df %>% 
  group_by(reference, neighbor) %>% 
  summarize(enrichment = sum(enrichment), 
            depletion = sum(depletion))



## AUC ---------------------------------------------------------------------

dat <- bind_rows(dats)

## calculate auc
dat_auc <- dat %>% 
  group_by(id, scale, reference, neighbor) %>% ## mean across perms
  summarise(Z = mean(Z)) %>% 
  group_by(id, reference, neighbor) %>% ## calculate auc
  summarise(auc = pracma::trapz(scale, Z)) 

## select which neighbor has the minimum auc for each reference
dat_auc %>% 
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
## [1] 39

## select which neighbor has the minimum auc for each reference
dat_auc %>% 
  mutate(pair = paste(reference, neighbor)) %>% ## calculate min and max auc
  group_by(id, reference) %>% 
  filter(auc == max(auc)) %>% 
  ungroup() %>% 
  mutate(right = case_when((reference == 'A') & (neighbor == 'A') ~ T,
                           (reference == 'B') & (neighbor == 'B') ~ T,
                           (reference == 'C') & (neighbor == 'C') ~ T,
                           (reference == 'D') & (neighbor == 'D') ~ T,
                           T ~ F)) %>%
  pull(right) %>% 
  sum()
## [1] 37