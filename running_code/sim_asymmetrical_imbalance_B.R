library(crawdad)
library(tidyverse)
ncores <- 7
set.seed(42)

# Create sim -------------------------------------------------------------

data(sim)

sim$celltypes <- 'A'

centers <- list(c(500, 500), c(1500, 500), c(1500, 1500), c(500, 1500))

c1idx <- which(sqrt((sim$x - centers[[1]][1])**2 + (sim$y - centers[[1]][2])**2) <  250)
c2idx <- which(sqrt((sim$x - centers[[2]][1])**2 + (sim$y - centers[[2]][2])**2) <  250)
c3idx <- which(sqrt((sim$x - centers[[3]][1])**2 + (sim$y - centers[[3]][2])**2) <  250)
c4idx <- which(sqrt((sim$x - centers[[4]][1])**2 + (sim$y - centers[[4]][2])**2) <  250)

c1idxb <- sample(c1idx, size = length(c1idx)/2)
c2idxb <- sample(c2idx, size = length(c2idx)/2)
c3idxb <- sample(c3idx, size = length(c3idx)/2)
c4idxb <- sample(c4idx, size = length(c4idx)/2)

sim[c1idxb, 'celltypes'] <- 'B'
sim[c2idxb, 'celltypes'] <- 'B'
sim[c3idxb, 'celltypes'] <- 'B'
sim[c4idxb, 'celltypes'] <- 'B'

c1idxac <- setdiff(c1idx, c1idxb)
c1idxc <- sample(c1idxac, length(c1idxac)/10)
c1idxa <- setdiff(c1idxac, c1idxc)

sim[c1idxc, 'celltypes'] <- 'C'
sim[c1idxa, 'celltypes'] <- 'A'

sim %>% ggplot() +
  geom_point(aes(x, y, color = celltypes)) +
  scale_color_manual(values = rainbow(3))

saveRDS(sim, 'running_code/processed_data/sim_asymmetrical_imbalance_B.RDS')



# Run CRAWDAD -------------------------------------------------------------

## convert to sf
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)

scales <- seq(100, 1000, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(sim,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 5.23 mins

## find trends, dist 50
results_50 <- crawdad::findTrends(sim,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE)
## Time was 0.73 mins
dat_50 <- crawdad::meltResultsList(results_50, withPerms = T)
saveRDS(dat_50, 'running_code/processed_data/dat_50_sim_asymmetrical_imbalance_B.RDS')

zsig <- correctZBonferroni(dat_50)
p <- vizColocDotplot(dat_50, reorder = TRUE, zSigThresh = zsig, 
                     zscoreLimit = zsig*2,
                     dotSizes = c(5, 25)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
p
