library(crawdad)
library(tidyverse)
ncores <- 7


# Process sim -------------------------------------------------------------

data(sim)

sim %>% ggplot() +
  geom_point(aes(x, y, color = celltypes))
sim <- sim %>% 
  filter(!((celltypes == 'C') & (x < 1000) & (y < 1000)),
         !((celltypes == 'C') & (x > 1000) & (y > 1000)))
sim %>% ggplot() +
  geom_point(aes(x, y, color = celltypes))

saveRDS(sim, 'running_code/processed_data/sim_asymmetrical_imbalance_B.RDS')



# Run CRAWDAD -------------------------------------------------------------

## convert to sf
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)

scales <- seq(100, 1000, by=50)

## generate background
shuffle.list <- crawdad:::makeShuffledCells(sim,
                                            scales = scales,
                                            perms = 10,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
## Time was 5.25 mins

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
