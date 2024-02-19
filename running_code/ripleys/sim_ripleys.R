library(crawdad)
library(tidyverse)
library(spatstat)

data(sim)

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=sim$x, y=sim$y, 
            window = owin(range(sim$x), range(sim$y)), 
            marks=factor(sim$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo,'B', x)
  ## isotropic-corrected minus theoretical
  test$iso - test$theo
}))
rownames(rkc) <- levels(ppo$marks)
# colnames(rkc) <- ppo$r

## plot
rkc.melt <- reshape2::melt(rkc)
colnames(rkc.melt) <- c('celltype', 'distance', 'score')
head(rkc.melt)

rkc.melt %>%  
  ggplot() +
  geom_line(aes(x=distance, y=score, color=celltype)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.6, linetype = "dotted") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc.melt$celltype)))) +
  labs(y = "isotropic-corrected minus theoretical score")

df <- rkc.melt
colnames(df) <- c('neighbor', 'distance', 'score')
df$reference <- 'B'
df$permutation <- 1
saveRDS(df, "running_code/processed_data/dat_sim_ripleys.RDS")