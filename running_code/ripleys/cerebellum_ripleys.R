library(crawdad)
library(tidyverse)
library(spatstat)

data(slide)

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=slide$x, y=slide$y, 
            window = owin(range(slide$x), range(slide$y)), 
            marks=factor(slide$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo,'Purkinje', x)
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
colnames(df) <- c('neighbor', 'scale', 'score')
df$reference <- 'Purkinje'
df$permutation <- 1
saveRDS(df, "running_code/processed_data/dat_cerebellum_ripleys.RDS")