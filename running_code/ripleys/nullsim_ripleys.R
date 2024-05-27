library(crawdad)
library(tidyverse)
library(spatstat)

nullsim <- read.csv(file = 'simulating_data/null_sim/df_nullsim_s1.csv', 
                row.names = 1)
head(nullsim)


# Homogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=nullsim$x, y=nullsim$y, 
            window = owin(range(nullsim$x), range(nullsim$y)), 
            marks=factor(nullsim$celltypes))
## run on all
reference_ct <- 'A'
## fixed scales
radii <- seq(0, 500, by=5)
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo, reference_ct, x, r = radii)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## isotropic-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  return(df_test)
}))

rkc %>%  
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "isotropic minus theoretical score")

df <- rkc
df$permutation <- 1
saveRDS(df, "running_code/processed_data/dat_nullsim_ripleys.RDS")


