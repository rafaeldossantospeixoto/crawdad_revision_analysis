library(crawdad)
library(tidyverse)
library(spatstat)

data(slide)


# Homogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=slide$x, y=slide$y, 
            window = owin(range(slide$x), range(slide$y)), 
            marks=factor(slide$celltypes))
## run on all
reference_ct <- 'Purkinje'
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo, reference_ct, x)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## border-corrected minus theoretical
  df_test$score <- test$border - test$theo
  return(df_test)
}))

rkc %>%  
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "border-corrected minus theoretical score")

df <- rkc
df$permutation <- 1
saveRDS(df, "running_code/processed_data/dat_cerebellum_ripleys.RDS")



# Inhomogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=slide$x, y=slide$y, 
            window = owin(range(slide$x), range(slide$y)), 
            marks=factor(slide$celltypes))
## run on all
reference_ct <- 'Purkinje'
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kinhom(X = ppo, i = reference_ct, j = x)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## border-corrected minus theoretical
  df_test$score <- test$border - test$theo
  return(df_test)
}))

rkc %>%
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "border-corrected minus theoretical score") +
  facet_wrap('neighbor')

df <- rkc
df$permutation <- 1
saveRDS(df, "running_code/processed_data/dat_slide_ripleys_inhom.RDS")