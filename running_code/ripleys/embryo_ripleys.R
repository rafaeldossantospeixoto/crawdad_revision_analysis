library(crawdad)
library(tidyverse)
library(spatstat)

data(seq)
table(seq$celltypes)
cells <- crawdad::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
vizClusters(cells)
vizClusters(cells, ofInterest = c('Endothelium', 'NMP', 
                                  'Sclerotome', 'Allantois', 
                                  'Anterior somitic tissues'), 
            alpha = 1)

# Endothelium -------------------------------------------------------------

reference_ct <- 'Endothelium'
radii <- seq(0, 1000, by=5)

## Homogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=seq$x, y=seq$y, 
            window = owin(range(seq$x), range(seq$y)), 
            marks=factor(seq$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo, reference_ct, x, r = radii)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  return(df_test)
}))

rkc %>%  
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "iso-corrected minus theoretical score")

df <- rkc
df$permutation <- 1
saveRDS(df, paste0("running_code/processed_data/dat_embryo_ripleys_",
                   sub(" ", "_", reference_ct),
                   ".RDS"))


## Inhomogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=seq$x, y=seq$y, 
            window = owin(range(seq$x), range(seq$y)), 
            marks=factor(seq$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross.inhom(X = ppo, i = reference_ct, j = x)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  return(df_test)
}))

rkc %>%  
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "iso-corrected minus theoretical score") +
  facet_wrap('neighbor')

rkc %>%  
  filter(neighbor != 'NMP') %>% 
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "iso-corrected minus theoretical score") +
  facet_wrap('neighbor')

df <- rkc
df$permutation <- 1
saveRDS(df, paste0("running_code/processed_data/dat_embryo_ripleys_inhom_",
                   sub(" ", "_", reference_ct),
                   ".RDS"))




# Intermediate mesoderm -----------------------------------------------------

reference_ct <- 'Intermediate mesoderm'
radii <- seq(0, 1000, by=5)

## Homogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=seq$x, y=seq$y, 
            window = owin(range(seq$x), range(seq$y)), 
            marks=factor(seq$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross(ppo, reference_ct, x, r = radii)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  return(df_test)
}))

rkc %>%  
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "iso-corrected minus theoretical score")

df <- rkc
df$permutation <- 1
saveRDS(df, paste0("running_code/processed_data/dat_embryo_ripleys_",
                   sub(" ", "_", reference_ct),
                   ".RDS"))


## Inhomogeneous -------------------------------------------------------------

## run ripley's L
## need to convert to these point process objects
ppo  <- ppp(x=seq$x, y=seq$y, 
            window = owin(range(seq$x), range(seq$y)), 
            marks=factor(seq$celltypes))
## run on all
rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
  test <- Kcross.inhom(X = ppo, i = reference_ct, j = x)
  df_test <- data.frame(reference = reference_ct,
                        neighbor = x,
                        radius = test$r)
  ## iso-corrected minus theoretical
  df_test$score <- test$iso - test$theo
  return(df_test)
}))

rkc %>%  
  filter(neighbor != 'NMP') %>% 
  ggplot() +
  geom_line(aes(x=radius, y=score, color=neighbor)) +
  ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3, linetype = "solid") +
  ggplot2::theme_classic() +
  ggplot2::scale_color_manual(values = rainbow(length(unique(rkc$neighbor)))) +
  labs(y = "iso-corrected minus theoretical score") +
  facet_wrap('neighbor')

df <- rkc
df$permutation <- 1
saveRDS(df, paste0("running_code/processed_data/dat_embryo_ripleys_inhom_",
                   sub(" ", "_", reference_ct),
                   ".RDS"))


