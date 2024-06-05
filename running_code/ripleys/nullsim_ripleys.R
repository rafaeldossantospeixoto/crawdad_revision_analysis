library(crawdad)
library(tidyverse)
library(spatstat)

nullsim <- read.csv(file = 'simulating_data/null_sim/df_nullsim_s1.csv', 
                row.names = 1)
head(nullsim)


# Seed 1 -------------------------------------------------------------

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


## AUC ---------------------------------------------------------------------

## calculate auc
dat_auc <- rkc %>% 
  group_by(reference, neighbor) %>% ## calculate auc
  summarise(auc = pracma::trapz(radius, score)) 





# Seeds 1 to 10 -----------------------------------------------------------

dfs <- read.csv('simulating_data/null_sim/dfs_nullsim.csv', row.names = 1)

## calculate K for all datasets
rkcs_list <- list()
for (i in 1:length(unique(dfs$id))) {
  
  nullsim <- dfs %>% 
    filter(id == i)
  
  ## need to convert to these point process objects
  ppo  <- ppp(x=nullsim$x, y=nullsim$y, 
              window = owin(range(nullsim$x), range(nullsim$y)), 
              marks = factor(nullsim$celltypes))
  
  radii <- seq(0, 500, by=5)
  
  ## calculate for all references
  rkc_references <- list()
  for (j in 1:4) {
    rkc <- do.call(rbind, lapply(levels(ppo$marks), function(x) {
      test <- Kcross(ppo, LETTERS[j], x, r = radii)
      df_test <- data.frame(reference = LETTERS[j],
                            neighbor = x,
                            radius = test$r)
      ## isotropic-corrected minus theoretical
      df_test$score <- test$iso - test$theo
      return(df_test)
    }))
    rkc_references[[j]] <- rkc
  }
  rkc <- bind_rows(rkc_references)
  
  rkc$id <- i
  
  rkcs_list[[i]] <- rkc
}

rkcs <- bind_rows(rkcs_list)



## AUC ---------------------------------------------------------------------

## calculate auc
rkc_auc <- rkcs %>% 
  group_by(id, reference, neighbor) %>% ## calculate auc
  summarise(auc = pracma::trapz(radius, score)) 

## select which neighbor has the minimum auc for each reference
rkc_auc %>% 
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
## [1] 31

## select which neighbor has the minimum auc for each reference
rkc_auc %>% 
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
## [1] 33


