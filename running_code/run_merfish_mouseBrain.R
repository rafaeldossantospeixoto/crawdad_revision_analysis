## This file is for running CRAWDAD on MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(here)
library(ggplot2)

# Run method --------------------------------------------------------------

## based on the tutorial https://jef.works/CRAWDAD/

slices <- seq(1,3)
replicates <- seq(1,3)

## define the scales to analyze the data
scales <- seq(100, 1000, by=100)

## define the number of cores for parallelization
# ncores <- parallel::detectCores() - 2
ncores <- 10

## whether to use original or cleaned celltypes ("original" or "cleaned")
ct_type <- "cleaned"
ct_type

for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Running CRAWDAD on ", dataset_name))
    
    ## load dataset
    spe <- readRDS(file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))
    
    ## convert dataframe to sf data.frame
    if (ct_type == "original") {
      # use original cell types
      cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)
    } else if (ct_type == "cleaned") {
      # use cleaned cell types
      cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype_cleaned)
    }
    
    ## visualize
    # crawdad::vizEachCluster(cells = cells,
    #                         coms = as.factor(cells$celltypes),
    #                         s = 2)
    
    ## shuffle cells to create null background
    shuffle.list <- crawdad:::makeShuffledCells(cells,
                                                scales = scales,
                                                perms = 3,
                                                ncores = ncores,
                                                seed = 1,
                                                verbose = TRUE)
    if (ct_type == "original") {
      saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
    } else if (ct_type == "cleaned") {
      saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells_ct_cleaned.RDS")))
    }
    ## calculate the zscore for the cell-type pairs at different scales
    ## error: Error in FUN(X[[i]], ...) : object 'neigh.cells' not found
    results <- crawdad::findTrends(cells = cells,
                                   dist = 100,
                                   shuffle.list = shuffle.list,
                                   ncores = ncores,
                                   verbose = TRUE,
                                   returnMeans = FALSE)
    if (ct_type == "original") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))
    } else if (ct_type == "cleaned") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    }
  }
}


# Plot --------------------------------------------------------------------

slices <- seq(1,3)
replicates <- seq(1,3)

## original cell types
for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Plotting CRAWDAD results for ", dataset_name))
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))
    
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    ## calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    ## summary visualization
    vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, dot.sizes = c(2, 8)) +
      theme(axis.text.x = element_text(angle = 35, h = 0)) +
      ggtitle(dataset_name)
    ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0(dataset_name, "_summary.pdf")), width = 12, height = 8, dpi = 300)
  }
}

## cleaned cell types
for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(paste0("Plotting CRAWDAD results for ", dataset_name))
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    ## calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    ## summary visualization
    vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, reorder = TRUE, dot.sizes = c(2, 8)) +
      theme(axis.text.x = element_text(angle = 35, h = 0)) +
      ggtitle(dataset_name)
    ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0(dataset_name, "_summary_ct_cleaned.pdf")), width = 12, height = 8, dpi = 300)
  }
}

## use results with cleaned cell-types to analyze similarities/differences across slices
## get scale and Z-score that first crossed the significant Z-score threshold
df <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)

    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    ## calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    
    ## create data.frame with the Z-scores and scales at the first scale the trend becomes significant
    ## get mean Z
    mean_dat <- dat %>% 
      group_by(neighbor, scale, reference) %>% 
      summarize(Z = mean(Z))
    ## get values before filtering
    max_scale <- max(dat$scale)
    u_cts <- unique(dat$reference)
    ## calculate sig z scores
    sig_dat <- mean_dat %>%
      filter(abs(Z) >= zsig) %>% 
      group_by(neighbor, reference) %>% 
      filter(scale == min(scale, na.rm = TRUE))
    
    return(data.frame(slice = slice, replicate = replicate, sig_dat))
  }))
}))

## summarize for each pair
## define cell-type pairs
ct_pairs <- expand.grid(unique(colData(spe)$celltype_cleaned), unique(colData(spe)$celltype_cleaned))
colnames(ct_pairs) <- c("neighbor", "reference")

## define slice pairs
slice_pairs <- combn(slices, 2)

df_diff <- do.call(rbind, lapply(seq(nrow(ct_pairs)), function(i) {
  ## subset to chosen neighbor and reference
  neighbor <- ct_pairs$neighbor[i]
  reference <- ct_pairs$reference[i]
  df_sub <- df[df$neighbor == neighbor & df$reference == reference,]
  
  ## summarize scale and Z-score by slide
  df_sub <- df[df$neighbor == neighbor & df$reference == reference,] %>%
    group_by(slice) %>%
    summarize(
      scale_mean = mean(scale),
      Z_mean = mean(Z)
    )
  
  ## compute difference in scale between each slide pair
  out <- do.call(rbind, lapply(seq(dim(pairs)[2]), function(j) {
    return(data.frame(slice_A = pairs[1,j], slice_B = pairs[2,j], neighbor = neighbor, reference = reference, diff_scale = df_sub$scale_mean[pairs[1,j]] - df_sub$scale_mean[pairs[2,j]], diff_scale_abs = abs(df_sub$scale_mean[pairs[1,j]] - df_sub$scale_mean[pairs[2,j]])))
  }))
  return(out)
}))

## plot difference in scale for each slide pair
## subtraction
for (i in seq(dim(pairs)[2])) {
  df_diff_sub <- df_diff[df_diff$slice_A == pairs[1,i] & df_diff$slice_B == pairs[2,i],]
  
  ggplot(df_diff_sub, aes(x = reference, y = neighbor, fill = diff_scale)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(position = "top") +
    labs(title = paste0("slice ", pairs[1,i], " - ", "slice ", pairs[2,i]),
         fill = "Difference\nin scale") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("slice", pairs[1,i], "-", "slice", pairs[2,i], "_diff_scale_ct_cleaned.pdf")), dpi = 300)
}

## absolute(subtraction)
for (i in seq(dim(pairs)[2])) {
  df_diff_sub <- df_diff[df_diff$slice_A == pairs[1,i] & df_diff$slice_B == pairs[2,i],]
  
  ggplot(df_diff_sub, aes(x = reference, y = neighbor, fill = diff_scale_abs)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(position = "top") +
    labs(title = paste0("|slice ", pairs[1,i], " - ", "slice ", pairs[2,i], "|"),
         fill = "Absolute\ndifference\nin scale") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("slice", pairs[1,i], "-", "slice", pairs[2,i], "_diff_scale_abs_ct_cleaned.pdf")), dpi = 300)
}

## select specific pair
pairs_selected <- rbind(
  data.frame(reference = "GABA-ergic Interneurons", neighbor = "Excitatory Neurons")
)
i <- 1

## validation with scale vs. Z curve
df_curve <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    return(data.frame(slice = slice, replicate = replicate, dat))
  }))
}))

df_plt <- df_curve %>%
  filter(reference == pairs_selected$reference[i]) %>%
  filter(neighbor == pairs_selected$neighbor[i])
  

## validation with single-cell visualization




neighbor <- "Astrocytes"
reference <- "Ependymal Cells"
df_sub <- df[df$neighbor == neighbor & df$reference == reference,]

pairs <- combn(slices, 2)

df_sub <- df[df$neighbor == neighbor & df$reference == reference,] %>%
  group_by(slice) %>%
  summarize(
    scale_mean = mean(scale),
    Z_mean = mean(Z)
  )

out <- do.call(rbind, lapply(seq(dim(pairs)[2]), function(i) {
  return(data.frame(slice_A = pairs[1,i], slice_B = pairs[2,i], neighbor = neighbor, reference = reference, diff_scale = df_sub$scale_mean[pairs[1,i]] - df_sub$scale_mean[pairs[2,i]]))
}))

df_sub <- df[df$neighbor == neighbor & df$reference == reference,] %>%
  select(c("slice", "scale"))

## statistical test? (anova, kruskal-wills)
shapiro.test(residuals(lm(scale ~ slice, data = df_sub)))

bartlett.test(scale ~ slice, data = df_sub)

anova_result <- aov(scale ~ slice, data = df_sub)
summary(anova_result)

# Perform post hoc test (Tukey's HSD)
tukey_result <- TukeyHSD(anova_result)
tukey_result

kruskal_result <- stats::kruskal.test(scale ~ slice, data = df_sub)
kruskal_result







dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
print(paste0("Plotting CRAWDAD results for ", dataset_name))

findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)

## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)

## create data.frame with the Z-scores and scales at the first scale the trend becomes significant
## get mean Z
mean_dat <- dat %>% 
  group_by(neighbor, scale, reference) %>% 
  summarize(Z = mean(Z))
## get values before filtering
max_scale <- max(dat$scale)
u_cts <- unique(dat$reference)
## calculate sig z scores
sig_dat <- mean_dat %>%
  filter(abs(Z) >= zsig) %>% 
  group_by(neighbor, reference) %>% 
  filter(scale == min(scale, na.rm = TRUE))

df <- data.frame(slice = slice, replicate = replicate, sig_dat)