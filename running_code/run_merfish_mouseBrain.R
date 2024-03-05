## This file is for running CRAWDAD on MERFISH mouse whole brain datasets

# Set up ------------------------------------------------------------------

source("running_code/functions.R")

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(here)
library(ggplot2)
library(gridExtra)
library(pracma)

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

n_dist <- 50

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
      cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype_merged)
    }
    
    ## visualize
    # crawdad::vizEachCluster(cells = cells,
    #                         coms = as.factor(cells$celltypes),
    #                         s = 2)
    
    # ## shuffle cells to create null background
    # shuffle.list <- crawdad:::makeShuffledCells(cells,
    #                                             scales = scales,
    #                                             perms = 3,
    #                                             ncores = ncores,
    #                                             seed = 1,
    #                                             verbose = TRUE)
    # if (ct_type == "original") {
    #   saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
    # } else if (ct_type == "cleaned") {
    #   saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells_ct_cleaned.RDS")))
    # }
    ## calculate the zscore for the cell-type pairs at different scales
    if (ct_type == "original") {
      shuffle.list <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
    } else if (ct_type == "cleaned") {
      shuffle.list <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells_ct_cleaned.RDS")))
    }
    results <- crawdad::findTrends(cells = cells,
                                   dist = n_dist,
                                   shuffle.list = shuffle.list,
                                   ncores = ncores,
                                   verbose = TRUE,
                                   returnMeans = FALSE)
    if (ct_type == "original") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends_dist_", n_dist, ".RDS")))
    } else if (ct_type == "cleaned") {
      saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned_dist_", n_dist, ".RDS")))
    }
  }
}


# Plot --------------------------------------------------------------------

slices <- seq(1,3)
replicates <- seq(1,3)

## set color values for each cell type for consistency
spe_s2r1 <- readRDS(file = here("running_code", "processed_data", paste0(paste0("merfish_mouseBrain_s", 2, "_r", 1), ".RDS")))
col_ct <- gg_color_hue(length(levels(spe_s2r1$celltype_merged)))
names(col_ct) <- levels(spe_s2r1$celltype_merged)

## Figure x (summary plots)
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
    print(paste0("Corrected Z-score threshold: ", zsig))
    ## summary visualization
    vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, reorder = TRUE, dot.sizes = c(2, 8)) +
      theme(axis.text.x = element_text(angle = 35, h = 0)) +
      ggtitle(dataset_name)
    ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0(dataset_name, "_summary_ct_cleaned.pdf")), width = 12, height = 8, dpi = 300)
  }
}

## use results with cleaned cell-types to analyze similarities/differences across slices

## Figure x (difference in scale after crossing the significant Z-score threshold)
# get scale and Z-score that first crossed the significant Z-score threshold
df <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)

    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    # calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    
    # create data.frame with the Z-scores and scales at the first scale the trend becomes significant
    # get mean Z
    mean_dat <- dat %>% 
      group_by(neighbor, scale, reference) %>% 
      summarize(Z = mean(Z))
    # get values before filtering
    max_scale <- max(dat$scale)
    u_cts <- unique(dat$reference)
    # calculate sig z scores
    sig_dat <- mean_dat %>%
      filter(abs(Z) >= zsig) %>% 
      group_by(neighbor, reference) %>% 
      filter(scale == min(scale, na.rm = TRUE))
    
    return(data.frame(slice = slice, replicate = replicate, sig_dat))
  }))
}))

# summarize for each pair
# define cell-type pairs
ct_pairs <- expand.grid(unique(colData(spe)$celltype_merged), unique(colData(spe)$celltype_merged))
colnames(ct_pairs) <- c("neighbor", "reference")

# define slice pairs
slice_pairs <- combn(slices, 2)

df_diff <- do.call(rbind, lapply(seq(nrow(ct_pairs)), function(i) {
  # subset to chosen neighbor and reference
  neighbor <- ct_pairs$neighbor[i]
  reference <- ct_pairs$reference[i]
  df_sub <- df[df$neighbor == neighbor & df$reference == reference,]
  
  # summarize scale and Z-score by slide
  df_sub <- df[df$neighbor == neighbor & df$reference == reference,] %>%
    group_by(slice) %>%
    summarize(
      scale_mean = mean(scale),
      Z_mean = mean(Z)
    )
  
  # compute difference in scale between each slide pair
  out <- do.call(rbind, lapply(seq(dim(slice_pairs)[2]), function(j) {
    return(data.frame(slice_A = slice_pairs[1,j], slice_B = slice_pairs[2,j], neighbor = neighbor, reference = reference, diff_scale = df_sub$scale_mean[slice_pairs[1,j]] - df_sub$scale_mean[slice_pairs[2,j]], diff_scale_abs = abs(df_sub$scale_mean[slice_pairs[1,j]] - df_sub$scale_mean[slice_pairs[2,j]])))
  }))
  return(out)
}))

# plot difference in scale for each slide pair
# subtraction
for (i in seq(dim(slice_pairs)[2])) {
  df_diff_sub <- df_diff[df_diff$slice_A == slice_pairs[1,i] & df_diff$slice_B == slice_pairs[2,i],]
  
  ggplot(df_diff_sub, aes(x = reference, y = neighbor, fill = diff_scale)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(position = "top") +
    labs(title = paste0("slice ", slice_pairs[1,i], " - ", "slice ", slice_pairs[2,i]),
         fill = "Difference\nin scale") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("slice", slice_pairs[1,i], "-", "slice", slice_pairs[2,i], "_diff_scale_ct_cleaned.pdf")), dpi = 300)
}

# absolute(subtraction)
for (i in seq(dim(slice_pairs)[2])) {
  df_diff_sub <- df_diff[df_diff$slice_A == slice_pairs[1,i] & df_diff$slice_B == slice_pairs[2,i],]
  
  ggplot(df_diff_sub, aes(x = reference, y = neighbor, fill = diff_scale_abs)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(position = "top") +
    labs(title = paste0("|slice ", slice_pairs[1,i], " - ", "slice ", slice_pairs[2,i], "|"),
         fill = "Absolute\ndifference\nin scale") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("slice", slice_pairs[1,i], "-", "slice", slice_pairs[2,i], "_diff_scale_abs_ct_cleaned.pdf")), dpi = 300)
}

## Figure x (difference in area under the scale vs. Z-score curves)
# compute AUC for each scale vs. Z-score curves from each dataset 
df_auc <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    
    # format data
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    # define pairs
    pairs <- distinct(dat[c("neighbor", "reference")])
    
    # compute AUC for each pair
    out <- do.call(rbind, lapply(seq(nrow(pairs)), function(i) {
      out2 <- dat %>%
        filter(reference == pairs$reference[i]) %>%
        filter(neighbor == pairs$neighbor[i]) %>%
        group_by(neighbor, scale, reference) %>%
        summarise(mean = mean(Z),
                  sd = sd(Z))
      return(data.frame(slice = slice, replicate = replicate, reference = pairs$reference[i], neighbor = pairs$neighbor[i], auc = trapz(out2$scale, out2$mean)))
    }))
  }))
}))
# save df_auc
saveRDS(df_auc, file = here("running_code", "outputs", "merfish_mouseBrain_diff_auc.RDS"))
# load df_auc
df_auc <- readRDS(file = here("running_code", "outputs", "merfish_mouseBrain_diff_auc.RDS"))

# compare AUC for each celltype for each pair of datasets
# create a data.frame of slice and replicate ids
df_ids <- distinct(df_auc[,c("slice", "replicate")])
# define dataset pairs
dataset_pairs <- combn(nrow(df_ids), 2)
# define cell type pairs
celltype_pairs <- distinct(df_auc[c("neighbor", "reference")])

# plot for each pair
for (i in seq(ncol(dataset_pairs))) {
  # get slice and replicate ids
  dataset1 <- df_ids[dataset_pairs[1,i],]
  dataset2 <- df_ids[dataset_pairs[2,i],]
  dataset1_name <- paste0("s", dataset1$slice, "r", dataset1$replicate)
  dataset2_name <- paste0("s", dataset2$slice, "r", dataset2$replicate)
  print(paste0(dataset1_name, " vs. ", dataset2_name))
  
  # subset to selected datasets
  df_auc_sub <- df_auc %>%
    filter((slice == dataset1$slice & replicate == dataset1$replicate) | (slice == dataset2$slice & replicate == dataset2$replicate))
  
  # compute difference (absolute value) in auc for each cell type pair between two datasets
  df_plt <- do.call(rbind, lapply(seq(nrow(celltype_pairs)), function(j) {
    out <- df_auc_sub %>%
      filter(reference == celltype_pairs$reference[j]) %>%
      filter(neighbor == celltype_pairs$neighbor[j])
    return(data.frame(reference = celltype_pairs$reference[j], neighbor = celltype_pairs$neighbor[j], diff_auc_abs = abs(out$auc[1] - out$auc[2])))
  }))
  
  # plot
  ggplot(df_plt, aes(x = reference, y = neighbor, fill = diff_auc_abs)) +
    coord_fixed() +
    geom_tile() +
    scale_fill_viridis_c() +
    scale_x_discrete(position = "top") +
    labs(title = paste0(dataset1_name, " vs. ", dataset2_name),
         fill = "Absolute\ndifference\nin AUC") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 35, h = 0))
  
  # save
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("diff_auc_abs_ct_cleaned_", dataset1_name, " vs. ", dataset2_name, ".pdf")), dpi = 300)
}

# combine everything into one plot to standardize the scale
df_plt <- do.call(rbind, lapply(seq(ncol(dataset_pairs)), function(i) {
  # get slice and replicate ids
  dataset1 <- df_ids[dataset_pairs[1,i],]
  dataset2 <- df_ids[dataset_pairs[2,i],]
  dataset1_name <- paste0("s", dataset1$slice, "r", dataset1$replicate)
  dataset2_name <- paste0("s", dataset2$slice, "r", dataset2$replicate)
  print(paste0(dataset1_name, " vs. ", dataset2_name))
  
  # subset to selected datasets
  df_auc_sub <- df_auc %>%
    filter((slice == dataset1$slice & replicate == dataset1$replicate) | (slice == dataset2$slice & replicate == dataset2$replicate))
  
  # compute difference (absolute value) in auc for each cell type pair between two datasets
  temp <- do.call(rbind, lapply(seq(nrow(celltype_pairs)), function(j) {
    out <- df_auc_sub %>%
      filter(reference == celltype_pairs$reference[j]) %>%
      filter(neighbor == celltype_pairs$neighbor[j])
    return(data.frame(dataset1 = dataset1_name, dataset2 = dataset2_name, reference = celltype_pairs$reference[j], neighbor = celltype_pairs$neighbor[j], diff_auc_abs = abs(out$auc[1] - out$auc[2])))
  }))
}))
ggplot(df_plt, aes(x = reference, y = neighbor, fill = diff_auc_abs)) +
  facet_grid(dataset1 ~ dataset2) +
  coord_fixed() +
  geom_tile() +
  scale_fill_viridis_c() +
  # scale_x_discrete(position = "top") +
  labs(fill = "Absolute\ndifference\nin AUC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, h = 1))
ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("diff_auc_abs_ct_cleaned_all_ctpairs.pdf")), width = 20, height = 20, dpi = 300)

## select interesting pairs
# distribution of absolute difference in AUC
hist(df_plt$diff_auc_abs)
# select top
threshold <- 6000
pairs_selected <- distinct(df_plt[df_plt$diff_auc_abs >= threshold,c("reference", "neighbor")])
pairs_selected <- pairs_selected[complete.cases(pairs_selected),] %>%
  mutate(reference = as.character(reference),
         neighbor = as.character(neighbor))
pairs_selected

# # checking NaN values
# ## s1r3 vs. s2s1
# i <- 16
# i <- 116
# slice <- 1
# replicate <- 3
# reference <- "Ependymal Cells"
# neighbor <- "Medium Spiny Neurons"
# # compute Z score for each scale
# df_curve <- do.call(rbind, lapply(slices, function(slice) {
#   out <- do.call(rbind, lapply(replicates, function(replicate) {
#     dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
#     
#     findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
#     dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
#     
#     return(data.frame(slice = slice, replicate = replicate, dat))
#   }))
# }))
# # filter to selected reference, neighbor pair
# # df_plt <- df_curve %>%
# #   filter(reference == reference) %>%
# #   filter(neighbor ==  neighbor)
# df_plt <- df_curve[(df_curve$neighbor == neighbor) & (df_curve$reference == reference),]
# # compute error bar
# df_plt <- df_plt %>%
#   group_by(slice, replicate, neighbor, scale, reference) %>%
#   summarise(mean = mean(Z),
#             sd = sd(Z)) %>%
#   mutate(slice = factor(slice),
#          replicate = factor(replicate))
# 
# # plot
# ggplot(df_plt) +
#   geom_point(aes(x = scale, y = mean, col = slice, shape = replicate)) +
#   geom_path(aes(x = scale, y = mean, col = slice, shape = replicate)) +
#   # geom_errorbar(aes(x = scale, ymin = mean-sd, ymax = mean+sd, col = slice), width = 10) +
#   # geom_hline(yintercept = zsig_mean, color = "black", linetype = "dotted") +
#   # geom_hline(yintercept = -zsig_mean, color = "black", linetype = "dotted") +
#   labs(title = paste0("Reference: ", reference, ", Neighbor: ", neighbor),
#        x = "scale (µm)",
#        y = "Z") +
#   theme_bw()

## Further analysis of cell type colocalizations that differ across datasets
# choose pair based on a threshold
diff_threshold <- 200
# df_diff_sub <- na.omit(df_diff[df_diff$diff_scale >= diff_threshold,])
df_diff_sub <- na.omit(df_diff[df_diff$diff_scale_abs >= diff_threshold,])
df_diff_sub
pairs_selected <- distinct(df_diff_sub[c("neighbor", "reference")])

## Figure x (validation with scale vs. Z curve)
# compute Z score for each scale
df_curve <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    return(data.frame(slice = slice, replicate = replicate, dat))
  }))
}))
# compute multiple test corrected Z-score threshold
df_zsig <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    
    findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends_ct_cleaned.RDS")))
    dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
    
    ## calculate the zscore for the multiple-test correction
    ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
    psig <- 0.05/ntests
    zsig <- round(qnorm(psig/2, lower.tail = F), 2)
    
    return(data.frame(slice = slice, replicate = replicate, num_ct = length(unique(dat$reference)), zsig = zsig))
  }))
}))
zsig_mean <- mean(df_zsig$zsig)

# loop through pair
for (i in seq(nrow(pairs_selected))) {
  # filter to selected reference, neighbor pair
  df_plt <- df_curve %>%
    filter(reference == pairs_selected$reference[i]) %>%
    filter(neighbor == pairs_selected$neighbor[i])
  # compute error bar
  df_plt <- df_plt %>%
    group_by(slice, replicate, neighbor, scale, reference) %>%
    summarise(mean = mean(Z),
              sd = sd(Z)) %>%
    mutate(slice = factor(slice),
           replicate = factor(replicate))
  # plot
  ggplot(df_plt) +
    geom_point(aes(x = scale, y = mean, col = slice, shape = replicate)) +
    geom_path(aes(x = scale, y = mean, col = slice, shape = replicate)) +
    # geom_errorbar(aes(x = scale, ymin = mean-sd, ymax = mean+sd, col = slice), width = 10) +
    geom_hline(yintercept = zsig_mean, color = "black", linetype = "dotted") +
    geom_hline(yintercept = -zsig_mean, color = "black", linetype = "dotted") +
    labs(title = paste0("Reference: ", pairs_selected$reference[i], ", Neighbor: ", pairs_selected$neighbor[i]),
         x = "scale (µm)",
         y = "Z") +
    theme_bw()
  # save
  reference_label <- gsub("/", "-", pairs_selected$reference[i])
  neighbor_label <- gsub("/", "-", pairs_selected$neighbor[i])
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("merfish_mouseBrain_scale_vs_Zscore_reference_", reference_label, "_neighbor_", neighbor_label, ".png")))
}

## Figure x (validation with single-cell visualization)
for (i in seq(nrow(pairs_selected))) {
  print(paste0("Reference: ", pairs_selected$reference[i], ", Neighbor: ", pairs_selected$neighbor[i]))
  pairs <- c(pairs_selected$reference[i], pairs_selected$neighbor[i])
  
  # get spatial coords for all datasets
  df_plt <- do.call(rbind, lapply(slices, function(slice) {
    out <- do.call(rbind, lapply(replicates, function(replicate) {
      dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
      
      # load dataset
      spe <- readRDS(file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))

      return(data.frame(slice = slice, replicate = replicate, spatialCoords(spe), celltype = spe$celltype_merged))
    }))
  }))
  
  # subset to a specific pair
  df_plt_sub <- df_plt[df_plt$celltype %in% pairs,]
  
  # plot
  ggplot() +
    facet_grid(replicate ~ slice, labeller = label_both) +
    coord_fixed() +
    geom_point(data = df_plt, aes(x = x, y = y), color = "lightgray", size = 0.5, stroke = 0) +
    geom_point(data = df_plt_sub, aes(x = x, y = y, col = celltype), size = 0.5, stroke = 0) +
    scale_color_manual(values = col_ct[pairs]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0("Reference: ", pairs_selected$reference[i], ", Neighbor: ", pairs_selected$neighbor[i]),
         x = "x (um)",
         y = "y (um)",
         col = "Cell types") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
    )
  
  # save
  reference_label <- gsub("/", "-", pairs_selected$reference[i])
  neighbor_label <- gsub("/", "-", pairs_selected$neighbor[i])
  ggsave(filename = here("running_code", "plots", "merfish_mouseBrain", paste0("merfish_mouseBrain_singlecell_reference_", reference_label, "_neighbor_", neighbor_label, ".png")), width = 10, height = 8, dpi = 300)
}

## Figure x (cell type proportion)
# compute cell type proportion
df_ct_comparison <- do.call(rbind, lapply(slices, function(slice) {
  out <- do.call(rbind, lapply(replicates, function(replicate) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(dataset_name)
    
    ## load dataset
    spe <- readRDS(file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))
    
    ## compute the cell type frequency
    count <- table(spe$celltype_merged)
    proportion <- prop.table(table(spe$celltype_merged))
    
    ## combine count and prop
    comb <- rownames_to_column(data.frame(cbind(count, proportion)), var = "celltype")
    
    ## format output
    temp <- data.frame(slice = slice, replicate = replicate, sample = paste0("s", slice, "_r", replicate), comb)
    
    return(temp)
  }))
}))

# plot
# count
ggplot(df_ct_comparison, aes(x = sample, y = count, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()

# proportion
ggplot(df_ct_comparison, aes(x = sample, y = proportion, fill = celltype)) +
  geom_bar(position = "stack", stat = "identity") +
  theme_bw()

# which cell types are missing from slice 1
spe_s1 <- readRDS(file = here("running_code", "processed_data", paste0(paste0("merfish_mouseBrain_s", 1, "_r", 1), ".RDS")))
spe_s2 <- readRDS(file = here("running_code", "processed_data", paste0(paste0("merfish_mouseBrain_s", 2, "_r", 2), ".RDS")))
levels(spe_s1$celltype_merged)
levels(spe_s2$celltype_merged)
