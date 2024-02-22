## This file is for running CRAWDAD on 10X Genomics Xenium human breast cancer dataset

# Set up ------------------------------------------------------------------

source("running_code/functions.R")

library(SpatialExperiment)
library(Matrix)
library(crawdad)
library(tidyverse)
library(ggplot2)
library(ggrastr)
library(here)

dataset_name <- "xenium_humanBreastCancer"

# Load dataset ------------------------------------------------------------

spe <- readRDS(file = "running_code/processed_data/xenium_humanBreastCancer_preprocessed.RDS")

df <- data.frame(spatialCoords(spe), colData(spe))
ggplot(df, aes(x = x, y = y, col = celltype)) +
  geom_point(size = 0.5, stroke = 0) +
  theme_classic()

ct_labels <- colData(spe)$celltype

# Run method --------------------------------------------------------------

## based on the tutorial https://jef.works/CRAWDAD/

## convert dataframe to spatial points (SP)
cells <- crawdad::toSF(pos = data.frame(spatialCoords(spe)), celltypes = colData(spe)$celltype)

## visualize
# crawdad::vizEachCluster(cells = cells,
#                         coms = as.factor(cells$celltypes),
#                         s = 2)

## define the scales to analyze the data
scales <- seq(100, 1000, by=100)
## shuffle cells to create null background
ncores <- parallel::detectCores() - 2
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perms = 3,
                                            ncores = ncores,
                                            seed = 1,
                                            verbose = TRUE)
saveRDS(shuffle.list, here("running_code", "outputs", paste0(dataset_name, "_makeShuffledCells.RDS")))
## calculate the zscore for the cell-type pairs at different scales
## error: Error in FUN(X[[i]], ...) : object 'neigh.cells' not found
results <- crawdad::findTrends(cells = cells,
                               dist = 100,
                               shuffle.list = shuffle.list,
                               ncores = ncores,
                               verbose = TRUE,
                               returnMeans = FALSE)
saveRDS(results, here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))


# Plot --------------------------------------------------------------------

col_ct <- gg_color_hue(length(levels(ct_labels)))
names(col_ct) <- levels(ct_labels)

## Figure x (singecell visualization)
df <- data.frame(spatialCoords(spe), celltype = colData(spe)$celltype)
ggplot(df, aes(x = x, y = y, col = celltype)) +
  rasterise(geom_point(size = 0.5, stroke = 0.01), dpi = 300) +
  scale_color_manual(values = col_ct) +
  guides(col = guide_legend(override.aes = list(size = 3))) +
  labs(x = "x (um)",
       y = "y (um)",
       col = "Cell types") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
  )
ggsave(filename = here("running_code", "plots", dataset_name, paste0(dataset_name, "_singlecell.pdf")), dpi = 300)

## Figure x (CRAWDAD summary plots)
findTrends_results <- readRDS(here("running_code", "outputs", paste0(dataset_name, "_findTrends.RDS")))

dat <- crawdad::meltResultsList(findTrends_results, withPerms = TRUE)
## calculate the zscore for the multiple-test correction
ntests <- length(unique(dat$reference)) * length(unique(dat$reference))
psig <- 0.05/ntests
zsig <- round(qnorm(psig/2, lower.tail = F), 2)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig, reorder = TRUE, dot.sizes = c(2, 10)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 35, h = 0))
ggsave(filename = here("running_code", "plots", dataset_name, paste0(dataset_name, "_summary.pdf")), dpi = 300)

## remove unlabeled
dat_filtered <- dat[!(dat$neighbor == "Unlabeled" | dat$reference == "Unlabeled"),]
vizColocDotplot(dat_filtered, zsig.thresh = zsig, zscore.limit = 2*zsig, reorder = TRUE, dot.sizes = c(2, 10)) +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 35, h = 0))
ggsave(filename = here("running_code", "plots", dataset_name, paste0(dataset_name, "_summary_unlabeled_removed.pdf")), dpi = 300)

## Figure x (colocalizations at single-cell)

# manually identified using a CRAWDAD summary plot (reordered, unlabeled cells removed)
niches <- list(
  factor(c("DCIS_1", "DCIS_2", "Myoepi_ACTA2+", "Myoepi_KRT15+"), levels = levels(ct_labels)),
  factor(c("DCIS_1", "DCIS_2", "Myoepi_ACTA2+", "Myoepi_KRT15+", "Invasive_Tumor", "Prolif_Invasive_Tumor"), levels = levels(ct_labels)),
  factor(c("Invasive_Tumor", "Prolif_Invasive_Tumor"), levels = levels(ct_labels)),
  factor(c("B_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "Endothelial", "IRF7+_DCs", "LAMP3+_DCs", "Macrophages_1", "Macrophages_2", "Mast_Cells", "Perivascular-Like", "Stromal", "Stromal_&_T_Cell_Hybrid", "T_Cell_&_Tumor_Hybrid"), levels = levels(ct_labels)),
  factor(c("B_Cells", "CD4+_T_Cells", "CD8+_T_Cells", "Endothelial", "IRF7+_DCs", "LAMP3+_DCs", "Macrophages_1", "Macrophages_2", "Mast_Cells", "Perivascular-Like", "Stromal"), levels = levels(ct_labels))
)

for (celltypes in niches) {
  ## subset by clusters
  spe_sub <- spe[,spe$celltype %in% celltypes]
  
  ## plot
  df <- data.frame(spatialCoords(spe))
  df_sub <- data.frame(spatialCoords(spe_sub), celltype = colData(spe_sub)$celltype)
  ggplot(df, aes(x = x, y = y)) +
    coord_fixed() +
    rasterise(geom_point(color = "lightgray", size = 0.5, stroke = 0), dpi = 300) +
    rasterise(geom_point(data = df_sub, aes(x = x, y = y, col = celltype), size = 0.5, stroke = 0), dpi = 300) +
    scale_color_manual(values = col_ct[celltypes]) +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    labs(title = paste0("Cell types: ", paste(celltypes, collapse = ", ")),
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
  ggsave(filename = here("running_code", "plots", dataset_name, paste0(dataset_name, "_singlecell_niche_celltypes_", paste(celltypes, collapse = "_"), ".pdf")), dpi = 300)
}
