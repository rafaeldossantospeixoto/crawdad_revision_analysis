## This file is for loading and preprocessing MERFISH mouse whole brain datasets

## cell-type annotations are obtained from STalign (https://zenodo.org/records/8384019)

# Set up ------------------------------------------------------------------

source("running_code/functions.R")

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(here)

par(mfrow=c(1,1))

# Pipeline ------------------------------------------------------------

## load cell type annotations
# ct_labels <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mouseBrain/STalign_celltypeannotations_merfishslices.csv.gz', row.names = 1)
ct_labels <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mouseBrain/STalign_celltypeannotations_merfishslices_v2.csv.gz', row.names = 1)

slices <- seq(1,3)
replicates <- seq(1,3)

for (slice in slices) {
  for (replicate in replicates) {
    dataset_name <- paste0("merfish_mouseBrain_s", slice, "_r", replicate)
    print(dataset_name)
    
# Load dataset ------------------------------------------------------------
    
    ## structure ids (mainly used to make sure cells have matching cell IDs across different objects)
    structure_labels <- read.csv(paste0('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mouseBrain/STalign_S', slice, 'R', replicate, '_with_structure_id_name.csv.gz'), row.names = 1)
    
    ## subset cell-type labels
    ct_labels_sub <- ct_labels[rownames(ct_labels) %in% rownames(structure_labels), , drop = FALSE]
    
    ## load feature x gene counts matrix
    gexp <- read.csv(paste0('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mouseBrain/datasets_mouse_brain_map_BrainReceptorShowcase_Slice', slice, '_Replicate', replicate, '_cell_by_gene_S', slice, 'R', replicate, '.csv.gz'), row.names = 1)
    
    ## convert to sparse matrix
    gexp <- as(t(gexp), "CsparseMatrix")
    
    ## load meta data
    meta <- read.csv(paste0('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/MERFISH_mouseBrain/datasets_mouse_brain_map_BrainReceptorShowcase_Slice', slice, '_Replicate', replicate, '_cell_metadata_S', slice, 'R', replicate, '.csv.gz'), row.names = 1)
    pos <- meta[, c('center_x', 'center_y')]
    
    colnames(pos) <- c("x", "y")
    colnames(gexp) <- rownames(meta) <- rownames(pos) <- rownames(ct_labels_sub) <- paste0("cell-", rownames(structure_labels))
    
    ## flip Y coordinates
    pos[,2] <- -pos[,2]
    ## make all coordinates positive
    pos[,1] <- pos[,1] - min(pos[,1])
    pos[,2] <- pos[,2] - min(pos[,2])
    
    ## plot
    plot(pos, pch=".")
    
    ## density of feature x observation matrix
    calculateDensity(gexp)
    
# preprocessing -----------------------------------------------------------

    ## filter genes
    bad_genes <- rownames(gexp)[grepl('Blank', rownames(gexp))]
    gexp <- gexp[!(rownames(gexp) %in% bad_genes),]
    
    ## filter cells
    par(mfrow=c(2,3))
    hist(log10(colSums(gexp)+1))
    hist(log10(colSums(gexp>0)+1))
    hist(log10(colSums(gexp>1)+1))
    
    good_cells <- colnames(gexp)[colSums(gexp>2) > 0]
    gexp <- gexp[,colnames(gexp) %in% good_cells]
    hist(log10(colSums(gexp)+1))
    hist(log10(colSums(gexp>0)+1))
    hist(log10(colSums(gexp>1)+1))
    
    pos <- pos[rownames(pos) %in% good_cells,]
    
    par(mfrow=c(1,1))
    plot(pos, pch=".")
    
    calculateDensity(gexp)
    
    ## normalization
    
    vol <- meta$volume
    names(vol) <- rownames(meta)
    vol <- vol[good_cells]
    
    hist(log10(vol))
    plot(vol, colSums(gexp), pch = ".")
    
    par(mfrow=c(2,1))
    hist(log10(colSums(gexp)+1))
    hist(log10(colSums(gexp/vol*mean(vol))+1))
    
    gexp_lognorm <- as(log10(gexp/vol*mean(vol) + 1), "CsparseMatrix")
    
    calculateDensity(gexp_lognorm)

# format into SpatialExperiment class -------------------------------------
    coldata <- ct_labels_sub[rownames(ct_labels_sub) %in% good_cells, c("celltype", "celltype_merged"), drop = FALSE]
    coldata$celltype <- as.factor(coldata$celltype)
    coldata$celltype_merged <- as.factor(coldata$celltype_merged)
    
    spe <- SpatialExperiment::SpatialExperiment(
      assays = list(counts = gexp, lognorm = gexp_lognorm),
      spatialCoords = as.matrix(pos),
      colData = coldata
    )
    
    saveRDS(spe, file = here("running_code", "processed_data", paste0(dataset_name, ".RDS")))
  }
}
