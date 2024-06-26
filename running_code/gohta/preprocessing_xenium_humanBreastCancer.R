## This file is for loading and preprocessing 10X Genomics Xenium human breast cancer dataset

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(rhdf5)
library(dplyr)

par(mfrow=c(1,1))


# Functions ---------------------------------------------------------------

calculateDensity <- function(matrix.array) {
  sum(matrix.array != 0)/(dim(matrix.array)[1] * dim(matrix.array)[2])
}

# Load dataset ------------------------------------------------------------

## load feature x observation count matrix
file <- '~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_humanBreastCancer/cell_feature_matrix.h5'
rhdf5::h5ls(file)
data <- rhdf5::h5read(file, 'matrix')
rhdf5::h5closeAll()
names(data)

counts <- data$data
indices <- data$indices
indptr <- data$indptr
shp <- data$shape
features <- data$features
barcodes <- data$barcodes

## convert to a sparse matrix
gexp <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                   x = as.numeric(x = counts[]), dims = shp[], repr = "T")
rownames(gexp) <- features$name
colnames(gexp) <- barcodes
class(gexp)
gexp <- as(gexp, "CsparseMatrix")
dim(gexp)

## load metadata

meta <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_humanBreastCancer/cells.csv.gz')
head(meta)
pos <- cbind(meta$x_centroid, meta$y_centroid)
colnames(pos) <- c("x", "y")
rownames(pos) <- meta$cell_id
plot(meta$cell_area, meta$total_counts, pch=".")
plot(meta$nucleus_area, meta$total_counts, pch=".")
plot(meta$cell_area, meta$nucleus_area, pch=".")

## load cell types

# # graph based clustering
clusters <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_humanBreastCancer/analysis/clustering/gene_expression_graphclust/clusters.csv')
rownames(clusters) <- paste0("cell-", clusters$Barcode)

# supervised cell type
celltype <- read.csv('~/Library/CloudStorage/OneDrive-JohnsHopkins/JEFworks Gohta Aihara/Data/Xenium_humanBreastCancer/Cell_Barcode_Type_Matrices_R1_Supervised.csv')

colnames(gexp) <- rownames(meta) <- rownames(pos) <- rownames(celltype) <- paste0("cell-", celltype$Barcode)

## flip x-axis
pos[,1] <- -pos[,1]
## make all coordinates positive
pos[,1] <- pos[,1] - min(pos[,1])
pos[,2] <- pos[,2] - min(pos[,2])

## plot
plot(pos, pch=".")

df <- data.frame(pos, celltype = celltype$Cluster)
plt <- ggplot(df, aes(x = x, y = y, col = celltype)) +
  coord_fixed() +
  geom_point(size = 0.5, stroke = 0) +
  theme_bw()
plt
plotly::ggplotly(plt)

## density of feature x observation matrix
calculateDensity(gexp)

# Preprocessing -----------------------------------------------------------

## filter genes
bad_genes <- c(rownames(gexp)[grepl('BLANK', rownames(gexp))], rownames(gexp)[grepl('NegControl', rownames(gexp))])
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
cell_area <- meta$cell_area
nuc_area <- meta$nucleus_area
names(cell_area) <- names(nuc_area) <- rownames(meta)
cell_area <- cell_area[good_cells]
nuc_area <- nuc_area[good_cells]

hist(cell_area)
hist(nuc_area)

par(mfrow=c(3,1))
hist(log10(colSums(gexp)+1))
hist(log10(colSums(gexp/cell_area*mean(cell_area))+1))
hist(log10(colSums(gexp/nuc_area*mean(nuc_area))+1))

## normalize by nucleus area
gexp_lognorm <- as(log10(gexp/nuc_area*mean(nuc_area) + 1), "CsparseMatrix")

calculateDensity(gexp_lognorm)

# Format into SpatialExperiment class -------------------------------------

coldata <- data.frame(
  cluster = clusters[rownames(clusters) %in% good_cells,"Cluster"],
  celltype = celltype[rownames(celltype) %in% good_cells,"Cluster"]
) %>%
  mutate(
    cluster = as.factor(cluster),
    celltype = as.factor(celltype)
  )
rownames(coldata) <- names(good_cells)

spe <- SpatialExperiment::SpatialExperiment(
  assays = list(counts = gexp, lognorm = gexp_lognorm),
  spatialCoords = as.matrix(pos),
  colData = coldata
)

saveRDS(spe, file = "running_code/processed_data/xenium_humanBreastCancer_preprocessed.RDS")