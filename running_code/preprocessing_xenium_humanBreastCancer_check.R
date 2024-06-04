# check if the analyzed dataset is indeed 10X Genomics Xenium breast cancer dataset in situ sample 1, replicate 1

# Set up ------------------------------------------------------------------

library(SpatialExperiment)
library(Matrix)
library(ggplot2)
library(rhdf5)
library(dplyr)

par(mfrow=c(1,1))

# Load dataset (analyzed) ------------------------------------------------------------

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

# Load dataset (in situ sample 1, replicate 1 just downloaded) ------------------------------------------------------------

## load feature x observation count matrix
file <- '~/Downloads/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs/cell_feature_matrix.h5'
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
gexp2 <- sparseMatrix(i = indices[] + 1, p = indptr[], 
                     x = as.numeric(x = counts[]), dims = shp[], repr = "T")

# check -------------------------------------------------------------------

identical(gexp, gexp2)
