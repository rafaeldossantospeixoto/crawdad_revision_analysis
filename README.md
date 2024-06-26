# crawdad_revision_analysis

This repository contains code scripts to reproduce analyses and figures for the review process of our *bioRxiv* preprint [Peixoto R. dos S. et al. (2023), "Characterizing cell-type spatial relationships across length scales in spatially resolved omics data", *bioRxiv*](https://doi.org/10.1101/2023.10.05.560733), which evaluates the utility of an R package CRAWDAD.

## CRAWDAD

`CRAWDAD` is a statistical framework that uses labeled spatial omics data to identify the colocalization or separation of cell types at different scales. `CRAWDAD` identifies regions where multiple cells colocalize, the scale of such colocalization, and also subsets the cell types based on their proximity to others. Therefore, `CRAWDAD` is a powerful tool for tissue characterization and comparison.

`CRAWDAD` is available on [GitHub](https://github.com/JEFworks-Lab/CRAWDAD)

## Contents

Preprocessing and analyses done by Aihara G. is organized into the [running_code/gohta](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/tree/main/running_code/gohta) directory. See below for the descriptions of each files.

-   [functions.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/functions.R): This file contains functions used in other analyses files. The `calculateDensity` computes the density of an input matrix. The `gg_color_hue` divides a standard ggplot color palette into an input integer value.
-   [preprocessing_merfish_mouseBrain.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/preprocessing_merfish_mouseBrain.R): This file is for loading and preprocessing the [MERFISH mouse brain datasets](https://info.vizgen.com/mouse-brain-map) into `SpatialExperiment` objects (`SpatialExperiment` objects are not used as an input for `CRAWDAD`, but it is used as an useful object to store both cell x,y coordinates and cell type annotations).
-   [preprocessing_xenium_humanBreastCancer.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/preprocessing_xenium_humanBreastCancer.R): This file is for loading and preprocessing the [10X Genomics Xenium human breast cancer (in situ sample 1, replicate 1) dataset](https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast) into a `SpatialExperiment` object.
-   [preprocessing_xenium_humanBreastCancer_check.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/preprocessing_xenium_humanBreastCancer_check.R): This file ensures that the analyzed 10X Genomics human breast cancer dataset is in situ sample 1, replicate 1.
-   [run_merfish_mouseBrain.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/run_merfish_mouseBrain.R): This file applies `CRAWDAD` to the preprocessed MERFISH mouse brain datasets and creates plots for analyzing the differences in the cell type spatial relationships across replicates or slices.
-   [run_xeniumhumanBreastCancer.R](https://github.com/rafaeldossantospeixoto/crawdad_revision_analysis/blob/main/running_code/gohta/run_xeniumhumanBreastCancer.R): This file applies `CRAWDAD` to the preprocessed 10X Genomics Xenium human breast cancer dataset and creates plots for analyzing the cell type spatial relationships within the tumor tissue.
