library(crawdad)

#' Calculate AUC for each cell type pair from a dat file
#' 
#' @description
#' Calculate the AUC values for each cell type pair from the dat variable of a
#' sample.
#' 
#' @param dat `findTrends()` list of data.frames; the information about the 
#' scale, Z-score, reference and the neighbor cell. The input data.frame should 
#' be the results list from `findTrends()` that has been melted into a 
#' data.frame using `meltResultsList()`.
#' @param sharedPairs boolean; if true, will remove the AUC values from the 
#' cell-type pairs that are not available in all samples (default: TRUE).
#' 
#' @return data.frame; the AUC values for the cell-type pairs in the samples
#' 
#' @export
calculateAUC <- function(datList, sharedPairs = TRUE) {
  
  auc_list <- list()
  
  ## calculate
  for (dat in datList) {
    sample_id <- dat$id[1]
    
    ## define pairs
    pairs <- dplyr::distinct(dat[c("neighbor", "reference")])
    
    ## compute AUC for each pair
    auc_sample <- do.call(rbind, lapply(seq(nrow(pairs)), function(i) {
      ## process dat
      summarized_dat <- dat %>%
        dplyr::filter(reference == pairs$reference[i]) %>%
        dplyr::filter(neighbor == pairs$neighbor[i]) %>%
        dplyr::group_by(neighbor, scale, reference) %>%
        dplyr::summarise(mean = mean(Z),
                         sd = sd(Z))
      
      ## calculate auc
      return(data.frame(id = sample_id,
                        reference = pairs$reference[i], 
                        neighbor = pairs$neighbor[i], 
                        auc = pracma::trapz(summarized_dat$scale, summarized_dat$mean)))
    }))
    
    auc_list[[sample_id]] <- auc_sample
  }
  
  ## bind auc from all samples
  auc_df <- dplyr::bind_rows(auc_list) 
  
  ## process auc
  ## create ct pair colunm
  auc_df <- auc_df %>% 
    dplyr::mutate(pair = paste0(paste0(reference, ' - ', neighbor))) 
  
  ## if true, will remove cell type pairs that are not present in all samples
  if (sharedPairs) {
    ## calculate the ct pairs that are shared across all samples
    shared_pairs <- auc_df %>% 
      dplyr::group_by(pair) %>%
      dplyr::summarize(n_unique = dplyr::n_distinct(id)) %>% 
      dplyr::filter(n_unique == dplyr::n_distinct(auc_df$id)) %>% 
      dplyr::pull(pair)
    
    ## filter data based on shared pairs
    auc_df <- auc_df %>% 
      dplyr::filter(pair %in% shared_pairs)
  }
  
  return(auc_df)
}



#' Visualize the samples using their AUC and PCA
#' 
#' @description
#' Uses the AUC values calculated from each cell-type pair in each sample to 
#' represent the samples in the PCA reduced dimension space.
#' 
#' @param aucSamples data.frame; the AUC values for each cell-type pair in each
#' sample as calculated by `calculateAUC()`
#' 
#' @return gglot2 plot; the PCA visualization of the samples
#' 
#' @export
#' 
vizPCASamples <- function(aucSamples) {
  sample_ids <- unique(aucSamples$id)
  
  ## create matrix
  auc_mtx <- aucSamples %>% 
    dplyr::select(c(pair, id, auc)) %>% 
    tidyr::pivot_wider(names_from = id, values_from = auc) %>% 
    dplyr::select(!pair) %>% 
    tidyr::drop_na() %>% ## drop nas
    as.matrix() %>% 
    t()
  
  ## normalize
  auc_mtx <- scale(auc_mtx)
  apply(auc_mtx, 2, mean)
  
  ## calculate pca
  pca <- prcomp(auc_mtx)
  pcs <- pca$x[, 1:2] %>% 
    as.data.frame() %>% 
    dplyr::mutate(id = rownames(pca$x)) %>% 
    dplyr::mutate(patient = 
                    sapply(rownames(pca$x), 
                           FUN = function(x) stringr::str_split(x, '_')[[1]][2]))
  
  ## plot
  pcs %>% 
    ggplot2::ggplot() + 
    ggplot2::geom_point(ggplot2::aes(x = PC1, y = PC2, color = id)) + 
    ggplot2::scale_color_manual(values = rainbow(length(sample_ids))) +
    ggplot2::theme_bw() +
    ggplot2::coord_equal()
}



#' Visualize the variance of each cell-type relationship in the samples
#' 
#' @description
#' Visualize the variance of the AUC values for each cell-type pair in all 
#' samples.
#' 
#' @param aucSamples data.frame; the AUC values for each cell-type pair in each
#' sample as calculated by `calculateAUC()`
#' 
#' @return gglot2 plot; the dot plot visualization of the variance in the 
#' samples
#' 
#' @export
#' 
vizVarianceSamples <- function(aucSamples) {
  
  ## all samples
  aucSamples %>% 
    dplyr::group_by(reference, neighbor) %>%
    dplyr::filter(reference != 'indistinct',
           neighbor != 'indistinct') %>% 
    dplyr::summarize(variance = (var(auc))) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = reference, y = neighbor, 
                                     size = variance), color = '#006437') +
    ggplot2::scale_radius(range = c(1, 10)) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(position = 'top') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 0.5, 
                                                       hjust=0)) +
    ggplot2::coord_equal()
}







# Testing -----------------------------------------------------------------

## load data
data(pkhl)
pkhl <- crawdad::toSF(pos = pkhl[,c("x", "y")], celltypes = pkhl$celltypes)
dat_pkhl <- readRDS('running_code/processed_data/spleen/dat_pkhl_50.RDS')
dat_pkhl$id <- 'pkhl'

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
xxcd <- crawdad::toSF(pos = xxcd[,c("x", "y")], celltypes = xxcd$celltypes)
dat_xxcd <- readRDS('running_code/processed_data/spleen/dat_xxcd_50.RDS')
dat_xxcd$id <- 'xxcd'

fsld <- read.csv2(file = '../CRAWDAD/data/spleen/FSLD.meta.csv.gz', row.names = 1)
fsld <- crawdad::toSF(pos = fsld, celltypes = fsld$celltypes)
dat_fsld <- readRDS('running_code/processed_data/spleen/dat_fsld_50.RDS')
dat_fsld$id <- 'fsld'

## calculate auc
auc_samples <- calculateAUC(list(dat_pkhl, dat_xxcd, dat_fsld))
vizPCASamples(auc_samples)
vizVarianceSamples(auc_samples)
