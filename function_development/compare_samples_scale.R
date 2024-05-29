# Functions ---------------------------------------------------------------

#' Define Relationship Type
#' 
#' @description
#' Define the type of relationship based on the Z score (enrichment or 
#' depletion).
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame 
#' using `meltResultsList()`.
#' @param zSigThresh numeric; the Z score significance threshold (default: 1.96).
#' 
define_relationship_type <- function(dat, zSigThresh){
  dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference, id) %>% 
    dplyr::summarize(Z = mean(Z)) %>% 
    dplyr::filter(abs(Z) >= zSigThresh) %>% 
    dplyr::group_by(neighbor, reference, id) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE)) %>% 
    dplyr::mutate(relationship = dplyr::case_when(Z > 0 ~ 'enrichment',
                                                  Z < 0 ~ 'depletion',
                                                  T ~ 'other')) %>% 
    dplyr::mutate(enrichment = (relationship == 'enrichment'),
                  depletion = (relationship == 'depletion'))
  return(dat)
}



vizMutualRelationships <- function(dats, zSigThresh = 1.96, dotSizes = c(6,31)) {
  
  ## create columns and filter shared pairs
  df <- dplyr::bind_rows(dats) %>% 
    create_pair_colum() %>% 
    filter_shared_pairs()
  
  ## filter and classify significant relationships
  df <- define_relationship_type(df, zSigThresh) %>% 
    group_by(reference, neighbor, relationship) %>% 
    summarize(Z = mean(Z),
              scale = mean(scale))
  
  ## separate enrichment and depletion
  enr_df <- df %>% 
    dplyr::filter(relationship == 'enrichment')
  dep_df <- df %>% 
    dplyr::filter(relationship == 'depletion')
  
  ## plotting parameters
  ## all cell types
  all_cts <- sort(unique(c(as.character(df$reference), 
                           as.character(df$neighbor))))
  ## scale sizes
  lsizes <- sort(unique(df$scale))
  legend_sizes <- c(lsizes[1],
                    round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                    lsizes[length(lsizes)])
  
  ## plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = enr_df,
                        ggplot2::aes(x=reference, y=neighbor, size=scale), 
                        color='red', alpha = 0.5) + 
    ggplot2::geom_point(data = dep_df,
                        ggplot2::aes(x=reference, y=neighbor, size=scale), 
                        color='blue', alpha = 0.5) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = dotSizes) + 
    ggplot2::theme_bw() +
    ggplot2::scale_x_discrete(limits = all_cts, position = 'top') +
    ggplot2::scale_y_discrete(limits = all_cts) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, 
                                                       vjust = 0.5, 
                                                       hjust=0)) +
    labs(size = 'mean scale') +
    ggplot2::coord_equal()
  
  return(p)
}



# Testing -----------------------------------------------------------------

library(crawdad)

## S1R*
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r1')
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r2_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r2')
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r3')
dats <- list(dat1, dat2, dat3)

zsig <- correctZBonferroni(dplyr::bind_rows(dats))
vizMutualRelationships(dats, zSigThresh =  zsig, dotSizes = c(5, 15))

## S*R1
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE)
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s2_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE)
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s3_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE)

zsig <- correctZBonferroni(dplyr::bind_rows(dats))
vizMutualRelationships(dats, zSigThresh =  zsig, dotSizes = c(5, 15))





# Paper figures -----------------------------------------------------------


## S1R* --------------------------------------------------------------------

## S1R*
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r1')
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r2_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r2')
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r3_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r3')
dats <- list(dat1, dat2, dat3)

zsig <- correctZBonferroni(dplyr::bind_rows(dats))
p <- vizMutualRelationships(dats, zSigThresh =  zsig, dotSizes = c(1, 7))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_mean_scale_s1.pdf'),
    height = 5, width = 7)
p
dev.off()


## S*R1 --------------------------------------------------------------------

## S*R1
dat1 <- readRDS('running_code/outputs/merfish_mouseBrain_s1_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's1r1')
dat2 <- readRDS('running_code/outputs/merfish_mouseBrain_s2_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's2r1')
dat3 <- readRDS('running_code/outputs/merfish_mouseBrain_s3_r1_findTrends_ct_cleaned_dist_50.RDS') %>% 
  crawdad::meltResultsList(withPerms = TRUE) %>% 
  dplyr::mutate(id = 's3r1')
dats <- list(dat1, dat2, dat3)

zsig <- correctZBonferroni(dplyr::bind_rows(dats))
p <- vizMutualRelationships(dats, zSigThresh =  zsig, dotSizes = c(1, 7))
p
pdf(paste0('function_development/comparing_samples/paper_figures/',
           'merfish_brains_mean_scale_r1.pdf'),
    height = 5, width = 7)
p
dev.off()




# Deprecated --------------------------------------------------------------

## Dotplots comparing ------------------------------------------------------



#' @param colors character vector; colors for the gradient heatmap (low, mid, 
#' high). NAs will be colored in light gray.

vizDiffZscores <- function(dat1, dat2, 
                          scale.thresh, 
                          zscore.limit = NULL,  
                          reorder = FALSE,
                          only.significant = FALSE,
                          colors = c("blue", "white", "red"),
                          dot.size = 30){
  
  ## check if cell types are the same
  if (!identical(sort(unique(dat1$reference)), sort(unique(dat2$reference)))) {
    stop("the number of reference cell types should be the same in both datasets")
  }
  if (!identical(sort(unique(dat1$neighbor)), sort(unique(dat2$neighbor)))) {
    stop("the number of neighbor cell types should be the same in both datasets")
  }
  
  ## use the selectSigDat function to select the z scores and scales that are 
  ## above a threshold. The difference now is that the threshold is not the 
  ## zscore of singificance, but it is one defined by the user
  sig_dat1 <- selectZscoresDat(dat1, scale.thresh, zscore.limit)
  sig_dat2 <- selectZscoresDat(dat2, scale.thresh, zscore.limit)
  
  sig_join <- full_join(sig_dat1, sig_dat2, by = c("neighbor", "reference"),
                        suffix = c("1", "2"))
  
  ## did not need to do this, but there are nans in the dat variable, so it is 
  ## best to check before
  ## we are considering the z score of the non significant as 0
  sig_dif <- sig_join %>% 
    dplyr::mutate(sig_change = 
                    case_when(( (is.na(Z1) & (! is.na(Z2))) ) ~ 'became significant',
                              ( (! is.na(Z1) & (is.na(Z2))) ) ~ 'became not significant',
                              TRUE ~ 'no change in sigificance') ) %>% 
    dplyr::mutate(z_dif = 
                    case_when(sig_change == 'became significant' ~ Z2, 
                              sig_change == 'not significant' ~ Z1,
                              TRUE ~ Z2 - Z1))
  
  ## reorder based on dat1
  if (reorder) {
    ## get values from dat1
    u_cts <- unique(dat1$reference)
    max_scale <- max(dat1$scale)
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat1, by = c('reference', 'neighbor'))
    df_all <- df_all %>% 
      dplyr::mutate(Z = dplyr::coalesce(Z, 0),
                    scale = dplyr::coalesce(scale, max_scale))
    ## create matrix
    sig_mat <- reshape::cast(df_all, neighbor~reference, value='Z')
    rownames(sig_mat) <- sig_mat[,1]
    sig_mat <- sig_mat[,-1]
    sig_mat[is.na(sig_mat)] <- 0
    ## cluster
    hc <- hclust(dist(sig_mat))
    sig_dif$neighbor <- factor(sig_dif$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dif$reference <- factor(sig_dif$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## plot
  p <- sig_dif %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=reference, y=neighbor, 
                                     color=z_dif), size = dot.size) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    ggplot2::scale_color_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "lightgray"
    ) + 
    ggplot2::scale_x_discrete(position = "top")  + 
    ggplot2::theme_bw()
  
  ## show all values (even those that do not reach significance)
  if (!only.significant) {
    all_cts <- sort(unique(dat$reference))
    alpha_cts <- sort(unique(dat$reference))
    if (reorder) {
      ct_order <- rownames(sig_mat)[hc$order]
      all_cts <- c(ct_order, setdiff(alpha_cts, ct_order))
    }
    p <- p +
      scale_x_discrete(limits = all_cts, position = 'top') +
      scale_y_discrete(limits = all_cts) 
  }
  
  return(p)
}



#' Plot Co-localization Dotplot
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the co-localization at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' @param zsig.thresh numeric; the Z score significance threshold (default: 1.96).
#' @param zscore.limit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zscore.limit will be represented as 
#' zscore.limit, scores below -zscore.limit will be represented as -zscore.limit (default: 3).
#' 
selectZscoresDat <- function(dat, scale.thresh, zscore.limit){
  ## create data.frame with the Z-scores and scales at the first scale
  ## the trend becomes significant
  ## get mean Z
  mean_dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference) %>% 
    dplyr::summarize(Z = mean(Z))
  ## calculate z scores that match scale
  sig_dat <- mean_dat %>%
    dplyr::filter(scale == scale.thresh)
  
  ## only cap if parameter is passed
  ## user could choose not to cap when comparing samples
  if (!is.null(zscore.limit)) {
    ## limit the z-score for the gradient in the figure to look better
    sig_dat$Z[sig_dat$Z > zscore.limit] <- zscore.limit
    sig_dat$Z[sig_dat$Z < -zscore.limit] <- -zscore.limit
  }
  
  return(sig_dat)
}






vizDiffScales <- function(dat1, dat2, 
                          zscore.thresh = 1.96, 
                          pvalue = NULL,  
                          zscore.limit = NULL,
                          reorder = FALSE,
                          only.significant = FALSE,
                          colors = c("blue", "white", "red"),
                          dot.size = 30){
  
  ## check if cell types are the same
  if (!identical(sort(unique(dat1$reference)), sort(unique(dat2$reference)))) {
    stop("the number of reference cell types should be the same in both datasets")
  }
  if (!identical(sort(unique(dat1$neighbor)), sort(unique(dat2$neighbor)))) {
    stop("the number of neighbor cell types should be the same in both datasets")
  }
  
  ## change pvalues to zscores
  if (!is.null(pvalue)) {
    zscore.thresh = round(qnorm(pvalue/2, lower.tail = F), 2)
  }
  
  ## use the selectSigDat function to select the z scores and scales that are 
  ## above a threshold. The difference now is that the threshold is not the 
  ## zscore of singificance, but it is one defined by the user
  sig_dat1 <- selectScalesDat(dat1, zscore.thresh, zscore.limit)
  sig_dat2 <- selectScalesDat(dat2, zscore.thresh, zscore.limit)
  
  sig_join <- full_join(sig_dat1, sig_dat2, by = c("neighbor", "reference"),
                        suffix = c("1", "2"))
  
  sig_dif <- sig_join %>% 
    dplyr::mutate(sig_change = 
                    if ( (is.na(Z1) & (! is.na(Z2))) ) {'became significant'} 
                  else if ( (! is.na(Z1) & (is.na(Z2))) ) {'became not significant'}
                  else {'no change in sigificance'}) %>% 
    dplyr::mutate(scale_dif = 
                    if ( sig_change == 'no change in sigificance' ) {scale2 - scale1} 
                  ## we are considering the scale of the non significant as 0
                  else if ( sig_change == 'not significant' ) { scale1 } 
                  else if ( sig_change == 'became significant' ) { scale2 } ) %>% 
    dplyr::mutate(z_dif = 
                    if ( sig_change == 'no change in sigificance' ) {Z2 - Z1} 
                  ## we are considering the zscore of the non significant as 0
                  else if ( sig_change == 'not significant' ) { Z1 } 
                  else if ( sig_change == 'became significant' ) { Z2 } )
  
  ## reorder based on dat1
  if (reorder) {
    ## get values from dat1
    u_cts <- unique(dat1$reference)
    max_scale <- max(dat1$scale)
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat1, by = c('reference', 'neighbor'))
    df_all <- df_all %>% 
      dplyr::mutate(Z = dplyr::coalesce(Z, 0),
                    scale = dplyr::coalesce(scale, max_scale))
    ## create matrix
    sig_mat <- reshape::cast(df_all, neighbor~reference, value='Z')
    rownames(sig_mat) <- sig_mat[,1]
    sig_mat <- sig_mat[,-1]
    sig_mat[is.na(sig_mat)] <- 0
    ## cluster
    hc <- hclust(dist(sig_mat))
    sig_dif$neighbor <- factor(sig_dif$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dif$reference <- factor(sig_dif$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## plot
  p <- sig_dif %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=reference, y=neighbor, 
                                     color=scale_dif), size = dot.size) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    ggplot2::scale_color_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "lightgray"
    ) + 
    ggplot2::scale_x_discrete(position = "top")  + 
    ggplot2::theme_bw()
  
  ## show all values (even those that do not reach significance)
  if (!only.significant) {
    all_cts <- sort(unique(dat$reference))
    alpha_cts <- sort(unique(dat$reference))
    if (reorder) {
      ct_order <- rownames(sig_mat)[hc$order]
      all_cts <- c(ct_order, setdiff(alpha_cts, ct_order))
    }
    p <- p +
      scale_x_discrete(limits = all_cts, position = 'top') +
      scale_y_discrete(limits = all_cts) 
  }
  
  return(p)
}



#' Plot Co-localization Dotplot
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the co-localization at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' @param zsig.thresh numeric; the Z score significance threshold (default: 1.96).
#' @param zscore.limit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zscore.limit will be represented as 
#' zscore.limit, scores below -zscore.limit will be represented as -zscore.limit (default: 3).
#' 
selectScalesDat <- function(dat, zsig.thresh, zscore.limit = NULL){
  ## create data.frame with the Z-scores and scales at the first scale
  ## the trend becomes significant
  ## get mean Z
  mean_dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference) %>% 
    dplyr::summarize(Z = mean(Z))
  ## calculate sig z scores
  sig_dat <- mean_dat %>%
    filter(case_when(zsig.thresh >= 0 ~ Z >= zsig.thresh,
                     TRUE ~ Z <= zsig.thresh)) %>% 
    dplyr::group_by(neighbor, reference) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE))
  
  ## only cap if parameter is passed
  ## user could choose not to cap when comparing samples
  if (!is.null(zscore.limit)) {
    ## limit the z-score for the gradient in the figure to look better
    sig_dat$Z[sig_dat$Z > zscore.limit] <- zscore.limit
    sig_dat$Z[sig_dat$Z < -zscore.limit] <- -zscore.limit
  }
  
  return(sig_dat)
}




## Draft -------------------------------------------------------------------

## The idea was to put the difference of Z scores and scales in the same plot.
## The difference in z score would be hue from blue to red and the difference
## in scales would be just a circle with the original scale and a filled one
## with the actual value. This was too much information for a plot, so we chose
## to represent it by fixing Z or scale.

#' Plot Co-localization Dotplot
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the co-localization at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' @param zsig.thresh numeric; the Z score significance threshold (default: 1.96).
#' @param psig.tresh numeric; the two-sided P value significance threshold. It 
#' can be used in place of the zsig.thresh parameter. If no value is provided, 
#' the zsig.thresh will be used.
#' @param zscore.limit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zscore.limit will be represented as 
#' zscore.limit, scores below -zscore.limit will be represented as -zscore.limit (default: 3).
#' @param reorder boolean; if TRUE, reorder the cell types by clustering on the 
#' z-score. If false, orders in alphabetical order (default: FALSE).
#' @param only.significant boolean; plot only cell types with significant relationships (default: TRUE).
#' @param colors character vector; colors for the gradient heatmap (low, mid, high).
#' @param title character; plot title (default: NULL).
#' 
#' @examples 
#' \dontrun{
#' data(sim)
#' cells <- toSF(pos = sim[,c("x", "y")], celltypes = sim$celltypes)
#' shuffle.list <- makeShuffledCells(cells, scales = c(150, 250, 500, 750, 1000, 1500, 2000), ncores = 2)
#' results <- findTrends(cells, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' dat <- meltResultsList(results)
#' vizColocDotplot(dat)
#' }
#' 
#' @export
vizRelDiff <- function(dat1, dat2, 
                       zsig.thresh1 = 1.96, zsig.thresh2 = 1.96, 
                       psig.tresh1 = NULL, psig.tresh2 = NULL,
                       zscore.limit = 3,  reorder = FALSE,
                       only.significant = FALSE,
                       colors = c("blue", "white", "red")){
  
  ## check if cell types are the same
  if (!identical(sort(unique(dat1$reference)), sort(unique(dat2$reference)))) {
    stop("the number of reference cell types should be the same in both datasets")
  }
  if (!identical(sort(unique(dat1$neighbor)), sort(unique(dat2$neighbor)))) {
    stop("the number of neighbor cell types should be the same in both datasets")
  }
  
  ## change pvalues to zscores
  if (!is.null(psig.tresh1)) {
    zsig.thresh1 = round(qnorm(psig.tresh/2, lower.tail = F), 2)
  }
  if (!is.null(psig.tresh2)) {
    zsig.thresh2 = round(qnorm(psig.tresh/2, lower.tail = F), 2)
  }
  
  ## get sig data
  sig_dat1 <- selectSigDat(dat1, zsig.thresh1, zscore.limit)
  sig_dat2 <- selectSigDat(dat2, zsig.thresh2, zscore.limit)
  
  sig_dat1$sample <- '1'
  sig_dat2$sample <- '2'
  
  sig_join <- full_join(sig_dat1, sig_dat2, by = c("neighbor", "reference"),
                        suffix = c("1", "2"))
  
  sig_dif <- sig_join %>% 
    dplyr::mutate(sig_change = 
                    if ( (is.na(Z1) & (! is.na(Z2))) ) {'became significant'} 
                  else if ( (! is.na(Z1) & (is.na(Z2))) ) {'became not significant'}
                  else {'no change in sigificance'}) %>% 
    dplyr::mutate(scale_dif = 
                    if ( sig_change == 'no change in sigificance' ) {scale2 - scale1} 
                  ## we are considering the scale of the non significant as 0
                  else if ( sig_change == 'not significant' ) { scale1 } 
                  else if ( sig_change == 'became significant' ) { scale2 } ) %>% 
    dplyr::mutate(z_dif = 
                    if ( sig_change == 'no change in sigificance' ) {Z2 - Z1} 
                  ## we are considering the zscore of the non significant as 0
                  else if ( sig_change == 'not significant' ) { Z1 } 
                  else if ( sig_change == 'became significant' ) { Z2 } )
  
  ## reorder based on dat1
  if (reorder) {
    ## get values from dat1
    u_cts <- unique(dat1$reference)
    max_scale <- max(dat1$scale)
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat1, by = c('reference', 'neighbor'))
    df_all <- df_all %>% 
      dplyr::mutate(Z = dplyr::coalesce(Z, 0),
                    scale = dplyr::coalesce(scale, max_scale))
    ## create matrix
    sig_mat <- reshape::cast(df_all, neighbor~reference, value='Z')
    rownames(sig_mat) <- sig_mat[,1]
    sig_mat <- sig_mat[,-1]
    sig_mat[is.na(sig_mat)] <- 0
    ## cluster
    hc <- hclust(dist(sig_mat))
    sig_dif$neighbor <- factor(sig_dif$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dif$reference <- factor(sig_dif$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## scale sizes
  lsizes <- sort(unique(sig_dif$scale2))
  legend_sizes <- c(lsizes[1],
                    round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                    lsizes[length(lsizes)])
  
  ## plot
  sig_dif %>% 
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x=reference, y=neighbor, 
                                     size=scale1), shape=1) + 
    ggplot2::geom_point(ggplot2::aes(x=reference, y=neighbor, 
                                     color=z_dif, size=scale2)) + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    ggplot2::scale_colour_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "#eeeeee"
    ) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = c(1, 11)) + 
    ggplot2::scale_x_discrete(position = "top")  + 
    ggplot2::theme_bw()
  
  
  ## plot figure
  p <- sig_dat %>% 
    ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, 
                                 color=Z, size=scale)) +
    ggplot2::geom_point() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    ggplot2::scale_colour_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "#eeeeee"
    ) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = c(1, 11)) + 
    ggplot2::scale_x_discrete(position = "top")  + 
    ggplot2::theme_bw()
  
  if (!only.significant) {
    all_cts <- sort(unique(dat$reference))
    alpha_cts <- sort(unique(dat$reference))
    if (reorder) {
      ct_order <- rownames(sig_mat)[hc$order]
      all_cts <- c(ct_order, setdiff(alpha_cts, ct_order))
    }
    p <- p +
      scale_x_discrete(limits = all_cts, position = 'top') +
      scale_y_discrete(limits = all_cts) 
  }
}


#' Plot Co-localization Dotplot
#' 
#' @description This function takes the `findTrends()` melted data frame and 
#' plots the scale and the Z-score in which the trend first crossed the 
#' significance line (Z = 1.96). The Z-score was capped between -3 and 3. Since 
#' the co-localization at smaller scales are more important than those at 
#' greater ones, we plotted the inverse of the scale so smaller ones would 
#' correspond to larger dots.
#' 
#' @param dat `findTrends()` data.frame; the information about the scale, 
#' Z-score, reference and the neighbor cell. The input data.frame should be the 
#' results list from `findTrends()` that has been melted into a data.frame using `meltResultsList()`.
#' @param zsig.thresh numeric; the Z score significance threshold (default: 1.96).
#' @param zscore.limit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zscore.limit will be represented as 
#' zscore.limit, scores below -zscore.limit will be represented as -zscore.limit (default: 3).
#' 
selectSigDat <- function(dat, zsig.thresh, zscore.limit = NULL){
  ## create data.frame with the Z-scores and scales at the first scale
  ## the trend becomes significant
  ## get mean Z
  mean_dat <- dat %>% 
    dplyr::group_by(neighbor, scale, reference) %>% 
    dplyr::summarize(Z = mean(Z))
  ## calculate sig z scores
  sig_dat <- mean_dat %>%
    dplyr::filter(abs(Z) >= zsig.thresh) %>% 
    dplyr::group_by(neighbor, reference) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE))
  
  ## only cap if parameter is passed
  ## user could choose not to cap when comparing samples
  if (!is.null(zscore.limit)) {
    ## limit the z-score for the gradient in the figure to look better
    sig_dat$Z[sig_dat$Z > zscore.limit] <- zscore.limit
    sig_dat$Z[sig_dat$Z < -zscore.limit] <- -zscore.limit
  }
  
  return(sig_dat)
}




## Compare mutual relationships --------------------------------------------


### Draft -------------------------------------------------------------------

library(tidyverse)
library(crawdad)

dat_vhck <- readRDS('running_code/processed_data/thymus/dat_vhck_50.RDS')
dat_ktjk <- readRDS('running_code/processed_data/thymus/dat_ktjk_50.RDS')

ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')
zsig <- correctZBonferroni(dat_vhck)

## Dotplot vhck
dat_vhck %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct') %>% 
  vizColocDotplot(zsigThresh = zsig, zscoreLimit = zsig*2, 
                  reorder = TRUE, mutual = T, dotSizes = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')

## Dotplot ktjk
dat_ktjk %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct') %>% 
  vizColocDotplot(zsigThresh = zsig, zscoreLimit = zsig*2, 
                  reorder = TRUE, mutual = T, dotSizes = c(2, 14)) +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='bottom',
        axis.text.x = element_text(angle = 45, h = 0),
        legend.box = 'vertical')

## define relationships
dat_ktjk <- dat_ktjk %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct') %>% 
  dplyr::group_by(neighbor, scale, reference) %>% 
  dplyr::summarize(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE)) %>% 
  dplyr::mutate(relationship = case_when(Z > 0 ~ 'enrichment',
                                         Z < 0 ~ 'depletion',
                                         T ~ 'other'))
dat_vhck <- dat_vhck %>% 
  filter(neighbor != 'indistinct') %>% 
  filter(reference != 'indistinct') %>% 
  dplyr::group_by(neighbor, scale, reference) %>% 
  dplyr::summarize(Z = mean(Z)) %>% 
  dplyr::filter(abs(Z) >= zsig) %>% 
  dplyr::group_by(neighbor, reference) %>% 
  dplyr::filter(scale == min(scale, na.rm = TRUE)) %>% 
  dplyr::mutate(relationship = case_when(Z > 0 ~ 'enrichment',
                                         Z < 0 ~ 'depletion',
                                         T ~ 'other'))


#### Enrichment --------------------------------------------------------------

## join
samples <- c('_ktjk', '_vhck')
merged_dat <- full_join(dat_ktjk, dat_vhck, by = c('neighbor', 'reference'), 
          suffix = samples) %>% 
  mutate(mutual = (relationship_ktjk == 'enrichment') & 
           (relationship_vhck == 'enrichment')) %>% 
  mutate(mean_scale = (scale_ktjk + scale_vhck)/2)

## scale sizes
lsizes <- sort(unique(merged_dat$scale))
legend_sizes <- c(lsizes[1],
                  round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                  lsizes[length(lsizes)])
## order
ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')
## plot
merged_dat %>% 
  filter(pct == T) %>% 
  ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, size=mean_scale)) +
  ggplot2::geom_point(color='darkgreen') + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     vjust = 0.5, 
                                                     hjust=1)) +
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = legend_sizes,
                        range = c(5, 15)) + 
  ggplot2::scale_x_discrete(position = "top") + 
  ggplot2::theme_bw() +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



#### Depletion ---------------------------------------------------------------

## join
samples <- c('_ktjk', '_vhck')
merged_dat <- full_join(dat_ktjk, dat_vhck, by = c('neighbor', 'reference'), 
                        suffix = samples) %>% 
  mutate(mutual = (relationship_ktjk == 'depletion') & 
           (relationship_vhck == 'depletion')) %>% 
  mutate(mean_scale = (scale_ktjk + scale_vhck)/2)

## scale sizes
lsizes <- sort(unique(merged_dat$scale))
legend_sizes <- c(lsizes[1],
                  round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                  lsizes[length(lsizes)])
## order
ct_order <- readRDS('running_code/processed_data/ct_order_thymus.RDS')
## plot
merged_dat %>% 
  filter(mutual == T) %>% 
  ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, size=mean_scale)) +
  ggplot2::geom_point(color='yellow') + 
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                     vjust = 0.5, 
                                                     hjust=1)) +
  ggplot2::scale_radius(trans = 'reverse',
                        breaks = legend_sizes,
                        range = c(5, 15)) + 
  ggplot2::scale_x_discrete(position = "top") + 
  ggplot2::theme_bw() +
  scale_x_discrete(limits = ct_order, position = 'top') +
  scale_y_discrete(limits = ct_order, position = 'right') +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))


