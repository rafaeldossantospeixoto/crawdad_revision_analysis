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
                  else if ( sig_change == 'not significant' ) { scale1 } ## we are considering the scale of the non significant as 0
                  else if ( sig_change == 'became significant' ) { scale2 } ) %>% 
    dplyr::mutate(z_dif = 
                    if ( sig_change == 'no change in sigificance' ) {Z2 - Z1} 
                  else if ( sig_change == 'not significant' ) { Z1 } ## we are considering the zscore of the non significant as 0
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
  ## get values before filtering
  max_scale <- max(dat$scale)
  u_cts <- unique(dat$reference)
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
