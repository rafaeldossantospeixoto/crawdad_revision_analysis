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
#' results list from `findTrends()` that has been melted into a data.frame 
#' using `meltResultsList()`.
#' @param zsigThresh numeric; the Z score significance threshold (default: 1.96).
#' @param psigThresh numeric; the two-sided P value significance threshold. It 
#' can be used in place of the zsigThresh parameter. If no value is provided, 
#' the zsigThresh will be used.
#' @param zscoreLimit numeric; limit the Z-score to look better in the graph 
#' scale gradient. Z-score values above zscoreLimit will be represented as 
#' zscoreLimit, scores below -zscoreLimit will be represented as -zscoreLimit
#' (default: NULL).
#' @param reorder boolean; if TRUE, reorder the cell types by clustering on the 
#' z-score. If false, orders in alphabetical order (default: FALSE).
#' @param mutual boolean; highligh relationships that are mutual between cell 
#' type pairs (default: TRUE).
#' @param onlySignificant boolean; plot only cell types with significant 
#' relationships (default: TRUE).
#' @param colors character vector; colors for the gradient heatmap (low, mid, high) 
#' (default: c("blue", "white", "red")).
#' @param dotSizes numeric vector; minimum and maximum size of the dot 
#' (default: c(6,31)). 
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
vizColocDotplot <- function(dat, zsigThresh = 1.96, psigThresh = NULL,
                            zscoreLimit = NULL,  reorder = FALSE,
                            mutual = FALSE, onlySignificant = FALSE,
                            colors = c("blue", "white", "red"),
                            dotSizes = c(6,31),
                            title = NULL){
  if (!is.null(psigThresh)) {
    zsigThresh = round(qnorm(psigThresh/2, lower.tail = F), 2)
  }
  
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
    dplyr::filter(abs(Z) >= zsigThresh) %>% 
    dplyr::group_by(neighbor, reference) %>% 
    dplyr::filter(scale == min(scale, na.rm = TRUE))
  
  ## limit the z-score for the gradient in the figure to look better
  if (!is.null(zscoreLimit)) {
    sig_dat$Z[sig_dat$Z > zscoreLimit] <- zscoreLimit
    sig_dat$Z[sig_dat$Z < -zscoreLimit] <- -zscoreLimit
  }
  
  ## scale sizes
  lsizes <- sort(unique(sig_dat$scale))
  legend_sizes <- c(lsizes[1],
                    round(mean(c(lsizes[1], lsizes[length(lsizes)]))),
                    lsizes[length(lsizes)])
  ## reorder
  if (reorder) {
    ## merge with all cts
    comb_cts <- tidyr::expand_grid(u_cts, u_cts)
    colnames(comb_cts) <- c('reference', 'neighbor')
    df_all <- dplyr::left_join(comb_cts, sig_dat, by = c('reference', 'neighbor'))
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
    sig_dat$neighbor <- factor(sig_dat$neighbor, 
                               levels=rownames(sig_mat)[hc$order])
    sig_dat$reference <- factor(sig_dat$reference, 
                                levels=colnames(sig_mat)[hc$order])
  }
  
  ## highligh mutual
  if (mutual) {
    # ## collect all pairs
    # ref_cts <- as.character(sig_dat$reference)
    # ngb_cts <- as.character(sig_dat$neighbor)
    # ## get Z scores and derive type of relationship to compare if it is the same
    # zscores <- sig_dat$Z
    df_pairs <- sig_dat %>% 
      select(reference, neighbor, Z) %>% 
      ## calculate the type of relationship
      mutate(type = case_when(Z > 0 ~ 'enrichment',
                              Z < 0 ~ 'depletion',
                              T ~ NA)) %>% 
      mutate(pair = paste(sort(c(gsub(" ", "", reference), 
                                 gsub(" ", "", neighbor))), collapse = '_'))
    df_same_type <- df_pairs %>% 
      group_by(pair) %>% 
      ## check if the type is not different for each ref of the pair and
      ## check if there are two relationships by checking distinct references
      summarise(same_type = (n_distinct(type) != 2) & 
                  (n_distinct(reference) == 2))
    ## merge to reorder
    df_pairs <- df_pairs %>% 
      left_join(df_same_type, by = 'pair')
    ## check if pairs are duplicate
    mutual_same_relationships <- df_pairs$same_type
    sig_dat$mutual <- mutual_same_relationships
  }
  
  ## plot figure
  p <- sig_dat %>% 
    ggplot2::ggplot(ggplot2::aes(x=reference, y=neighbor, 
                                 color=Z, size=scale)) +
    ggplot2::geom_point() + 
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, 
                                                       vjust = 0.5, 
                                                       hjust=1)) +
    {if(mutual)ggplot2::geom_point(data = ~dplyr::filter(.x, mutual == T),
                                   ggplot2::aes(x=reference, y=neighbor),
                                   shape = 8, color = 'gold', size = 2*dotSizes[1]/3)} + 
    ggplot2::scale_colour_gradient2(
      low = colors[1],
      mid = colors[2],
      high = colors[3],
      na.value = "lightgray"
    ) + 
    ggplot2::scale_radius(trans = 'reverse',
                          breaks = legend_sizes,
                          range = dotSizes) + 
    ggplot2::scale_x_discrete(position = "top") + 
    ggplot2::theme_bw()
  
  if (!onlySignificant) {
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



# Run ---------------------------------------------------------------------

library(crawdad)
library(tidyverse)

dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')
zsig <- correctZBonferroni(dat_50)

vizColocDotplot(dat_50, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = T,
                     mutual = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))
# p
# pdf('function_development/embryo_dotplot_crawdad.pdf',
#     height = 7.5, width = 9)
# p
# dev.off()

## vizColocDotplot will be deprecated
vizRelationships(dat_50, zSigThresh = zsig, zScoreLimit = 2*zsig, reorder = T,
                 symmetrical = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

