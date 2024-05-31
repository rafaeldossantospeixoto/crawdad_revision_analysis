## funcitons that need to be copied to the CRAWDAD package



# From compare_samples_scale ----------------------------------------------

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




# From compare_samples_auc ------------------------------------------------


create_id_column <- function(df){
  df_ids <- df %>%
    dplyr::mutate(id = paste0('s', slice, 'r', replicate))
  return(df_ids)
}
create_pair_column <- function(df){
  df_pairs <- df %>%
    dplyr::mutate(pair = paste0(paste0(reference, ' - ', neighbor)))
  return(df_pairs)
}
## filter data based on common values
filter_shared_pairs <- function(df_pairs){
  shared_pairs <- df_pairs %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarize(n_unique = dplyr::n_distinct(id)) %>% 
    dplyr::filter(n_unique == dplyr::n_distinct(df_pairs$id)) %>% 
    dplyr::pull(pair)
  
  filtered_df <- df_pairs %>% 
    filter(pair %in% shared_pairs)
  
  return(filtered_df)
}

