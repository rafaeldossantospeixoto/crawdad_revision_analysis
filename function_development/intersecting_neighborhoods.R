library(crawdad)
library(tidyverse)


# Embryo ------------------------------------------------------------------

## Without dups ---------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_seq_50.RDS')
zsig <- correctZBonferroni(dat_50)


## vizColocDotplot
vizColocDotplot(dat_50, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = T,
                mutual = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

## vizSigTrends
## plot significant trends
all_cts <- unique(dat_50$reference)
all_cts <- c('Presomitic mesoderm')
for (ref_ct in all_cts) {
  ## spat
  p1 <- vizClusters(seq, ofInterest = c(ref_ct), alpha = 1) + 
    theme_void() +
    theme(legend.position = "none")
  
  
  ## trends
  sig_cts <- dat_50 %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(abs(Z) > zsig) %>% 
    pull(neighbor) %>% 
    unique() %>% 
    as.character()
  not_sig_cts <- all_cts[! all_cts %in% sig_cts]
  
  filtered_dat <- dat_50 %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(neighbor %in% sig_cts)
  # mutate(selected_neighbor = fct_other(neighbor, keep = sig_cts)) %>% 
  p2 <- filtered_dat %>% ggplot() +
    geom_line(aes(scale, Z, color = neighbor)) +
    scale_color_manual(values = rainbow(length(sig_cts))) +
    geom_hline(yintercept = zsig, linetype = 2) + 
    geom_hline(yintercept = -zsig, linetype = 2) +
    xlim(c(50, 1250)) + 
    ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                                 label = neighbor), 
                             data = filtered_dat %>%
                               filter(scale == max(scale)),
                             hjust = 'right', box.padding = 0.5,
                             nudge_x = .5, nudge_y = .5,
                             segment.alpha = .5,
                             xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
    labs(title = ref_ct)
  
  p <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1/3, 2/3))
  
  print(p2)
  # png(paste0('running_code/exploring_embryo/trends_ref_', 
  #            gsub('[ /]', '_', ref_ct), '.png'),
  #     height = 700, width = 1400)
  # print(p)
  # dev.off()
}



## With dups ---------------------------------------------------

ncores <- 7

data(seq)
seq %>% ggplot() +
  geom_point(aes(x, y, color = celltypes), size = .5) + 
  facet_wrap('celltypes')

## convert to sf
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                      celltypes = seq$celltypes)

# scales <- seq(100, 1000, by=50)
# ## generate background
# shuffle.list <- crawdad:::makeShuffledCells(seq,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
saveRDS(shuffle.list, 'running_code/processed_data/shufflelist_seq.RDS')
## No need to edit the function, just set removeDups to FALSE.
## find trends, dist 50
# results_50_dups <- crawdad::findTrends(seq,
#                                   dist = 50,
#                                   shuffle.list = shuffle.list,
#                                   ncores = ncores,
#                                   verbose = TRUE,
#                                   returnMeans = FALSE,
#                                   removeDups = FALSE)
# ## Time was 17.62 mins
# dat_50_dups <- crawdad::meltResultsList(results_50_dups, withPerms = T)
saveRDS(dat_50_dups, 'running_code/processed_data/dat_seq_50_dups.RDS')
dat_50_dups <- readRDS('running_code/processed_data/dat_seq_50_dups.RDS')


## vizColocDotplot
zsig <- correctZBonferroni(dat_50_dups)
vizColocDotplot(dat_50_dups, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = T,
                mutual = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

## vizSigTrends
## plot significant trends
all_cts <- unique(dat_50_dups$reference)
all_cts <- c('Definitive endoderm', 'Cranial mesoderm', 'Endothelium')
for (ref_ct in all_cts) {
  ## spat
  # p1 <- vizClusters(seq, ofInterest = c(ref_ct), alpha = 1) + 
  #   theme_void() +
  #   theme(legend.position = "none")
  
  
  ## trends
  sig_cts <- dat_50_dups %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(abs(Z) > zsig) %>% 
    pull(neighbor) %>% 
    unique() %>% 
    as.character()
  not_sig_cts <- all_cts[! all_cts %in% sig_cts]
  
  filtered_dat <- dat_50_dups %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    # filter(neighbor %in% sig_cts) %>% 
    filter(neighbor %in% all_cts) %>% 
    filter(neighbor != ref_ct) #%>% ## remove reference cell type
    # filter(scale <= 500) ## remove scale greater than 5000
  # mutate(selected_neighbor = fct_other(neighbor, keep = sig_cts)) %>% 
  p2 <- filtered_dat %>% 
    ggplot() +
    geom_line(aes(scale, Z, color = neighbor)) +
    scale_color_manual(values = rainbow(length(sig_cts))) +
    geom_hline(yintercept = zsig, linetype = 2) + 
    geom_hline(yintercept = -zsig, linetype = 2) +
    xlim(c(50, 1250)) + 
    ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                                 label = neighbor), 
                             data = filtered_dat %>%
                               filter(scale == max(scale)),
                             hjust = 'right', box.padding = 0.5,
                             nudge_x = .5, nudge_y = .5,
                             segment.alpha = .5,
                             xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
    labs(title = ref_ct)
  
  # p <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1/3, 2/3))
  
  print(p2)
  # png(paste0('running_code/exploring_embryo/trends_ref_', 
  #            gsub('[ /]', '_', ref_ct), '.png'),
  #     height = 700, width = 1400)
  # print(p)
  # dev.off()
}




# Sim ------------------------------------------------------------------


## Without dups ------------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_sim_50.RDS')
zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = F,
                mutual = T, dotSizes = c(20, 50))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

ref_ct <- 'A'
filtered_dat <- dat_50 %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == ref_ct) 
  # filter(!(neighbor %in% ref_ct))
filtered_dat %>%   
  ggplot() +
  geom_line(aes(scale, Z, color = neighbor)) +
  scale_color_manual(values = rainbow(length(unique(filtered_dat$neighbor)))) +
  geom_hline(yintercept = zsig, linetype = 2) + 
  geom_hline(yintercept = -zsig, linetype = 2) +
  xlim(c(50, 1250)) + 
  ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                               label = neighbor), 
                           data = filtered_dat %>%
                             filter(scale == max(scale)),
                           hjust = 'right', box.padding = 0.5,
                           nudge_x = .5, nudge_y = .5,
                           segment.alpha = .5,
                           xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
  labs(title = ref_ct)

## With dups ---------------------------------------------------------------

data(sim)
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)
scales <- seq(100, 1000, by=50)
# shuffle.list <- crawdad:::makeShuffledCells(sim,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
saveRDS(shuffle.list, 'running_code/processed_data/shufflelist_sim.RDS')
## find trends, dist 50
# results_50_dups <- crawdad::findTrends(sim,
#                                        dist = 50,
#                                        shuffle.list = shuffle.list,
#                                        ncores = ncores,
#                                        verbose = TRUE,
#                                        returnMeans = FALSE,
#                                        removeDups = FALSE)
# dat_50_dups <- crawdad::meltResultsList(results_50_dups, withPerms = T)
saveRDS(dat_50_dups, 'running_code/processed_data/dat_sim_50_dups.RDS')

## vizColocDotplot
zsig <- correctZBonferroni(dat_50_dups)
vizColocDotplot(dat_50_dups, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = F,
                mutual = T, dotSizes = c(20, 50))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

ref_ct <- 'A'
filtered_dat <- dat_50_dups %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == ref_ct) 
# filter(!(neighbor %in% ref_ct))
filtered_dat %>%   
  ggplot() +
  geom_line(aes(scale, Z, color = neighbor)) +
  scale_color_manual(values = rainbow(length(unique(filtered_dat$neighbor)))) +
  geom_hline(yintercept = zsig, linetype = 2) + 
  geom_hline(yintercept = -zsig, linetype = 2) +
  xlim(c(50, 1250)) + 
  ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                               label = neighbor), 
                           data = filtered_dat %>%
                             filter(scale == max(scale)),
                           hjust = 'right', box.padding = 0.5,
                           nudge_x = .5, nudge_y = .5,
                           segment.alpha = .5,
                           xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
  labs(title = ref_ct)



## With dups and ref -------------------------------------------------------

data(sim)
sim <- crawdad:::toSF(pos = sim[,c("x", "y")],
                      celltypes = sim$celltypes)
scales <- seq(100, 1000, by=50)
shuffle.list <- readRDS('running_code/processed_data/shufflelist_sim.RDS')
## find trends, dist 50
results_50_dups_ref <- findTrendsRef(sim,
                                       dist = 50,
                                       shuffle.list = shuffle.list,
                                       ncores = ncores,
                                       verbose = TRUE,
                                       returnMeans = FALSE,
                                       removeDups = FALSE)
dat_50_dups_ref <- crawdad::meltResultsList(results_50_dups_ref, withPerms = T)
saveRDS(dat_50_dups_ref, 'running_code/processed_data/dat_sim_50_dups_ref.RDS')

## vizColocDotplot
zsig <- correctZBonferroni(dat_50_dups_ref)
zsig <- 1.96
vizColocDotplot(dat_50_dups_ref, zsigThresh = zsig, zscoreLimit = 2*zsig, reorder = F,
                mutual = T, dotSizes = c(20, 50))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

ref_ct <- 'B'
filtered_dat <- dat_50_dups_ref %>% 
  group_by(reference, neighbor, scale) %>% 
  summarize(Z = mean(Z)) %>% 
  filter(reference == ref_ct) 
# filter(!(neighbor %in% ref_ct))
filtered_dat %>%   
  ggplot() +
  geom_line(aes(scale, Z, color = neighbor)) +
  scale_color_manual(values = rainbow(length(unique(filtered_dat$neighbor)))) +
  geom_hline(yintercept = zsig, linetype = 2) + 
  geom_hline(yintercept = -zsig, linetype = 2) +
  xlim(c(50, 1250)) + 
  ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                               label = neighbor), 
                           data = filtered_dat %>%
                             filter(scale == max(scale)),
                           hjust = 'right', box.padding = 0.5,
                           nudge_x = .5, nudge_y = .5,
                           segment.alpha = .5,
                           xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
  labs(title = ref_ct)




# Cerebellum --------------------------------------------------------------


## Without dups ------------------------------------------------------------

dat_50 <- readRDS('running_code/processed_data/dat_slide_50.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50, reorder = TRUE, zsigThresh = zsig, 
                zscoreLimit = zsig*2, mutual = T,
                dotSizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))



## Without dups ------------------------------------------------------------

ncores <- 7
data(slide)
slide <- crawdad:::toSF(pos = slide[,c("x", "y")],
                        celltypes = slide$celltypes)
scales <- seq(100, 1000, by=100)
# shuffle.list <- crawdad:::makeShuffledCells(slide,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
saveRDS(shuffle.list, 'running_code/processed_data/shufflelist_slide.RDS')
## find trends, dist 50
# results_50_dups <- crawdad::findTrends(slide,
#                                   dist = 50,
#                                   shuffle.list = shuffle.list,
#                                   ncores = ncores,
#                                   verbose = TRUE,
#                                   returnMeans = FALSE, 
#                                   removeDups = FALSE)
# dat_50_dups <- crawdad::meltResultsList(results_50_dups, withPerms = T)
saveRDS(dat_50_dups, 'running_code/processed_data/dat_slide_50_dups.RDS')

vizColocDotplot(dat_50_dups, reorder = TRUE, zsigThresh = zsig, 
                zscoreLimit = zsig*2, mutual = T,
                dotSizes = c(2, 14)) +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))




# CRAWDAD functions -------------------------------------------------------

## No need to edit the function, just set removeDups to FALSE.

#' Compute significant different between real and randomly shuffled cell neighbor proportions
#' 
#' @description for a given reference cell type, computes the Z scores for being 
#' colocalized or separated from each query cell type. If `returnMeans = TRUE`, 
#' then the result will be a data.frame where each row is a scale, each column 
#' is a query cell type, and each value is the Z score. If 
#' `returnMeans = FALSE`, then the result will be a list of data.frames, where 
#' each data.frame is a scale,the rows are permutations, the columns are query 
#' cell types, and each value is a Z score.
#'
#' @param cells sf object of all the cells
#' @param randomcellslist list of lists of randomly shuffled cell type labels 
#' produced from `makeShuffledCells`
#' @param trueNeighCells Simple feature collection of real cells for a given 
#' reference cell type, with geometries of a given dist (from sf::st_buffer)
#' @param cellBuffer Simple feature collection of the neighbor cells that are 
#' within "dist" of the ref cells (from sf::intersection)
#' @param ncores number of cores for parallelization (default 1)
#' @param removeDups remove duplicate neighbor cells to prevent them from being 
#' counted multiple times and inflate the Z scores (default: TRUE)
#' @param returnMeans if multiple permutations, return the mean Z score across 
#' the permutations in each scale with respect to each neighbor cell type 
#' (default: TRUE) 
#'
#'@export
evaluateSignificance <- function(cells,
                                 randomcellslist,
                                 trueNeighCells,
                                 cellBuffer,
                                 ncores = 1,
                                 removeDups = TRUE,
                                 returnMeans = TRUE){
  
  allcells <- cells
  trueNeighCells <- trueNeighCells
  cellBuffer <- cellBuffer
  
  ## if true, will return a data.frame
  if(returnMeans){
    
    ## for each scale:
    results <- do.call(rbind, BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the mean Z score across permutations for the given scale
      return(colMeans(scores))
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
    
    ## otherwise, returns a list
  } else {
    
    ## for each scale:
    results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        y1 <- table(trueNeighCells$celltypes)
        y2 <- table(bufferrandomcells$celltypes)
        n1 <- length(trueNeighCells$celltypes)
        n2 <- length(bufferrandomcells$celltypes)
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the data.frame of Z scores for each permutation (row)
      return(scores)
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
    
    ## will be list of data.frame in this case
    ## each data.frame is a scale
    names(results) <- names(randomcellslist)
    
  }
  
  rm(allcells)
  rm(trueNeighCells)
  rm(cellBuffer)
  gc(verbose = FALSE, reset = TRUE)
  
  return(results)
}





#' Compute trends of cell type colocalization for each cell type combination across specified scales
#'
#' @description Trends are based on significant differences in cell type proportions 
#' between the real and randomly shuffled datasets.
#' Cell type proportions are with respect to the different cell types that are 
#' neighboring the cells of a given reference cell type within a certain defined distance.
#' This is done at difference scales, where a scale is whether the cell type 
#' labels are shuffled locally or globally.
#' Trends are essentially built from significance values. The significance test 
#' basically asks if two cell types are localized or separated by assessing if 
#' the proportion of the neighboring cell type is significantly greater, or less 
#' than, random chance.
#'
#' @param cells sf object, with celltypes features and point geometries
#' @param dist numeric distance to define neighbor cells with respect to each 
#' reference cell (default: 100)
#' @param ncores number of cores for parallelization (default 1)
#' @param shuffle.list a list of cell type labels shuffled at different scales 
#' (output from `makeShuffledCells()`)
#' @param subset.list a subset list (output from `selectSubsets()`). Required if 
#' computing trends for subsets (default NULL)
#' @param verbose Boolean for verbosity (default TRUE)
#' @param removeDups remove duplicate neighbor cells to prevent them from being 
#' counted multiple times and inflate the Z scores (default: TRUE)
#' @param returnMeans if multiple permutations, return the mean Z score across 
#' the permutations in each scale with respect to each neighbor cell type (default: TRUE)
#'
#' @return A list that contains a dataframe for each reference cell type, where 
#' the dataframe contains the significance values for each neighbor cell type at each scale
#' 
#' @examples
#' \dontrun{
#' data(sim)
#' shuffle.list <- makeShuffledCells(sim, scales = c(50, 100, 200, 300, 400, 500))
#' findTrends(sim, dist = 100, shuffle.list = shuffle.list, ncores = 2)
#' }
#' 
#' @export
findTrends <- function(cells,
                       dist = 100,
                       ncores = 1,
                       shuffle.list,
                       subset.list = NULL,
                       verbose = TRUE,
                       removeDups = TRUE,
                       returnMeans = TRUE){
  
  if(!is.list(shuffle.list)){
    stop("`shuffle.list` is not a list. You can make this using `makeShuffledCells()`")
  }
  
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  sf::st_agr(cells) <- "constant"
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  if(is.null(subset.list)){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    celltypes <- factor(cells$celltypes)
    
    d <- dist
    results.all <- lapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      # use this to assess neighbors within "d" um of each cell
      ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans)
      return(results)
    }) 
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
  }
  
  if(!is.null(subset.list)){
    
    ## load in the subset file if it exists, or make it and probably be a good
    ## idea to save it, too
    if(!is.list(subset.list)){
      stop(paste0("`subset.list` is not a list. You can build this using: `binomialTestMatrix()` then `selectSubsets()`"))
    }
    
    combo_ids <- names(subset.list)
    d <- dist
    
    if(verbose){
      message("Calculating trends for each subset in `subset.list` with respect to the cell types in `cells$celltypes`")
    }
    
    ## initialize list
    results.all <- list()
    
    ## for each subset of cells, evaluate significance against each reference cell type across scales
    for(i in combo_ids){
      
      if(verbose){
        message(i)
      }
      
      ## get area around the subset cells to identify neighbors
      ref.buffer <- sf::st_buffer(cells[subset.list[[i]], ], d)
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        ## remove only the original cells, as the duplicates will have 
        ## .something in the name
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificance(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans)
      results.all[[i]] <- results
      
      rm(ref.buffer)
      rm(neigh.cells)
      rm(results)
      gc(verbose = FALSE, reset = TRUE)
    }
    
  }
  
  ## return results
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  return(results.all)
  
}



# Testing functions -------------------------------------------------------

## Params ------------------------------------------------------------------

ncores <- 7
data(seq)
cells <- crawdad:::toSF(pos = seq[,c("x", "y")], celltypes = seq$celltypes)
dist = 50
verbose = TRUE
returnMeans = FALSE
removeDups = FALSE
shuffle.list <- readRDS('running_code/processed_data/shufflelist_seq.RDS')


## findTrends ----------------------------------------------------------------

sf::st_agr(cells) <- "constant"
celltypes <- factor(cells$celltypes)
d <- dist

ct <- 'Definitive endoderm'
vizClusters(cells, ofInterest = ct)
results.all <- lapply(levels(celltypes), function(ct) {
  
  # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
  # use this to assess neighbors within "d" um of each cell
  ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
  ## union polygons to avoid memory overflow, too slow
  # ref.buffer_union <- sf::st_union(ref.buffer)
  # get the different types of neighbor cells that are within "d" of the ref cells
  neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
  
  ## remove duplicate neighbor cells to prevent them from being counted multiple times
  ## and inflate the Z scores
  if(removeDups){
    ## need to remove self cells else have trivial enrichment when d~0
    self.cells <- cells[cells$celltypes == ct,]
    ## remove only the original cells, as the duplicates will have 
    ## .something in the name
    neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
    
    # message("number of neighbor cells before: ", nrow(neigh.cells))
    #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
    ## hack to accommodate self cells that are neighbors of another self cell
    neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
    # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
  }
  
  ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
  ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
  ## I think a bottle neck originally was waiting for certain cell types to finish
  ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
  ## I could also split up the permutations, but then each cell type for each scale is done one by one
  ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one

  ## not run
  results <- evaluateSignificance(cells = cells,
                                  randomcellslist = shuffle.list,
                                  trueNeighCells = neigh.cells,
                                  cellBuffer = ref.buffer,
                                  ncores = ncores,
                                  removeDups = removeDups,
                                  returnMeans = returnMeans)
  return(results)
}) 
names(results.all) <- levels(celltypes)


## evaluateSignificance ----------------------------------------------------

randomcellslist = shuffle.list
trueNeighCells = neigh.cells
cellBuffer = ref.buffer

allcells <- cells
trueNeighCells <- trueNeighCells
cellBuffer <- cellBuffer

cellsAtRes <- randomcellslist[['250']]
## for each scale:
results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
  
  ## Iterate through each permutation of a given scale and 
  ## produce the scores for each neighbor cell type. 
  ## Scores for a given permutation added as a row to a dataframe.
  ## Take the column mean of the scores for each neighbor cell type across permutations.
  ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
  ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
  
  randomcellslabels <- cellsAtRes[['3']]
  ## for each res
  scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
    
    randomcells <- allcells
    randomcells$celltypes <- as.factor(randomcellslabels)
    sf::st_agr(randomcells) <- "constant"
    
    ## should it be 19416? it's the same as the total number of cells
    bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
    
    ## should the self be removed from here too?
    ## remove duplicate neighbor cells to prevent them from being counted multiple times
    ## and inflate the Z scores
    if(removeDups){
      # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
      bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
      # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
    }
    
    ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
    n_ref <- table(cells$celltypes)[[ct]]
    y1 <- table(trueNeighCells$celltypes)
    y2 <- table(bufferrandomcells$celltypes)
    n1 <- length(trueNeighCells$celltypes)
    n2 <- length(bufferrandomcells$celltypes)
    p1 <- y1/n1
    p2 <- y2/n2
    p <- (y1+y2)/(n1+n2)
    Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
    
    rm(bufferrandomcells)
    rm(randomcells)
    gc(verbose = FALSE, reset = TRUE)
    
    return(Z)
  }))
  
  ## returning the data.frame of Z scores for each permutation (row)
  return(scores)
  
}, BPPARAM=BiocParallel::SnowParam(workers=ncores))

## will be list of data.frame in this case
## each data.frame is a scale
names(results) <- names(randomcellslist)

rm(allcells)
rm(trueNeighCells)
rm(cellBuffer)
gc(verbose = FALSE, reset = TRUE)

return(results)
results.all



# Divide by n ref cells -------------------------------------------------------

evaluateSignificanceRef <- function(cells,
                                 randomcellslist,
                                 trueNeighCells,
                                 cellBuffer,
                                 ncores = 1,
                                 removeDups = TRUE,
                                 returnMeans = TRUE,
                                 ct){
  
  allcells <- cells
  trueNeighCells <- trueNeighCells
  cellBuffer <- cellBuffer
  
  ## if true, will return a data.frame
  if(returnMeans){
    
    ## for each scale:
    results <- do.call(rbind, BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        n_ref <- table(cells$celltypes)[[ct]]
        # print(n_ref)
        y1 <- table(trueNeighCells$celltypes) /n_ref
        y2 <- table(bufferrandomcells$celltypes) /n_ref
        n1 <- length(trueNeighCells$celltypes) /n_ref
        n2 <- length(bufferrandomcells$celltypes) /n_ref
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the mean Z score across permutations for the given scale
      return(colMeans(scores))
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores)))
    
    ## otherwise, returns a list
  } else {
    
    ## for each scale:
    results <- BiocParallel::bplapply(randomcellslist, function(cellsAtRes){
      
      ## Iterate through each permutation of a given scale and 
      ## produce the scores for each neighbor cell type. 
      ## Scores for a given permutation added as a row to a dataframe.
      ## Take the column mean of the scores for each neighbor cell type across permutations.
      ## The resulting vector of score means is returned and appended as a row in the `results` data,frame,
      ## where each row is a scale, and contains Z scores for each neighbor cell type (the columns)
      
      scores <- do.call(rbind, lapply(cellsAtRes, function(randomcellslabels){
        
        randomcells <- allcells
        randomcells$celltypes <- as.factor(randomcellslabels)
        sf::st_agr(randomcells) <- "constant"
        
        bufferrandomcells <- sf::st_intersection(randomcells, cellBuffer$geometry)
        
        ## remove duplicate neighbor cells to prevent them from being counted multiple times
        ## and inflate the Z scores
        if(removeDups){
          # message("number of permuted neighbor cells before: ", nrow(bufferrandomcells))
          bufferrandomcells <- bufferrandomcells[intersect(rownames(bufferrandomcells), rownames(randomcells)),]
          # message("number of permuted neighbor cells after removing dups: ", nrow(bufferrandomcells))
        }
        
        ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
        n_ref <- table(cells$celltypes)[[ct]]
        # print(n_ref)
        y1 <- table(trueNeighCells$celltypes) /n_ref
        y2 <- table(bufferrandomcells$celltypes) /n_ref
        n1 <- length(trueNeighCells$celltypes) /n_ref
        n2 <- length(bufferrandomcells$celltypes) /n_ref
        p1 <- y1/n1
        p2 <- y2/n2
        p <- (y1+y2)/(n1+n2)
        Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
        
        rm(bufferrandomcells)
        rm(randomcells)
        gc(verbose = FALSE, reset = TRUE)
        
        return(Z)
      }))
      
      ## returning the data.frame of Z scores for each permutation (row)
      return(scores)
      
    }, BPPARAM=BiocParallel::SnowParam(workers=ncores))
    
    ## will be list of data.frame in this case
    ## each data.frame is a scale
    names(results) <- names(randomcellslist)
    
  }
  
  rm(allcells)
  rm(trueNeighCells)
  rm(cellBuffer)
  gc(verbose = FALSE, reset = TRUE)
  
  return(results)
}


findTrendsRef <- function(cells,
                       dist = 100,
                       ncores = 1,
                       shuffle.list,
                       subset.list = NULL,
                       verbose = TRUE,
                       removeDups = TRUE,
                       returnMeans = TRUE){
  
  if(!is.list(shuffle.list)){
    stop("`shuffle.list` is not a list. You can make this using `makeShuffledCells()`")
  }
  
  if( !any(class(cells) == "sf") ){
    stop("`cells` needs to be an `sf` object. You can make this using `toSF()`")
  }
  
  if( !any(grepl("celltypes", colnames(cells))) ){
    stop("`cells` needs a column named `celltypes`. You can make this using `toSF()`")
  }
  
  if(verbose){
    start_time <- Sys.time()
  }
  
  sf::st_agr(cells) <- "constant"
  
  if(verbose){
    message("Evaluating significance for each cell type")
    message("using neighbor distance of ", dist)
  }
  
  ## Evaluate significance (pairwise)
  if(is.null(subset.list)){
    
    if(verbose){
      message("Calculating for pairwise combinations")
    }
    
    celltypes <- factor(cells$celltypes)
    
    d <- dist
    results.all <- lapply(levels(celltypes), function(ct) {
      
      if(verbose){
        message(ct)
      }
      
      # get polygon geometries of reference cells of "celltype" up to defined distance "dist"
      # use this to assess neighbors within "d" um of each cell
      ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct,], d) 
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificanceRef(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans,
                                      ct = ct)
      return(results)
    }) 
    names(results.all) <- levels(celltypes)
    
    ## Evaluate significance (cell type subsets)
  }
  
  if(!is.null(subset.list)){
    
    ## load in the subset file if it exists, or make it and probably be a good
    ## idea to save it, too
    if(!is.list(subset.list)){
      stop(paste0("`subset.list` is not a list. You can build this using: `binomialTestMatrix()` then `selectSubsets()`"))
    }
    
    combo_ids <- names(subset.list)
    d <- dist
    
    if(verbose){
      message("Calculating trends for each subset in `subset.list` with respect to the cell types in `cells$celltypes`")
    }
    
    ## initialize list
    results.all <- list()
    
    ## for each subset of cells, evaluate significance against each reference cell type across scales
    for(i in combo_ids){
      
      if(verbose){
        message(i)
      }
      
      ## get area around the subset cells to identify neighbors
      ref.buffer <- sf::st_buffer(cells[subset.list[[i]], ], d)
      ## union polygons to avoid memory overflow, too slow
      # ref.buffer_union <- sf::st_union(ref.buffer)
      # get the different types of neighbor cells that are within "d" of the ref cells
      neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
      
      ## remove duplicate neighbor cells to prevent them from being counted multiple times
      ## and inflate the Z scores
      if(removeDups){
        ## need to remove self cells else have trivial enrichment when d~0
        self.cells <- cells[cells$celltypes == ct,]
        ## remove only the original cells, as the duplicates will have 
        ## .something in the name
        neigh.cells <- neigh.cells[setdiff(rownames(neigh.cells), rownames(self.cells)),]
        
        # message("number of neighbor cells before: ", nrow(neigh.cells))
        #neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), rownames(cells)),]
        ## hack to accommodate self cells that are neighbors of another self cell
        neigh.cells <- neigh.cells[intersect(rownames(neigh.cells), c(rownames(cells), paste0(rownames(self.cells), '.1'))),]
        # message("number of neighbor cells after removing dups: ", nrow(neigh.cells))
      }
      
      ## evaluate significance https://online.stat.psu.edu/stat415/lesson/9/9.4
      ## chose to shuffle the scales in parallel, but in each scale, the perms done linearly
      ## I think a bottle neck originally was waiting for certain cell types to finish
      ## so if I split up the scales, might be able to get through each cell type faster and speed up entire process?
      ## I could also split up the permutations, but then each cell type for each scale is done one by one
      ## I could do cell types in parallel, but then for each cell type need to go through each res and each perm one by one
      results <- evaluateSignificanceRef(cells = cells,
                                      randomcellslist = shuffle.list,
                                      trueNeighCells = neigh.cells,
                                      cellBuffer = ref.buffer,
                                      ncores = ncores,
                                      removeDups = removeDups,
                                      returnMeans = returnMeans,
                                      ct = ct)
      results.all[[i]] <- results
      
      rm(ref.buffer)
      rm(neigh.cells)
      rm(results)
      gc(verbose = FALSE, reset = TRUE)
    }
    
  }
  
  ## return results
  if(verbose){
    total_t <- round(difftime(Sys.time(), start_time, units="mins"), 2)
    message(sprintf("Time was %s mins", total_t))
  }
  
  return(results.all)
  
}





## With dups ---------------------------------------------------

ncores <- 7

data(seq)
seq <- crawdad:::toSF(pos = seq[,c("x", "y")],
                      celltypes = seq$celltypes)

scales <- seq(100, 1000, by=50)
shuffle.list <- readRDS('running_code/processed_data/shufflelist_seq.RDS')
## No need to edit the function, just set removeDups to FALSE.
# ## find trends, dist 50
# results_50_dups_ref <- findTrendsRef(seq,
#                                   dist = 50,
#                                   shuffle.list = shuffle.list,
#                                   ncores = ncores,
#                                   verbose = TRUE,
#                                   returnMeans = FALSE,
#                                   removeDups = FALSE)
# ## Time was 12.03 mins
# dat_50_dups_ref <- crawdad::meltResultsList(results_50_dups_ref, withPerms = T)
# saveRDS(dat_50_dups_ref, 'running_code/processed_data/dat_seq_50_dups_ref.RDS')
dat_50_dups_ref <- readRDS('running_code/processed_data/dat_seq_50_dups_ref.RDS')


## vizColocDotplot
zsig <- correctZBonferroni(dat_50_dups_ref)
zsig <- 1.96
vizColocDotplot(dat_50_dups_ref, zsigThresh = zsig, zscoreLimit = 2*zsig, 
                reorder = F, mutual = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))


## vizSigTrends
## plot significant trends
all_cts <- unique(dat_50_dups_ref$reference)
all_cts <- c('Presomitic mesoderm')
sig_cts <- c('Dermomyotome', 'Spinal cord', 'Intermediate mesoderm')
for (ref_ct in all_cts) {
  ## spat
  # p1 <- vizClusters(seq, ofInterest = c(ref_ct), alpha = 1) + 
  #   theme_void() +
  #   theme(legend.position = "none")
  
  
  ## trends
  # sig_cts <- dat_50_dups_ref %>% 
  #   group_by(reference, neighbor, scale) %>% 
  #   summarize(Z = mean(Z)) %>% 
  #   filter(reference == ref_ct) %>% 
  #   filter(abs(Z) > zsig) %>% 
  #   pull(neighbor) %>% 
  #   unique() %>% 
  #   as.character()
  # not_sig_cts <- all_cts[! all_cts %in% sig_cts]
  
  filtered_dat <- dat_50_dups_ref %>% 
    group_by(reference, neighbor, scale) %>% 
    summarize(Z = mean(Z)) %>% 
    filter(reference == ref_ct) %>% 
    filter(neighbor %in% sig_cts) %>% 
    # filter(neighbor %in% all_cts) %>% 
    filter(neighbor != ref_ct) #%>% ## remove reference cell type
  # filter(scale <= 500) ## remove scale greater than 5000
  # mutate(selected_neighbor = fct_other(neighbor, keep = sig_cts)) %>% 
  p2 <- filtered_dat %>% 
    ggplot() +
    geom_line(aes(scale, Z, color = neighbor)) +
    scale_color_manual(values = rainbow(length(sig_cts))) +
    geom_hline(yintercept = zsig, linetype = 2) + 
    geom_hline(yintercept = -zsig, linetype = 2) +
    xlim(c(50, 1250)) + 
    ggrepel::geom_text_repel(aes(x=scale, y = Z, colour = neighbor, 
                                 label = neighbor), 
                             data = filtered_dat %>%
                               filter(scale == max(scale)),
                             hjust = 'right', box.padding = 0.5,
                             nudge_x = .5, nudge_y = .5,
                             segment.alpha = .5,
                             xlim = c(-Inf, 1250), ylim = c(-Inf, Inf)) + 
    labs(title = ref_ct)
  
  # p <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1/3, 2/3))
  
  print(p2)
  # png(paste0('running_code/exploring_embryo/trends_ref_', 
  #            gsub('[ /]', '_', ref_ct), '.png'),
  #     height = 700, width = 1400)
  # print(p)
  # dev.off()
}




# Null Sim ----------------------------------------------------------------

library(crawdad)
library(tidyverse)
ncores <- 7

cells <- readRDS('simulating_data/null_sim/cells_nullsim_s0.RDS')
vizClusters(cells)

## reduced to 500
scales <- seq(100, 500, by=50)

## With dups ---------------------------------------------------

# ## generate background
# shuffle.list <- crawdad:::makeShuffledCells(cells,
#                                             scales = scales,
#                                             perms = 10,
#                                             ncores = ncores,
#                                             seed = 1,
#                                             verbose = TRUE)
# saveRDS(shuffle.list, 'running_code/processed_data/shufflelist_nullsim.RDS')
# # No need to edit the function, just set removeDups to FALSE.
# # find trends, dist 50
# results_50_dups <- crawdad::findTrends(cells,
#                                   dist = 50,
#                                   shuffle.list = shuffle.list,
#                                   ncores = ncores,
#                                   verbose = TRUE,
#                                   returnMeans = FALSE,
#                                   removeDups = FALSE)
# ## Time was 17.62 mins
# dat_50_dups <- crawdad::meltResultsList(results_50_dups, withPerms = T)
# saveRDS(dat_50_dups, 'running_code/processed_data/dat_nullsim_50_dups.RDS')
dat_50_dups <- readRDS('running_code/processed_data/dat_nullsim_50_dups.RDS')

zsig <- correctZBonferroni(dat_50)
vizColocDotplot(dat_50_dups, reorder = F, 
                zsigThresh = zsig, zscoreLimit = zsig*2,
                dotSizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1))



## With dups and corrected -------------------------------------------------

shuffle.list <- readRDS('running_code/processed_data/shufflelist_nullsim.RDS')
## No need to edit the function, just set removeDups to FALSE.
## find trends, dist 50
results_50_dups_ref <- findTrendsRef(cells,
                                  dist = 50,
                                  shuffle.list = shuffle.list,
                                  ncores = ncores,
                                  verbose = TRUE,
                                  returnMeans = FALSE,
                                  removeDups = FALSE)
## Time was 12.03 mins
dat_50_dups_ref <- crawdad::meltResultsList(results_50_dups_ref, withPerms = T)
saveRDS(dat_50_dups_ref, 'running_code/processed_data/dat_nullsim_50_dups_ref.RDS')
dat_50_dups_ref <- readRDS('running_code/processed_data/dat_nullsim_50_dups_ref.RDS')


## vizColocDotplot
zsig <- correctZBonferroni(dat_50_dups_ref)
zsig <- 1.96
## will not reach significance
vizColocDotplot(dat_50_dups_ref, zscoreLimit = 2*zsig, 
                reorder = F, mutual = T, dotSizes = c(3,10))  +
  theme(legend.position='right',
        axis.text.x = element_text(angle = 45, h = 0))

vizColocDotplot(dat_50_dups_ref, reorder = F, 
                zscoreLimit = zsig*2,
                dotSizes = c(5, 20)) +
  theme(legend.position='bottom',
        axis.text.x = element_text(vjust = 1))
