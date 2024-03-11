## simulate true positives and negatives

## Trevor's recommendation, simulate two spatially autocorrelated but independent cell-types
## https://hastie.su.domains/Papers/biodiversity/Biodiversity.pdf

## RandomFields not available anymore
## write own code / copy from their source


#results <- sapply(1:10, function(i) {
#set.seed(i)
#print(i)

# Set parameters
set.seed(10)
N <- 2000  # Number of locations
nugget_variance <- 0.1  # Nugget variance (τ^2)
range_parameter <- 0.5  # Range parameter (κ)
smoothness_parameter <- 0.3  # Smoothness parameter (φ)

# Define locations in [0, 1] × [0, 1]
locations <- cbind(runif(N), runif(N)) 
colnames(locations) <- c('x', 'y')

# Function to calculate Matern covariance
matern_covariance <- function(h, range, smoothness) {
  nu <- smoothness
  term1 <- 2^(1 - nu)/gamma(nu)
  term2 <- (sqrt(2 * nu) * h / range)
  term3 <- besselK(term2, nu)
  return(term1 * (term2^nu) * term3)
}

# Function to calculate the covariance matrix
calculate_covariance_matrix <- function(locations, range, smoothness) {
  N <- nrow(locations)
  cov_matrix <- matrix(0, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in 1:N) {
      h <- sqrt(sum((locations[i,] - locations[j,])^2))
      cov_matrix[i, j] <- matern_covariance(h, range, smoothness)
    }
  }
  return(cov_matrix)
}

# Generate covariance matrix for the Gaussian random field
cov_matrix <- calculate_covariance_matrix(locations, range_parameter, smoothness_parameter)
diag(cov_matrix) <- 1

# Generate realizations of the Gaussian random field
W1 <- MASS::mvrnorm(1, rep(0, N), cov_matrix)
W2 <- MASS::mvrnorm(1, rep(0, N), cov_matrix)

# Generate independent and identically distributed errors
Z1 <- rnorm(N, mean = 0, sd = sqrt(nugget_variance))
Z2 <- rnorm(N, mean = 0, sd = sqrt(nugget_variance))

# Generate realizations of the process Xsi
Xs1 <- W1 + Z1
Xs2 <- W2 + Z2

# Print or plot the results as needed
print(Xs1)
print(Xs2) 

# Set names
rownames(locations) <- names(Xs1) <- names(Xs2) <- paste0('cell', 1:N)
MERINGUE::plotEmbedding(locations, col=Xs1, cex=1)
MERINGUE::plotEmbedding(locations, col=Xs2, cex=1)

# Confirm Rob and Trevor's original concerns with pixel data
plot(Xs1, Xs2)
cor(Xs1, Xs2) ## indeed not 0
cor.test(Xs1, Xs2) ## indeed significant

# versus just uniform
Ys1 <- rnorm(N)
Ys2 <- rnorm(N)
names(Ys1) <- names(Ys2) <- names(Xs1)
MERINGUE::plotEmbedding(locations, col=Ys1, cex=1)
MERINGUE::plotEmbedding(locations, col=Ys2, cex=1)
plot(Ys1, Ys2)
cor(Ys1, Ys2) ## indeed not 0
cor.test(Ys1, Ys2) ## indeed significant

# Split into two groups
set.seed(0)
vi <- sample(1:N, N/2)
vi2 <- setdiff(1:N, vi)
par(mfrow=c(1,2), mar=rep(2,4))
MERINGUE::plotEmbedding(locations[vi,], col=Xs1[vi], cex=1)
MERINGUE::plotEmbedding(locations[vi2,], col=Xs2[vi2], cex=1)

######################### For us, we have categorical variables 
celltypes <- rep('X', N)
names(celltypes) <- names(Xs2)
celltypes[vi[Xs1[vi] > 0]] <- 'A'
celltypes[vi[Xs1[vi] < 0]] <- 'B'
celltypes[vi2[Xs2[vi2] > 0]] <- 'C'
celltypes[vi2[Xs2[vi2] < 0]] <- 'D'
table(celltypes)

par(mfrow=c(2,2))
test <- celltypes %in% c('A')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('B')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('C')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('D')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)

# Plot
par(mfrow=c(1,1), mar=rep(2,4))
MERINGUE::plotEmbedding(locations, groups=celltypes, cex=1, show.legend = TRUE)

par(mfrow=c(3,2))
## A and B separated, C and D separated
test <- celltypes %in% c('A', 'B')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('C', 'D')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)

## but no spatial relationship with each other
test <- celltypes %in% c('A', 'C')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('B', 'C')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('A', 'D')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)
test <- celltypes %in% c('B', 'D')
MERINGUE::plotEmbedding(locations[test,], groups=celltypes, cex=1)

## all should be autocorrelated with self


######################### Run crawdad
library(crawdad)
library(tidyverse)
## convert dataframe to spatial points (SP)
## convert to more conventional micron ranges
cells <- crawdad::toSF(pos = data.frame(locations)*1000, celltypes = as.factor(celltypes))
plot(cells, pch=16, axes=TRUE)
1000/100 
grid(nx=10, col='black')

## define the scales to analyze the data
scales <- seq(50, 1000, by=100)
scales
## shuffle cells to create null background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perm = 10,
                                            ncores = 10,
                                            seed = 1,
                                            verbose = TRUE)
## calculate the zscore for the cell-type pairs at different scales
ref.buffer <- sf::st_buffer(cells[cells$celltypes == 'B', ], 30)
neigh.cells <- sf::st_intersection(cells, ref.buffer$geometry)
rand.cells <- cells; rand.cells$celltypes <- shuffle.list$`1000`$`1`
rand.neigh.cells <-  sf::st_intersection(rand.cells, ref.buffer$geometry)
plot(neigh.cells, pch=16, axes=TRUE)
plot(rand.neigh.cells, pch=16, axes=TRUE)
table(neigh.cells$celltypes)/table(rand.neigh.cells$celltypes)

## plot variogram
ct <- 'A'
Xs <- as.numeric(celltypes == ct)
variog <- geoR::variog(data = Xs, coords = data.frame(locations)*1000, max.dist = 100, option = "bin", messages = FALSE)
plot(variog, type="l", main='variogram')

## what if we choose an inappropriate dist? 
## -> if too large, looks like global so no significant results
## seems more liable to false negatives
## -> if too small, only see self by definition (Jean believes this to be a bug that leads to trivial error)
results <- crawdad::findTrends(cells,
                               dist = 30, 
                               shuffle.list = shuffle.list,
                               ncores = 10,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
range(dat$Z)
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))

# dat %>% 
#   filter(reference == 'A') %>% 
#   filter(neighbor == 'B') %>% 
#   vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
# dat %>% 
#   filter(reference == 'B') %>% 
#   filter(neighbor == 'A') %>% 
#   vizTrends(lines = TRUE, withPerms = TRUE, sig.thresh = zsig)
# 
# dat.diag <- dat %>% 
#   filter(reference == neighbor)
# table(dat.diag$Z > zsig)
# 
# foo <- dat %>% 
#   filter(reference == 'A') %>% 
#   filter(neighbor == 'B') 
# foo
# table(foo$Z > zsig)


### not sure how to deal with all scales frankly...should it be sum or any?

## true positives for separation? (how to distinguish sep from colocalized)
bar <- dat %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'B') 
ab <- sum(bar$Z < -zsig)
bar <- dat %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'A') 
ba <- sum(bar$Z < -zsig)
bar <- dat %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'D') 
cd <- sum(bar$Z < -zsig)
bar <- dat %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'C') 
dc <- sum(bar$Z < -zsig)
tn <- sum(c(ab, ba, cd, dc))
tn

## true positives for colocalization
bar <- dat %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'A') 
aa <- sum(bar$Z > zsig)
bar <- dat %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'B') 
bb <- sum(bar$Z > zsig)
bar <- dat %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'C') 
cc <- sum(bar$Z > zsig)
bar <- dat %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'D') 
dd <- sum(bar$Z > zsig)
tp <- sum(c(aa, bb, cc, dd))
tp

## false positives
bar <- dat %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'D') 
ad <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'A') 
da <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'D') 
bd <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'D') %>% 
  filter(neighbor == 'B') 
db <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'A') %>% 
  filter(neighbor == 'C') 
ac <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'A') 
ca <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'B') %>% 
  filter(neighbor == 'C') 
bc <- sum(abs(bar$Z) > zsig)
bar <- dat %>% 
  filter(reference == 'C') %>% 
  filter(neighbor == 'B') 
cb <- sum(abs(bar$Z) > zsig)

fp <- sum(c(ad, da, bd, db, ac, ca, bc, cb))
fp

## Jean: this is slightly wrong
fdr = fp / (fp + tp)
fdr

tnr = tn/ (tn + fp)
tnr

#})

## what is the distribution of p-values?
#hist(pnorm(dat$Z, mean = 0, sd = 1, lower.tail = TRUE))

############## Try ripley's?
library(spatstat.core)
## convert to more conventional micron ranges
pp <- as.ppp(locations*1000, c(range(locations[,1]*1000), range(locations[,2]*1000)))
pp <- pp %mark% as.factor(celltypes)
par(mfrow=c(1,1))
plot(pp)

sapply(levels(as.factor(celltypes)), function(i) {
  sapply(levels(as.factor(celltypes)), function(j) {
    test <- Lcross.inhom(pp, i, j)
    plot(test)
    
    ## observed - theoretical across distances
    plot(x = test$r, y = test$iso - test$theo, type="l", main=paste0(i, ':', j))
  })
})


data(pkhl)
plot(pkhl[,1:2], pch=".")

pp <- as.ppp(pkhl[,1:2])
pp <- pp %mark% as.factor(celltypes)
par(mfrow=c(1,1))
plot(pp)

library(spatstat.explore)
Gdot