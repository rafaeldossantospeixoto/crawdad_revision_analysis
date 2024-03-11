## simulate true positives and negatives
results <- do.call(rbind, lapply(1:10, function(s) {

print(s)

# Set parameters
set.seed(s)

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

rownames(locations) <- names(Xs1) <- names(Xs2) <- paste0('cell', 1:N)

# Split cells into two groups
set.seed(0)
vi <- sample(1:N, N/2)
vi2 <- setdiff(1:N, vi)

# Convert to categorical cell-types
celltypes <- rep('X', N)
names(celltypes) <- names(Xs2)
celltypes[vi[Xs1[vi] > 0]] <- 'A'
celltypes[vi[Xs1[vi] < 0]] <- 'B'
celltypes[vi2[Xs2[vi2] > 0]] <- 'C'
celltypes[vi2[Xs2[vi2] < 0]] <- 'D'
table(celltypes)

# Plot
par(mfrow=c(1,1), mar=rep(2,4))
MERINGUE::plotEmbedding(locations[names(Xs1),], col=Xs1, cex=1, alpha=0.5, s=1, v=1)
MERINGUE::plotEmbedding(locations[names(Xs2),], col=Xs2, cex=1, alpha=0.5, s=1, v=1)
MERINGUE::plotEmbedding(locations, groups=celltypes, cex=1, show.legend = TRUE, alpha=0.5, s=1, v=1)

library(crawdad)
library(tidyverse)
## convert dataframe to spatial points (SP)
## convert to more conventional micron ranges
scalefactor <- 1000
cells <- crawdad::toSF(pos = data.frame(locations)*scalefactor, celltypes = as.factor(celltypes))
plot(cells, pch=16, axes=TRUE)
plot(locations*1000)
grid <- sf::st_make_grid(cells, cellsize = 100, square=TRUE)
plot(grid, add=TRUE)

## define the scales to analyze the data
scales <- seq(100, 500, by=50)
scales
## shuffle cells to create null background
shuffle.list <- crawdad:::makeShuffledCells(cells,
                                            scales = scales,
                                            perm = 1,
                                            ncores = 10,
                                            seed = 1,
                                            verbose = TRUE)

## plot variogram
ct <- 'A'
Xs <- as.numeric(celltypes == ct)
variog <- geoR::variog(data = Xs, coords = data.frame(locations)*scalefactor, max.dist = 200, option = "bin", messages = FALSE)
plot(variog, type="l", main='variogram')

############# visualizing d
#d = smoothness_parameter*scalefactor/10*5
d = 300
ct = 'B'
ref.buffer <- sf::st_buffer(cells[cells$celltypes == ct, ], d)

## need to visualize buffer and maybe also quantify number of cells in each buffer?
foo <- celltypes == ct; names(foo) <- names(celltypes)
MERINGUE::plotEmbedding(data.frame(locations)*scalefactor, groups=foo, s=1, v=1)
plot(ref.buffer$geometry, add=TRUE)

## what if we choose an inappropriate dist? 
results <- findTrends(cells,
                               dist = d, 
                               shuffle.list = shuffle.list,
                               ncores = 10,
                               verbose = TRUE,
                               returnMeans = FALSE)
dat <- crawdad::meltResultsList(results, withPerms = TRUE)
#dat$Z[is.nan(dat$Z)] <- 0
## calculate the zscore for the multiple-test correction
zsig <- correctZBonferroni(dat)
## summary visualization
vizColocDotplot(dat, zsig.thresh = zsig, zscore.limit = 2*zsig) +
  theme(axis.text.x = element_text(angle = 35, h = 0))


################################ quantify performance
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

fpr = fp / (fp + tp) ## type 1 error rate
fpr

tnr = tn/ (tn + fp) ## specificity
tnr

tpr = tp/ (tp + fp) ## sensitivity
tpr

out <- data.frame(fpr=fpr, tnr=tnr, tpr=tpr, tp=tp, fp=fp, tn=tn)
print(out)

return(out)

}))

print(results)

## seems like we are very liable to false negatives but not false positives
## difficult to evaluate false negatives because we expect negatives at certain scales due to autocorrelation anyway

## is performance a function of d?

## d = 30
#   fpr tnr tpr tp fp tn
#1    0   1   1  9  0 11
#2    0   1   1 11  0 10
#3    0   1   1 11  0  9
#4    0   1   1 10  0  9
#5    0   1   1 12  0  9
#6    0   1   1  6  0 13
#7    0   1   1  3  0  1
#8    0   1   1  3  0  8
#9    0   1   1 10  0 20
#10   0   1   1 10  0 16

## d = 0.01 (super small, should only contain self?)
#   fpr tnr tpr tp fp tn
#1   NA  NA  NA NA NA NA
#2   NA  NA  NA NA NA NA
#3   NA  NA  NA NA NA NA
#4   NA  NA  NA NA NA NA
#5   NA  NA  NA NA NA NA
#6   NA  NA  NA NA NA NA
#7   NA  NA  NA NA NA NA
#8   NA  NA  NA NA NA NA
#9   NA  NA  NA NA NA NA
#10  NA  NA  NA NA NA NA

## d = 300 (way too big, so neighbors would be all cells, should also lead to false negatives)
#   fpr tnr tpr tp fp tn
#1  NaN NaN NaN  0  0  0
#2  NaN NaN NaN  0  0  0
#3  NaN NaN NaN  0  0  0
#4  NaN NaN NaN  0  0  0
#5  NaN NaN NaN  0  0  0
#6  NaN NaN NaN  0  0  0
#7  NaN NaN NaN  0  0  0
#8  NaN NaN NaN  0  0  0
#9  NaN NaN NaN  0  0  0
#10 NaN NaN NaN  0  0  0

############## behaves as expected






