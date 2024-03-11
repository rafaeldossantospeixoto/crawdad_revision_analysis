library(crawdad)
library(tidyverse)


# Seed 1:10 ---------------------------------------------------------------

## simulate true positives and negatives
dfs <- lapply(1:10, function(s) {
  
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
  
  ## convert dataframe to spatial points (SF)
  ## convert to more conventional micron ranges
  scalefactor <- 1000
  cells <- crawdad::toSF(pos = data.frame(locations)*scalefactor, celltypes = as.factor(celltypes))
  plot(cells, pch=16, axes=TRUE)
  plot(locations*1000)
  grid <- sf::st_make_grid(cells, cellsize = 100, square=TRUE)
  plot(grid, add=TRUE)
  
  return(cells)
})

saveRDS(dfs, 'simulating_data/null_sim/cells_nullsim.RDS')


# Seed 0 ------------------------------------------------------------------

s <- 0

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

## convert dataframe to spatial points (SF)
## convert to more conventional micron ranges
scalefactor <- 1000
cells <- crawdad::toSF(pos = data.frame(locations)*scalefactor, celltypes = as.factor(celltypes))


## Save data as csv --------------------------------------------------------

vizClusters(cells, pointSize = 2)

cells_df <- data.frame(locations)*scalefactor
cells_df$celltypes = as.factor(celltypes)
head(cells_df)

write.csv(cells_df, file = 'simulating_data/null_sim/df_nullsim_s0.csv')
saveRDS(cells, 'simulating_data/null_sim/cells_nullsim_s0.RDS')
