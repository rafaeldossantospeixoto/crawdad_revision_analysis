## Stores all functions used to complete analysis that are not in the SEraster package.

calculateDensity <- function(matrix.array) {
  sum(matrix.array != 0)/(dim(matrix.array)[1] * dim(matrix.array)[2])
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}