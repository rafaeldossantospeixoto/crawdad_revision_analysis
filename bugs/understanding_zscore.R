
calculateZScore <- function(y1, y2, n1, n2) {
  p1 <- y1/n1
  p2 <- y2/n2
  p <- (y1+y2)/(n1+n2)
  Z <- (p1-p2)/sqrt(p*(1-p)*(1/n1+1/n2))
  return(Z)
}



# Individual cells --------------------------------------------------------

## star cells with intersection
calculateZScore(5,6,11,11)
# [1] -0.4264014

## star cells individually
calculateZScore(3,3,6,6)
# [1] 0
calculateZScore(2,3,5,5)
# [1] -0.6324555


## square cells with intersection
calculateZScore(4,1,11,11)
# [1] 1.526241

## square cells individually
calculateZScore(2,1,6,6)
# [1] 0.6666667
calculateZScore(2,0,5,5)
# [1] 1.581139

## square cells individually, changing remaning cell
calculateZScore(2,0,6,6)
# [1] 1.549193
calculateZScore(2,1,5,5)
# [1] 0.6900656




# Normalizing by ref ------------------------------------------------------

calculateZScore(5, 6, 11, 11)
# [1] -0.4264014

calculateZScore(5/2, 6/2, 11/2, 11/2)
# [1] -0.3015113