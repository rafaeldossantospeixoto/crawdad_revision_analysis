library(crawdad)
library(tidyverse)

setwd("/Users/rafaeldossantospeixoto/Library/CloudStorage/OneDrive-JohnsHopkins/jefworks/crawdad/repos/CRAWDAD")

# convert pkhl coordinates ------------------------------------------------

data(pkhl)
paste(min(pkhl$x), max(pkhl$x), min(pkhl$y), max(pkhl$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/1393ff7729eff13cd2f4d5d58af96fdb
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

pkhl$x <- pkhl$x * x_res
pkhl$y <- pkhl$y * y_res
paste(min(pkhl$x), max(pkhl$x), min(pkhl$y), max(pkhl$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

save(pkhl, file = "data/pkhl.rda")



# convert xxcd ------------------------------------------------------------

xxcd <- read.csv2(file = paste0("data/spleen/XXCD.meta.csv.gz"), row.names = 1)
xxcd <- xxcd[,c("x", "y", "celltypes_folBcombined")]
xxcd <- xxcd %>%
  dplyr::mutate_at(vars(x, y), as.numeric)
colnames(xxcd) <- c('x', 'y', 'celltypes')
paste(min(xxcd$x), max(xxcd$x), min(xxcd$y), max(xxcd$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/d73a864f5ad6bdffdf148f43423e2a01
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

xxcd$x <- xxcd$x * x_res
xxcd$y <- xxcd$y * y_res
paste(min(xxcd$x), max(xxcd$x), min(xxcd$y), max(xxcd$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

write.csv2(xxcd, file = paste0("data/spleen/XXCD.meta.csv.gz"), row.names = TRUE)



# convert fsld ------------------------------------------------------------

fsld <- read.csv2(file = paste0("data/spleen/FSLD.meta.csv.gz"), row.names = 1)
fsld <- fsld[,c("x", "y", "celltypes_folBcombined")]
fsld <- fsld %>%
  dplyr::mutate_at(vars(x, y), as.numeric)
colnames(fsld) <- c('x', 'y', 'celltypes')
paste(min(fsld$x), max(fsld$x), min(fsld$y), max(fsld$y))
## [1] "0 12095 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/2b57ac6fea7aa96d26d2685a297a4e7a
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

fsld$x <- fsld$x * x_res
fsld$y <- fsld$y * y_res
paste(min(fsld$x), max(fsld$x), min(fsld$y), max(fsld$y))
## [1] "0 4564.699519789 -3423.4302888802 0"

write.csv2(fsld, file = paste0("data/spleen/FSLD.meta.csv.gz"), row.names = TRUE)



# convert pbvn ------------------------------------------------------------

pbvn <- read.csv2(file = paste0("data/spleen/PBVN.meta.csv.gz"), row.names = 1)
pbvn <- pbvn[,c("x", "y", "celltypes_folBcombined")]
pbvn <- pbvn %>%
  dplyr::mutate_at(vars(x, y), as.numeric)
colnames(pbvn) <- c('x', 'y', 'celltypes')
paste(min(pbvn$x), max(pbvn$x), min(pbvn$y), max(pbvn$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/3de2804ec3219228758dba8382d33423
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

pbvn$x <- pbvn$x * x_res
pbvn$y <- pbvn$y * y_res
paste(min(pbvn$x), max(pbvn$x), min(pbvn$y), max(pbvn$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

write.csv2(pbvn, file = paste0("data/spleen/PBVN.meta.csv.gz"), row.names = TRUE)



# convert ksfb ------------------------------------------------------------

ksfb <- read.csv2(file = paste0("data/spleen/KSFB.meta.csv.gz"), row.names = 1)
ksfb <- ksfb[,c("x", "y", "celltypes_folBcombined")]
ksfb <- ksfb %>%
  dplyr::mutate_at(vars(x, y), as.numeric)
colnames(ksfb) <- c('x', 'y', 'celltypes')
paste(min(ksfb$x), max(ksfb$x), min(ksfb$y), max(ksfb$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/00d1a3623dac388773bc7780fcb42797
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

ksfb$x <- ksfb$x * x_res
ksfb$y <- ksfb$y * y_res
paste(min(ksfb$x), max(ksfb$x), min(ksfb$y), max(ksfb$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

write.csv2(ksfb, file = paste0("data/spleen/KSFB.meta.csv.gz"), row.names = TRUE)



# Conver ngpl -------------------------------------------------------------

ngpl <- read.csv2(file = paste0("data/spleen/NGPL.meta.csv.gz"), row.names = 1)
ngpl <- ngpl[,c("x", "y", "celltypes_folBcombined")]
ngpl <- ngpl %>%
  dplyr::mutate_at(vars(x, y), as.numeric)
colnames(ngpl) <- c('x', 'y', 'celltypes')
paste(min(ngpl$x), max(ngpl$x), min(ngpl$y), max(ngpl$y))
## [1] "0 9407 -9071 0"

## https://portal.hubmapconsortium.org/browse/dataset/654418415bed5ecb9596b17a0320a2c6
x_res <- 377.4038462 / 1000 ## to micrometers
y_res <- 377.4038462 / 1000 ## to micrometers

ngpl$x <- ngpl$x * x_res
ngpl$y <- ngpl$y * y_res
paste(min(ngpl$x), max(ngpl$x), min(ngpl$y), max(ngpl$y))
## [1] "0 3550.2379812034 -3423.4302888802 0"

write.csv2(ngpl, file = paste0("data/spleen/NGPL.meta.csv.gz"), row.names = TRUE)
