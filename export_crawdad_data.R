library(crawdad)
library(tidyverse)

data(sim)
write.csv(sim, 'running_code/squidpy/exported_data/sim.csv')

data(slide)
write.csv(slide, 'running_code/squidpy/exported_data/slide.csv')

data(seq)
write.csv(seq, 'running_code/squidpy/exported_data/seq.csv')

data(pkhl)
write.csv(pkhl, 'running_code/squidpy/exported_data/spleen/pkhl.csv')

spleen <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
head(spleen)
range(spleen$x)
write.csv(spleen, 'running_code/squidpy/exported_data/spleen/xxcd.csv')

spleen <- read.csv2(file = '../CRAWDAD/data/spleen/FSLD.meta.csv.gz', row.names = 1)
head(spleen)
range(spleen$x)
write.csv(spleen, 'running_code/squidpy/exported_data/spleen/fsld.csv')

spleen <- read.csv2(file = '../CRAWDAD/data/spleen/KSFB.meta.csv.gz', row.names = 1)
head(spleen)
range(spleen$x)
write.csv(spleen, 'running_code/squidpy/exported_data/spleen/ksfb.csv')

spleen <- read.csv2(file = '../CRAWDAD/data/spleen/NGPL.meta.csv.gz', row.names = 1)
head(spleen)
range(spleen$x)
write.csv(spleen, 'running_code/squidpy/exported_data/spleen/ngpl.csv')

spleen <- read.csv2(file = '../CRAWDAD/data/spleen/PBVN.meta.csv.gz', row.names = 1)
head(spleen)
range(spleen$x)
write.csv(spleen, 'running_code/squidpy/exported_data/spleen/pbvn.csv')