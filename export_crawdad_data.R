library(crawdad)
library(tidyverse)

data(sim)
write.csv(sim, 'running_code/squidpy/exported_data/sim.csv')

data(slide)
write.csv(slide, 'running_code/squidpy/exported_data/slide.csv')

data(seq)
write.csv(seq, 'running_code/squidpy/exported_data/seq.csv')