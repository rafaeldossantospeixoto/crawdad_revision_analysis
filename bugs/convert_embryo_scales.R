library(crawdad)
library(tidyverse)


# seqFISH ----------------------------------------------------------------

data('seq')
max(seq$x) - min(seq$x) # 5069.662
max(seq$y) - min(seq$y) # 6984.168
## it should be aprox 1000 microns, not 5000 or 6000
ggplot(seq, aes(x, -y)) + geom_point(size = .1)