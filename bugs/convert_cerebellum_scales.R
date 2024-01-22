library(crawdad)
library(tidyverse)


# slideSeq ----------------------------------------------------------------

data('slide')
max(slide$x) - min(slide$x) # 4628.55
max(slide$y) - min(slide$y) # 3800.465
## it should be aprox 2904 and 2425 microns, not 4628.55 or 3800.465
ggplot(slide, aes(x, -y)) + geom_point(size = .1)