library(crawdad)
library(dplyr)
library(sf)

# Spleen ------------------------------------------------------------------

## pkhl ------------------------------------------------------------------

data('pkhl')
head(pkhl)
range(pkhl$x)
range(pkhl$y)
write.csv(pkhl, file = 'data_publication/pkhl.csv')

## xxcd ------------------------------------------------------------------

xxcd <- read.csv2(file = '../CRAWDAD/data/spleen/XXCD.meta.csv.gz', row.names = 1)
head(xxcd)
range(xxcd$x)
range(xxcd$y)
write.csv(xxcd, file = 'data_publication/xxcd.csv')

## fsld --------------------------------------------------------------------

fsld <- read.csv2(file = '../CRAWDAD/data/spleen/FSLD.meta.csv.gz', row.names = 1)
head(fsld)
range(fsld$x)
range(fsld$y)
write.csv(fsld, file = 'data_publication/fsld.csv')

## pbvn --------------------------------------------------------------------

pbvn <- read.csv2(file = '../CRAWDAD/data/spleen/PBVN.meta.csv.gz', row.names = 1)
head(pbvn)
range(pbvn$x)
range(pbvn$y)
write.csv(pbvn, file = 'data_publication/pbvn.csv')

## ksfb --------------------------------------------------------------------

ksfb <- read.csv2(file = '../CRAWDAD/data/spleen/KSFB.meta.csv.gz', row.names = 1)
head(ksfb)
range(ksfb$x)
range(ksfb$y)
write.csv(ksfb, file = 'data_publication/ksfb.csv')

## ngpl --------------------------------------------------------------------

ngpl <- read.csv2(file = '../CRAWDAD/data/spleen/NGPL.meta.csv.gz', row.names = 1)
head(ngpl)
range(ngpl$x)
range(ngpl$y)
write.csv(ngpl, file = 'data_publication/ngpl.csv')




# Sim ---------------------------------------------------------------------

data(sim)
head(sim)
range(sim$x)
range(sim$y)
write.csv(sim, file = 'data_publication/sim.csv')





# Ext sim -----------------------------------------------------------------

data(sim)
simbr <- sim %>% 
  mutate(celltypes = 'E',
         x = x + 2000,
         y = y) 
simtr <- sim %>% 
  mutate(celltypes = 'E',
         x = x + 2000,
         y = y + 2000) 
simtl <- sim %>% 
  mutate(celltypes = 'E',
         x = x,
         y = y + 2000) 
ext_sim <- bind_rows(sim, simbr, simtr, simtl)
head(ext_sim)
range(ext_sim$x)
range(ext_sim$y)
write.csv(ext_sim, file = 'data_publication/ext_sim.csv')



# Null sim ----------------------------------------------------------------

dfs <- readRDS('simulating_data/null_sim/cells_nullsim.RDS')
head(dfs)

for (i in 1:length(dfs)){
  dfi <- dfs[[i]]
  coords <- st_coordinates(dfi)
  dfi$x <- coords[, 1]
  dfi$y <- coords[, 2]
  dfi <- st_drop_geometry(dfi)
  head(dfi)
  write.csv(dfi, file = paste0('data_publication/null_sim_',i,'.csv'))
}

## df 
df1 <- readRDS("simulating_data/null_sim/seed_1/df1.RDS")
df2 <- readRDS("simulating_data/null_sim/seed_1/df2.RDS")
df_viz <- rbind(df1, df2)
head(df_viz)
write.csv(df_viz, file = 'data_publication/null_sim_visualization.csv')


