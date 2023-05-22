# download and load packages for reading and managing data
# install.packages('tidyverse')
library(tidyverse)

# read data and get an overview
meta <-read_tsv('rawdata/Metadata_SOB_share.tsv')
meta
meta %>% 
  group_by(site_code) %>% 
  summarise(n())

# finding and incorporating additional data 
library(raster)

# can get precipition and temperature by querying worldclim within R
# see https://www.worldclim.org/data/worldclim21.html for more info
prec <- getData('worldclim', var='prec', res=2.5)
tavg <- getData('worldclim', var='tmean', res=2.5)

# example of loading and merging multiple geotiffs for solar irradiation
files <- list.files('wc2.1_2.5m_srad', full.names=TRUE)
srad <- brick(files[1])
for(i in 2:length(files)) {
  temp <- brick(files[i])
  srad <- merge(srad, temp)
}

## other databases
# aridity and PET: https://cgiarcsi.community/data/global-aridity-and-pet-database/
# GPP and NPP: https://modis.gsfc.nasa.gov/data/dataprod/mod17.php
# harmonised world soil database: https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v12/en/

# extract values at sampled coordinates and add to metadata table
meta$prec <- apply(extract(prec, meta[, c('long', 'lat')]), 1, sum)
meta$tmean <- apply(extract(tavg, meta[, c('long', 'lat')]), 1, mean)
meta$tmean <- meta$tmean/10  # because temperature stored in way to avoid decimal
meta$srad <- apply(extract(srad, meta[, c('long', 'lat')]), 1, max)
summary(meta)

# read in species-sample table
otu <- read_delim('rawdata/OTU_table_SOB.csv', delim=';')
# rows are species, columns are samples, need to transpose
otu %>% 
  column_to_rownames('OTU_ID') %>% 
  t() -> mat

# are there many zero values? check this for each species
summary(apply(mat, 2, function(x)sum(x==0)))

# is there variable sampling intensity? check this for each sample
summary(rowSums(mat))

## calculate some alpha diversity metrics

# most of this uses the vegan library
library(vegan)

# create dataframe with three different diversity metrics: 
divers <- data.frame(SampleID=rownames(mat), 
                     rich=specnumber(mat), # richness
                     shan=diversity(mat, index='shannon'), # shannon diversity
                     chao1=apply(mat, 1, fossil::chao1)) # rarefied richness using the Chao1 estimator
summary(divers)


