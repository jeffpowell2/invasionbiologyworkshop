
# load libraries
library(tidyverse)
# library(phyloseq)
library(sf)
library(tmap)
library(terra)
library(bioclim)


# read data
meta <- read_tsv('rawdata/Metadata_SOB.tsv', na=c('', 'NA', '#N/A'))
soil <- read_tsv('rawdata/Soil_CNP.txt')
tax <- read_delim('rawdata/Taxa_table_SOB.csv', delim=';')
otu <- read_delim('rawdata/OTU_table_SOB.csv', delim=';')


# check data
table(tax$Phylum)
any(is.na(tax$Phylum))
summary(soil)

meta %>% 
  group_by(site_code, site_host) %>% 
  summarise(n = n())
# PL = plantation
# UN = native forest
# INV = native forest with invasive pines


### maps
geo <- meta  %>% 
  filter(!is.na(lat), site_code %in% c('INV', 'PL', 'UN')) %>% 
  st_as_sf(coords=c('long', 'lat'), crs=4326, remove=FALSE) %>%  # EPSG:4326 = WGS84
  st_transform(crs=3577)  # EPSG:3577 = Australian Albers


## static map using `ggplot2`
library(rnaturalearth)
library(rnaturalearthdata)
AU <- ne_countries(country='australia', returnclass='sf')
AU <- ne_states(country='australia', returnclass='sf')
st_geometry(AU) 

# change the CRS to Australian Albers system to match our data 
AU <- st_transform(AU, 3577)  # EPSG:3577 = Australian Albers
ggplot(AU) + 
  geom_sf() + 
  geom_sf(data=geo, aes(colour=site_code))


## interactive map using `tmap`
tmap_mode('view')
tm_basemap(c('OpenStreetMap', 'OpenTopoMap')) +
  tm_shape(geo) + 
  tm_dots(col='site_code', size=0.1, alpha=1, title='site_code', 
          palette=c(INV='red', UN='grey', PL='blue'), popup.vars='SampleID')



### write output
write_tsv(meta %>% 
            select(SampleID, site_name, site_host, site_code, state, lat, long, elev, date), 
          'deriveddata/Metadata_SOB_share.tsv')
write_tsv(soil %>% 
            select(SampleID, percN, percC, P), 
          'deriveddata/Soil_CNP_share.txt')
