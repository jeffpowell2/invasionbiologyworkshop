###### Monday morning ######

# download and load packages for reading and managing data
install.packages('tidyverse')
library(tidyverse)

# read data and get an overview
meta <-read_tsv('rawdata/Metadata_SOB_share.tsv')
meta
meta %>% 
  group_by(site_code) %>% 
  summarise(n())

#### mapping

# install and load sf, for working with geographical data
install.packages('sf')
library(sf)

# identify columns with geographical data
# name them for the coords argument in x, y order
# need to tell it what coordinate reference system to use (here WGS84=4326)
summary(meta)
geo <- meta %>% 
  st_as_sf(coords=c('long', 'lat'), crs=4326, remove=FALSE)
st_geometry(geo)
bbox(filter(geo, site_code=='INV'))


# get map layer data needed for plotting
install.packages('rnaturalearth')
install.packages('rnaturalearthdata')
install.packages('remotes')
remotes::install_github("ropensci/rnaturalearthhires")
library(rnaturalearth)
library(rnaturalearthdata)

# get map layer to plot the points onto (choose one of the following two lines)
au <- ne_countries(country='Australia', returnclass='sf') # use this for just country
au <- ne_states(country='Australia', returnclass='sf') # use this for the state boundaries too
st_geometry(au)

# plot points on the map
ggplot(au) + 
  geom_sf() + 
  geom_sf(data=geo, aes(colour=site_code))

# produce an interactive map to be able to zoom in
install.packages('tmap')
library(tmap)
tmap_mode('view')
tm_basemap(c('OpenStreetMap', 'OpenTopoMap')) + 
  tm_shape(geo) + 
  tm_dots(col='site_code', size=0.1)


