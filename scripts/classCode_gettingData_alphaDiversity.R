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
meta$srad <- extract(srad, meta[, c('long', 'lat')])
# below doesn't work, need to fix reading multiple geotiffs for srad
# meta$srad <- apply(extract(srad, meta[, c('long', 'lat')]), 1, max)
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
dim(divers)
dim(meta)

# join our table with diversity metrics to the metadata table, for further analyses
dat <- left_join(meta, divers) %>% 
  mutate(site_code = fct_relevel(site_code, 'UN')) # set native forest as intercept

# look at the relationship between precipitation and richness, and see if this differs among site types
ggplot(dat, aes(y=rich, x=prec, colour=site_code)) + 
  geom_point() + # scatter plot
  geom_smooth(method='lm', formula='y~poly(x, 2)') + # show relationship, linear and quadratic
  facet_wrap(~site_code, nrow=1) # separate panel for each site_code

# do same for average temperature, but only linear relationship
ggplot(dat, aes(y=rich, x=tmean, colour=state)) + 
  geom_point() + 
  geom_smooth(method='lm') + 
  facet_wrap(~site_code, nrow=1)

# fit a linear model for temperature
m1 <- lm(log10(rich) ~ tmean * site_code, data=dat)

# check model diagnostic plots
library(car)
qqPlot(m1)
residualPlot(m1)

# Anova table, what terms are significant
Anova(m1)
# evaluate model estimates for each level of site_code
summary(m1)
# plot model predictions visually
library(ggeffects)
ggpredict(m1, c('tmean', 'site_code')) %>% plot()

# option #1 for accounting for nonindependence of sample sites
# mixed effects model
library(lme4)
m1e <- lmer(logrich ~ tmean * site_code + (1|site_name), data=dat)
# check diagnostic plot
plot(m1e)
# produce ANOVA table
Anova(m1e, test='F')
# check coefficients, including random effects
summary(m1e)
# visualise model predictions
ggpredict(m1e, c('tmean', 'site_code'), back.transform = FALSE) %>% plot()

# option #2 for accounting for nonindependence of samples sites
# spatial regression
library(spatialreg)
library(spdep)
# up to 15 samples collected per site
dat %>% 
  group_by(site_name) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))
# create sample neighbourhoods for fitting spatial regression model
nb <- knn2nb(knearneigh(as.matrix(dat[, c('long', 'lat')]), k=15, longlat=TRUE))
lw <- nb2listw(nb)
# fit model and inspect parameter estimates 
m1s <- lagsarlm(m1, data=dat, listw=lw)
summary(m1s)


