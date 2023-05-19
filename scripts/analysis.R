
# load libraries
library(tidyverse)
library(sf)
library(tmap)
library(raster)


# read data
meta <- read_tsv('deriveddata/Metadata_SOB_share.tsv', na=c('', 'NA', '#N/A'))
soil <- read_csv('deriveddata/hwsdvars.csv')
tax <- read_delim('rawdata/Taxa_table_SOB.csv', delim=';')
otu <- read_delim('rawdata/OTU_table_SOB.csv', delim=';')






## get climate data
library(raster)
summary(meta$long)
summary(meta$lat)
# https://www.worldclim.org/data/worldclim21.html
tmax <- getData('worldclim', var='tmax', res=2.5)
# tmax <- getData('worldclim', var='tmax', res=0.5, lon=mean(meta$long), lat=mean(meta$lat))
tmin <- getData('worldclim', var='tmin', res=2.5)
prec <- getData('worldclim', var='prec', res=2.5)

extract(prec, meta[1:5, c('long', 'lat')])
meta$tmax <- apply(extract(tmax, meta[, c('long', 'lat')]), 1, max)
meta$tmin <- apply(extract(tmin, meta[, c('long', 'lat')]), 1, min)
meta$prec <- apply(extract(prec, meta[, c('long', 'lat')]), 1, sum)
summary(meta)
meta$tmax <- meta$tmax/10


## get soil data
devtools::install_github("stineb/rhwsd")
library(rhwsd)
# https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v12/en/
curl::curl_download('https://www.fao.org/fileadmin/user_upload/soils/HWSD%20Viewer/HWSD_RASTER.zip', 
                    destfile = 'hwsd/HWSD_RASTER.zip')
unzip(zipfile = 'hwsd/HWSD_RASTER.zip', exdir = 'hwsd/HWSD_RASTER')
con <- get_hwsd_con()
# ans <- get_hwsd_siteset(x=meta %>% rename(lon=long), con=con, hwsd.bil="hwsd/HWSD_RASTER/hwsd.bil")
ans <- get_hwsd_siteset(x=meta[1:3, ] %>% rename(lon=long), con=con, hwsd.bil="hwsd/HWSD_RASTER/hwsd.bil")
ans
names(unnest(ans, cols=data))
unnest(ans, cols=data) %>%
  ungroup() %>%
  # # added the next line to account for multiple lines per site being added
  filter(!duplicated(SampleID)) %>%
  dplyr::select(SampleID, T_CLAY, T_SILT, T_SAND, T_PH_H2O, T_OC, T_BULK_DENSITY) -> soil
# write_csv(soil, 'deriveddata/hwsdvars.csv')
meta <- left_join(meta, soil)


summary(meta)

library(vegan)
otu %>% 
  column_to_rownames('OTU_ID') %>% 
  t() -> mat
summary(rowSums(mat))
sort(rowSums(mat))[1:20]
rarecurve(mat, step=500)
matr <- rrarefy(mat, 5000)
summary(rowSums(matr))
matr <- matr[rowSums(matr) >= 5000, ]
summary(rowSums(matr))
divers <- data.frame(SampleID=rownames(mat), rich=specnumber(mat), 
                     shan=diversity(mat, index='shannon'), 
                     chao1=apply(mat, 1, fossil::chao1))
divers_r <- data.frame(SampleID=rownames(matr), rich_r=specnumber(matr), 
                     shan_r=diversity(matr, index='shannon'), 
                     chao1_r=apply(matr, 1, fossil::chao1))
divers <- full_join(divers, divers_r)

dat <- left_join(meta, divers) %>% 
  mutate(site_code=fct_relevel(site_code, 'UN', 'PL', 'INV'))


library(car)
library(lme4)
m1 <- lm(chao1_r ~ site_code + tmax + prec + T_CLAY, data=dat)
m1 <- lm(chao1_r ~ site_code*tmax + site_code*prec + site_code*T_CLAY, data=dat)
Anova(m1)
summary(m1)
residualPlot(m1)
ggpredict(m1, c('tmax', 'site_code')) %>% plot()

m1 <- lmer(chao1_r ~ site_code*tmax + site_code*prec + site_code*T_CLAY + (1|site_name), data=dat)
Anova(m1)
summary(m1)
plot(m1)


dat %>% 
  group_by(site_name) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

library(spdep)
library(spatialreg)
nb <- knn2nb(knearneigh(as.matrix(dat[, c('long', 'lat')]), k=15, longlat=TRUE))
lw <- nb2listw(nb, zero.policy=NULL)
m1s <- lagsarlm(m1, data=dat, listw=lw)
# m1s <- errorsarlm(m1, data=dat, listw=lw)
summary(m1s)

coef(m1); coef(m1s)

library(ggeffects)
m1.pred <- ggpredict(m1, 'site_code')
m1.pred
names(m1.pred)
ggplot(m1.pred, aes(x=x, y=predicted)) + 
  geom_linerange(aes(ymin=conf.low, ymax=conf.high), linewidth=1) + 
  geom_point(size=2) + 
  geom_jitter(data=dat, 
              mapping=aes(x=site_code, y=chao1_r), 
              size=0.5, alpha=0.2, width=0.1)

m1.pred <- ggpredict(m1, c('tmax', 'site_code'))
ggplot(m1.pred, aes(x=x, y=predicted, colour=group, fill=group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high, colour=NULL), alpha=0.2) +
  geom_line(linewidth=1) +
  geom_point(data=dat,
             mapping=aes(x=tmax, y=chao1_r, colour=site_code),
             size=0.5, alpha=0.2, inherit.aes=FALSE) +
  theme_classic()


adonis
betadisper
indicators

funguild and ectos

