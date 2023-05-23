# download and load packages for reading and managing data
# install.packages('tidyverse')
library(tidyverse)

# read data files
# metadata (note where you saved the files, you might need to change the path)
meta <- read_tsv('deriveddata/meta_with_worldclim.tsv')
soil <- read_csv('rawdata/hwsdvars.csv')
# community table
otu <- read_delim('rawdata/OTU_table_SOB.csv', delim=';')
# funguild output
ecm <- read_tsv('rawdata/funguild_ecm.tsv')
path <- read_tsv('rawdata/funguild_path.tsv')
sap <- read_tsv('rawdata/funguild_sap.tsv')
# taxonomy table
tax <- read_delim('rawdata/Taxa_table_SOB.csv', delim=';') 

# add available guild information to taxonomy table
tax %>% 
  mutate(across(Kingdom:Species, ~gsub('^[a-z]\\_\\_', '', .)), 
         guild = case_when(Genus %in% ecm$genus ~ 'ectomycorrhizal', 
                           Genus %in% path$genus ~ 'pathogenic', 
                           Genus %in% sap$genus ~ 'saprotrophic')) -> tax

# transpose community table, to sample-taxon, for vegan functions
otu %>% 
  column_to_rownames('OTU_ID') %>% 
  t() -> mat


library(vegan)
temp <- rarecurve(mat, step=1000, tidy=TRUE)
left_join(meta, temp %>% rename(SampleID=Site)) %>% 
  ggplot(aes(x=Sample, y=Species, colour=site_code, group=SampleID)) + 
  geom_line() + 
  facet_wrap(~site_code)
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:25]
temp<-rrarefy(mat, 5000)
matr <- temp[rowSums(temp)==5000, ]
dim(matr)

ecm_otus <- filter(tax, guild=='ectomycorrhizal')$OTU_ID
left_join(meta, mat %>% 
            as.data.frame() %>% 
            dplyr::select(ecm_otus) %>% 
            rownames_to_column('SampleID')) -> dat_ecm

library(indicspecies)
multipatt(dat_ecm %>% select(starts_with('OTU_')), 
          dat_ecm$site_code) -> res
summary(res)

mat_ecm <- dat_ecm %>% select(starts_with('OTU_'))
adonis2(mat_ecm ~ site_code, data=dat_ecm, add=TRUE)

rda1 <- rda(mat_ecm ~ site_code, data=dat_ecm)
rda1
plot(rda1)
anova(rda1)

mat_ecm1 <- decostand(mat_ecm, method='pa')
cap1 <- capscale(mat_ecm1 ~ site_code, data=dat_ecm, distance='bray', add=TRUE)
cap1
anova(cap1)
plot(cap1)
ordiplot(cap1)
summary(cap1)

scrs_spp <- scores(cap1, display='species')
scrs_site <- scores(cap1, display='sites') # TIDY
head(scrs_site)

cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2, colour=site_code)) + 
  geom_point() -> p
plotly::ggplotly(p)

head(scores(cap1, display='lc'))

head(scrs_spp[order(scrs_spp[, 'CAP2'], decreasing=TRUE), ]) -> inv_otus
tax %>% 
  filter(OTU_ID %in% rownames(inv_otus))
