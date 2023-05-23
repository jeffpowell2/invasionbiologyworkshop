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


# load vegan library for betadiversity analyses
library(vegan)

# assess variation in sampling effort, plotting sample effort curves
temp <- rarecurve(mat, step=1000, tidy=TRUE)
left_join(meta, temp %>% rename(SampleID=Site)) %>% 
  ggplot(aes(x=Sample, y=Species, colour=site_code, group=SampleID)) + 
  geom_line() + 
  facet_wrap(~site_code)

# check variation in sample effort, looking for break points
hist(log10(rowSums(mat)))
sort(rowSums(mat))[1:25]

# rarefy the community matrix, using an arbitrary cut-off
temp<-rrarefy(mat, 5000)
matr <- temp[rowSums(temp)==5000, ]
dim(matr)

# prepare our data for betadiversity analyses:
# 1 - focus analysis on ectomycorrhizal fungal taxa (create vector of OTU ids)
# 2 - join the community table with our metadata table (need to also modify the community table so that it can be joined)
ecm_otus <- filter(tax, guild=='ectomycorrhizal')$OTU_ID
left_join(meta, mat %>% 
            as.data.frame() %>% 
            dplyr::select(ecm_otus) %>% # just ecto OTUs
            rownames_to_column('SampleID')) -> dat_ecm

# first analysis - indicator species analysis 
# identify OTUs that are overrepresented in samples coming from one or more groups of site_code
library(indicspecies)
multipatt(dat_ecm %>% select(starts_with('OTU_')), # first argument is the community table, select only those columns
          dat_ecm$site_code) -> res
summary(res)

# to visualise differences in taxonomic composition, using output from multipatt()
# first prepare the data in the res object so that it can be joined with taxonomic info
# output is only those otus significant for one or more site_codes
out <- res[['sign']] %>% 
  filter(p.value <= 0.05) %>% 
  rownames_to_column('OTU_ID') %>% 
  pivot_longer(cols=starts_with('s.'), names_to='group', values_to='value') %>% 
  filter(value==1) %>% 
  mutate(site_code = gsub('^s.', '', group))
# then join with the taxonomy table, then the relevant community data, 
# and reorder the otu levels by decreasing abundance
out <- left_join(out, tax) %>% 
  left_join(dat_ecm %>% 
              select(SampleID, site_code, starts_with('OTU_')) %>% 
              pivot_longer(cols=starts_with('OTU_'), names_to='OTU_ID', 
                           values_to='count')) %>% 
  mutate(OTU_ID = fct_reorder(OTU_ID, count, max), 
         site_code = as_factor(site_code))
# finally produce the barplot
out %>% 
  ggplot(aes(x=site_code, y=count, fill=Order, text=OTU_ID)) + # text aesthetic is for the ggplotly visualisation below
  geom_bar(stat='identity', position=position_fill()) + 
  scale_x_discrete(drop=FALSE) + 
  scale_fill_brewer(palette='Set3') + 
  scale_y_continuous(labels = scales::percent) + 
  labs(y='Percentage') + 
  theme_bw() -> p1
p1
# use this next lines to interactively get OTU IDs
plotly::ggplotly(p1)


# next analysis - permanova
# extract the community table, save as a new object
mat_ecm <- dat_ecm %>% select(starts_with('OTU_'))

# run three permanovas, each with a different distance index / raw data input
adonis2(mat_ecm ~ site_name, data=dat_ecm, distance='bray', add=TRUE)
adonis2(mat_ecm ~ site_name, data=dat_ecm, distance='bray', binary=TRUE, add=TRUE)
adonis2(mat_ecm ~ site_name, data=dat_ecm, distance='robust.aitchison', add=TRUE)

# for more distance/similarity indices, look at:
?dist
?vegdist
?labdsv::dsvdis

# for ways to standardise/transform the data prior to analysis, look at:
?decostand

# a redundancy analysis (constrained PCA) - ordination result is poor, don't interpret this
rda1 <- rda(mat_ecm ~ site_code, data=dat_ecm)
plot(rda1)

# a constrained analysis of principal coordinates using bray-curtis distances - ordination result is better, still some skew on the top of the y-axis
cap1 <- capscale(mat_ecm ~ site_code, data=dat_ecm, distance='bray', add=TRUE)
plot(cap1)

# try the same but with presence-absence instead of raw counts - still some skew but better
cap1 <- capscale(decostand(mat_ecm, method='pa') ~ 
                   site_code, data=dat_ecm, distance='bray', add=TRUE)
plot(cap1)

# try a different distance index - result is quite good
cap1 <- capscale(mat_ecm ~ site_code, data=dat_ecm, distance='robust.aitchison', add=TRUE)
plot(cap1)
cap1 # summary of inertia
cap1$CCA$eig/cap1$tot.chi # proportion of variation associated with each axis
anova(cap1) # statistical significance of the constrait

# produce a nice plot
# first extract scores from the resulting object and subset out different types of scores
scrs <- scores(cap1, tidy=TRUE)
scrs_spp <- scrs %>% filter(score=='species')
scrs_site <- scrs %>% filter(score=='sites')
scrs_cent <- scrs %>% filter(score=='centroids')

# first plot - site scores along with centroids for each group
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=CAP1, y=CAP2, colour=site_code)) + 
  geom_point(size=0.5, alpha=0.5) + 
  geom_point(data=scrs_cent %>% 
               rename(site_code=label) %>% 
               mutate(site_code=gsub('site_code', '', site_code)), size=2) + 
  theme_bw() + 
  theme(legend.position='top') -> p1
# second plot - species scores for those loaded heavily along at least one axis
# this one still needs some work
ggplot(scrs_spp %>% 
         filter(abs(CAP1) > 0.5 | abs(CAP2) > 0.5), aes(x=CAP1, y=CAP2, label=label)) + 
  geom_text() + 
  xlim(c(min(scrs_site[, 'CAP1']), max(scrs_site[, 'CAP1']))) + 
  ylim(c(min(scrs_site[, 'CAP2']), max(scrs_site[, 'CAP2']))) + 
  theme_bw() -> p2
# plot side-by-side using the patchwork package
library(patchwork)
p1 + p2


# top ten otus associated with the top of the ordination, presumably 'INV' samples
scrs_spp %>% 
  arrange(desc(CAP2)) %>% 
  head(10) -> inv_otus
# taxonomic information for those otus
tax %>% 
  filter(OTU_ID %in% inv_otus$label)







# still tidying - plotting turnover in space
pco1 <- capscale(mat_ecm1 ~ 1, data=dat_ecm, distance='bray', add=TRUE)
scrs_site <- scores(pco1, display='sites') # TIDY
cbind(dat_ecm, scrs_site) %>% 
  ggplot(aes(x=long, y=lat, colour=MDS1)) + 
  geom_point()
