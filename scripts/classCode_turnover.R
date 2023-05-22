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


