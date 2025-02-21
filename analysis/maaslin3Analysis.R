library(maaslin3)
library(tidyverse)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
metadata = metadata%>%
  column_to_rownames('Sample')%>%
  as.data.frame()
phylumLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_phylum-MATRIX.txt')
phylumLevel = phylumLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()


maaslin3(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = phylumLevel,
  input_metadata = metadata,
  output = "plots/phylum/maaslin3Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration','Family,Date_of_Isolation'),
  reference = c("Status,Unaffected"))