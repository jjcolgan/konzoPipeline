library(Maaslin2)
library(tidyverse)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
metadata = metadata%>%
  column_to_rownames('Sample')
phylumLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_phylum-MATRIX.txt')
phylumLevel = phylumLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = phylumLevel,
  input_metadata = metadata,
  output = "plots/phylum/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration','Family,Date_of_Isolation'),
  reference = c ("Status,Unaffected"))

classLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_class-MATRIX.txt')
classLevel = classLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = classLevel,
  input_metadata = metadata,
  output = "plots/class/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration','Family,Date_of_Isolation'),
  reference = c ("Status,Unaffected"))

orderLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_order-MATRIX.txt')
orderLevel = orderLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = orderLevel,
  input_metadata = metadata,
  output = "plots/order/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration'),
  reference = c ("Status,Unaffected"))

familyLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_family-MATRIX.txt')
familyLevel = familyLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plots/family/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration,Date_of_Isolation','Family'),
  reference = c ("Status,Unaffected"))

genusLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
genusLevel = genusLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = genusLevel,
  input_metadata = metadata,
  output = "plots/genus/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration,Date_of_Isolation','Family'),
  reference = c ("Status,Unaffected"))

speciesLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')
speciesLevel = speciesLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 5,
  min_prevalence =0.05136986,
  input_data = speciesLevel,
  input_metadata = metadata,
  output = "plots/species/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration,Date_of_Isolation','Family'),
  reference = c ("Status,Unaffected"))
