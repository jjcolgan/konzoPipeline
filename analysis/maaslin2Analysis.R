library(Maaslin2)
library(tidyverse)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
metadata = metadata%>%
  column_to_rownames('Sample')
phylumLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_phylum-MATRIX.txt')
phylumLevel = phylumLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = phylumLevel,
  input_metadata = metadata,
  output = "plots/phylum/maaslin2Res",
  fixed_effects = c("Sex", "Status", 'Location'),
  reference =c(("Status,Unaffected"), ('Location,Masi-Manimba')))

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = phylumLevel,
  input_metadata = metadata,
  output = "plots/phylum/maaslin2ResLocation",
  fixed_effects = c("Sex",'Location'),
  reference =c('Location,Masi-Manimba'))

classLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_class-MATRIX.txt')
classLevel = classLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
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
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = orderLevel,
  input_metadata = metadata,
  output = "plots/order/maaslin2Res",
  fixed_effects = c("Status"),
  random_effects = c('Date_Of_Isolation'),
  reference = ("Status,Unaffected"))

familyLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_family-MATRIX.txt')
familyLevel = familyLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plots/family/maaslin2Res",
  fixed_effects = c("Status"),
  random_effects = c('Concentration','Date_of_Isolation', 'Family'),
  reference = c ("Status,Unaffected"))

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plots/family/maaslin2ResLocation",
  fixed_effects = c("Location"),
  random_effects = c('Concentration','Date_of_Isolation','Family'),
  reference = c ("Location,Masi-Manimba"))

genusLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
genusLevel = genusLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = genusLevel,
  input_metadata = metadata,
  output = "plots/genus/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Concentration','Date_of_Isolation','Family'),
  reference = c ("Status,Unaffected"))

speciesLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')
speciesLevel = speciesLevel%>%
  column_to_rownames('taxon')

phylumFit = Maaslin2(
  min_abundance = 10,
  min_prevalence =0.1,
  input_data = speciesLevel,
  input_metadata = metadata,
  output = "plots/species/maaslin2Res",
  fixed_effects = c("Sex", "Status"),
  random_effects = c('Age_T0','Concentration','Date_of_Isolation','Family'),
  reference = c ("Status,Unaffected"))
