library(maaslin3)
library(tidyverse)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
readCounts = read_tsv('read_counts.tsv')

familiesPass = metadata %>%
  group_by(Family)%>%
  summarise('nFamily'= n())%>%
  filter(nFamily > 1)%>%
  .$Family

metadata = metadata %>%
  filter(Family %in% familiesPass)

metadata = metadata%>%
  left_join(readCounts, by ='Seq_ID')%>%
  column_to_rownames('Sample')%>%
  as.data.frame()

phylumLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_phylum-MATRIX.txt')
phylumLevel = phylumLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()


maaslin3(
  input_data = phylumLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/phylum/maaslin3Res",
  formula = '~ Status +(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = phylumLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/phylum/maaslin3ResInteraction",
  formula = '~ Status+(1|Family)+Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

classLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_class-MATRIX.txt')
classLevel = classLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()

maaslin3(
  input_data = classLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/class/maaslin3Res",
  formula = '~ Status +(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = classLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/class/maaslin3ResInteraction",
  formula = '~ Status+(1|Family)+Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

orderLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_order-MATRIX.txt')
orderLevel = orderLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = orderLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/order/maaslin3Res",
  formula = '~ Status +(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = orderLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/order/maaslin3ResInteraction",
  formula = '~ Status+(1|Family)+Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

familyLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_family-MATRIX.txt')
familyLevel = familyLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plotsKahemba/family/maaslin3Res",
  formula = '~ Status +(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/family/maaslin3ResInteraction",
  formula = '~ Status+(1|Family)+Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

fit=maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plotsKahemba/family/maaslin3ResFamily",
  formula = '~ Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

fit$fit_data_abundance

 maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/family/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

genusLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
genusLevel = genusLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = genusLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/genus/maaslin3Res",
  formula = '~ Status+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = genusLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/genus/maaslin3ResInteraction",
  formula = '~ Status+(1|Family)+Status*(1|Family)+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = genusLevel,
  input_metadata = metadata,
  min_prevalence = .1,
  output = "plotsKahemba/genus/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

speciesLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')
speciesLevel = speciesLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = speciesLevel,
  cores = 1,
  min_prevalence = .1,
  input_metadata = metadata,
  output = "plotsKahemba/species/maaslin3Res",
  formula = '~ Status+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = speciesLevel,
  input_metadata = metadata,
  output = "plotsKahemba/species/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

