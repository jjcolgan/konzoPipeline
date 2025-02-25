library(maaslin3)
library(tidyverse)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
readCounts = read_tsv('read_counts.tsv')
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
  output = "plots/phylum/maaslin3Res",
  formula = '~ Status + Sex+readCounts',
  reference = c("Status,Unaffected"))

classLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_class-MATRIX.txt')
classLevel = classLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = classLevel,
  input_metadata = metadata,
  output = "plots/class/maaslin3Res",
  formula = '~ Status+Sex+readCounts',
  reference = c("Status,Unaffected"))

orderLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_order-MATRIX.txt')
orderLevel = orderLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = orderLevel,
  input_metadata = metadata,
  output = "plots/order/maaslin3Res",
  formula = '~ Status+Sex+readCounts',
  reference = c("Status,Unaffected"))

familyLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_family-MATRIX.txt')
familyLevel = familyLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plots/family/maaslin3Res",
  formula = '~ Status+Sex+readCounts',
  reference = c("Status,Unaffected"))
maaslin3(
  input_data = familyLevel,
  input_metadata = metadata,
  output = "plots/family/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

genusLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
genusLevel = genusLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = genusLevel,
  input_metadata = metadata,
  output = "plots/genus/maaslin3Res",
  formula = '~ Status+Sex+readCounts',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = genusLevel,
  input_metadata = metadata,
  output = "plots/genus/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

speciesLevel = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')
speciesLevel = speciesLevel%>%
  column_to_rownames('taxon')%>%
  as.data.frame()
maaslin3(
  input_data = speciesLevel,
  cores = 4,
  input_metadata = metadata,
  output = "plots/species/maaslin3Res",
  formula = '~ Status+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

maaslin3(
  input_data = speciesLevel,
  input_metadata = metadata,
  output = "plots/species/maaslin3ResStage",
  formula = '~ Stage+readCounts+(1|Family)',
  reference = c("Status,Unaffected"))

