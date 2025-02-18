library(tidyverse)
library(ggpubr)
'10% or 15% is probably the best for species level'
tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
meta = meta %>%
  rename('sample'= Sample)

taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')

taxaDataLong=taxaData %>%
  pivot_longer(cols = c(2:293),
               names_to = 'sample',
               values_to = 'coverage')

fullDataLong=taxaDataLong %>%
  left_join(meta, by = 'sample')

fullDataAbundanceFiltered=fullDataLong %>%
  filter(coverage >= 3)


groupCounts=meta %>%
  group_by(Status) %>%
  summarise('groupCounts' = n())

#54 taxa passing konzo at 25% prevelance
# 79 at 20%
#107 at 15%
#176 at 10%
#289 at 5%
taxaPassingPrevelanceFilterKonzo=fullDataAbundanceFiltered %>%
  filter(Status == 'Konzo') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 95*.15)%>%
  .$taxon
length(taxaPassingPrevelanceFilterKonzo)


#53 taxa passing unaffected at 25% prevelance
#77 at 20%
#115 at 15%
#161 at 10%
#276 at 5%
taxaPassingPrevelanceFilterUnaffected=fullDataAbundanceFiltered %>%
  filter(Status == 'Unaffected') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 93*.15)%>%
  .$taxon
length(taxaPassingPrevelanceFilterUnaffected)

#2 taxa passing 25%
#80 at 20%
#127 at 15%
#161 at 10%
#289 at 5%
taxaPassingPrevelanceFilterMasiManimba=fullDataAbundanceFiltered %>%
  filter(Status == 'Masi-Manimba') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 77*.15)%>%
  .$taxon
length(taxaPassingPrevelanceFilterMasiManimba)

#2 taxa passing 25%
#56 at 20%
#79 at 15%
#123 at 10%
#255 at 5%
taxaPassingPrevelanceFilterKinshasa=fullDataAbundanceFiltered %>%
  filter(Status == 'Kinshasa') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 35*.15) %>%
.$taxon
length(taxaPassingPrevelanceFilterKinshasa)

#81 total taxa passing at 25% prevelance
#135 at 20%
#187 at 15%
# 252 at 10%
#430 at 5%

taxaPassingPrevelanceFilter=union(taxaPassingPrevelanceFilterKinshasa,
      taxaPassingPrevelanceFilterKonzo)%>%
  union(taxaPassingPrevelanceFilterUnaffected)%>%
  union( taxaPassingPrevelanceFilterMasiManimba)

taxaDataTss= taxaData%>%
  filter(taxon %in% taxaPassingPrevelanceFilter)%>%
  column_to_rownames('taxon')%>%
  tss_normalize()

meta=taxaData%>%
  column_to_rownames('taxon')%>%
  t()%>%
  rowSums()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  rename('librarySize' = c(2))%>%
  left_join(meta, by = "sample")

speciesPCAOut=taxaDataTss %>%
  t()%>%
  prcomp(center = T,
         scale. = T)

# 16% of bugs retained at 15% prevelance
# 23% of bugs retained at 10% prevelance
# 31% of bugs retained at 5% prevelance

length(taxaPassingPrevelanceFilter)/length(taxaData$taxon)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status))+
  geom_point(alpha = .5)+
  stat_ellipse()

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point(alpha = .5)+
  stat_ellipse(level = .95)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC3,
             y = PC4,
             col = Status))+
  geom_point(alpha = .5)+
  stat_ellipse(level = .95)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC4,
             y = PC5,
             col = Status))+
  geom_point(alpha = .5)+
  stat_ellipse(level = .95)








