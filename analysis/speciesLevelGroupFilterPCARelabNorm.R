library(tidyverse)
library(ggpubr)
library(vegan)
'10% or 15% is probably the best for species level'
'25% looks better when using bray curtis.'
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
  filter(prevelance >= 95*.20)%>%
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
  filter(prevelance >= 93*.20)%>%
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
  filter(prevelance >= 77*.20)%>%
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
  filter(prevelance >= 35*.20) %>%
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
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status,
             group = Family))+
  geom_point(alpha = .5)+
  geom_line()

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

'Horseshoe with jaccard at 25% prevelance, not ideal'
'Robust aitchison looks weird, need to read about it. Maybe filtering
is not as nessicary'
'I like chao, I believe this is placing weight on low abundance taxa.
Significant at p < .001 for status. Low R2 however. Not sure this is surpising'
'Family explains the most variation at .82, p < .001'
distMatBray=taxaDataTss %>%
  t()%>%
  vegdist(method = 'chao')

adonis2(distMatBray~Status, data = column_to_rownames(meta,'sample'))
adonis2(distMatBray~Location, data = column_to_rownames(meta,'sample'))
adonis2(distMatBray~Sex, data = column_to_rownames(meta,'sample'))
#adonis2(distMatBray~Family, data = column_to_rownames(meta,'sample'))


pcoaOut=cmdscale(distMatBray, k = 4 , eig = T)
totalvar = sum(pcoaOut$eig)

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  stat_ellipse()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  stat_ellipse()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Status,
             group = Family))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  geom_line()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))+
  stat_ellipse()

'Sex is still insignficant with bray, status and location remain significant, with smaller R2'
distMatBray=taxaDataTss %>%
  t()%>%
  vegdist(method = 'bray')

adonis2(distMatBray~Status, data = column_to_rownames(meta,'sample'))
adonis2(distMatBray~Location, data = column_to_rownames(meta,'sample'))
adonis2(distMatBray~Sex, data = column_to_rownames(meta,'sample'))
adonis2(distMatBray~Family, data = column_to_rownames(meta,'sample'))


pcoaOut=cmdscale(distMatBray, k = 4 , eig = T)
totalvar = sum(pcoaOut$eig)

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  stat_ellipse()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Status,
             group = Family))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  geom_line()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))+
  stat_ellipse()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))+
  stat_ellipse()








