library(tidyverse)
library(ggpubr)
library(vegan)

tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
meta$tage = factor(meta$Stage, levels = c(NA, 0,1,2,3))
meta = meta %>%
  rename('sample'= Sample)%>%
  filter(is.na(Stage)!= T)

taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')

taxaData %>%
  column_to_rownames('taxon')%>%
  colSums()%>%
  view()
taxaData %>%
  column_to_rownames('taxon')%>%
  rowSums()%>%
  view()

taxaDataLong=taxaData %>%
  pivot_longer(cols = c(2:293),
               names_to = 'sample',
               values_to = 'coverage')%>%
  filter(sample %in% meta$sample)

fullDataLong=taxaDataLong %>%
  left_join(meta, by = 'sample')

fullDataAbundanceFiltered=fullDataLong %>%
  filter(coverage >= 3)

groupCounts=meta %>%
  group_by(Status) %>%
  summarise('groupCounts' = n())
#54 passing at 25%
#78 passing at 20$
#107 at 15%
taxaPassingPrevelanceFilterKonzo=fullDataAbundanceFiltered %>%
  filter(Status == 'Konzo') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 95*.1)%>%
  .$taxon
length(taxaPassingPrevelanceFilterKonzo)

#50 taxa passing at 25%, this is not the same as when filtering by Status??
#76 passing at 20%
#115 at 15%
taxaPassingPrevelanceFilterUnaffected=fullDataAbundanceFiltered %>%
  filter(Status == 'Unaffected') %>%
  group_by(taxon)%>%
  summarise('prevelance' = n())%>%
  filter(prevelance >= 93*.1)%>%
  .$taxon
length(taxaPassingPrevelanceFilterUnaffected)



taxaPassingPrevelanceFilter=union(taxaPassingPrevelanceFilterUnaffected,
                                  taxaPassingPrevelanceFilterKonzo)

#4% of taxa pass this filter(25%)
#8 percent pass at 20%
length(taxaPassingPrevelanceFilter)/length(taxaData$taxon)

taxaDataTss= taxaData%>%
  filter(taxon %in% taxaPassingPrevelanceFilter)%>%
  column_to_rownames('taxon')%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  filter(sample %in% meta$sample)%>%
  column_to_rownames('sample')%>%
  t()%>%
  tss_normalize()

speciesPCAOut=taxaDataTss %>%
  t()%>%
  prcomp(center = T,
         scale. = T)



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
             col = Status))+
  geom_point(alpha = .5)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(Stage)))+
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
  stat_ellipse()

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status,
             group = Family))+
  geom_point(alpha = .5)+
  geom_line(alpha = .5)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point(alpha = .5)

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = as.factor(Stage)))+
  geom_point(alpha = .5)+
  stat_ellipse()

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC3,
             y = PC4,
             col = Status))+
  geom_point(alpha = .5)+
  stat_ellipse()

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC3,
             y = PC4,
             col = as.factor(Stage)))+
  geom_point(alpha = .5)+
  stat_ellipse()

speciesPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC3,
             y = PC4,
             col = Status,
             group = Family))+
  geom_point(alpha = .5)+
  geom_line(alpha = .5,
            col = 'grey')

distMatBray=taxaDataTss %>%
  t()%>%
  vegdist(method = 'jaccard')

adonisMeta = meta %>%
  filter(sample %in% colnames(taxaDataTss))%>%
  column_to_rownames('sample')
'Nothing significant using jaccard'
adonis2(distMatBray~Status, data = adonisMeta)
adonis2(distMatBray~Sex, data =  adonisMeta)
#adonis2(distMatBray~Family,  data = adonisMeta)
adonis2(distMatBray~Stage,  data = adonisMeta)


pcoaOut=cmdscale(distMatBray, k = 4 , eig = T)
totalvar = sum(pcoaOut$eig)

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = as.factor(Stage)))+
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
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))

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
             col = Stage))+
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
             col = Stage))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))+
  stat_ellipse()

distMatBray=taxaDataTss %>%
  t()%>%
  vegdist(method = 'bray')

adonisMeta = meta %>%
  filter(sample %in% colnames(taxaDataTss))%>%
  column_to_rownames('sample')
'Nothing significant using bray-curtis'
adonis2(distMatBray~Status, data = adonisMeta)
adonis2(distMatBray~Sex, data =  adonisMeta)
#adonis2(distMatBray~Family,  data = adonisMeta)
adonis2(distMatBray~Stage,  data = adonisMeta)


pcoaOut=cmdscale(distMatBray, k = 4 , eig = T)
totalvar = sum(pcoaOut$eig)

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = as.factor(Stage)))+
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
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))

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
             col = Stage))+
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
             col = Stage))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))+
  stat_ellipse()