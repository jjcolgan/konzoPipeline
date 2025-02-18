library(vegan)
library(tidyverse)
library(ggExtra)

tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

'Sep driven by pc 2 (location)'
taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
taxaDataFiltered <- taxaData[rowSums(taxaData >= 5) >= 3, ]
taxaDataTss= taxaDataFiltered%>%
  column_to_rownames('taxon')%>%
  tss_normalize()
taxaDataTss = taxaDataTss %>%
  rownames_to_column('taxon')

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')

distMatBray=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  vegdist()

pcoaOut=cmdscale(distMatBray, k = 4 , eig = T)
totalvar = sum(pcoaOut$eig)

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V1,
             y = V2))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))
pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))+
  stat_ellipse()

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V1,
             y = V2,
             col = Disease))+
  geom_point()+
  labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
       y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V2,
             y = V3,
             col = Disease))+
  geom_point()+
  labs(x = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
       y = paste0('PC3 - ', pcoaOut$eig[3]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))

pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Location))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))
dir.create('plots')
p=pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Status))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))


pcoaOut$points%>%
  as.data.frame()%>%
  rownames_to_column('Sample')%>%
  left_join(meta, by = 'Sample')%>%
  ggplot(aes(x = V3,
             y = V4,
             col = Disease))+
  geom_point()+
  labs(x = paste0('PC3 - ', pcoaOut$eig[3]/totalvar),
       y = paste0('PC4 - ', pcoaOut$eig[4]/totalvar))
