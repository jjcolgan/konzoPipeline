library(tidyverse)
library(vegan)
tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}
'Looks like region is driving some level of segration, otherwise not seeing much.
Better seperation on region when using the scaled PCAs. Mainly driven
by PC1, pattern is much less apparent when looking at pc2 - 3. Scaled PC2 / 3
seems to show greater variance in composition of konzo relative to unaffected'

taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_order-MATRIX.txt')
taxaDataFiltered <- taxaData[rowSums(taxaData >= 5) >= 3, ]
taxaDataTss= taxaDataFiltered%>%
  column_to_rownames('taxon')%>%
  tss_normalize()
taxaDataTss = taxaDataTss %>%
  rownames_to_column('taxon')

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
meta = meta %>%
  rename('sample'= Sample)

taxaData[c(2:293)]%>%
  colSums()%>%
  hist()

meta=taxaData%>%
  column_to_rownames('taxon')%>%
  t()%>%
  rowSums()%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  rename('librarySize' = c(2))%>%
  left_join(meta, by = "sample")

orderPCAOutNoScaling=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  prcomp(center = T)

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Location))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = librarySize))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Disease))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Sex))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')



orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Age_T0))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(Stage)))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Concentration))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Location))+
  geom_point()+stat_ellipse()


orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = librarySize))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Disease))+
  geom_point()+stat_ellipse()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Sex))+
  geom_point()

;orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Age_T0))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Stage))+
  geom_point()

orderPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Concentration))+
  geom_point()

orderPCAOut=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  prcomp(center = T,
         scale = T)


orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Location))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = librarySize))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Disease))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Sex))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Age_T0))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Stage))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Concentration))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Location))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = librarySize))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Disease))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Sex))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Age_T0))+
  geom_point()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = as.factor(Stage)))+
  geom_point()+stat_ellipse()

orderPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Concentration))+
  geom_point()