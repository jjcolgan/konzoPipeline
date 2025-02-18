library(tidyverse)
tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

'Unscaled pc1 v 2 shows seperation on region being driven by pc 2.
There ellipsessuggest that there may be seperation in konzo v non konzo.
The larger variance in previous plots at the order and family level are not really
apparant for Konzo.'

taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_genus-MATRIX.txt')
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

genusPCAOutNoScaling=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  prcomp(center = T)

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Location))+
  geom_point()+stat_ellipse(level = .9)

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = librarySize))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status))+
  geom_point()+stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Disease))+
  geom_point()+stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Sex))+
  geom_point()+stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')



genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Age_T0))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(Stage)))+
  geom_point()+
  stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Concentration))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Location))+
  geom_point()


genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = librarySize))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point()+stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Disease))+
  geom_point()+stat_ellipse()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Sex))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Age_T0))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Stage))+
  geom_point()

genusPCAOutNoScaling$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Concentration))+
  geom_point()

genusPCAOut=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  prcomp(center = T,
         scale = T)

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Location))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = librarySize))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Status))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Disease))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Sex))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')


genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Age_T0))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = as.factor(Stage)))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC1,
             y = PC2,
             col = Concentration))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Location))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = librarySize))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Status))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Disease))+
  geom_point()+stat_ellipse()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Sex))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Family))+
  geom_point()+
  theme(legend.position = 'none')

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Age_T0))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Stage))+
  geom_point()

genusPCAOut$x%>%
  as.data.frame()%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = PC2,
             y = PC3,
             col = Concentration))+
  geom_point()