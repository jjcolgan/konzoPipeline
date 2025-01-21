library(tidyverse)
library(vegan)
library(ggpubr)
library(ape)
library(ggExtra)
library(Maaslin2)

tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_species-MATRIX.txt')
taxaDataFiltered <- taxaData[rowSums(taxaData >= 5) >= 3, ]
taxaDataTss= taxaDataFiltered%>%
  column_to_rownames('taxon')%>%
  tss_normalize()
taxaDataTss = taxaDataTss %>%
  rownames_to_column('taxon')

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
meta = meta %>%
  rename('sample'= Sample)

shannonOut= taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t() %>%
  diversity(index = 'shannon')%>%
  as.data.frame()%>%
  rename('shannon' = c(1))%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')

shannonOut %>%
  ggplot(aes(x = Status,
             y = shannon))+
  geom_violin()+stat_compare_means(comparisons = list(c('Kinshasa', 'Konzo'),
                                                      c('Kinshasa', 'Masi-Manimba'),
                                                      c('Kinshasa', 'Unaffected'),
                                                      c('Konzo', 'Masi-Manimba'),
                                                      c('Konzo', 'Unaffected'),
                                                      c('Masi-Manimba', 'Unaffected')))

simpsonOut= taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t() %>%
  diversity(index = 'invsimpson')%>%
  as.data.frame()%>%
  rename('simpson' = c(1))%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')

simpsonOut %>%
  ggplot(aes(x = Status,
             y = simpson))+
  geom_violin()+stat_compare_means(comparisons = list(c('Kinshasa', 'Konzo'),
                                                      c('Kinshasa', 'Masi-Manimba'),
                                                      c('Kinshasa', 'Unaffected'),
                                                      c('Konzo', 'Masi-Manimba'),
                                                      c('Konzo', 'Unaffected'),
                                                      c('Masi-Manimba', 'Unaffected')))


scores=taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t() %>%
  prcomp(center = TRUE)

scores$x %>%
  as.data.frame() %>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample') %>%
  ggplot(aes(x = PC1,
             col = Status,
             y = PC2))+
  geom_point()

brayDist<-taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  vegdist()
pcoaRes=wcmdscale(brayDist)

adonisData= taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()

adonisMeta = meta %>%
  filter(sample %in% rownames(adonisData))

adonis2(brayDist~Status, data = adonisMeta)

pcoaRes %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = Status,
             y = V2))+
geom_point()

ggMarginal(p = brayCoa,
           type = 'box',
           groupColour = T)

brayDist<-taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()%>%
  vegdist(method = 'jaccard')
pcoaRes=wcmdscale(brayDist,eig = T,add=T)
adonisData= taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()

adonisMeta = meta %>%
  filter(sample %in% rownames(adonisData))

adonis2(brayDist~Status, data = adonisMeta)
(pcoaRes$eig / sum(pcoaRes$eig))*100
pcoaRes$points %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(meta, by = 'sample')%>%
  ggplot(aes(x = Dim1,
             col = Status,
             y = Dim2))+
  geom_point()+
  labs(x = 'PCo1 = 9.73%',
       y = 'PCo2 = 2.044%')

ruralMeta=meta %>%
  filter(Location != 'Kinshasa')
ruralTaxa=taxaData %>%
  column_to_rownames('taxon')%>%
  t() %>%
  as.data.frame()%>%
  rownames_to_column('sample') %>%
  filter(sample %in% ruralMeta$sample)%>%
  column_to_rownames('sample')%>%
  t()%>%
  as.data.frame()

ruralTaxaFiltered <- ruralTaxa[rowSums(ruralTaxa >= 5) >= 3, ]
ruralTss = ruralTaxaFiltered %>%
  tss_normalize()

shannonOut= ruralTss %>%t()%>%
  diversity(index = 'shannon')%>%
  as.data.frame()%>%
  rename('shannon' = c(1))%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')

shannonOut %>%
  ggplot(aes(x = Status,
             y = shannon))+
  geom_violin()+stat_compare_means(comparisons = list(
                                                      c('Konzo', 'Masi-Manimba'),
                                                      c('Konzo', 'Unaffected'),
                                                      c('Masi-Manimba', 'Unaffected')))

simpsonOut= ruralTss %>%t()%>%
  diversity(index = 'invsimp')%>%
  as.data.frame()%>%
  rename('simpson' = c(1))%>%
  rownames_to_column('sample')%>%
  left_join(meta, by = 'sample')

simpsonOut %>%
  ggplot(aes(x = Status,
             y = simpson))+
  geom_violin()+stat_compare_means(comparisons = list(
  c('Konzo', 'Masi-Manimba'),
  c('Konzo', 'Unaffected'),
  c('Masi-Manimba', 'Unaffected')))

brayDist<-ruralTss %>%
  t()%>%
  vegdist()
pcoaRes=wcmdscale(brayDist)


adonisMeta = meta %>%
  filter(sample %in% colnames(ruralTss))

adonis2(brayDist~Status, data = adonisMeta)
adonis2(brayDist~Disease, data = adonisMeta)

pcoaRes %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(meta, by = 'sample')%>%
  ggplot(aes(x = V1,
             col = Disease,
             y = V2))+
  geom_point()

brayDist<-ruralTss %>%
  t()%>%
  vegdist(method = 'jaccard')
pcoaRes=wcmdscale(brayDist,eig = T,add=T)
adonisData= taxaDataTss %>%
  column_to_rownames('taxon')%>%
  t()

adonisMeta = meta %>%
  filter(sample %in% rownames(adonisData))

adonis2(brayDist~Status, data = adonisMeta)
(pcoaRes$eig / sum(pcoaRes$eig))*100
pcoaRes$points %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  merge(meta, by = 'sample')%>%
  ggplot(aes(x = Dim1,
             col = Status,
             y = Dim2))+
  geom_point()+
  labs(x = 'PCo1 = 9.73%',
       y = 'PCo2 = 2.044%')
maaslinRuralMeta = ruralMeta %>%
  column_to_rownames('sample')
Maaslin2(input_data = ruralTaxa,
         output = 'firstPassRural',
         fixed_effects = 'Status',
         reference = c('Status,Masi-Manimba'),
         input_metadata = maaslinRuralMeta,
)