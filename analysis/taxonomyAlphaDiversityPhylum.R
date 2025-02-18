library(tidyverse,
        vegan)
library(ggpubr)
tss_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}
taxaData = read_tsv('08_TAXONOMY/taxonomyResults-t_phylum-MATRIX.txt')
taxaDataFiltered <- taxaData[rowSums(taxaData >= 5) >= 3, ]
taxaDataTss= taxaDataFiltered%>%
  column_to_rownames('taxon')%>%
  tss_normalize()
taxaDataTss = taxaDataTss %>%
  rownames_to_column('taxon')

meta = read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')

shannon=taxaDataTss%>%
  column_to_rownames('taxon')%>%
  t()%>%
  vegan::diversity(index = c('shannon'))%>%
  as.data.frame()%>%
  rename('shannon '=c(1))


invSimpson=taxaDataTss%>%
  column_to_rownames('taxon')%>%
  t()%>%
  vegan::diversity(index = c('invsimpson'))%>%
  as.data.frame()%>%
  rename('invSimpson '=c(1))

phylumAlpha = cbind(shannon,
                    invSimpson)%>%
  rownames_to_column('Sample')%>%
  left_join(meta,
            by = 'Sample')

phylumAlpha%>%
  ggplot(aes(x = Status,
             y = phylumAlpha$`invSimpson`))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('Kinshasa', 'Konzo'),
                                        c('Kinshasa', 'Masi-Manimba'),
                                        c('Kinshasa', 'Unaffected'),
                                        c('Konzo', 'Masi-Manimba'),
                                        c('Konzo', 'Unaffected'),
                                        c('Masi-Manimba', 'Unaffected')))

phylumAlpha%>%
  ggplot(aes(x = Disease,
             y = phylumAlpha$`invSimpson`))+
  geom_violin()+
  stat_compare_means()

phylumAlpha%>%
  ggplot(aes(x = Location,
             y = phylumAlpha$`invSimpson`))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('Kahemba', 'Kinshasa'),
                                        c('Kahemba', 'Masi-Manimba'),
                                        c('Kinshasa', 'Masi-Manimba')))

phylumAlpha%>%
  ggplot(aes(x = Status,
             y = phylumAlpha$`shannon`))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('Kinshasa', 'Konzo'),
                                        c('Kinshasa', 'Masi-Manimba'),
                                        c('Kinshasa', 'Unaffected'),
                                        c('Konzo', 'Masi-Manimba'),
                                        c('Konzo', 'Unaffected'),
                                        c('Masi-Manimba', 'Unaffected')))

phylumAlpha%>%
  ggplot(aes(x = Disease,
             y = phylumAlpha$`shannon`))+
  geom_violin()+
  stat_compare_means()

phylumAlpha%>%
  ggplot(aes(x = Location,
             y = phylumAlpha$shannon))+
  geom_violin()+
  stat_compare_means(comparisons = list(c('Kahemba', 'Kinshasa'),
                                        c('Kahemba', 'Masi-Manimba'),
                                        c('Kinshasa', 'Masi-Manimba')))