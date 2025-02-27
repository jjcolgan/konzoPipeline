library(tidyverse)
library(vegan)
filterBasic = function(prevalance,abundance, unfilteredTaxTab){
  filteredTaxTab <- unfilteredTaxTab[rowSums(unfilteredTaxTab >= abundance) >= prevalance, ]
  return(filteredTaxTab)
}

#Relab normalize, accepts a filtered taxonomic data frame, returns a relab norm
relab_normalize <- function(df) {
  # Ensure all values are numeric
  # Calculate the column sums (total abundance per sample)
  col_sums <- colSums(df, na.rm = TRUE)

  # Perform TSS normalization: divide each value by the sum of its column
  df_normalized <- sweep(df, 2, col_sums, FUN = "/")

  return(df_normalized)
}

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
levels = c('phylum', 'class', 'order', 'family', 'genus', 'species')
metadata=metadata %>%
  filter(Location =='Kahemba')%>%
  filter(is.na(Age_T0 ) != T)

for (taxonomy in levels){
  file = paste0('08_TAXONOMY/taxonomyResults-t_',taxonomy,'-MATRIX.txt')
  input = read_tsv(file)
  input <- input %>%
    column_to_rownames('taxon')%>%
    select(any_of(metadata$Sample))
  'need to account for taxon column'
  nsamples = ncol(input)
  p = .1 *nsamples
  a = 5
  filteredAndNormTaxaTab=filterBasic(unfilteredTaxTab = input,abundance = a, prevalance = p)%>%
    relab_normalize()*100
  adonisIn=filteredAndNormTaxaTab %>%
    t()
  adonisMeta = metadata %>%
    filter(Sample %in% row.names(adonisIn))
  print(taxonomy)

  adonis2(adonisIn ~ Disease * Sex, data = adonisMeta, by = 'terms')%>%
    print()
#use this model, it accouts for disease and family. Since the interaction is still significant, that
  #means that is capturing variation which cannot be attributibale to family alone
  adonis2(adonisIn ~ Disease + Family+Disease * Family, data = adonisMeta, by = 'terms')%>%
    print()

  adonis2(adonisIn ~ Disease * Family, data = adonisMeta, by = 'terms')%>%
    print()

  adonis2(adonisIn ~ Disease * Age_T0, data = adonisMeta, by = 'terms')%>%
    print()

}

konzoFamilies=metadata %>%
  filter(Location == 'Kahemba')%>%
  filter(Disease == 'Konzo')%>%
  select(Family)

unaffectedFamilies=metadata %>%
  filter(Location == 'Kahemba')%>%
  filter(Disease == 'Unaffected')%>%
  select(Family)

metadata %>% filter(Family %in% konzoFamilies$Family)

metadata %>%
  filter(Location == 'Kahemba')