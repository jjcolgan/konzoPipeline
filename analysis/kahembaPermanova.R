library(aPCoA)
library(tidyverse)
library(vegan)


bray = function(longFull, comparison, meta){
  adonisIn=longFull %>%
    pivot_wider(id_cols = 'Sample', names_from = 'taxon', values_from = 'relativeAbundance')%>%
    column_to_rownames('Sample')

  distMat=vegdist(method = 'bray', adonisIn)

  pcoaOut=cmdscale(distMat, k = 4 , eig = T)
  totalvar = sum(pcoaOut$eig)
  adonisMeta = meta%>%
    filter(Sample %in% rownames(adonisIn))%>%
    column_to_rownames('Sample')
  adonisRes=adonis2(as.formula(paste('adonisIn ~', comparison)),
                    data = adonisMeta)

  p=pcoaOut$points%>%
    as.data.frame()%>%
    rownames_to_column('Sample')%>%
    left_join(metadata, by = 'Sample')%>%
    ggplot(aes(x = V1,
               group = Family,
               y = V2,
               col = .data[[comparison]]))+
    geom_point()+
    geom_line(alpha = .25,
              col = 'black' )+
    labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
         y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
         title = paste("Bray PCoA", comparison),
         caption = paste('P = ', adonisRes$`Pr(>F)`))
  #ggsave(filename = paste0(path,'/','bray', comparison,'.pdf'), plot = p,width = 8, height = 6, units = 'in')
  return(p)
}

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
#filter for families in both konzo and unaffected
  metadata$Family = as.factor(metadata$Family)

  konzoFamilies =metadata %>%
    filter(Location == 'Kahemba')%>%
    group_by(Family)%>%
    summarise('nFamilyMembers'= n())%>%
    filter(nFamilyMembers > 1)%>%
    filter(Family != '98')%>%
    .$Family


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

  filteredAndNormTaxaTabLong = filteredAndNormTaxaTab %>%
    rownames_to_column('taxon')%>%
    pivot_longer(cols = c(2:183), names_to = 'Sample', values_to = 'relativeAbundance')%>%
    left_join(metadata, by = 'Sample')

  adonisIn=filteredAndNormTaxaTab %>%
    t()

  adonisMeta = metadata %>%
    filter(Sample %in% row.names(adonisIn),
           Family %in% konzoFamilies)%>%
    column_to_rownames('Sample')

  adonisIn=adonisIn %>%
    as.data.frame()%>%
    rownames_to_column('Sample')%>%
    filter(Sample %in% rownames(adonisMeta))%>%
    column_to_rownames('Sample')

  adonisMeta%>%
    group_by(Family)%>%
    summarise(n())%>%
    filter(`n()` < 2)

  bray=vegdist(as.matrix(adonisIn), method = 'bray')
  results=aPCoA(bray ~ Family, data = adonisMeta,maincov = Disease, drawCenter =F, drawEllipse = T)
  #adonis2(adonisIn ~ Disease * Sex, data = adonisMeta, by = 'terms')%>%
    #print()
#use this model, it accouts for disease and family. Since the interaction is still significant, that
  #means that is capturing variation which cannot be attributibale to family alone
  adonis2(adonisIn ~ Disease + Family+Disease * Family, data = adonisMeta, by = 'terms')%>%
    print()
  adonisMeta %>%
    row
  bd=betadisper( d = bray, group = adonisMeta$Disease)
  permutest(bd)%>%
    print()

  #adonis2(adonisIn ~ Disease * Family, data = adonisMeta, by = 'terms')%>%
    #print()

  #adonis2(adonisIn ~ Disease * Age_T0, data = adonisMeta, by = 'terms')%>%
    #print()

  results$plotMatrix
  #p=bray(longFull = filteredAndNormTaxaTabLong, comparison = 'Disease', meta = adonisMeta)
  ##ggsave(plot = p, filename = paste0('/Users/johnjamescolgan/Library/CloudStorage/Box-Box/automatedAssemblyPipeline/plotsKahemba/pcoasFamilyMembersConnected/brayCurtis', taxonomy, '.pdf'))

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