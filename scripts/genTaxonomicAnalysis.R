library(vegan)
library(ggpubr)
library(tidyverse)

'Function to filter the SCGs in a manner similar to 16S studies
Accepts an numceric for the min abundance in each sample, and a
integer for the number of samples this taxon needs to be present to pass.
Also needs a dataframe from anvio to filter(duh)'
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

# wrapper for the previous functions, just call one function and get
#processed dataframe, reduce the error ig. p is prevelance (see filter basic)
# a is abundance (see filter basic)
filterAndNorm <- function(rawTaxTab, p,a){
  rawTaxTab = rawTaxTab %>%
    column_to_rownames('taxon')
  filteredTab = filterBasic(rawTaxTab, prevalance = p, abundance = a)
  normalizedTab = relab_normalize(filteredTab)*100
  return(normalizedTab)
}
#just going to be easier to keep the metadata and taxonomic data in a single place1
pivotMergeMeta = function(normalizedTaxtab, metadata){
  colCount= ncol(normalizedTaxtab)+1
  fulldataLong=normalizedTaxtab %>%
    rownames_to_column('taxon')%>%
    pivot_longer(cols = c(2:colCount),
                 names_to = 'Sample',
                 values_to = 'relativeAbundance')%>%
    left_join(metadata, by = 'Sample')
  return(fulldataLong)
}
#Only adding labels for what is what if the fuckin things are
#manangable
stackedAreaChart = function(longFull, facet, path){
  taxCount=longFull%>%
    select(taxon)%>%
    distinct()%>%
    nrow()
  if (taxCount < 50){


    p=ggplot(data = longFull,
           aes(x = Sample,
               y = relativeAbundance,
               fill = taxon))+
      geom_col()+
      facet_wrap(as.formula(paste0('~',facet)), scales = 'free_x', nrow = 1)+
      scale_x_discrete(breaks = NULL)+
      theme(legend.text = element_text(size = 6))
  }
  else{
    p=ggplot(data = longFull,
           aes(x = Sample,
               y = relativeAbundance,
               fill = taxon))+
      geom_col()+
      facet_wrap(as.formula(paste0('~', facet)), scales = 'free_x', nrow = 1)+
      scale_x_discrete(breaks = NULL)+
      theme(legend.position = 'none')
  }
  ggsave(filename = paste0(path,'/','stackedArea', facet,'.pdf'), plot = p,width = 8, height = 6, units = 'in')
  return(p)
}

shannon = function(longFull, comparison, meta,path){
  comparisons <- combn(unique(meta[[comparison]]), 2, simplify = FALSE)
  p=longFull %>%
    pivot_wider(id_cols = 'Sample', names_from = 'taxon', values_from = 'relativeAbundance')%>%
    column_to_rownames('Sample')%>%
    diversity(index = 'shannon')%>%
    as.data.frame()%>%
    rename('shannon'=c(1))%>%
    rownames_to_column('Sample')%>%
    left_join(meta, by = 'Sample')%>%
    ggplot(aes(x = .data[[comparison]],
               y = shannon))+
    geom_boxplot()+
    stat_compare_means(comparisons = comparisons)+
    labs(title = paste("Shannon", comparison))
  ggsave(filename = paste0(path,'/','shannon', comparison,'.pdf'), plot = p, width = 8, height = 6, units = 'in')
  return(p)
}

simpson = function(longFull, comparison, meta,path){
  comparisons <- combn(unique(meta[[comparison]]), 2, simplify = FALSE)
  p=longFull %>%
    pivot_wider(id_cols = 'Sample', names_from = 'taxon', values_from = 'relativeAbundance')%>%
    column_to_rownames('Sample')%>%
    diversity(index = 'simpson')%>%
    as.data.frame()%>%
    rename('simpson'=c(1))%>%
    rownames_to_column('Sample')%>%
    left_join(meta, by = 'Sample')%>%
    ggplot(aes(x = .data[[comparison]],
               y = simpson))+
    geom_boxplot()+
    stat_compare_means(comparisons = comparisons)+
    labs(title =(paste('Simpson diversity by ', comparison)))
  ggsave(filename = paste0(path,'/','simpson', comparison,'.pdf'), plot = p,width = 8, height = 6, units = 'in')
  return(p)
}

bray = function(longFull, comparison, meta,path){
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
               y = V2,
               col = .data[[comparison]]))+
    geom_point()+
    labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
         y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
         title = paste("Bray PCoA", comparison),
         caption = paste('P = ', adonisRes$`Pr(>F)`))
  ggsave(filename = paste0(path,'/','bray', comparison,'.pdf'), plot = p,width = 8, height = 6, units = 'in')
  return(p)
}

jaccard = function(longFull, comparison, meta,path){
  adonisIn=longFull %>%
    pivot_wider(id_cols = 'Sample', names_from = 'taxon', values_from = 'relativeAbundance')%>%
    column_to_rownames('Sample')

  distMat=vegdist(method = 'jaccard', adonisIn)

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
               y = V2,
               col = .data[[comparison]]))+
    geom_point()+
    labs(x = paste0('PC1 - ', pcoaOut$eig[1]/totalvar),
         y = paste0('PC2 - ', pcoaOut$eig[2]/totalvar),
         title = paste("Jaccard PCoA", comparison),
         caption = paste('P = ', adonisRes$`Pr(>F)`))
  ggsave(filename = paste0(path,'/','jaccard', comparison,'.pdf'), plot = p,width = 8, height = 6, units = 'in')
  return(p)
}
wrapper = function(taxonomy, metadata){
  file = paste0('08_TAXONOMY/taxonomyResults-t_',taxonomy,'-MATRIX.txt')
  input=read_tsv(file)
  path = paste0('plots/',taxonomy)
  dir.create(path)
  filterdAndNormalizedTaxTab=filterAndNorm(rawTaxTab =input, a = 5, p = 15)
  longfilteredAndNormalizedTaxa=pivotMergeMeta(filterdAndNormalizedTaxTab, metadata = metadata)
  stackedAreaChart(longFull=longfilteredAndNormalizedTaxa, facet = 'Location', path =path)
  stackedAreaChart(longFull=longfilteredAndNormalizedTaxa, facet = 'Status', path =path)
  stackedAreaChart(longFull=longfilteredAndNormalizedTaxa, facet = 'Disease', path =path)
  shannon(longFull=longfilteredAndNormalizedTaxa, comparison='Location', meta=metadata,path=path)
  shannon(longFull=longfilteredAndNormalizedTaxa, comparison='Status', meta=metadata,path=path)
  shannon(longFull=longfilteredAndNormalizedTaxa, comparison='Disease', meta=metadata,path=path)
  simpson(longFull=longfilteredAndNormalizedTaxa, comparison='Location', meta=metadata,path=path)
  simpson(longFull=longfilteredAndNormalizedTaxa, comparison='Status', meta=metadata,path=path)
  simpson(longFull=longfilteredAndNormalizedTaxa, comparison='Disease', meta=metadata,path=path)
  col_names <- colnames(metadata)[!colnames(metadata) %in% exclude_columns]
  for (col_name in col_names) {
    # Call the jaccard function with the current column name as comparison
    jaccard(longFull = longfilteredAndNormalizedTaxa, comparison = col_name, meta = metadata, path = path)
  }
  for (col_name in col_names) {
    # Call the jaccard function with the current column name as comparison
    bray(longFull = longfilteredAndNormalizedTaxa, comparison = col_name, meta = metadata, path = path)
  }
}

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')
exclude_columns <- c("Sample", "Number", "Total Elution Volume", "Sample_ID", "Seq_ID", "Stage","Age_T0", "Discordant")
levels = c('phylum', 'class', 'order', 'family', 'genus', 'species')

for (level in levels){
  wrapper(level, metadata=metadata)
}

