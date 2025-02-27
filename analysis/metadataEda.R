library(tidyverse)
library(ggpubr)

metadata=read_csv('/Users/johnjamescolgan/Downloads/Konzo_Metagenomics_2021_Meta (2).csv')

metadata %>%
  ggplot(aes(x = Age_T0,
         fill = Location))+
  geom_density()

metadata %>%
  ggplot(aes(x = Age_T0,
             fill = Status,
             alpha = .25))+
  geom_density()

summary(metadata$Age_T0)

metadata %>%
  ggplot(aes(x = Sex,
             fill = Status,
             alpha = .25))+
  geom_density()

metadata %>%
  ggplot(aes(x = Sex))+
  geom_histogram(bins = 2, stat = 'count')+
  facet_wrap(~Status)