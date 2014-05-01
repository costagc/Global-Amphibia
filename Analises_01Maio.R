setwd('~/Documents/Manuscripts/ProjetoMarilia')

tabela <- read.csv('~/Dropbox/allamphibians IUCN_July/Final/tabela_preds_final.csv')
tabela.cell <- tabela[,-c(1,8)]
rm(tabela)
tabela <- read.csv('~/Dropbox/allamphibians IUCN_July/Final/tabela_species_final.csv')
tabela.species <- tabela[,2:30]
head(tabela.species)
