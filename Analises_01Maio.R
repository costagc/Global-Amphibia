setwd('~/Documents/Manuscripts/ProjetoMarilia')
library(raster)
library(rgdal)
###Load in Tables already done##########
tabela <- read.csv('~/Dropbox/allamphibians IUCN_July/Final/tabela_preds_final.csv')
tabela.cell <- tabela[,-c(1,8)]
rm(tabela)
tabela <- read.csv('~/Dropbox/allamphibians IUCN_July/Final/tabela_species_final.csv')
tabela.species <- tabela[,2:30]
head(tabela.species)
head(tabela.cell)
######Load in table with larvae development info#######
tabela.larvae <- read.csv('~/Dropbox/allamphibians IUCN_July/2951 spp.csv')
head(tabela.larvae)
tabela.larvae<- tabela.larvae[,c(1,12)]
names(tabela.larvae)<- c('Species','Larvae.Site')
tabela.larvae <- subset(tabela.larvae, !is.na(Larvae.Site))

####Load Maps and produce SiteXspp table#####
labels <- paste(tabela.larvae[,1],'.grd',sep='')
labels <- gsub('_',' ',labels)
maps <- stack(paste('~/Documents/Manuscripts/ProjetoMarilia/mapas raster/',labels,sep=''))
names(maps) <- gsub('.grd',' ',labels)
sp.pres <- as.data.frame(maps)
xy.values <- xyFromCell(maps[[1]], 1:ncell(maps[[1]]))
siteXspp <- cbind(xy.values,sp.pres)
head(siteXspp[,2870:2873])
siteXspp$rich <- rowSums(siteXspp[,3:2872],na.rm = T)
siteXspp <- subset(siteXspp, rich != 0)

#####test by plotting richness####
source('~/hubiC/Documents/Coisas do R/functions/DFtoRaster.R')
test<- DFtoRaster(siteXspp,d=6)
plot(test)
######Calculate richness of aquatic larvae and proportion of aquatic and terrestrial larvae#####
aqua <- subset(tabela.larvae, Larvae.Site == 4)
sp.aqua <- aqua$Species
sp.aqua <- gsub('_','.',sp.aqua)
sp.aqua2 <- subset(siteXspp, select=c(sp.aqua))
aqua.rich <- rowSums(sp.aqua2, na.rm = T)
prop.aqua <- aqua.rich/siteXspp$rich
prop.terr <- (siteXspp$rich - aqua.rich)/siteXspp$rich
siteXspp <- siteXspp[,c(1,2,2873,2874,2875,3:2872)]
siteXspp <- cbind(siteXspp,aqua.rich,prop.aqua)
siteXspp <- cbind(siteXspp,prop.terr)
siteXspp <- siteXspp[,c(1:5,2876,6:2875)]
head(siteXspp[,1:10])
write.csv(siteXspp,'siteXspp.csv',row.names=F)

#####Hidro analysis########
#####Open all Hidro categories#######
library(rgdal)
floodedforest <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/floodedforest/hdr.adf")
intermitent <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/intermitent/hdr.adf")
lake <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/lake/hdr.adf")
marsh <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/marsh/hdr.adf")
peatland <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/peatland/hdr.adf")
reservoir <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/reservoir/hdr.adf")
riveronly <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/riveronly/hdr.adf")
wet0_25 <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/wet0_25/hdr.adf")
wet25_50 <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/wet25_50/hdr.adf")
wet50_100 <-raster("~/Documents/Manuscripts/ProjetoMarilia/rastersArcGIS/wet50_100/hdr.adf")
#####Stack them#######
hidro <- stack(floodedforest,intermitent,lake,marsh,peatland,reservoir,riveronly,wet0_25,wet25_50,wet50_100)
names(hidro) <- c('floodedforest','intermitent','lake','marsh','peatland','reservoir','riveronly','wet0_25','wet25_50','wet50_100')
######Extract number of 1km cells within each 1 dgr cell for each hidro category###########
alt <- raster('~/Documents/Dados/Bioclimatic Variables/alt_30s_bil/alt.bil')
e <- extent(min(siteXspp$x),max(siteXspp$x),min(siteXspp$y),max(siteXspp$y))
alt.crop <- crop(alt,e)
alt.1dgr <- aggregate(alt.crop,120)
plot(alt.1gdr)
zone <- alt.1dgr
values(zone) <- 1:43596
zone <- disaggregate(zone,120)

hidro.crop <- crop(hidro,alt.crop)
hidro.crop2 <- resample(hidro.crop,alt.crop)

hidro.sum <- zonal(hidro.crop2, zone, fun = 'sum')
xy.zone <- xyFromCell(alt.1dgr,1:ncell(alt.1dgr))
xy.zone <- xy.zone + 0.5
hidro.table <- as.data.frame(cbind(xy.zone, hidro.sum))
head(hidro.table)
dim(hidro.table)

#####merging everything#####
head(siteXspp[,1:6])
head(tabela.cell[,c(1,2,7:27,33)])
head(hidro.table[,-c(3)])

table.cell.int <- merge(siteXspp[,1:6],tabela.cell[,c(1,2,7:27,33)])
table.cell.final <- merge(table.cell.int,hidro.table[,-c(3)])
head(table.cell.final)
write.csv(table.cell.final, '~/Dropbox/allamphibians IUCN_July/Final/table_cell_May14.csv',row.names=F)

#############Species database###########
head(subset(tabela.species, Species == 'Rhinella marina'))
head(tabela.species[,1:20])
variables= stack(alt.crop,hidro.crop2)

table.int <-as.data.frame(matrix(nrow=length(tabela.species[,1]),ncol=13))
colnames(table.int) <- c('alt_mean','alt_sd','floodedforest','intermitent','lake','marsh','peatland','reservoir','riveronly','wet0_25','wet25_50','wet50_100','dist_area')
head(table.int)

poly <- readOGR(dsn='mapas subset',layer='amph_phylo') ## ler shapefile
a <- unique(poly$BINOMIAL)

for (i in 1:length(tabela.species[,1])){
  maps <- poly[poly$BINOMIAL == a[i],]
  table.int[i,1] <- mean(extract(variables[[1]],maps,fun=mean,na.rm=T),na.rm=T)
  table.int[i,2] <- mean(extract(variables[[1]],maps,fun=sd,na.rm=T),na.rm=T)
  table.int[i,3] <- sum(extract(variables[[2]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,4] <- sum(extract(variables[[3]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,5] <- sum(extract(variables[[4]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,6] <- sum(extract(variables[[5]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,7] <- sum(extract(variables[[6]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,8] <- sum(extract(variables[[7]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,9] <- sum(extract(variables[[8]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,10] <- sum(extract(variables[[9]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,11] <- sum(extract(variables[[10]],maps,fun=sum,na.rm=T),na.rm=T)
  table.int[i,12] <- sum(extract(variables[[11]],maps,fun=sum,na.rm=T),na.rm=T)
  summary(maps)
  getClass("Polygon")
  ur.area<-sum(sapply(slot(maps, "polygons"), function(x) sapply(slot(x, "Polygons"), slot, "area")))
  table.int[i,13] <- ur.area
}
table.final <- cbind(tabela.species[,1],table.int)
colnames(table.final)[1] <- 'Species'
head(table.final)
write.csv(table.final, '~/Dropbox/allamphibians IUCN_July/Final/table_species_May14.csv',row.names=F)
