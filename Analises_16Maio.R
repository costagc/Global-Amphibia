setwd('~/Documents/Manuscripts/ProjetoMarilia')
library(raster)
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

bio <- stack(list.files('~/Documents/Dados/Bioclimatic Variables/bioclim',pattern='.bil',full.names=T))
bio1dgr <- aggregate(bio,6)
alt <- raster('~/Documents/Dados/Bioclimatic Variables/alt_30s_bil/alt.bil')
alt.1dgr <- aggregate(alt,120)

hidro <- stack(floodedforest,intermitent,lake,marsh,peatland,reservoir,riveronly,wet0_25,wet25_50,wet50_100)
names(hidro) <- c('floodedforest','intermitent','lake','marsh','peatland','reservoir','riveronly','wet0_25','wet25_50','wet50_100')
rm(floodedforest,intermitent,lake,marsh,peatland,reservoir,riveronly,wet0_25,wet25_50,wet50_100)

hidro1dgr <- aggregate(hidro,fact=120,fun=sum)
hidro1dgr.crop <- crop(hidro1dgr,alt.1dgr)
hidro1dgr.crop.res <- resample(hidro1dgr.crop,alt.1dgr)
cafs = raster('~/Documents/Manuscripts/ProjetoMarilia/CAFS/cafs.grd')
variables <- stack(alt.1dgr,bio1dgr,hidro1dgr.crop.res,cafs)

tabela <- read.csv('~/Dropbox/allamphibians IUCN_July/Final/tabela_species_final.csv')

poligonos <- readOGR(dsn='/Users/gabriel/Documents/Dados/Dados Distribuição/Distribuição IUCN/AMPHANURA',layer='Amphibians_Anura')
##especies <- sort(unique(poligonos@data$BINOMIAL))###
especies <- as.vector(tabela[,2])

table.int <-as.data.frame(matrix(nrow=length(especies),ncol=33))
colnames(table.int) <- c('alt_mean','alt_sd',names(bio1dgr),names(hidro1dgr.crop.res),'cafs','area')
head(table.int)

require(doSNOW)
cl<-makeCluster(5,"SOCK")
registerDoSNOW(cl)
for (i in 1:5593){
  maps <- poligonos[poligonos$BINOMIAL == especies[i],]
  table.int[i,1] <- mean(extract(variables[[1]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,2] <- mean(extract(variables[[1]],maps,fun=sd,na.rm=T,small=T),na.rm=T)
  table.int[i,3] <- mean(extract(variables[[2]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,4] <- mean(extract(variables[[3]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,5] <- mean(extract(variables[[4]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,6] <- mean(extract(variables[[5]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,7] <- mean(extract(variables[[6]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,8] <- mean(extract(variables[[7]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,9] <- mean(extract(variables[[8]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,10] <- mean(extract(variables[[9]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,11] <- mean(extract(variables[[10]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,12] <- mean(extract(variables[[11]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,13] <- mean(extract(variables[[12]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,14] <- mean(extract(variables[[13]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,15] <- mean(extract(variables[[14]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,16] <- mean(extract(variables[[15]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,17] <- mean(extract(variables[[16]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,18] <- mean(extract(variables[[17]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,19] <- mean(extract(variables[[18]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,20] <- mean(extract(variables[[19]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,21] <- mean(extract(variables[[20]],maps,fun=mean,na.rm=T,small=T),na.rm=T)
  table.int[i,22] <- sum(extract(variables[[21]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,23] <- sum(extract(variables[[22]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,24] <- sum(extract(variables[[23]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,25] <- sum(extract(variables[[24]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,26] <- sum(extract(variables[[25]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,27] <- sum(extract(variables[[26]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,28] <- sum(extract(variables[[27]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,29] <- sum(extract(variables[[28]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,30] <- sum(extract(variables[[29]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,31] <- sum(extract(variables[[30]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,32] <- sum(extract(variables[[31]],maps,fun=sum,na.rm=T,small=T),na.rm=T)
  table.int[i,33] <- gArea(maps)
}
stopCluster(cl)
table.final <- cbind(especies,table.int)
write.csv(table.final, '~/Dropbox/allamphibians IUCN_July/Final/table_species_17May14.csv',row.names=F)

require(rgeos)
x <- matrix(nrow=5593,ncol=1)
cl<-makeCluster(5,"SOCK")
registerDoSNOW(cl)
for (i in 1:5593){
  maps <- poligonos[poligonos$BINOMIAL == especies[i],]
  area <- gArea(maps)
  x[i,1]<-area
}
stopCluster(cl)

