# The purpose of this code is to run the SSTTextract and create correlation plots for to resulting data

##### Step 0: add in libraries #####
pacman::p_load(devtools)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,elevatr,raster,devtools,rgdal,readxl,plot.matrix,grid,corrplot)
library('RColorBrewer')
#######

##### Step 1: Run function with all data of interest ####
#define common data
setwd("/groups/bse5304g/Elyce")
soilsreadin = read.csv("Soils/MoreSoils_LOC.csv") #soils data is not publically availible
outletcoords = matrix(c(-73.24896,44.1981),1,2) #LOC
#outletcoords = matrix(c(-73.34746, 44.08466),1,2) #DC -73.34647, 44.07354 for 1m
POIcoords = matrix(cbind(soilsreadin$Longitude[!is.na(soilsreadin$Latitude)],soilsreadin$Latitude[!is.na(soilsreadin$Latitude)]),length(which(!is.na(soilsreadin$Latitude))),2)

#1as USGS
rast = raster("DEMs/LOC_1as_projandtrim.tif")
D81as = SSTTextract(rast,outletcoords,POIcoords,"D8",10)
Dinf1as = SSTTextract(rast,outletcoords,POIcoords,"Dinf",10)

#1/3as USGS
rast = raster("DEMs/LOC_13as_projandtrim.tif")
D813as = SSTTextract(rast,outletcoords,POIcoords,"D8",10)
Dinf13as = SSTTextract(rast,outletcoords,POIcoords,"Dinf",10)

#GDEM
rast = raster("DEMs/LOC_GDEM_projandtrim.tif")
D8GDEM = SSTTextract(rast,outletcoords,POIcoords,"D8",10)
DinfGDEM = SSTTextract(rast,outletcoords,POIcoords,"Dinf",10)

#SRTM
rast = raster("DEMs/LOC_STRM_projandtrim.tif")
D8SRTM = SSTTextract(rast,outletcoords,POIcoords,"D8",10)
DinfSRTM = SSTTextract(rast,outletcoords,POIcoords,"Dinf",10)


#### These were saved locally using TuaDEM and exporting relevant rasters >
#1m USGS
D81m = D8SRTM; Dinf1m = D8SRTM
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
POI = SpatialPoints(POIcoords,proj4string = crs_ll)
POI = spTransform(POI,CRSobj = crs_utm)
setwd("/groups/bse5304g/test/1mD8/LOCSpecific")
sca = raster("scaclp.tif")
D81m$sca = raster::extract(sca,POI)
slp = raster("slpclp.tif")
D81m$slp = raster::extract(slp,POI)
TIC = raster("TICclp.tif")
D81m$TIC = raster::extract(TIC,POI)
TIV = raster("TIVclp.tif")
D81m$TIV = raster::extract(TIV,POI)

setwd("/groups/bse5304g/test/1mDinf/LOCSpecific")
sca = raster("scaclp.tif")
Dinf1m$sca = raster::extract(sca,POI)
slp = raster("slpclp.tif")
Dinf1m$slp = raster::extract(slp,POI)
TIC = raster("TICclp.tif")
Dinf1m$TIC = raster::extract(TIC,POI)
TIV = raster("TIVclp.tif")
Dinf1m$TIV = raster::extract(TIV,POI)


#LIDAR
D8LIDAR = D8SRTM; DinfLIDAR = D8SRTM
setwd("/groups/bse5304g/test/LIDARD8/LOCSpecific")
sca = raster("scaclp.tif")
D8LIDAR$sca = raster::extract(sca,POI)
slp = raster("slpclp.tif")
D8LIDAR$slp = raster::extract(slp,POI)
TIC = raster("TICclp.tif")
D8LIDAR$TIC = raster::extract(TIC,POI)
TIV = raster("TIVclp.tif")
D8LIDAR$TIV = raster::extract(TIV,POI)

setwd("/groups/bse5304g/test/LIDARDinf/LOCSpecfic")
sca = raster("scaclp.tif")
DinfLIDAR$sca = raster::extract(sca,POI)
slp = raster("slpclp.tif")
DinfLIDAR$slp = raster::extract(slp,POI)
TIC = raster("TICclp.tif")
DinfLIDAR$TIC = raster::extract(TIC,POI)
TIV = raster("TIVclp.tif")
DinfLIDAR$TIV = raster::extract(TIV,POI)


#### These were saved locally using TuaDEM and exporting relevant rasters <
########

##### Step 2: Create correlation matrix #####
soils = soilsreadin
cor4plot = data.frame(slp_LIDARdinf = matrix(NA,length(names(soils))-3,1),lnsca_LIDARdinf = NA,tiv_LIDARdinf = NA,tic_LIDARdinf = NA,
                         slp_LIDARd8 = NA,lnsca_LIDARd8 = NA,tiv_LIDARd8 = NA,tic_LIDARd8 = NA,
                         slp_1mdinf = matrix(NA,length(names(soils))-3,1),lnsca_1mdinf = NA,tiv_1mdinf = NA,tic_1mdinf = NA,
                         slp_1md8 = NA,lnsca_1md8 = NA,tiv_1md8 = NA,tic_1md8 = NA,
                         slp_1_3asdinf = NA,lnsca_1_3asdinf = NA,tiv_1_3asdinf = NA,tic_1_3asdinf = NA,
                         slp_1_3asd8 = NA,lnsca_1_3asd8 = NA,tiv_1_3asd8 = NA,tic_1_3asd8 = NA,
                         slp_1asdinf = NA,lnsca_1asdinf = NA,tiv_1asdinf = NA,tic_1asdinf = NA,
                         slp_1asd8 = NA,lnsca_1asd8 = NA,tiv_1asd8 = NA,tic_1asd8 = NA,
                         slp_SRTMdinf = NA,lnsca_SRTMdinf = NA,tiv_SRTMdinf = NA,tic_SRTMdinf = NA,
                         slp_SRTMd8 = NA,lnsca_SRTMd8 = NA,tiv_SRTMd8 = NA,tic_SRTMd8 = NA,
                         slp_GDEMdinf = NA,lnsca_GDEMdinf = NA,tiv_GDEMdinf = NA,tic_GDEMdinf = NA,
                         slp_GDEMd8 = NA,lnsca_GDEMd8 = NA,tiv_GDEMd8 = NA,tic_GDEMd8 = NA)
rownames(cor4plot) = names(soils)[4:length(names(soils))]
dfgets = c("DinfLIDAR","D8LIDAR","Dinf1m","D81m","Dinf13as","D813as","Dinf1as","D81as","DinfSRTM","D8SRTM","DinfGDEM","D8GDEM")

for(i in 1:length(rownames(cor4plot))){
  for(j in 1:length(dfgets)){
    df = get(dfgets[j])
    cor4plot[i,(j-1)*4+1] = cor(soils[,i+3],df$slp,use = "complete.obs")
    cor4plot[i,(j-1)*4+2] = cor(soils[,i+3],df$sca,use = "complete.obs")
    cor4plot[i,(j-1)*4+3] = cor(soils[,i+3],df$TIV,use = "complete.obs")
    cor4plot[i,(j-1)*4+4] = cor(soils[,i+3],df$TIC,use = "complete.obs")
  }
}


#visualize matrix
cor4plot[,which(is.na(cor4plot[1,]))] = 0
library('RColorBrewer')
corrplot(as.matrix(cor4plot),tl.col = "black",mar=c(0,1,1,1),tl.cex=.75,main = "Little Otter Creek Correlations", col = c(brewer.pal(5,"Reds")[5:2],"white",brewer.pal(5,"Blues")[2:5]),cl.length=3)

