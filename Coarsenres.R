#The purpose of this code is understand how the correlations of
# sand, clay, awc and OM, relate to DEM resolution


##### Step 0: add in libraries #####
pacman::p_load(devtools)
if (!require("pacman")) install.packages("pacman")
pacman::p_load(rgdal,elevatr,raster,devtools,rgdal,readxl,plot.matrix,grid,corrplot)
#######

#### Step 1: Load in soil data ####
setwd("/groups/bse5304g/Elyce")
soils = read.csv("RelevantSoils.csv") #soils data not publically availible
#####

#### Step 2: Create dataframe ####
coarsestresslp = 715 #aggregation steps
newstepres = 45 #aggregation steps
coarsestressca = 43 #aggregation steps
scaby = 1 #aggregation steps
slpby = 5 #aggregation steps
scad8 = data.frame(matrix(NA,length(soils$SampleNumber),length(seq(1,coarsestressca,scaby))))
names(scad8) = seq(1,coarsestressca,scaby); scadinf = scad8
slpd8 = data.frame(matrix(NA,length(soils$SampleNumber),length(c(seq(1,coarsestressca,scaby),seq(newstepres,coarsestresslp,slpby)))))
names(slpd8) = c(seq(1,coarsestressca,scaby),seq(newstepres,coarsestresslp,slpby)); slpdinf = slpd8
#######

#### Step 3: Load in DEM ####
demlidar = raster("DEMs/LIDAR_projandtrim.tif")
####

#### Step 4: Calculate cor for 0.7m ####
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
POIcoords = matrix(cbind(soils$Longitude[!is.na(soils$Latitude)],soils$Latitude[!is.na(soils$Latitude)]),length(which(!is.na(soils$Latitude))),2)
POI = SpatialPoints(POIcoords,proj4string = crs_ll)
POI = spTransform(POI,CRSobj = crs_utm)
slpscad8 = data.frame(slp = matrix(NA,length(soils$SampleNumber),1),sca=NA); slpscadinf=slpscad8

#LIDAR
setwd("/groups/bse5304g/test")
slp = raster("LIDARD8/slp.tif")
slpd8$`1` = raster::extract(slp,POI)
sca = raster("LIDARD8/sca.tif")
scad8$`1` = raster::extract(sca,POI)

slp = raster("LIDARDinf/slp.tif")
slpdinf$`1` = raster::extract(slp,POI)
sca = raster("LIDARDinf/sca.tif")
scadinf$`1` = raster::extract(sca,POI)
#####

#### Step 5a: Calculate scas and slps for the remaining resolutions ####
for(i in 21:length(names(scad8))){
  dem = aggregate(demlidar,as.numeric(names(scad8)[i]),fun = 'mean')
  slpscadinf = SSextract(dem,POI,"Dinf")
  slpscad8 = SSextract(dem,POI,"D8")
  scad8[,i] = slpscad8$sca
  scadinf[,i] = slpscadinf$sca
  slpd8[,i] = slpscad8$slp
  slpdinf[,i] = slpscadinf$slp
  print(i)
  save.image(paste0("/groups/bse5304g/Elyce/WorkspaceDump/LIDAR/sca/demcorseningwksp",i,".RData"))
}
#######

#### Step 5a: Calculate scas and slps for the remaining resolutions ####
for(i in (length(names(scad8))+1):length(names(slpd8))){
  dem = aggregate(demlidar,as.numeric(names(slpd8)[i]),fun = 'mean')
  slpscadinf = SSextract(dem,POI,"Dinf")
  slpscad8 = SSextract(dem,POI,"D8")
  slpd8[,i] = slpscad8$slp
  slpdinf[,i] = slpscadinf$slp
  print(i)
  save.image(paste0("/groups/bse5304g/Elyce/WorkspaceDump/LIDAR/slp/demcorseningwksp",i,".RData"))
}
#######
