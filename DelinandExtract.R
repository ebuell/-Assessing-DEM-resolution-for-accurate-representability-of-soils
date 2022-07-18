# This is a function file that will delineate a watershed and extract spatial values of interest
#   The purpose of this is to be able to simplify our analysis by not having to deal with any pushing
#   and pulling of large spatial rasters

#Function inputs
# -DEM: DEM shall be inputted as projected into UTM and clipped as needed
# -Outlet coordinates: a 2x1 matrix
# -Coordinates of interest: a 2xn matrix (for this paper these would be the soils locations that we are interested in)
# -An indication of either Dinf or D8 for the delineation
# -Number of TICs desired

#Function steps:
# (1) Delineate watershed
# (2) Create rasters for spatial data
# (3) Extract spatial data from rasters at coordinates of interest

#Function outputs
# -Output extracted data
#   Output will be a 6xn matrix (where n is the number of points of interest there are)
#   Output matrix columns: Lat  Long  Slp SCA TIV TIC

#Things to automate for a real user
# MPICH and TauDEM??
# Get UTM zone of their dem

SSTTextract = function(dem_utm,outletcoords,POIcoods,d8ordinf,numofTICs){

### 0b: Load in necessary libraries ###
if (!require("pacman")) install.packages("pacman")
pacman::p_load(httr,parallel,EcoHydRology,GSODR,curl,elevatr,raster,rgdal,shapefiles,sp,rgeos,reticulate,sf,classInt)
######

#### Step 1 #####
### 1a: set working directory ###
#save original working directory
wd = getwd()
setwd("/groups/bse5304g/test/deleterepo")
ncores = detectCores()
file.remove(list.files())

### 1b: Create spatial points for outlets and POIs ###
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
outlet = SpatialPoints(outletcoords,proj4string = crs_ll)
outlet = spTransform(outlet,CRSobj = crs_utm)
POI = SpatialPoints(POIcoords,proj4string = crs_ll)
POI = spTransform(POI,CRSobj = crs_utm)

### 1c: remove pits from DEM ###
writeRaster(dem_utm,"logan",format = "GTiff",overwrite = TRUE)
system("mpiexec -n 8 pitremove -z logan.tif -fel loganfel.tif")
dem_utm_pitrm=raster("loganfel.tif")

### 1d: slope and contributing area calculation ###
if(d8ordinf=="D8"){
  #D8 flow directions
  system(paste0("mpiexec -n ",ncores," d8flowdir -p loganp.tif -sd8 logans.tif -fel loganfel.tif"),show.output.on.console=F,invisible=F)
  p=raster("loganp.tif")
  slp=raster("logans.tif")
  # Contributing area
  system(paste0("mpiexec -n ",ncores," aread8 -p loganp.tif -ad8 logana.tif"))
  sca=raster("logana.tif")
}
if(d8ordinf=="Dinf"){
  # DInf flow directions
  system(paste0("mpiexec -n ",ncores," dinfflowdir -ang loganpdinf.tif -slp logansdinf.tif -fel loganfel.tif"),show.output.on.console=F,invisible=F)
  p=raster("loganpdinf.tif")
  slp=raster("logansdinf.tif")
  # Contributing area
  system(paste0("mpiexec -n ",ncores," areadinf -ang loganpdinf.tif -sca loganadinf.tif"))
  sca=raster("loganadinf.tif")
  system(paste0("mpiexec -n ",ncores," d8flowdir -p loganp.tif -sd8 logans.tif -fel loganfel.tif"),show.output.on.console=F,invisible=F)
  p=raster("loganp.tif")
}

### 1e: delineate watershed ###
shapefile(outlet, "approxoutlet.shp",overwrite = TRUE)
system(paste0("mpiexec -n ",ncores," aread8 -p loganp.tif -o approxoutlet.shp -ad8 loganssa.tif"))
ssa=raster("loganssa.tif")
yn = readline(paste("Outlet given yeilds a",round(length(values(ssa)[!is.na(values(ssa))])*res(dem_utm)[1]^2),"units^2 watershed.\n If you agree with this value enter y, else enter n"))
######

### 1f: if outlet point does not create a large enough watershed ###
n = 1
while(!(yn=="y"|yn=="Y")){
###change outlet to a larger sca ###
buf = buffer(outlet,n*res(dem_utm)[1])
multthresh = mask(sca,buf); 
multthresh = trim(multthresh,values=NA)
pnts = SpatialPoints(multthresh)
pnts = SpatialPointsDataFrame(pnts@coords,data.frame(SCA = (values(multthresh))))
outlet2 =SpatialPoints(matrix(pnts@coords[which.max(pnts$SCA),],1,2))

###delineate watershed ###
shapefile(outlet2, "approxoutlet.shp",overwrite = TRUE)
system(paste0("mpiexec -n ",ncores," aread8 -p loganp.tif -o approxoutlet.shp -ad8 loganssa.tif"))
ssa=raster("loganssa.tif")
yn = readline(paste("Outlet given yeilds a",round(length(values(ssa)[!is.na(values(ssa))])*res(dem_utm)[1]^2),"units^2 watershed.\n If you agree with this value enter y, else enter n"))
n = n+1
}
########


#### Step 2 #####
### 2a: calculate TIV and clip to watershed boundary ###
sca = mask(sca,ssa)
#sca = trim(sca,values=NA)
slp = mask(slp,ssa)
#slp = trim(slp,values=NA)
TIV = log((sca+1)/(slp+0.00001))

### 2c: calculate TIC ###
brks.qt = classIntervals(values(TIV), n = numofTICs, style = "quantile")$brks 
TIC = cut(TIV, breaks=brks.qt, include.lowest = T, right=T)
#####


#### Step 3 #####
### 3a: create data frame to pass back through function ###
extractedspatial = data.frame(slp = matrix(NA,length(POIcoords[,1]),1),sca=NA,TIV=NA,TIC=NA)

### 3b: extract spatial data ###
extractedspatial$slp = raster::extract(slp,POI)
extractedspatial$sca = raster::extract(sca,POI)
extractedspatial$TIV = raster::extract(TIV,POI)
extractedspatial$TIC = raster::extract(TIC,POI)

### 3c: set the wd back to its orignal wd ###
setwd(wd)

### 3d: pass data frame back to user ###
return(extractedspatial)

}
#####