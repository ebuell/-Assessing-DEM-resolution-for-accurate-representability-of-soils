# This is a function file is to extract spatial values: slope and SCA

#Function inputs
# -DEM: DEM shall be inputted as projected into UTM and clipped as needed
# -Coordinates of interest: a 2xn matrix (for this paper these would be the soils locations that we are interested in)
# -An indication of either Dinf or D8 for the delineation

#Function steps:
# (1) Create rasters for spatial data
# (2) Extract spatial data from rasters at coordinates of interest

#Function outputs
# -Output extracted data
#   Output will be a 4xn matrix (where n is the number of points of interest there are)
#   Output matrix columns: Lat  Long  Slp SCA

#Things to automate for a real user
# MPICH and TauDEM??
# Get UTM zone of their dem

SSextract = function(dem_utm,POI,d8ordinf){

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

### 1c: remove pits from DEM ###
writeRaster(dem_utm,"logan",format = "GTiff",overwrite = TRUE)
system("mpiexec -n 8 pitremove -z logan.tif -fel loganfel.tif")
dem_utm_pitrm=raster("loganfel.tif")

### 1d: slope and contributing area calculation ###
if(d8ordinf=="D8"){
  #D8 flow directions
  system(paste0("mpiexec -n ",ncores," d8flowdir -p loganp.tif -sd8 logans.tif -fel loganfel.tif"),show.output.on.console=F,invisible=F)
  slp=raster("logans.tif")
  # Contributing area
  system(paste0("mpiexec -n ",ncores," aread8 -p loganp.tif -ad8 logana.tif"))
  sca=raster("logana.tif")
}
if(d8ordinf=="Dinf"){
  # DInf flow directions
  system(paste0("mpiexec -n ",ncores," dinfflowdir -ang loganpdinf.tif -slp logansdinf.tif -fel loganfel.tif"),show.output.on.console=F,invisible=F)
  slp=raster("logansdinf.tif")
  # Contributing area
  system(paste0("mpiexec -n ",ncores," areadinf -ang loganpdinf.tif -sca loganadinf.tif"))
  sca=raster("loganadinf.tif")
}




#### Step 2 #####
### 2a: create data frame to pass back through function ###
extractedspatial = data.frame(slp = matrix(NA,length(POI@coords[,1]),1),sca=NA)

### 2b: extract spatial data ###
extractedspatial$slp = raster::extract(slp,POI)
extractedspatial$sca = raster::extract(sca,POI)

### 3c: set the wd back to its orignal wd ###
setwd(wd)

### 3d: pass data frame back to user ###
return(extractedspatial)

}
#####