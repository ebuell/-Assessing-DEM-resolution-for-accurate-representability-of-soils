Integrating multiple Digital Elevation Models into soil characteristic distribution [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6857407.svg)](https://doi.org/10.5281/zenodo.6857407)
=================

This repository contains all the data, spatial files, and code used to execute the analysis in Integrating multiple Digital Elevation Models into 
soil characteristic distribution. Contributing authors: Elyce Buell, Roja Kaveh Garna, Sabrina Mehzabin, Binyam Asfaw, Louise Koepele, Daniel R. Fuka, 
Amy S. Collick, William Auchincloss, and Zachary M. Easton


If you have any questions regarding this publication please contact Elyce Buell (<enb46@cornell.edu> or <ebuell@vt.edu>) or Roja Kaveh Garna (rojakaveh@vt.edu).

## Links
See the following links for more information on  `R` and `RStudio` download and installation:

- An introduction to `R`: <https://cran.r-project.org/doc/manuals/r-release/R-intro.pdf>
- `R` download: <https://www.r-project.org/>
- `RStudio` download: <https://www.rstudio.com/>

There is also a cloud-based `RStudio` sever at the following location:

- Cloud-based `RStudio` server: <https://rstudio.cloud/>

See the following likes for information regarding `DEM` data

- GDEM DEM download <https://lpdaac.usgs.gov/products/astgtmv003/>
- SRTM DEM download <https://ui.adsabs.harvard.edu/abs/2000EOSTr..81..583F%2F/abstract>
- USGS DEM download <https://apps.nationalmap.gov/downloader/#/>
- Aerial LiDAR DEM download <https://maps.vcgi.vermont.gov/LidarFinder/>

## Description
This repository contains R codes, excel files and DEM files for a project named *Integrating multiple Digital Elevation Models into 
soil characteristic distribution*. This study is conducted to propose a new method of distributing soils using multiple DEMs as inputs. Fifty-nine soil samples are collected by the Vermont Association of Conservative Districts and analyzes for a broad range of soil characteristics including sand & clay content, organic matter, and predicted water capacity in western Vermont. Six DEMs (GDEM, SRTM, USGS 1m, ⅓ and 1as, and aerial LiDAR) are analyzed for spatial differences between derived properties (slope, Specific Catchment Area (SCA), and TI). Using multivariate regression, these soil properties are predicted, and a framework for soil map distribution is proposed.

 ## Quick start

### R packages that need to be installed:
•   raster
•   rgdal
•   shapefiles
•   RColorBrewer
•   corrplot
•   rgeos
•   elevatr
•   devtools
•   plot.matrix
•   grid
•   httr
•   parallel
•   EcoHydRology
•   GSODR
•   curl
•   sp
•   sf
•   reticulate
•   classInt
•   ggridges
•   ggplot2
•   viridis
•   hrbrthemes
•   dplyr
•   tidyr
•   forcats
•   gridExtra


	if (!require("pacman")) install.packages("pacman")
     pacman::p_load(raster,rgdal,shapefiles,RColorBrewer,corrplot,rgeos,elevatr,devtools,plot.matrix,grid,httr,parallel,EcoHydRology,GSODR,curl,sp,sf,reticulate,classInt,ggridges,ggplot2,viridis,hrbrthemes,dplyr,tidyr,forcats,gridExtra)
	
### Additional R functionality required: TauDEM
Download and install per instuctions found here <https://hydrology.usu.edu/taudem/taudem5/>

## Running R code
Please note that soils data was provided by the Vermont Association of Conservative Districts (USDA NRCS, 2019). The have requested that soils data not be published 
publically with this code. Please reach out the the contacts listed above for more information. Because this information has been omitted, published codes cannot be run in entirity.

### DelinandExtract.R
Function that delinates full watershed using TauDEM: available at: https://github.com/dtarb/TauDEM. Both D8 and Dinf flow direction algorithms are options for delineation and extration (discussed in Tarboton DG. 1997; DOI: 10.1029/96WR03137). This is a function file that will delineate a watershed and extract spatial values of interest. The purpose of this is to be able to simplify our analysis by not having to deal with any pushing and pulling of large spatial rasters.

Function inputs:
-DEM: DEM shall be inputted as projected into UTM and clipped as needed
-Outlet coordinates: a 2x1 matrix
-Coordinates of interest: a 2xn matrix (for this paper these would be the soils locations that we are interested in)
-An indication of either Dinf or D8 for the delineation (either "Dinf" or "D8")
-Number of TICs desired
				
Function outputs
-Output extracted data
  Output will be a 6xn matrix (where n is the number of points of interest there are)
  Output matrix columns: Lat  Long  Slp SCA TIV TIC

Example: output = SSTTExtract(DEM,outletcoords,POI,"D8",10)
Please note that code asks you the anticipated size of watershed to see if you agree with the delineation. If you do not agree, the outlet is moved to the maxima raster SCA for all the rasters surrounding the original raster. This process occurs recursively until the user agrees with the rough size of the watershed. Please note that you will also need to download a DEM of choice to run this function
Please note the code below:
	
	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/DelinandExtract.R","DelinandExtract.R")
	file.edit("DelinandExtract.R")

### ExtractSlpSCA.R
Function does the processes described in DelinandExtract.R but does not run delineation. The purpose of this code is to extract slope and SCA for loactions of interest.

Function inputs
-DEM: DEM shall be inputted as projected into UTM and clipped as needed
-Coordinates of interest: a 2xn matrix (for this paper these would be the soils locations that we are interested in)
-An indication of either Dinf or D8 for the delineation

Function outputs
-Output extracted data
 Output will be a 4xn matrix (where n is the number of points of interest there are)
 Output matrix columns: Lat  Long  Slp SCA
				
Example: output = SSExtract(DEM,POI,"D8")
Please note that you will also need to download a DEM of choice to run this function.
Please note the code below:

	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/ExtractSlpSCA.R","ExtractSlpSCA.R")
	file.edit("ExtractSlpSCA.R")

### corrplotfiggenTI.R
The purpose of this code is to run the SSTTextract and create correlation plots for to resulting data. Must run DelinandExtract.R before running this code. 

	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/DelinandExtract.R","DelinandExtract.R")
	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/corrplotfiggenTI.R","corrplotfiggenTI.R")
	file.edit("DelinandExtract.R")
	file.edit("corrplotfiggenTI.R")

### corrplotfiggen.R		
The purpose of this code is to run the SSextract and create correlation plots for to resulting data. This code also produces the multiple regression summary figure. Must run ExtractSlpSCA.R before running this code.

	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/DelinandExtract.R","ExtractSlpSCA.R")
	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/corrplotfiggenTI.R","corrplotfiggen.R")
	file.edit("ExtractSlpSCA.R")
	file.edit("corrplotfiggen.R")

### soilsdistandcomparessurgo.R
This code generates the figures the compare distributed soil maps (using the multiple regressions resulting from corrplotfiggen) to SSURGO distributed soils. This code also compares (via bootstrapping) the informed multiple regressions to using the finest resolution for soil distribution.

	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/soildistandcomparessurgo.R","soildistandcomparessurgo.R")
	file.edit("soildistandcomparessurgo.R")

### Coarsenres.R
This code recursively coarsens the LiDAR DEM and extracts slopes and SCAs from each coarsened DEM. Currently this is setup for LiDAR only, though both 1m and 1/3as have also been coarsened for the analysis found in the corresponding paper. The workspace is saved after each coarsening because this process was time intense and ran unsupervised for several days. Workspace saving was done in order to avoid any data losses due to crashed code.
	
	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/blob/main/Coarsenres.R","Coarsenres.R")
	file.edit("Coarsenres.R")
	
### bestres.R
Data extracted from recursively coarsened DEMs generated from Coarsenres and is saved locally. This code generates the figures seen in the paper.
	
	download.file("https://github.com/ebuell/-Assessing-DEM-resolution-for-accurate-representability-of-soils/edit/main/bestres.R","bestres.R")
	file.edit("bestres.R")


## Excel summary:

### Soils Data
Soils data was provided by the Vermont Association of Conservative Districts (USDA NRCS, 2019). They have requested that soils data not be published 
publically with this code. Please reach out the the contacts listed above for more information. 
The general format of the relevant .csv file is (by column):
Latitude	
Longitude	
Sand content (% by mass)	
Clay content (% by mass)	
Predicted awc (cm/cm)	
Organic matter content (% by mass)

## Spatial data files

### 1as_projandtrim.tif
DEM raster for 1as sourced from USGS https://apps.nationalmap.gov/downloader/#/

### GDEM_projandtrim.tif
DEM raster for GDEM sourced from https://doi.org/10.5067/ASTER/ASTGTM.003

### STRM_projandtrim.tif
DEM raster for SRTM sourced from https://doi.org/10.1029/EO081i048p00583

### Remaining DEMs
The remaining DEMs used in this project are too large to include in this publication (per GitHub requirements).
LiDAR data can be found here: https://maps.vcgi.vermont.gov/LidarFinder/;
1m, and 1/3as can be found here: https://apps.nationalmap.gov/downloader/#/

Please note that a great deal of raster merging and clipping was done for the LiDAR data to make the spatial data small enough for the computer (64 cores)
to handle. This code is not published but can be shared upon request. Because of the size of many of the rasters, many of the resulting rasters are saved
locally once raster manipulation is completed. As a result, the code is not stand alone because of the intensity of the computation needed to make the 
code stand alone.

# License
Please see the LICENSE.md file for license information.
