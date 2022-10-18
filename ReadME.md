Integrating multiple Digital Elevation Models into soil characteristic distribution [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6857407.svg)](https://doi.org/10.5281/zenodo.6857407)
=================

This repository contains all the data, spatial files, and code used to execute the analysis in Integrating multiple Digital Elevation Models into 
soil characteristic distribution. Contributing authors: Elyce Buell, Roja Kaveh Garna, Sabrina Mehzabin, Binyam Asfaw, Louise Koepele, Daniel R. Fuka, 
Amy S. Collick, William Auchincloss, and Zachary M. Eastona


If you have any questions regarding this publication please contact Elyce Buell (<enb46@cornell.edu> or <ebuell@vt.edu>).

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

### Fig8to10
Figure generation for figures 8-10. Must run runDelin_taudem.R and save rasters locally in order for figure generation to work

	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/Fig8to10.R","Fig8to10.R")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/PhysicalProperties_datainpaper.xlsx","PhysicalProperties_datainpaper.xlsx")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/DEMderivedSpatialAtt.xlsx","DEMderivedSpatialAtt.xlsx")
	file.edit("Fig8to10.R")

### Figsupp
Supplemental figure generation
	
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/Figsupp.R","Figsupp.R")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/1arcsec_arcproj_bilin.tif","1arcsec_arcproj_bilin.tif")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/1_3arcsec_arcproj_bilin.tif","1_3arcsec_arcproj_bilin.tif")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/2010LIDAR_arcproj_bilin.tif","2010LIDAR_arcproj_bilin.tif")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/2018LIDAR_arcproj_bilin.tif","2018LIDAR_arcproj_bilin.tif")
	download.file("https://raw.githubusercontent.com/ebuell/Evaluating-Digital-Elevation-Models-and-Topographic-Indices-for-Geomorphic-Landscape-Representaton/main/MonitoringPoint.shp","MonitoringPoint.shp")
	file.edit("Figsupp.R")


## Excels (excel and sheet names) summary:

### PhysicalProperties_datainpaper.xlsx
Summarizes physical properties of soil cores taken at 36 locations in southwest virginia

#### PhysicalResults	
Final texture, organic matter, and horizon depth information found here. All categories included in this excel are used in the paper analysis. 
Details on direct measurements and QAQC can be found in KSatMeas, OrgMatter, and SoilTexture.

#### WSView 		      
Locations of soil cores (numbers correspond to WSView ID)

### PhysicalProperties_alldata.xlsx

#### PhysicalResults
Summarizes physical properties of soil cores taken at 36 locations in southwest virginia. Final texture, organic matter, and horizon depth information found here. Not all categories included in this excel are used in the paper analysis. Details on direct measurements and QAQC can be found in KSatMeas, OrgMatter, and SoilTexture.
	
#### WSView
Locations of soil cores (numbers correspond to WSView ID)

### OrganicMatter.xlsx

#### OM
Summarizes OM for all horizons. Loss on ignition tests are conducted on each horizon for all soil cores. Soil samples ranged from 50g to 250g. Samples were dried for 2 days in a Fisher Scientific Isotherm oven set to 105 degC and are subsequently weighed. These samples are incinerated and before and after weights are compared to yield percent organic matter by weight. The samples are incinerated at 425 degC. The length of incineration was dictated by the incinerator used. For the Barnstead International Model 1730 12 sample incinerator, samples are incinerated for a minimum of 6 hours. For the Thermolyne Furnace Type 10500 1 sample incinerator, samples are incinerated for a minimum of 4 hours. A 4 hour incineration in the Barnstead International Model 1730 was shown to have magnitude errors of up to 1% organic matter and thus 2 additional hours were added on to the incineration time.

#### WSView 
Locations of soil cores (numbers correspond to WSView ID)
	
#### OrgMatter
Direct weight measurements and calculations for organic matter

#### QAQC_NewSoil
New soil subsamples from the same core/horizon are analyzed for organic matter

#### QAQC_NewSoilAnalysis
Analysis to see if the errors seen in these QAQC runs is tied to horizon thickness or the amount of time that the soil has been left unfrozen in between runs
	
#### QAQC_SampleRerun
Already incinerated samples were rerun in order to check that the current methodology was suffcient accurate within reason.

### SoilTexture.xlxs

#### TextureFinal
Final texture (sand and clay) values decided by R code. For more information on texture descision making for this calculation is in the R code

#### WSView
Locations of soil cores (numbers correspond to WSView ID)

#### SoilGravel
Measured gravel content using a 2mm sieve

#### TextureSummary
Relevant measurement values for calculating soil texture

#### QAQCTextureSummary
Summarizes QAQC and orgininal runs; this includes two types of QAQC: reruns (i.e. shake the soil column again with same soil) and new soil from soil column

### HydrometerMeas.xlsx
Direct hydrometer measurements following the method described in ASTM D4222. In order to directly calculate the hydrometer corrects (see ASTM D4222), several corrections need to be measured. This ultimately did not pan out and texture was calculated using envalysis package.
	
#### TextureSummary
Relevant measurement values for calculating soil texture

#### QAQCTextureSummary
Relevant measurement values for calculating QAQC soil texture

#### WSView
Locations of soil cores (numbers correspond to WSView ID)	

#### SoilTextureA
All measurements for hydrometer trials for A horizon

#### SoilTextureBa
All measurements for hydrometer trials for Ba horizon
	
#### SoilTextureBt
All measurements for hydrometer trials for Bt horizon
	
#### TempCorrMeas
Measured hydrometer correction values in distilled water at a variety of temperatures. A very poor relationship between temperature and hydrometer correction is found implying the poor reliability of the thermometers being used in the hydrometer method.
	
#### CorrectionSummary
Direct measurement corrections.

### KSatMeas.xlsx

#### KSAT
Ksat values and lats and longs

#### WSView
Locations of soil cores (numbers correspond to WSView ID)

#### SatHydCondTestx
Sheets with direct field readings for the Eijkelkamp Soil & Water double-ring infiltrometer test following (Bouwer, 1986;  DOI: 10.2136/sssabookser5.1.2ed.c32)

### DEMderivedSpatialAtt.xlsx
All numbers in this document are derived using TauDEM: available at: <https://github.com/dtarb/TauDEM>. Both D8 and Dinf flow direction algorithms (discussed in Tarboton DG. 1997; DOI: 10.1029/96WR03137) are used as well as the four DEMs (LIDAR 2018, LIDAR 2010, 1/3as, and 1as).
	
#### SCA
Speicfic catchment area (expressed in m)

#### Slope
Slope (expressed in length/length)

#### TIValue
TIV = ln(SCA/slope)
	
#### TIClass
TIC bins the TIVs into 10 equally sized categories ranging from 1-10 (where the highest TIVs correspond to TIC = 10 etc.)

#### WSView
Locations of soil cores (numbers correspond to WSView ID)

#### Aspect
Aspects range from 1-8. Refer to TauDEM documentation to understand these values.

### MeasSpatProp.xlsx
#### SpatialProp
Measured spatial properties (slope and aspect); a phone app clinometer (<https://play.google.com/store/apps/details?id=com.plaincode.clinometer&hl=en_US&gl=US>) and phone compass are used to measure these values. The measurement was taken with the assistance of a 4-wheel vehicle; the vehicle was oriented perpendicularly to the slope and the slope that the utility vehicle is experiancing is recorded. The distance between the two front wheels of the utility vehicle is 1.2m. The aspect is measured by pointing the compass perpendicularly and downslope of how the utility vehicle is oriented for the slope measurement.

#### WSView 
Locations of soil cores (numbers correspond to WSView ID)

## Spatial data files

### Construction_exclusion
A box that surround highway 460 throughout the entire watershed so analysis can exclude this area

### Monitoring_point
Outlet location

### pathsunderoverpass
A path that is created to mimic storm drainage culverts that run below highway 460. Only used for LIDAR analysis; USGS DEM did not have any highway interfereance.

### SoilSampleLocations
Soil sampling locations

### 1as_arcproj_bilin
1 as DEM bilinear projection executed in ArcMap

### 1_3as_arcproj_bilin
1/3 as DEM bilinear projection executed in ArcMap

### 2010LIDAR_arcproj_bilin
2010LIDAR DEM bilinear projection executed in ArcMap

### 2018LIDAR_arcproj_bilin
2018LIDAR DEM bilinear projection executed in ArcMap; dataset too big to upload to github, please find via google drive here: <https://drive.google.com/file/d/1fc_7hbiMM8hfNFi38inOW2mjg_jeq-K6/view?usp=sharing>

# License
Please see the LICENSE.md file for license information.
