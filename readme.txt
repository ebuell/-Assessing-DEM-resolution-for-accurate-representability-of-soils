If you have any questions regarding this publication please contact 
Elyce Buell (enb46@cornell.edu or ebuell@vt.edu) or Roja Kaveh Garna (rojakaveh@vt.edu)

Soils data:
--------------------------------------------------------------------------------------
Soils data was provided by the Vermont Association of Conservative Districts (USDA NRCS, 2019. USDA Announces Water Quality Study in Vermont’s Lake Champlain Basin | NRCS [WWW Document]. URL https://www.nrcs.usda.gov/wps/portal/nrcs/detail/national/newsroom/releases/?cid=NRCSEPRD1495437 (accessed 4.8.22)).
They have requested that soils data not be published publically with this code. Please reach out the the contacts listed above for more information. 
The general format of the relevant .csv file is (by column):
Latitude	
Longitude	
Sand content (% by mass)	
Clay content (% by mass)	
Predicted awc (cm/cm)	
Organic matter content (% by mass)
--------------------------------------------------------------------------------------


R Code
---------------------------------------------------------------------------------------
DelinandExtract		-	Function that delinates full watershed using TauDEM: available at: https://github.com/dtarb/TauDEM. Both D8 and 
				Dinf flow direction algorithms are options for delineation and extration 
				(discussed in Tarboton DG. 1997; DOI: 10.1029/96WR03137). 
				This is a function file that will delineate a watershed and extract spatial values of interest
				   The purpose of this is to be able to simplify our analysis by not having to deal with any pushing
				   and pulling of large spatial rasters
				
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
				Please note that code asks you the anticipated size of watershed to see if you agree with the delineation. If you do noy 
				agree, the outlet is moved to the maxima raster SCA for all the rasters surrounding the original raster. This process 
				occurs recursively until the user agrees with the rough size of the watershed.

ExtractSlpSCA		-	Function does the processes described above but does not run delineation. The purpose of this code is to extract slope
				and SCA for loactions of interest
				 This is a function file is to extract spatial values: slope and SCA

				Function inputs
				 -DEM: DEM shall be inputted as projected into UTM and clipped as needed
				 -Coordinates of interest: a 2xn matrix (for this paper these would be the soils locations that we are interested in)
				 -An indication of either Dinf or D8 for the delineation

				Function outputs
				 -Output extracted data
				   Output will be a 4xn matrix (where n is the number of points of interest there are)
				   Output matrix columns: Lat  Long  Slp SCA
				
				Example: output = SSExtract(DEM,POI,"D8")

corrplotfiggenTI	-	The purpose of this code is to run the SSTTextract and create correlation plots for to resulting data

corrplotfiggen		-	The purpose of this code is to run the SSextract and create correlation plots for to resulting data
				This code also produces the multiple regression summary figure

soilsdistandcomparessurgo -	This code generates the figures the compare distributed soil maps (using the multiple regressions resulting from 
				corrplotfiggen) to SSURGO distributed soils. This code also compares (via bootstrapping) the informed multiple 
				regressions to blindly using the finest resolution for soil distribution.

Coarsenres		-	This code recursively coarsens the LiDAR DEM and extracts slopes and SCAs from each coarsened DEM. Currently this is setup
				for LiDAR only, though both 1m and 1/3as have also been coarsened for the analysis found in the corresponding paper. 
				The workspace is saved after each coarsening because this process was time intense and ran unsupervised for several days.
				Workspace saving was done in order to avoid any data losses due to crashed code.

bestres			-	Data extracted from recursively coarsened DEMs generated from Coarsenres and is saved locally. This code generates the figures
				seen in the paper.
--------------------------------------------------------------------------------------

Spatial data
--------------------------------------------------------------------------------------
1as_projandtrim.tif	-	DEM raster for 1as sourced from USGS https://apps.nationalmap.gov/downloader/#/
GDEM_projandtrim.tif	-	DEM raster for GDEM sourced from https://doi.org/10.5067/ASTER/ASTGTM.003
STRM_projandtrim.tif	-	DEM raster for SRTM sourced from https://doi.org/10.1029/EO081i048p00583

The remaining DEMs used in this project are too large to include in this publication (per GitHub requirements).
LiDAR data can be found here: https://maps.vcgi.vermont.gov/LidarFinder/
1m, and 1/3as can be found here: https://apps.nationalmap.gov/downloader/#/

Please note that a great deal of raster merging and clipping was done for the LiDAR data to make the spatial data small enough for the computer (64 cores)
to handle. This code is not published but can be shared upon request. Because of the size of many of the rasters, many of the resulting rasters are saved
locally once raster manipulation is completed. As a result, the code is not stand alone because of the intensity of the computation needed to make the 
code stand alone.
