# The purpose of this code is to run the SSTTextract and create correlation plots for to resulting data

##### Step 0: add in libraries #####

if (!require("pacman")) install.packages("pacman")
pacman::p_load(devtools)
pacman::p_load(rgdal,elevatr,raster,devtools,rgdal,readxl,plot.matrix,grid,corrplot)
old_path <- Sys.getenv("PATH")
Sys.setenv(PATH = paste(old_path, "/groups/bse5304g/Elyce/TauDEM/bin", sep = ":"))
Sys.unsetenv("http_proxy"); Sys.unsetenv("https_proxy")
#######

##### Step 1: Run function with all data of interest ####
#define common data
setwd("/groups/bse5304g/Elyce")
soilsreadin = read.csv("Soils/MoreSoils.csv") #soils data is not publically availible
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
POI = SpatialPoints(matrix(cbind(soilsreadin$Longitude[!is.na(soilsreadin$Latitude)],soilsreadin$Latitude[!is.na(soilsreadin$Latitude)]),length(which(!is.na(soilsreadin$Latitude))),2),proj4string = crs_ll)

#1as USGS
rast = raster("DEMs/1as_projandtrim.tif")
D81as = SSextract(rast,POI,"D8")
Dinf1as = SSextract(rast,POI,"Dinf")

#1/3as USGS
rast = raster("DEMs/13as_projandtrim.tif")
D813as = SSextract(rast,POI,"D8")
Dinf13as = SSextract(rast,POI,"Dinf")

#GDEM
rast = raster("DEMs/GDEM_projandtrim.tif")
D8GDEM = SSextract(rast,POI,"D8")
DinfGDEM = SSextract(rast,POI,"Dinf")

#SRTM
rast = raster("DEMs/STRM_projandtrim.tif")
D8SRTM = SSextract(rast,POI,"D8")
DinfSRTM = SSextract(rast,POI,"Dinf")

#### These were saved locally using TuaDEM and exporting relevant rasters >
#1m USGS
D81m = D8SRTM; Dinf1m = D8SRTM
setwd("/groups/bse5304g/test/1mD8")
sca = raster("sca.tif")
D81m$sca = raster::extract(sca,POI)
slp = raster("slp.tif")
D81m$slp = raster::extract(slp,POI)

setwd("/groups/bse5304g/test/1mDinf")
sca = raster("sca.tif")
Dinf1m$sca = raster::extract(sca,POI)
slp = raster("slp.tif")
Dinf1m$slp = raster::extract(slp,POI)

#LIDAR
D8LIDAR = D8SRTM; DinfLIDAR = D8SRTM
setwd("/groups/bse5304g/test/LIDARD8")
sca = raster("sca.tif")
D8LIDAR$sca = raster::extract(sca,POI)
slp = raster("slp.tif")
D8LIDAR$slp = raster::extract(slp,POI)

setwd("/groups/bse5304g/test/LIDARDinf")
sca = raster("sca.tif")
DinfLIDAR$sca = raster::extract(sca,POI)
slp = raster("slp.tif")
DinfLIDAR$slp = raster::extract(slp,POI)
#### These were saved locally using TuaDEM and exporting relevant rasters <

setwd("/groups/bse5304g/Elyce")
########

##### Step 2: Create correlation matrix #####
#clean up soils data
soils = soilsreadin
cor4plot = data.frame(slp_LIDARdinf = matrix(NA,length(names(soils))-3,1),lnsca_LIDARdinf = NA,slp_LIDARd8 = NA,lnsca_LIDARd8 = NA,
                      slp_1mdinf = NA,lnsca_1mdinf = NA,slp_1md8 = NA,lnsca_1md8 = NA,
                      slp_1_3asdinf = NA,lnsca_1_3asdinf = NA,slp_1_3asd8 = NA,lnsca_1_3asd8 = NA,
                      slp_1asdinf = NA,lnsca_1asdinf = NA,slp_1asd8 = NA,lnsca_1asd8 = NA,
                      slp_SRTMdinf = NA,lnsca_SRTMdinf = NA,slp_SRTMd8 = NA,lnsca_SRTMd8 = NA,
                      slp_GDEMdinf = NA,lnsca_GDEMdinf = NA,slp_GDEMd8 = NA,lnsca_GDEMd8 = NA)
rownames(cor4plot) = names(soils)[4:length(names(soils))]
dfgets = c("DinfLIDAR","D8LIDAR","Dinf1m","D81m","Dinf13as","D813as","Dinf1as","D81as","DinfSRTM","D8SRTM","DinfGDEM","D8GDEM")

for(i in 1:length(rownames(cor4plot))){
  for(j in 1:length(dfgets)){
    df = get(dfgets[j])
    cor4plot[i,(j-1)*2+1] = cor(soils[,i+3],df$slp,use = "complete.obs")
    cor4plot[i,(j-1)*2+2] = cor(soils[,i+3],log(df$sca),use = "complete.obs")
  }
}


#cor plot for slope and sca for DC and LOC
corrplot(as.matrix(cor4plot),cl.pos = 'n',tl.col = "black",mar=c(4,2,4,2),main = "Little Otter Creek and Dead Creek")


##### Step 3: Multiple (2) regression #####
regressioncoef = matrix(NA,length(4:length(names(soils))),2); rownames(regressioncoef) = names(soils)[4:length(names(soils))]
slpsca = data.frame(matrix(NA,length(which(!is.na(D8GDEM$slp))),length(dfgets)*2))

for(i in 1:length(dfgets)){
  df = get(dfgets[i])
  slpsca[,i] = df$slp[!is.na(df$slp)];  names(slpsca)[i] = paste0("Slope_",dfgets[i])
  slpsca[,i+length(dfgets)] = df$sca[!is.na(df$sca)];  names(slpsca)[i+length(dfgets)] = paste0("ln(SCA)_",dfgets[i])
  
}

cor4plot2 = cor4plot
cor4plot2$step2 = NA
iscolin = data.frame(matrix(TRUE,length(dfgets)*2,length(dfgets)*2))
names(cor4plot2)[seq(1,length(dfgets)*2,2)] = names(slpsca)[1:length(dfgets)]
names(cor4plot2)[seq(2,length(dfgets)*2,2)] = names(slpsca)[(length(dfgets)+1):(length(dfgets)*2)]
names(iscolin)[seq(1,length(dfgets)*2,2)] = names(slpsca)[1:length(dfgets)]
names(iscolin)[seq(2,length(dfgets)*2,2)] = names(slpsca)[(length(dfgets)+1):(length(dfgets)*2)]
iscolin[grepl("Slope",names(iscolin)),grepl("SCA",names(iscolin))] = FALSE
iscolin[grepl("SCA",names(iscolin)),grepl("Slope",names(iscolin))] = FALSE
slpsca[,(length(dfgets)+1):(length(dfgets)*2)] = log(slpsca[,(length(dfgets)+1):(length(dfgets)*2)])
pval = matrix(1,6,1); slpcoef = pval; lnscacoef = pval 
for(i in 4:length(names(soils))){
    response = soils[,i]
    regressioncoef[i-3,1] = names(which.max(abs(cor4plot2[i-3,!is.na(cor4plot2[i-3,])])))
    best = slpsca[names(slpsca)==regressioncoef[i-3,1]]
    #Take out any of the colinear inputs
    ind = which(!iscolin[which(names(best)==names(iscolin)),])
    hold = slpsca[,which(names(slpsca) %in% names(iscolin)[ind])]
    AIC = NA; AICbest = NA; AICone = AIC(lm(response~.,data=best))
    for(j in 1:length(names(hold))){
      linearmodel = lm(response~.,data = cbind(best,hold[,j]))
      AIC = AIC(linearmodel)
      if(j==1){AICbest=AIC}
      if(AIC<=AICbest){
        AICbest = AIC
        secondbest = hold[,j]; secondbestname = names(hold)[j]
      }
      print(AIC)
    }
    print(paste0("AIC of one: ",AICone))
    linearmodel = lm(response~.,data=cbind(best,secondbest))
    cor4plot2$step2[i-3] = cor(linearmodel$fitted.values,response)
    x = summary(linearmodel)
    pval[i-3] = 1-pf(x$fstatistic[1],x$fstatistic[2],x$fstatistic[3])
    slpcoef[i-3] = round(linearmodel$coefficients[2],2)
    lnscacoef[i-3] = round(linearmodel$coefficients[3],3)
    regressioncoef[i-3,2] = secondbestname
    
}
    
TFmatrix = matrix(NA,length(regressioncoef[,1]),(length(names(slpsca))))
rownames(TFmatrix) = rownames(regressioncoef)
colnames(TFmatrix) = names(slpsca)
rownames(regressioncoef)=NULL

for(i in 1:length(regressioncoef[,1])){
  TFmatrix[i,which(regressioncoef[i,1]==colnames(TFmatrix))] = "Best"
  TFmatrix[i,which(regressioncoef[i,2]==colnames(TFmatrix))] = "Second best"
}

names(cor4plot2)[25] = "Regression2"
library('RColorBrewer')
par(mfrow = c(1,2))
cor4plot5 = data.frame(cor4plot2[1:6,25])
rownames(cor4plot5) = c("Sand","Clay","AWC pred","Organic Matter","Total N","Mg")
colnames(cor4plot5) = "Multiple\nRegression"

TFnew = data.frame(Slope_LiDAR = TFmatrix[1:4,1],Slope_1m = TFmatrix[1:4,3],Slope_13as = TFmatrix[1:4,5],
                   Slope_1as = TFmatrix[1:4,7],Slope_SRTM = TFmatrix[1:4,9],Slope_GDEM = TFmatrix[1:4,11],
                   SCA_LiDAR = TFmatrix[1:4,13],SCA_1m = TFmatrix[1:4,15],SCA_13as = TFmatrix[1:4,17],
                   SCA_1as = TFmatrix[1:4,19],SCA_SRTM = TFmatrix[1:4,21],SCA_GDEM = TFmatrix[1:4,23])
for(i in 1:length(names(TFnew))){
  if(sum(which(!is.na(TFnew[,i])))!=0){TFnew[which(!is.na(TFnew[,i])),i] = 1}
  if(sum(which(is.na(TFnew[,i]))!=0)){TFnew[which(is.na(TFnew[,i])),i] = 0}
}
TFnew$Slope_1as[c(1,2,4)] = 2

layout(t(matrix(c(1:4))),widths = c(.95,2,2,1),heights = 1)
hold = matrix(0,4,1); colnames(hold) = "R Squared";rownames(hold) = c("Sand","Clay","Pred AWC","Org Matter")
par(mar = c(12,7,3,0))
plot(hold,col = 'white',las=2,main = "",xlab="",ylab="",key=NULL,cex.axis=1.45,axis.col=list(side=1))
text(1,1:4,round((cor4plot2$Regression2[1:4])^2,2),cex=1.5,col='black')

colorTF = c('white','black','grey40')
par(mar = c(12,1,3,1))
plot((as.matrix(TFnew[,1:6])),col = colorTF,las=2,main = "",xlab="",ylab="",key=NULL,cex.axis=1.45,axis.row = NULL,axis.col=list(side=1))
text(4,4,slpcoef[1],cex=1.25,col='white')
text(4,3,slpcoef[2],cex=1.25,col='white')
text(3,2,slpcoef[3],cex=1.25,col='white')
text(4,1,slpcoef[4],cex=1.25,col='white')

plot((as.matrix(TFnew[,7:12])),col = colorTF,las=2,main = "",xlab="",ylab="",key=NULL,cex.axis=1.45,axis.row = NULL,axis.col=list(side=1))
text(2,1:4,lnscacoef[1:4],cex=1.25,col='white')

par(mar = c(17.5,3,10,3))
legendmatrix = matrix(c(1,2))
rownames(legendmatrix) = c("Dinf","D8")
plot(t(legendmatrix),col = colorTF[2:3],las=2,main = "\nLegend",xlab="",ylab="",key = NULL,cex.axis=1.25,axis.row = NULL)

