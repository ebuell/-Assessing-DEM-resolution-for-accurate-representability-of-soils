library('sp')
library("hydroGOF")
library("raster")
library("RColorBrewer")
##### define training and boostrapping criteria #####
N = 10000
############

##### set up dataframes #####
#soils dataframe
setwd("/groups/bse5304g/Elyce")
soils = read.csv("Soils/RelevantSoils.csv")
soils$X = NULL
soils = soils[1:59,]
#get POIs
proj4_utm = paste0("+proj=utm +zone=", 18, " +datum=WGS84 +units=m +no_defs")
proj4_ll = "+proj=longlat"; crs_ll=CRS(proj4_ll); crs_utm=CRS(proj4_utm)
POI = SpatialPoints(matrix(cbind(soils$Longitude[!is.na(soils$Latitude)],soils$Latitude[!is.na(soils$Latitude)]),length(which(!is.na(soils$Latitude))),2),proj4string = crs_ll)
#POI = spTransform(POI,CRSobj = crs_utm)



#regression dataframes - only run for first iter
if(TRUE){
  DinfLIDAR = data.frame(slp = matrix(NA,length(soils$SampleNumber)),lnsca=NA)
  setwd("/groups/bse5304g/test/LIDARDinf")
  sca = raster("sca.tif")
  DinfLIDAR$lnsca = (raster::extract(sca,POI))
  slp = raster("slp.tif")
  DinfLIDAR$slp = raster::extract(slp,POI)
  D8LIDAR = data.frame(slp = matrix(NA,length(soils$SampleNumber)),lnsca=NA)
  setwd("/groups/bse5304g/test/LIDARD8")
  sca = raster("sca.tif")
  D8LIDAR$lnsca = (raster::extract(sca,POI))
  slp = raster("slp.tif")
  D8LIDAR$slp = raster::extract(slp,POI)
  BestSlopes = data.frame(D81as = matrix(NA,length(soils$SampleNumber)),Dinf13as = NA)
  setwd("/groups/bse5304g/Elyce")
  rast = raster("DEMs/1as_projandtrim.tif")
  D81as = SSextract(rast,POI,"D8")
  BestSlopes$D81as = D81as$slp
  rast = raster("DEMs/13as_projandtrim.tif")
  Dinf13as = SSextract(rast,POI,"Dinf")
  BestSlopes$Dinf13as = Dinf13as$slp
  BestlnSCAs = data.frame(Dinf1m = matrix(NA,length(soils$SampleNumber),1))
  setwd("/groups/bse5304g/test/1mDinf")
  sca = raster("sca.tif")
  BestlnSCAs$Dinf1m = log(raster::extract(sca,POI))
}
######

##### import SSURGO soils #####
if (!require("pacman")) install.packages("pacman")
pacman::p_load(raster,soilDB)
library('googlesheets4'); library('soilDB'); library('rgdal'); library('aqp')

q <- "SELECT mukey, muname
FROM mapunit
WHERE mukey IN (
SELECT * from SDA_Get_Mukey_from_intersection_with_WktWgs84('point(-73.118 44.136)')
)"

res <- SDA_query(q)
knitr::kable(res)
#################################
#####Ggetting mukey cokey and sand% for LOC 
###################################
myflowgage_id="04282650"
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2010-01-01",
                         end_date = "2020-01-01")
trunc((180+myflowgage$declon)/6+1)
proj4_ll = "+proj=longlat"
crs_ll=CRS(proj4_ll)
LittleOtterBsn=readOGR("/groups/bse5304g/Elyce/Figs/Roja/LOCwb.shp")
LittleOtterBsn_ll=spTransform(LittleOtterBsn,crs_ll)
mybbox=c(LittleOtterBsn_ll@bbox)
mysoil = mapunit_geom_by_ll_bbox(mybbox,source = 'sda') # Careful, this might error out.
DCwb=readOGR("/groups/bse5304g/Elyce/Figs/Roja/DCwb.shp")
DCwb=spTransform(DCwb,crs_ll)
mybbox=c(DCwb@bbox)
mysoil3 = mapunit_geom_by_ll_bbox(mybbox,source = 'sda') # Careful, this might error out.
mukey_statement = format_SQL_in_statement(unique(mysoil$mukey))
q_mu2co = paste("SELECT mukey,cokey FROM component WHERE mukey IN ", mukey_statement, sep="")
mu2co = SDA_query(q_mu2co)
cokey_statement = format_SQL_in_statement(unique(mu2co$cokey))

####percent sand from chorizon table
q_co2ch = paste("SELECT cokey,hzname,OM_r,awc_r,sandtotal_r,claytotal_r,hzname FROM chorizon WHERE cokey IN ", cokey_statement, sep="")
q.2 = paste("SELECT cokey,comppct_r FROM legend INNER JOIN mapunit ON mapunit.lkey = legend.lkey
       LEFT OUTER JOIN component ON component.mukey = mapunit.mukey WHERE cokey IN ", cokey_statement, sep="")
co2ch = SDA_query(q_co2ch)
co2dom = SDA_query(q.2)
mu2ch=merge(mu2co,co2ch)
mu2ch = merge(mu2ch,co2dom)
mu2ch = mu2ch[which((mu2ch$hzname=="H1")|mu2ch$hzname=="Ap"|mu2ch$hzname=="A"|mu2ch$hzname=="R"),]

mu2chdomH1 = data.frame(mukey = unique(mu2ch$mukey),cokey=NA,sand_r=NA,clay_r=NA,
                        OM_r=NA,AWC_r = NA,comppct_r=NA,horizon = NA)

for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(mu2ch$mukey==mu2chdomH1$mukey[i])
  ind3 = ind[which.max(mu2ch$comppct_r[ind])]
  mu2chdomH1$cokey[i]=mu2ch$cokey[ind3]
  mu2chdomH1$horizon[i]=mu2ch$hzname.1[ind3]
  mu2chdomH1$OM_r[i]=mu2ch$OM_r[ind3]
  mu2chdomH1$clay_r[i]=mu2ch$claytotal_r[ind3]
  mu2chdomH1$sand_r[i]=mu2ch$sandtotal_r[ind3]
  mu2chdomH1$AWC_r[i]=mu2ch$awc_r[ind3]
  mu2chdomH1$comppct_r[i] = mu2ch$comppct_r[ind3]
  if(mu2ch$hzname.1[i]=="R"){
    mu2chdomH1$OM_r[i]=0
    mu2chdomH1$clay_r[i]=0
    mu2chdomH1$sand_r[i]=0
    mu2chdomH1$AWC_r[i]=0
  }
}

######create SSURGO rasters
#generate OM map from SSURGO
setwd("/groups/bse5304g/test/bootstrapslopes")
LOCa = raster("LOCamukey.tif")
DCa = raster("DCamukey.tif")

#Sand named a;  
#Clay named b;        
#AWC names c;      
#OM named d;
LOCb = LOCa; LOCc = LOCa; LOCd = LOCa
DCb = DCa; DCc = DCa; DCd = DCa

#LOC
for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(values(LOCa)==mu2chdomH1$mukey[i])
  if(is.na(mu2chdomH1$sand_r[i])){
    values(LOCa)[ind] = 0 
    values(LOCb)[ind] = 0
    values(LOCc)[ind] = 0
    values(LOCd)[ind] = 0
  }else{
    values(LOCa)[ind] = mu2chdomH1$sand_r[i]
    values(LOCb)[ind] = mu2chdomH1$clay_r[i]
    values(LOCc)[ind] = mu2chdomH1$AWC_r[i]
    values(LOCd)[ind] = mu2chdomH1$OM_r[i]
    }
}
values(LOCa)[which(values(LOCa)>100)] = 0
values(LOCb)[which(values(LOCb)>100)] = 0
values(LOCc)[which(values(LOCc)>100)] = 0
values(LOCd)[which(values(LOCd)>100)] = 0


#DC
for(i in 1:length(mu2chdomH1$mukey)){
  ind = which(values(DCa)==mu2chdomH1$mukey[i])
  #if(is.na(mu2chdomH1$sand_r[i])){
   # values(LOCa)[ind] = 0 
   # values(LOCb)[ind] = 0
   # values(LOCc)[ind] = 0
   # values(LOCd)[ind] = 0
  #}else{
    values(DCa)[ind] = mu2chdomH1$sand_r[i]
    values(DCb)[ind] = mu2chdomH1$clay_r[i]
    values(DCc)[ind] = mu2chdomH1$AWC_r[i]
    values(DCd)[ind] = mu2chdomH1$OM_r[i]
  #}
}
values(DCa)[which(values(DCa)>100)] = 0
values(DCb)[which(values(DCb)>100)] = 0
values(DCc)[which(values(DCc)>100)] = 0
values(DCd)[which(values(DCd)>100)] = 0



##### Run bootstrap ####

######create relevant rasters for bootstrapping fig creation 
#import relevant rasters
setwd("/groups/bse5304g/test/bootstrapslopes")
slp2DC = raster("slp13asdinf1mDC.tif")
slp1DC = raster("slp1asd81mDC.tif")
slp2LOC = raster("slp13asdinf1mLOC.tif")
slp1LOC = raster("slp1asd81mLOC.tif")

setwd("/groups/bse5304g/test/1mDinf/DCspecific")
onemscaDC = log(raster("scaclip.tif"))
setwd("/groups/bse5304g/test/1mDinf/LOCSpecific")
onemscaLOC = log(raster("scaclp.tif"))


#Sand;  slp1: slp1asd8as1m,     onemsca: ln(sca1mdinf)
#Clay;  slp1: slp1asd8as1m,     onemsca: ln(sca1mdinf)
#AWC;   slp2: slp13asdinfas1m,  onemsca: ln(sca1mdinf)
#OM;    slp1: slp1asd8as1m,     onemsca: ln(sca1mdinf)

#seq for letters
letterseq = c("a","b","c","d")

#rast zlims
zlimits = rbind(c(0,100),
                c(0,100),
                c(0,.4),
                c(0,10))

colorforplots = brewer.pal(9,"YlGnBu")[2:9]
colorforplots2 = c(brewer.pal(9,"YlGnBu")[c(2,4,6,8:9)],brewer.pal(9,"Reds")[2:9])
colorforplots4 = c(brewer.pal(9,"YlGnBu")[c(2,6,9)],brewer.pal(9,"Reds")[1:9])
xlabtext = c("LOC Sand content",
             "LOC Clay content",
             "LOC Predicted AWC",
             "LOC Organic matter",
             "DC Sand content",
             "DC Clay content",
             "DC Predicted AWC",
             "DC Organic matter")
unitstext = c("%","%","cm/cm","%")

#stats holder
fitstats = data.frame(iter = 1:N, trainx = NA, trainsig = NA, testx = NA,testsig = NA,
                      finestresdinf = NA,finestresd8 = NA,
                      bestguesstrainxbar = NA, bestguesstrainsd = NA, bestguesstestxbar = NA, bestguesstestsd = NA,
                      bestguesscoefp = NA, bestguesscorrtrain=NA, bestguesscorrtest=NA)
forplottcoefvals = data.frame(iter = 1:N,sandc = NA, sandb1 = NA, sandb2 = NA,
                              clayc = NA, clayb1 = NA, clayb2 = NA,
                              awcc = NA, awcb1 = NA, awcb2 = NA,
                              OMc = NA, OMb1 = NA, OMb2 = NA)
forplottcoefvalsind = rbind(c(2:4),c(5:7),c(8:10),c(11:13))

spatial = data.frame(lnsca1mdinf = BestlnSCAs$Dinf1m, slp1asd8 = D81as$slp, slp13asdinf = Dinf13as$slp)
spatialind = rbind(c(2,1),
                   c(2,1),
                   c(3,1),
                   c(2,1))
rownames(spatialind) = c("sand","clay","awc","om")
colnames(spatialind) = c("Slope","lnSCA")
spatialind = data.frame(spatialind)

#set up soils df of interest
soils = data.frame(sand = soils$soil_texture_sand,clay = soils$soil_texture_clay,
                       awc = soils$pred_water_capacity, om = soils$organic_matter)



##### Run bootstrap ####
#set up bootstraph dfs


#generate soil dist LOC
trainnum = round(length(soilsdata)/2)

layout(matrix(rbind(c(1:5),c(6:10)),2,5),widths = c(1.1,1.2,1.2,.8,1.3),heights = c(.5,.5))
#bootstrap
for(j in 3:4){
  indnna = !is.na(soils[,j])
  soilsdata = soils[indnna,j]
  for(i in 1:N){
    #define subsamples
    subsample = sample(1:length(soilsdata),trainnum,replace=FALSE)
    trainingsample = soils[subsample,j]
    subtest = !(1:length(soilsdata) %in% subsample)
    testsample = soilsdata[subtest]
    
    #define statistics for subsamples
    fitstats$trainx[i] = mean(trainingsample)
    fitstats$testx[i] = mean(testsample)
    fitstats$trainsig[i] = sd(trainingsample)
    fitstats$testsig[i] = sd(testsample)
    
    #select data for bootstrap
    dataforbootstrap = data.frame(bestSCA = spatial[indnna,spatialind$lnSCA[j]],bestslp = spatial[indnna,spatialind$Slope[j]])
    
    
    #for finest res D8
    df = data.frame(slplid8 = D8LIDAR$slp[indnna],lnscalid8 = log(D8LIDAR$lnsca)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subsample,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subtest,])
    SEline = (testsample-testvals)^2
    SEybar = (testsample-mean(testsample,na.rm=TRUE))^2
    fitstats$finestresd8[i] = 1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE)
    
    #for finest res Dinf
    df = data.frame(slplidinf = DinfLIDAR$slp[indnna],lnscalidinf = log(DinfLIDAR$lnsca)[indnna])
    linearmodel = lm(trainingsample~.,data = df[subsample,])
    #predict values
    testvals = predict(linearmodel,newdata = df[subtest,])
    SEline = (testsample-testvals)^2
    SEybar = (testsample-mean(testsample,na.rm=TRUE))^2
    fitstats$finestresdinf[i] = 1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE)
    
    #for best guess
    linearmodelbest = lm(trainingsample~.,data = dataforbootstrap[subsample,])
    fitstats$bestguesstrainxbar[i] = mean(linearmodelbest$fitted.values)
    fitstats$bestguesstrainsd[i] = sd(linearmodelbest$fitted.values)
    fitstats$bestguessinterp[i] = summary(linearmodelbest)$coefficients[1,4]
    fitstats$bestguessp[i] = pf(summary(linearmodelbest)$fstatistic[1],summary(linearmodelbest)$fstatistic[2],summary(linearmodelbest)$fstatistic[3],lower.tail = FALSE)
    fitstats$bestguesscorrtrain[i] = sqrt(summary(linearmodelbest)$r.squared)
    #predict values
    testvals = predict(linearmodelbest,newdata = dataforbootstrap[subtest,])
    fitstats$bestguesstestxbar[i] = mean(testvals)
    fitstats$bestguesstestsd[i] = sd(testvals)
    SEline = (testsample-testvals)^2
    SEybar = (testsample-mean(testsample,na.rm=TRUE))^2
    fitstats$bestguesscorrtest[i] = 1-sum(SEline,na.rm=TRUE)/sum(SEybar,na.rm=TRUE)
    
    forplottcoefvals[i,forplottcoefvalsind[j,]] = linearmodelbest$coefficients
    
  }
  colorforplots3 = colorforplots
 
  fitstats$bestguesscorrtest[which(fitstats$bestguesscorrtest< -1)] = -1
  fitstats$finestresdinf[which(fitstats$finestresdinf< -1)] = -1
  fitstats$finestresd8[which(fitstats$finestresd8< -1)] = -1
  
  #plot densities 
  par(mar = c(4,4,4,1))
  plot(density(fitstats$bestguesscorrtest),col = 'darkgoldenrod3',xlab = "R Squared",ylab="Density",bty='n',main=paste0("\n(A",j,") ",gsub("LOC ","",xlabtext[j])),xlim = c(-1,1),lwd=3,ylim = c(0,4.5))
  lines(density(fitstats$finestresd8),col = 'coral3',lwd=3)
  lines(density(fitstats$finestresdinf),col = 'darkcyan',lwd=3)
  legend('topleft',legend = c("Multivariate Regression","Finest res D8 slp&sca","Finest res Dinf slp&sca"),
         col = c("darkgoldenrod3","coral3","darkcyan"),lwd=3,bty = 'n')
  par(xpd=TRUE); text(-1.22,-.61,">"); par(xpd=FALSE)
  
  
  #plot proposed soil dist
  par(mar=c(4,0,4,2))
  ylimit = c(4880050, 4904696); xlimit = c(638940, 656183)
  slope = slp1LOC; if(j==3){slope = slp2LOC}
  lnsca = onemscaLOC
  rast = mean(forplottcoefvals[,forplottcoefvalsind[j,][1]])+mean(forplottcoefvals[,forplottcoefvalsind[j,][2]])*lnsca+mean(forplottcoefvals[,forplottcoefvalsind[j,][3]])*slope
  values(rast)[which(values(rast)<0)] = 0
  plot(rast,main = paste0("\n (B",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  if(j==1){
    par(xpd=TRUE)
    text(647000,4909200,"Multivariate TIV Regression",cex=1.5)
    par(xpd=FALSE)
  }
  #plot ssurgo soil dist
  plot(rast,main = paste0("\n (C",j,") ",xlabtext[j]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  plot(get(paste0("LOC",letterseq[j])),add=TRUE,legend=FALSE,col=colorforplots3)
  
  if(j==1){
    par(xpd=TRUE)
    text(645700,4909200,"SSURGO Estimated",cex=1.5)
    par(xpd=FALSE)
  }
  
  #plot proposed soil dist
  par(mar=c(4,0,4,4))
  ylimit = c(4829000, 4884500); xlimit = c(630000, 638500)
  slope = slp1DC; if(j==3){slope = slp2DC}
  lnsca = onemscaDC
  rast = mean(forplottcoefvals[,forplottcoefvalsind[j,][1]])+mean(forplottcoefvals[,forplottcoefvalsind[j,][2]])*lnsca+mean(forplottcoefvals[,forplottcoefvalsind[j,][3]])*slope
  values(rast)[which(values(rast)<0)] = 0
  plot(rast,main = paste0("\n (D",j,") ",xlabtext[j+4]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  if(j==1){
    par(xpd=TRUE)
    text(634851,4887819,"Multivariate TIV Regression",cex=1.5)
    par(xpd=FALSE)
  }
  #plot ssurgo soil dist
  plot(rast,main = paste0("\n (E",j,") ",xlabtext[j+4]),bty='n',xaxt='n',yaxt='n',xlim=xlimit,ylim=ylimit,zlim=zlimits[j,],col=colorforplots3,legend=FALSE)
  plot(get(paste0("DC",letterseq[j])),add=TRUE,legend=FALSE,col=colorforplots3)
  
  par(xpd=TRUE)
  legend(x = 641851, y = 4879319,bty='n',legend = c(paste0(zlimits[j,2],unitstext[j]),"","","","","","",paste0(zlimits[j,1],unitstext[j])),fill = rev(colorforplots3),
         border = NA,y.intersp = 0.5,cex = 1, text.font = 1)
  par(xpd=FALSE)
  if(j==1){
    par(xpd=TRUE)
    text(632851,4887819,"SSURGO Estimated",cex=1.5)
    par(xpd=FALSE)
  }
}

####Measured Properties vs SSURGO ####
SP_ll=SpatialPoints(matrix(c(soils$Longitude,soils$Latitude), 
                           ncol = 2, byrow = FALSE),proj4string =crs_ll)
SP_utm = spTransform(SP_ll,crs_utm)
compareprops = data.frame(sand_meas = soils$soil_texture_sand,
                          clay_meas = soils$soil_texture_clay,
                          awc_meas = soils$pred_water_capacity,
                          om_meas = soils$organic_matter,
                          sand_ssurgo_r = extract(LOCa,SP_utm),
                          clay_ssurgo_r = extract(LOCb,SP_utm),
                          awc_ssurgo_r = extract(LOCc,SP_utm),
                          om_ssurgo_r = extract(LOCd,SP_utm))


for (i in 1:length(compareprops$sand_meas)){
  if(is.na(compareprops$sand_ssurgo_r[i])){
    compareprops$sand_ssurgo_r[i] = extract(DCa,SP_utm[i])
    compareprops$clay_ssurgo_r[i] = extract(DCb,SP_utm[i])
    compareprops$awc_ssurgo_r[i] = extract(DCc,SP_utm[i])
    compareprops$om_ssurgo_r[i] = extract(DCd,SP_utm[i])
  }
}


par(mfrow=c(2,2))
par(mar = c(4,4,4,1))

linm1 = lm(soils$soil_texture_sand~spatial$slp1asd8+spatial$lnsca1mdinf)
plot(compareprops$sand_meas,compareprops$sand_ssurgo_r,xlim = c(0,120),ylim = c(0,120),ylab = "Predicted sand content (%)",xlab = "Measured sand content (%)",main = "Predicted vs Measured\nsand content",pch=3,col = 'black',cex.main = 1.5)
points(compareprops$sand_meas,linm1$fitted.values,pch=16)
segments(0,0,100,100,lty="dotted")
SEline = (compareprops$sand_meas-compareprops$sand_ssurgo_r)^2
SEybar = (compareprops$sand_meas-mean(compareprops$sand_meas))^2
r21 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
SEline = (compareprops$sand_meas-linm1$fitted.values)^2
r22 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
legend(0,120,box.lty=0,legend = c(paste("SSURGO R2 =",r21),paste0("Multivariate regression R2=",r22)),pch=c(3,16),col = 'black',border='white')

linm1 = lm(soils$soil_texture_clay~spatial$slp1asd8+spatial$lnsca1mdinf)
plot(compareprops$clay_meas,compareprops$clay_ssurgo_r,xlim = c(0,87),ylim = c(0,87),ylab = "Predicted clay content (%)",xlab = "Measured clay content (%)",main = "Predicted vs Measured\nclay content",pch=3,col = 'black',cex.main = 1.5)
points(compareprops$clay_meas,linm1$fitted.values,pch=16)
segments(0,0,70,70,lty="dotted")
SEline = (compareprops$clay_meas-compareprops$clay_ssurgo_r)^2
SEybar = (compareprops$clay_meas-mean(compareprops$clay_meas))^2
r21 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
SEline = (compareprops$clay_meas-linm1$fitted.values)^2
r22 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
legend(0,87,box.lty=0,legend = c(paste("SSURGO R2 =",r21),paste0("Multivariate regression R2=",r22)),pch=c(3,16),col = 'black',border='white')

linm1 = lm(soils$pred_water_capacity~spatial$slp13asdinf+spatial$lnsca1mdinf)
plot(compareprops$awc_meas,compareprops$awc_ssurgo_r,xlim = c(.15,.32),ylim = c(.15,.32),ylab = "Predicted AWC (cm/cm)",xlab = "Measured AWC (cm/cm)",main = "Predicted vs Measured\nAWC",pch=3,col = 'black',cex.main = 1.5)
points(compareprops$awc_meas,linm1$fitted.values,pch=16)
segments(0,0,.3,.3,lty="dotted")
SEline = (compareprops$awc_meas-compareprops$awc_ssurgo_r)^2
SEybar = (compareprops$awc_meas-mean(compareprops$awc_meas))^2
r21 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
SEline = (compareprops$awc_meas-linm1$fitted.values)^2
r22 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
legend(.15,.32,box.lty=0,legend = c(paste("SSURGO R2 =",r21),paste0("Multivariate regression R2=",r22)),pch=c(3,16),col = 'black',border='white')

linm1 = lm(soils$organic_matter~spatial$slp1asd8+spatial$lnsca1mdinf)
plot(compareprops$om_meas,compareprops$om_ssurgo_r,xlim = c(0,15),ylim = c(0,15),ylab = "Predicted organic matter content (%)",xlab = "Measured organic matter content (%)",main = "Predicted vs Measured\norganic matter content",pch=3,col = 'black',cex.main = 1.5)
points(compareprops$om_meas,linm1$fitted.values,pch=16)
segments(0,0,12,12,lty="dotted")
SEline = (compareprops$om_meas-compareprops$om_ssurgo_r)^2
SEybar = (compareprops$om_meas-mean(compareprops$om_meas))^2
r21 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
SEline = (compareprops$om_meas-linm1$fitted.values)^2
r22 = round(1-sum(SEline,na.rm = T)/sum(SEybar,na.rm = T),2)
legend(0,15,box.lty=0,legend = c(paste("SSURGO R2 =",r21),paste0("Multivariate regression R2=",r22)),pch=c(3,16),col = 'black',border='white')

