#best res visualization

#the purpose of this code is to visualize the best resolution data

#### data is saved locally >
#load in slopes ad scas for 59 data points
setwd("/groups/bse5304g/Elyce/CreateVis/Recursive Coarsening")
slpd8_LI = read.csv("slpd8_LIDAR.csv")
slpdinf_LI = read.csv("slpdinf_LIDAR.csv")
scad8_LI = read.csv("scad8_LIDAR.csv")
scadinf_LI = read.csv("scadinf_LIDAR.csv")
for(i in 1:length(names(scad8_LI))){
  scad8_LI[,i] = scad8_LI[,i]*i*.7
}

slpd8_1m = read.csv("slpd8_1m.csv")
slpdinf_1m = read.csv("slpdinf_1m.csv")
scad8_1m = read.csv("scad8_1m.csv")
scadinf_1m = read.csv("scadinf_1m.csv")
scad8_1m = scad8_1m*c(1:length(names(scad8_LI)))
for(i in 1:length(names(scad8_1m))){
  scad8_1m[,i] = scad8_1m[,i]*i
}

slpd8_13as = read.csv("slpd8_13as.csv")
slpdinf_13as = read.csv("slpdinf_13as.csv")
scad8_13as = read.csv("scad8_13as.csv")
scadinf_13as = read.csv("scadinf_13as.csv")
#### data is saved locally <

#load in soils
setwd("/groups/bse5304g/Elyce")
soils = read.csv("Soils/MoreSoils.csv") #soils data not publically availible

#set up correlation matrix
setupmatrix = matrix(NA,length(names(scad8_LI))+length(names(scad8_1m))+length(names(scad8_13as)), (length(4:length(names(soils)))+2))
colnames(setupmatrix) = c("aggrgation","res",names(soils)[4:length(names(soils))])
cor4plotscad8 = data.frame(setupmatrix)
cor4plotscad8$aggrgation = gsub("X","",c(names(scad8_LI),names(scad8_1m),names(scadinf_13as)))

setupmatrix = matrix(NA,length(names(slpd8_LI))+length(names(slpd8_1m))+length(names(slpd8_13as)), (length(4:length(names(soils)))+2))
colnames(setupmatrix) = c("aggrgation","res",names(soils)[4:length(names(soils))])
cor4plotslpd8 = data.frame(setupmatrix)
cor4plotslpd8$aggrgation = gsub("X","",names(slpd8_13as))

scanum = 43; slpnum = 178
cor4plotscad8$res[1:scanum] = .7*as.numeric(cor4plotscad8$aggrgation[1:scanum])
cor4plotslpd8$res[1:slpnum] = .7*as.numeric(cor4plotslpd8$aggrgation[1:slpnum])
cor4plotscad8$res[(scanum+1):(scanum*2)] = 1*as.numeric(cor4plotscad8$aggrgation[(scanum+1):(scanum*2)])
cor4plotslpd8$res[(slpnum+1):(slpnum*2)] = 1*as.numeric(cor4plotslpd8$aggrgation[(slpnum+1):(slpnum*2)])
cor4plotscad8$res[(scanum*2+1):length(cor4plotscad8$res)] = 10.3*as.numeric(cor4plotscad8$aggrgation[(scanum*2+1):length(cor4plotscad8$res)])
cor4plotslpd8$res[(slpnum*2+1):length(cor4plotslpd8$res)] = 10.3*as.numeric(cor4plotslpd8$aggrgation[(slpnum*2+1):length(cor4plotslpd8$res)])
cor4plotscadinf = cor4plotscad8; cor4plotslpdinf = cor4plotslpd8

#calculate correlations
#sca
for(i in 1:scanum){
  for(j in 4:length(names(soils))){
    #d8
    x = log(scad8_LI[,i])
    y = soils[,j]
    cor4plotscad8[i,j-1] = cor(x,y)
    #dinf
    x = log(scadinf_LI[,i])
    y = soils[,j]
    cor4plotscadinf[i,j-1] = cor(x,y)
  }
}

#slp
for(i in 1:slpnum){
  for(j in 4:length(names(soils))){
    #d8
    x = slpd8_LI[,i]
    y = soils[,j]
    cor4plotslpd8[i,j-1] = cor(x,y)
    #dinf
    x = slpdinf_LI[,i]
    y = soils[,j]
    cor4plotslpdinf[i,j-1] = cor(x,y)
  }
}

#calculate correlations
#sca
for(i in (scanum+1):(scanum*2)){
  for(j in 4:length(names(soils))){
    #d8
    x = log(scad8_1m[,i-scanum])
    y = soils[,j]
    cor4plotscad8[i,j-1] = cor(x,y)
    #dinf
    x = log(scadinf_1m[,i-scanum])
    y = soils[,j]
    cor4plotscadinf[i,j-1] = cor(x,y)
  }
}

#slp
for(i in (slpnum+1):(slpnum*2)){
  for(j in 4:length(names(soils))){
    #d8
    x = slpd8_1m[,i-slpnum]
    y = soils[,j]
    cor4plotslpd8[i,j-1] = cor(x,y)
    #dinf
    x = slpdinf_1m[,i-slpnum]
    y = soils[,j]
    cor4plotslpdinf[i,j-1] = cor(x,y)
  }
}

#calculate correlations
#sca
for(i in (scanum*2+1):(length(cor4plotscad8$aggrgation))){
  for(j in 4:length(names(soils))){
    #d8
    x = log(scad8_13as[,i-scanum*2])
    y = soils[,j]
    cor4plotscad8[i,j-1] = cor(x,y)
    #dinf
    x = log(scadinf_13as[,i-scanum*2])
    y = soils[,j]
    cor4plotscadinf[i,j-1] = cor(x,y)
  }
}

#slp
for(i in (slpnum*2+1):(length(cor4plotslpd8$aggrgation))){
  for(j in 4:length(names(soils))){
    #d8
    x = slpd8_13as[,i-slpnum*2]
    y = soils[,j]
    cor4plotslpd8[i,j-1] = cor(x,y)
    #dinf
    x = slpdinf_13as[,i-slpnum*2]
    y = soils[,j]
    cor4plotslpdinf[i,j-1] = cor(x,y)
  }
}

#plot for slp dinf only
par(mfrow = c(1,1))
par(mar = c(4,4,4,1))
ylims=.7
plot(cor4plotslpdinf$res[1:(slpnum*2)],cor4plotslpdinf$soil_texture_sand[1:(slpnum*2)],main = "Slope Dinf",type='l',col = 'cyan4',lwd=1.5,ylim = c(-ylims,ylims),xlim = c(1,31),xlab = 'DEM Resolution (m)',ylab = 'Correlation')
points(cor4plotslpdinf$res[1:slpnum],cor4plotslpdinf$soil_texture_sand[1:slpnum],pch=8,cex=.75,col = 'cyan4')
points(cor4plotslpdinf$res[(slpnum+1):(slpnum*2)],cor4plotslpdinf$soil_texture_sand[(slpnum+1):(slpnum*2)],pch=2,cex=.75,col = 'cyan4')
points(cor4plotslpdinf$res[357:359],cor4plotslpdinf$soil_texture_sand[357:359],pch=12,cex=2,col = 'cyan4')
lines(cor4plotslpdinf$res[1:(slpnum*2)], cor4plotslpdinf$soil_texture_clay[1:(slpnum*2)],col = 'darkgoldenrod',lwd=2)
points(cor4plotslpdinf$res[1:slpnum],cor4plotslpdinf$soil_texture_clay[1:slpnum],pch=8,cex=.75,col = 'darkgoldenrod')
points(cor4plotslpdinf$res[(slpnum+1):(slpnum*2)],cor4plotslpdinf$soil_texture_clay[(slpnum+1):(slpnum*2)],pch=2,cex=.75,col = 'darkgoldenrod')
points(cor4plotslpdinf$res[357:359],cor4plotslpdinf$soil_texture_clay[357:359],pch=12,cex=2,col = 'darkgoldenrod')
lines(cor4plotslpdinf$res[1:(slpnum*2)], cor4plotslpdinf$pred_water_capacity[1:(slpnum*2)],col = 'darkblue',lwd=2)
points(cor4plotslpdinf$res[1:slpnum],cor4plotslpdinf$pred_water_capacity[1:slpnum],pch=8,cex=.75,col = 'darkblue')
points(cor4plotslpdinf$res[(slpnum+1):(slpnum*2)],cor4plotslpdinf$pred_water_capacity[(slpnum+1):(slpnum*2)],pch=2,cex=.75,col = 'darkblue')
points(cor4plotslpdinf$res[357:359],cor4plotslpdinf$pred_water_capacity[357:359],pch=12,cex=2,col = 'darkblue')
lines(cor4plotslpdinf$res[1:(slpnum*2)], cor4plotslpdinf$organic_matter[1:(slpnum*2)],col = 'darkred',lwd=2)
points(cor4plotslpdinf$res[1:slpnum],cor4plotslpdinf$organic_matter[1:slpnum],pch=8,cex=.75,col = 'darkred')
points(cor4plotslpdinf$res[(slpnum+1):(slpnum*2)],cor4plotslpdinf$organic_matter[(slpnum+1):(slpnum*2)],pch=2,cex=.5,col = 'darkred')
points(cor4plotslpdinf$res[357:359],cor4plotslpdinf$organic_matter[357:359],pch=12,cex=2,col = 'darkred')

points(30.9,cor4plot$slp_1asdinf[1],pch=5,col='cyan4',cex=4); points(30.9,cor4plot$slp_1asdinf[2],pch=5,col='darkgoldenrod',cex=4)
points(30.9,cor4plot$slp_1asdinf[3],pch=5,col='darkblue',cex=4); points(30.9,cor4plot$slp_1asdinf[2],pch=5,col='darkred',cex=4)
points(30.9,cor4plot$slp_STRMdinf[1],pch=4,col='cyan4',cex=4); points(30.9,cor4plot$slp_STRMdinf[2],pch=4,col='darkgoldenrod',cex=5)
points(30.9,cor4plot$slp_STRMdinf[3],pch=4,col='darkblue',cex=4); points(30.9,cor4plot$slp_STRMdinf[2],pch=4,col='darkred',cex=4)
points(30.9,cor4plot$slp_GDEMdinf[1],pch=1,col='cyan4',cex=4); points(30.9,cor4plot$slp_GDEMdinf[2],pch=1,col='darkgoldenrod',cex=4.5)
points(30.9,cor4plot$slp_GDEMdinf[3],pch=1,col='darkblue',cex=4); points(30.9,cor4plot$slp_GDEMdinf[2],pch=1,col='darkred',cex=4)
legend(24,0.3,c("Sand","Clay","AWC pred","Org Matter","LIDAR","1m","1/3as","1as","SRTM","GDEM"),col=c("cyan4","darkgoldenrod","darkblue","darkred","black","black","black","black","black","black"),lwd=c(2,2,2,2,NA,NA,NA,NA,NA,NA),pch=c(NA,NA,NA,NA,8,2,12,5,4,1),cex=.75)

#plot for SCA dinf only
par(xpd = FALSE)
par(mfrow = c(1,1))
par(mar = c(4,4,4,1))
ylims=.5
plot(cor4plotscadinf$res[1:(scanum)],cor4plotscadinf$soil_texture_sand[1:(scanum)],main = "ln(SCA) Dinf",type='l',col = 'cyan4',lwd=1.5,ylim = c(-ylims,ylims),xlim = c(1,31),xlab = 'DEM Resolution (m)',ylab = 'Correlation')
lines(cor4plotscadinf$res[(scanum+1):(scanum*2)], cor4plotscadinf$soil_texture_sand[(scanum+1):(scanum*2)],col = 'cyan4',lwd=2)
points(cor4plotscadinf$res[1:scanum],cor4plotscadinf$soil_texture_sand[1:scanum],pch=8,cex=.75,col = 'cyan4')
points(cor4plotscadinf$res[(scanum+1):(scanum*2)],cor4plotscadinf$soil_texture_sand[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'cyan4')
points(cor4plotscadinf$res[87:89],cor4plotscadinf$soil_texture_sand[87:89],pch=12,cex=2,col = 'cyan4')
lines(cor4plotscadinf$res[1:(scanum)], cor4plotscadinf$soil_texture_clay[1:(scanum)],col = 'darkgoldenrod',lwd=2)
lines(cor4plotscadinf$res[(scanum+1):(scanum*2)], cor4plotscadinf$soil_texture_clay[(scanum+1):(scanum*2)],col = 'darkgoldenrod',lwd=2)
points(cor4plotscadinf$res[1:scanum],cor4plotscadinf$soil_texture_clay[1:scanum],pch=8,cex=.75,col = 'darkgoldenrod')
points(cor4plotscadinf$res[(scanum+1):(scanum*2)],cor4plotscadinf$soil_texture_clay[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'darkgoldenrod')
points(cor4plotscadinf$res[87:89],cor4plotscadinf$soil_texture_clay[87:89],pch=12,cex=2,col = 'darkgoldenrod')
lines(cor4plotscadinf$res[1:(scanum)], cor4plotscadinf$pred_water_capacity[1:(scanum)],col = 'darkblue',lwd=2)
lines(cor4plotscadinf$res[(scanum+1):(scanum*2)], cor4plotscadinf$pred_water_capacity[(scanum+1):(scanum*2)],col = 'darkblue',lwd=2)
points(cor4plotscadinf$res[1:scanum*2],cor4plotscadinf$pred_water_capacity[1:scanum*2],pch=8,cex=.75,col = 'darkblue')
points(cor4plotscadinf$res[(scanum+1):(scanum*2)],cor4plotscadinf$pred_water_capacity[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'darkblue')
points(cor4plotscadinf$res[87:89],cor4plotscadinf$pred_water_capacity[87:89],pch=12,cex=2,col = 'darkblue')
lines(cor4plotscadinf$res[1:(scanum)], cor4plotscadinf$organic_matter[1:(scanum)],col = 'darkred',lwd=2)
lines(cor4plotscadinf$res[(scanum+1):(scanum*2)], cor4plotscadinf$organic_matter[(scanum+1):(scanum*2)],col = 'darkred',lwd=2)
points(cor4plotscadinf$res[1:scanum],cor4plotscadinf$organic_matter[1:scanum],pch=8,cex=.75,col = 'darkred')
points(cor4plotscadinf$res[(scanum+1):(scanum*2)],cor4plotscadinf$organic_matter[(scanum+1):(scanum*2)],pch=2,cex=.5,col = 'darkred')
points(cor4plotscadinf$res[87:89],cor4plotscadinf$organic_matter[87:89],pch=12,cex=2,col = 'darkred')
points(30.9,cor4plot$lnsca_1asdinf[1],pch=5,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_1asdinf[2],pch=5,col='darkgoldenrod',cex=4)
points(30.9,cor4plot$lnsca_1asdinf[3],pch=5,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_1asdinf[2],pch=5,col='darkred',cex=4)
points(30.9,cor4plot$lnsca_SRTMdinf[1],pch=4,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_SRTMdinf[2],pch=4,col='darkgoldenrod',cex=5)
points(30.9,cor4plot$lnsca_SRTMdinf[3],pch=4,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_SRTMdinf[2],pch=4,col='darkred',cex=4)
points(30.9,cor4plot$lnsca_GDEMdinf[1],pch=1,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_GDEMdinf[2],pch=1,col='darkgoldenrod',cex=4.5)
points(30.9,cor4plot$lmsca_GDEMdinf[3],pch=1,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_GDEMdinf[2],pch=1,col='darkred',cex=4)
legend(24,0.5,c("Sand","Clay","AWC pred","Org Matter"),col=c("cyan4","darkgoldenrod","darkblue","darkred"),lwd=2,cex=.75)
legend(24,-0.25,c("LIDAR","1m","1/3as","1as","SRTM","GDEM"),col="black",pch=c(8,2,12,5,4,1),cex=.75)

plot(cor4plotscad8$res[1:(scanum)],cor4plotscad8$soil_texture_sand[1:(scanum)],main = "ln(SCA) D8",type='l',col = 'cyan4',lwd=1.5,ylim = c(-ylims,ylims),xlim = c(1,31),xlab = 'DEM Resolution (m)',ylab = 'Correlation')
lines(cor4plotscad8$res[(scanum+1):(scanum*2)], cor4plotscad8$soil_texture_sand[(scanum+1):(scanum*2)],col = 'cyan4',lwd=2)
points(cor4plotscad8$res[1:scanum],cor4plotscad8$soil_texture_sand[1:scanum],pch=8,cex=.75,col = 'cyan4')
points(cor4plotscad8$res[(scanum+1):(scanum*2)],cor4plotscad8$soil_texture_sand[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'cyan4')
points(cor4plotscad8$res[87:89],cor4plotscad8$soil_texture_sand[87:89],pch=12,cex=2,col = 'cyan4')
lines(cor4plotscad8$res[1:(scanum)], cor4plotscad8$soil_texture_clay[1:(scanum)],col = 'darkgoldenrod',lwd=2)
lines(cor4plotscad8$res[(scanum+1):(scanum*2)], cor4plotscad8$soil_texture_clay[(scanum+1):(scanum*2)],col = 'darkgoldenrod',lwd=2)
points(cor4plotscad8$res[1:scanum],cor4plotscad8$soil_texture_clay[1:scanum],pch=8,cex=.75,col = 'darkgoldenrod')
points(cor4plotscad8$res[(scanum+1):(scanum*2)],cor4plotscad8$soil_texture_clay[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'darkgoldenrod')
points(cor4plotscad8$res[87:89],cor4plotscad8$soil_texture_clay[87:89],pch=12,cex=2,col = 'darkgoldenrod')
lines(cor4plotscad8$res[1:(scanum)], cor4plotscad8$pred_water_capacity[1:(scanum)],col = 'darkblue',lwd=2)
lines(cor4plotscad8$res[(scanum+1):(scanum*2)], cor4plotscad8$pred_water_capacity[(scanum+1):(scanum*2)],col = 'darkblue',lwd=2)
points(cor4plotscad8$res[1:scanum*2],cor4plotscad8$pred_water_capacity[1:scanum*2],pch=8,cex=.75,col = 'darkblue')
points(cor4plotscad8$res[(scanum+1):(scanum*2)],cor4plotscad8$pred_water_capacity[(scanum+1):(scanum*2)],pch=2,cex=.75,col = 'darkblue')
points(cor4plotscad8$res[87:89],cor4plotscad8$pred_water_capacity[87:89],pch=12,cex=2,col = 'darkblue')
lines(cor4plotscad8$res[1:(scanum)], cor4plotscad8$organic_matter[1:(scanum)],col = 'darkred',lwd=2)
lines(cor4plotscad8$res[(scanum+1):(scanum*2)], cor4plotscad8$organic_matter[(scanum+1):(scanum*2)],col = 'darkred',lwd=2)
points(cor4plotscad8$res[1:scanum],cor4plotscad8$organic_matter[1:scanum],pch=8,cex=.75,col = 'darkred')
points(cor4plotscad8$res[(scanum+1):(scanum*2)],cor4plotscad8$organic_matter[(scanum+1):(scanum*2)],pch=2,cex=.5,col = 'darkred')
points(cor4plotscad8$res[87:89],cor4plotscad8$organic_matter[87:89],pch=12,cex=2,col = 'darkred')
points(30.9,cor4plot$lnsca_1asd8[1],pch=5,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_1asd8[2],pch=5,col='darkgoldenrod',cex=4)
points(30.9,cor4plot$lnsca_1asd8[3],pch=5,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_1asd8[2],pch=5,col='darkred',cex=4)
points(30.9,cor4plot$lnsca_SRTMd8[1],pch=4,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_SRTMd8[2],pch=4,col='darkgoldenrod',cex=5)
points(30.9,cor4plot$lnsca_SRTMd8[3],pch=4,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_SRTMd8[2],pch=4,col='darkred',cex=4)
points(30.9,cor4plot$lnsca_GDEMd8[1],pch=1,col='cyan4',cex=4); points(30.9,cor4plot$lnsca_GDEMd8[2],pch=1,col='darkgoldenrod',cex=4.5)
points(30.9,cor4plot$lmsca_GDEMd8[3],pch=1,col='darkblue',cex=4); points(30.9,cor4plot$lnsca_GDEMd8[2],pch=1,col='darkred',cex=4)
legend(24,0.5,c("Sand","Clay","AWC pred","Org Matter"),col=c("cyan4","darkgoldenrod","darkblue","darkred"),lwd=2,cex=.75)
legend(24,-0.25,c("LIDAR","1m","1/3as","1as","SRTM","GDEM"),col="black",pch=c(8,2,12,5,4,1),cex=.75)




#### histograms for sca and slp
# library
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(forcats)
library(gridExtra)
##################################
slpdinfhist_LI = slpdinf_LI[,1:43]
names(slpdinfhist_LI) = paste0(.7*1:43,"m")
slpdinfhist_1m = slpdinf_1m[,1:30]
names(slpdinfhist_1m) = paste0(1*1:30,".0m")
names(slpdinfhist_1m) = paste(c("ab","ad","ag","ai","al","an","aq","as","au","ax",
                                "ba","bd","bf","bi","bk","bm","bp","br","bu","bw",
                                "ca","cc","ce","ch","cj","cm","co","cr","ct","cv"),
                              names(slpdinfhist_1m))
slpdinfhist_13as = slpdinf_13as[,1:3]
names(slpdinfhist_13as) = paste0(10.3*1:3,"m")
names(slpdinfhist_13as) = paste(c("ay","by","cx"),
                              names(slpdinfhist_13as))
names(slpdinfhist_LI)[10] = "7.0m"
names(slpdinfhist_LI)[20] = "14.0m"
names(slpdinfhist_LI)[30] = "21.0m"
names(slpdinfhist_LI)[40] = "28.0m"
names(slpdinfhist_LI) = paste(c("aa","ac","ae","af","ah","aj","ak","am","ao","ap","ar","at","av","aw","az",
                                "bb","bc","be","bg","bh","bj","bl","bn","bo","bq","bs","bt","bv","bx","bz",
                                "cb","cd","cf","cg","ci","ck","cl","cn","cp","cq","cs","cu","cw"),
                              names(slpdinfhist_LI))
for(i in 1:length(names(slpdinfhist_LI))){
  df=as.data.frame(slpdinfhist_LI[,i])
  colnames(df)="Slope"
  df$Resolution=names(slpdinfhist_LI)[i]
  df$orig = "LIDAR"
  assign(paste0("Slope_LI",i),df)
}
for(i in 1:length(names(slpdinfhist_1m))){
  df=as.data.frame(slpdinfhist_1m[,i])
  colnames(df)="Slope"
  df$Resolution=names(slpdinfhist_1m)[i]
  df$orig = "1m"
  assign(paste0("Slope_1m",i),df)
}
for(i in 1:length(names(slpdinfhist_13as))){
  df=as.data.frame(slpdinfhist_13as[,i])
  colnames(df)="Slope"
  df$Resolution=names(slpdinfhist_13as)[i]
  df$orig = "1/3as"
  assign(paste0("Slope_13as",i),df)
}

NewDfslp=do.call("rbind", list(Slope_LI1,Slope_1m1,Slope_LI2,Slope_1m2,Slope_LI3,Slope_LI4,Slope_1m3,Slope_LI5,Slope_1m4,Slope_LI6,Slope_LI7,Slope_1m5,
                               Slope_LI8,Slope_1m6,Slope_LI9,Slope_LI10,Slope_1m7,Slope_LI11,Slope_1m8,Slope_LI12,Slope_1m9,Slope_LI13,Slope_LI14,Slope_1m10,Slope_13as1,
                               Slope_LI15,Slope_1m11,Slope_LI16,Slope_LI17,Slope_1m12,Slope_LI18,Slope_1m13,Slope_LI19,Slope_LI20,Slope_1m14,Slope_LI21,Slope_1m15,Slope_LI22,
                               Slope_1m16,Slope_LI23,Slope_LI24,Slope_1m17,Slope_LI25,Slope_1m18,Slope_LI26,Slope_LI27,Slope_1m19,Slope_LI28,Slope_1m20,Slope_LI29,Slope_13as2,
                               Slope_LI30,Slope_1m21,Slope_LI31,Slope_1m22,Slope_LI32,Slope_1m23,Slope_LI33,Slope_LI34,Slope_1m24,Slope_LI35,Slope_1m25,Slope_LI36,Slope_LI37,
                               Slope_1m26,Slope_LI38,Slope_1m27,Slope_LI39,Slope_LI40,Slope_1m28,Slope_LI41,Slope_1m29,Slope_LI42,Slope_1m30,Slope_LI43,Slope_13as3))

names(NewDfslp)[3] = "Source DEM"
#NewDfslp = stack(list(NewDfslp))
p <- ggplot(NewDfslp, aes(Resolution, y=Slope,col=`Source DEM`)) + 
  geom_violin()+ylim(0,.4)+theme(axis.text.x = element_text(angle = 90))
p

##################################
scadinfhist_LI = log(scadinf_LI)[,1:15]
names(scadinfhist_LI) = paste0(.7*1:15,"m")
scadinfhist_1m = log(scadinf_1m)[,1:10]
names(scadinfhist_1m) = paste0(1*1:10,".0m")
scadinfhist_13as = log(scadinf_13as)[,1]
names(scadinfhist_LI)[10] = "7.0m"
names(scadinfhist_LI)[15] = "z. 10.5m"
names(scadinfhist_1m)[10] = "x. 10.0m"

for(i in 1:length(names(scadinfhist_LI))){
  df=as.data.frame(scadinfhist_LI[,i])
  colnames(df)="SCA"
  df$Resolution=names(scadinfhist_LI)[i]
  df$orig = "LIDAR"
  assign(paste0("SCA_LI",i),df)
}
for(i in 1:length(names(scadinfhist_1m))){
  df=as.data.frame(scadinfhist_1m[,i])
  colnames(df)="SCA"
  df$Resolution=names(scadinfhist_1m)[i]
  df$orig = "1m"
  assign(paste0("SCA_1m",i),df)
}
SCA_13as1=data.frame(SCA = unlist(scadinfhist_13as), Resolution = "y. 10.3m",orig = "13as")

NewDfsca=do.call("rbind", list(SCA_LI1,SCA_1m1,SCA_LI2,SCA_1m2,SCA_LI3,SCA_LI4,SCA_1m3,SCA_LI5,SCA_1m4,SCA_LI6,SCA_LI7,SCA_1m5,
                               SCA_LI8,SCA_1m6,SCA_LI9,SCA_LI10,SCA_1m7,SCA_LI11,SCA_1m8,SCA_LI12,SCA_1m9,
                               SCA_LI13,SCA_LI14,SCA_1m10,SCA_13as1,SCA_LI15))

names(NewDfsca)[3] = "Source DEM"
p <- ggplot(NewDfsca, aes(Resolution, y=SCA,col=`Source DEM`)) + ggtitle("ln(SCA) Dinf")+
  geom_violin()+theme(axis.text.x = element_text(angle = 90))
p

##################################