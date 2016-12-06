##################################################################################################
# devtools::install_github("hunzikp/velox")
######################################################################################################

library(rgdal)
library(raster)
library(gdalUtils)
library(tiff)
library(leaflet)
library(rgeos)
library(velox)
library(mapview)
library(RColorBrewer)
library(classInt)
library(mgcv)
library(XLConnect)

######################################################################################################
setwd("/home/alf/Scrivania/aaa_lavori/lav_new_parma")

source("aux_aster.r")

######################################################################################################
# load geo data 

confini_comune_ls=readRDS("confini_comune_landsat.rds")
parma_urbanizzato_ls=readRDS("parma_urbanizzato_ls.rds")
fabbricati_civilires_all=readRDS("fabbricati_civilires_all.rds")
fabbricati_civilires_centroid=readRDS("fabbricati_civilires_centroid.rds")

area_urban=over(fabbricati_civilires_centroid,parma_urbanizzato_ls)[,1]

consumo_suolo_2015=raster("parma_consumo_suolo_2015.tif")
consumo_suolo_2012=raster("parma_consumo_suolo_2012.tif")

consumo_suolo_2015_vx <- velox(consumo_suolo_2015)
consumo_suolo_2012_vx <- velox(consumo_suolo_2012)
fabbricati_civilires_centroid_200 <- gBuffer(fabbricati_civilires_centroid, width=100, byid=TRUE)
fabbricati_civilires_centroid_1ha <- gBuffer(fabbricati_civilires_centroid, width=57, byid=TRUE)

fabbricati_civilires_centroid_100 <- SpatialPolygonsDataFrame(fabbricati_civilires_centroid_200, data.frame(id=1:length(fabbricati_civilires_centroid_200)), FALSE)
fabbricati_civilires_centroid_1ha<- SpatialPolygonsDataFrame(fabbricati_civilires_centroid_1ha, data.frame(id=1:length(fabbricati_civilires_centroid_1ha)), FALSE)

cs100_fabbricati=consumo_suolo_2015_vx$extract(fabbricati_civilires_centroid_200,fun=function(x) sum(x,na.rm=T))
cs1ha_fabbricati=consumo_suolo_2015_vx$extract(fabbricati_civilires_centroid_1ha,fun=function(x) sum(x,na.rm=T))

####################################################################################################
# Calcolo lst con emissività da tetti ad immagine aster 90 m TIRS 13 e 14

setwd("/home/alf/Scrivania/aaa_lavori/lav_new_parma/night")

file_list <- list.files(pattern = 'AST_L1T_.*hdf$')
file_name=gsub(".hdf","",file_list)
file_name_night=gsub(".hdf","",file_list)

b13_night=list()
b14_night=list()
b13_parma_night=list()
b14_parma_night=list()
b13_parma_tb_night=list()
b14_parma_tb_night=list()
b13_parma_lst_night_E90=list()
b14_parma_lst_night_E90=list()
b13_parma_lst_night_E95=list()
b14_parma_lst_night_E95=list()


for ( i in 1:length(file_list)) {

work=read_metadata_ASTER(file_list[i])

b13_night[[i]]=readband_ASTER(file_list[i],index_k=4,utm_zone=work$utm_zone,geoextent=work$raster_dims_90m)
b14_night[[i]]=readband_ASTER(file_list[i],index_k=5,utm_zone=work$utm_zone,geoextent=work$raster_dims_90m)

b13_parma_night[[i]]=crop(b13_night[[i]],extent(confini_comune_ls))
b14_parma_night[[i]]=crop(b14_night[[i]],extent(confini_comune_ls))

b13_parma_tb_night[[i]]=calc_tbright(b13_parma_night[[i]],klist_ASTER$B13$k1,klist_ASTER$B13$k2,klist_ASTER$B13$unitconversion)
b14_parma_tb_night[[i]]=calc_tbright(b14_parma_night[[i]],klist_ASTER$B14$k1,klist_ASTER$B14$k2,klist_ASTER$B14$unitconversion)

b13_parma_lst_night_E90[[i]]=calc_lst(b13_parma_tb_night[[i]],emis=0.90,Lwave=klist_ASTER$B13$wavelength)
b14_parma_lst_night_E90[[i]]=calc_lst(b14_parma_tb_night[[i]],emis=0.90,Lwave=klist_ASTER$B14$wavelength)

b13_parma_lst_night_E95[[i]]=calc_lst(b13_parma_tb_night[[i]],emis=0.95,Lwave=klist_ASTER$B13$wavelength)
b14_parma_lst_night_E95[[i]]=calc_lst(b14_parma_tb_night[[i]],emis=0.95,Lwave=klist_ASTER$B14$wavelength)

writeRaster(mask(b13_parma_lst_night_E90[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b13_PR_lstnight_E90.tif"), overwrite=TRUE)
writeRaster(mask(b13_parma_lst_night_E95[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b13_PR_lstnight_E95.tif"), overwrite=TRUE)
writeRaster(mask(b14_parma_lst_night_E90[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b14_PR_lstnight_E90.tif"), overwrite=TRUE)
writeRaster(mask(b14_parma_lst_night_E95[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b14_PR_lstnight_E95.tif"), overwrite=TRUE)

}


###########################################################################################################################
setwd("/home/alf/Scrivania/aaa_lavori/lav_new_parma/day")

#source("ASTERL1T_DN2REF.R")


file_list <- list.files(pattern = 'AST_L1T_.*hdf$')

file_name=gsub(".hdf","",file_list)
file_name_day=gsub(".hdf","",file_list)

b13_day=list()
b14_day=list()
b13_parma_day=list()
b14_parma_day=list()
b13_parma_tb_day=list()
b14_parma_tb_day=list()
b13_parma_lst_day_E90=list()
b14_parma_lst_day_E90=list()
b13_parma_lst_day_E95=list()
b14_parma_lst_day_E95=list()
blue_rd_day=list()
blue_rf_day=list()
green_rd_day=list()
green_rf_day=list()
red_rd_day=list()
red_rf_day=list()

for ( i in 1:length(file_list)) {

work=read_metadata_ASTER(file_list[i])


b13_day[[i]]=readband_ASTER(file_list[i],index_k=4,utm_zone=work$utm_zone,geoextent=work$raster_dims_90m)
b14_day[[i]]=readband_ASTER(file_list[i],index_k=5,utm_zone=work$utm_zone,geoextent=work$raster_dims_90m)

b13_parma_day[[i]]=crop(b13_day[[i]],extent(confini_comune_ls))
b14_parma_day[[i]]=crop(b14_day[[i]],extent(confini_comune_ls))

b13_parma_tb_day[[i]]=calc_tbright(b13_parma_day[[i]],klist_ASTER$B13$k1,klist_ASTER$B13$k2,klist_ASTER$B13$unitconversion)
b14_parma_tb_day[[i]]=calc_tbright(b14_parma_day[[i]],klist_ASTER$B14$k1,klist_ASTER$B14$k2,klist_ASTER$B14$unitconversion)

b13_parma_lst_day_E90[[i]]=calc_lst(b13_parma_tb_day[[i]],emis=0.90,Lwave=klist_ASTER$B13$wavelength)
b14_parma_lst_day_E90[[i]]=calc_lst(b14_parma_tb_day[[i]],emis=0.90,Lwave=klist_ASTER$B14$wavelength)

b13_parma_lst_day_E95[[i]]=calc_lst(b13_parma_tb_day[[i]],emis=0.95,Lwave=klist_ASTER$B13$wavelength)
b14_parma_lst_day_E95[[i]]=calc_lst(b14_parma_tb_day[[i]],emis=0.95,Lwave=klist_ASTER$B14$wavelength)

writeRaster(mask(b13_parma_lst_day_E90[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b13_PR_lstday_E90.tif"), overwrite=TRUE)
writeRaster(mask(b13_parma_lst_day_E95[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b13_PR_lstday_E95.tif"), overwrite=TRUE)
writeRaster(mask(b14_parma_lst_day_E90[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b14_PR_lstday_E90.tif"), overwrite=TRUE)
writeRaster(mask(b14_parma_lst_day_E95[[i]],confini_comune_ls), filename=paste0(file_name[[i]],"_b14_PR_lstday_E95.tif"), overwrite=TRUE)

blue_rd_day[[i]]=raster(paste0(file_name[[i]],"_ImageData1_radiance.tif"))
blue_rf_day[[i]]=raster(paste0(file_name[[i]],"_ImageData1_reflectance.tif"))
green_rd_day[[i]]=raster(paste0(file_name[[i]],"_ImageData2_radiance.tif"))
green_rf_day[[i]]=raster(paste0(file_name[[i]],"_ImageData2_reflectance.tif"))
red_rd_day[[i]]=raster(paste0(file_name[[i]],"_ImageData3N_radiance.tif"))
red_rf_day[[i]]=raster(paste0(file_name[[i]],"_ImageData3N_reflectance.tif"))

}

################################################################

setwd("/home/alf/Scrivania/aaa_lavori/lav_new_parma")

lst_fabbricati_b14_night_E90=list()
lst_fabbricati_b13_night_E95=list()
lst_fabbricati_b14_night_E95=list()
lst_fabbricati_b13_night_E90=list()


lst_fabbricati_b13_day_E90=list()
lst_fabbricati_b14_day_E90=list()
lst_fabbricati_b13_day_E95=list()
lst_fabbricati_b14_day_E95=list()

ed_B_rd=list()
ed_B_rf=list()
ed_G_rd=list()
ed_G_rf=list()
ed_R_rd=list()
ed_R_rf=list()

for ( i in 1:length(file_name_day)) {
  
  lst_fabbricati_b13_day_E90[[i]]=raster::extract(b13_parma_lst_day_E90[[i]],fabbricati_civilires_centroid)
  
  lst_fabbricati_b14_day_E90[[i]]=raster::extract(b14_parma_lst_day_E90[[i]],fabbricati_civilires_centroid)
  
  lst_fabbricati_b13_day_E95[[i]]=raster::extract(b13_parma_lst_day_E95[[i]],fabbricati_civilires_centroid)
  
  lst_fabbricati_b14_day_E95[[i]]=raster::extract(b14_parma_lst_day_E95[[i]],fabbricati_civilires_centroid)
  
  ed_B_rd[[i]]=raster::extract(blue_rd_day[[i]],fabbricati_civilires_centroid)
  ed_B_rf[[i]]=raster::extract(blue_rf_day[[i]],fabbricati_civilires_centroid)
  ed_G_rd[[i]]=raster::extract(green_rd_day[[i]],fabbricati_civilires_centroid)
  ed_G_rf[[i]]=raster::extract(green_rf_day[[i]],fabbricati_civilires_centroid)
  ed_R_rd[[i]]=raster::extract(red_rd_day[[i]],fabbricati_civilires_centroid)
  ed_R_rf[[i]]=raster::extract(red_rf_day[[i]],fabbricati_civilires_centroid)
  
}



for ( i in 1:length(file_name_night)) {

lst_fabbricati_b13_night_E90[[i]]=raster::extract(b13_parma_lst_night_E90[[i]],fabbricati_civilires_centroid)

lst_fabbricati_b14_night_E90[[i]]=raster::extract(b14_parma_lst_night_E90[[i]],fabbricati_civilires_centroid)

lst_fabbricati_b13_night_E95[[i]]=raster::extract(b13_parma_lst_night_E95[[i]],fabbricati_civilires_centroid)

lst_fabbricati_b14_night_E95[[i]]=raster::extract(b14_parma_lst_night_E95[[i]],fabbricati_civilires_centroid)


}



for ( i in seq_along(file_name_day)) {


res_lst=data.frame(id=fabbricati_civilires_all$FABCD,
                   area_urban=area_urban,
                   cs100_fabbricati=cs100_fabbricati,
                   cs57_fabbricati=cs1ha_fabbricati,
                   cbind(lst_fabbricati_b13_day_E90[[i]],
                         lst_fabbricati_b13_day_E95[[i]],
                         lst_fabbricati_b14_day_E90[[i]],
                         lst_fabbricati_b14_day_E95[[i]],
                         ed_B_rd[[i]],
                         ed_B_rf[[i]],
                         ed_G_rd[[i]],
                         ed_G_rf[[i]],
                         ed_R_rd[[i]],
                         ed_R_rf[[i]]))

names(res_lst)=c("id","area_urban","cs100_fabbricati","cs57_fabbricati","B13_90","B13_95","B14_90","B14_95","rad_B","rf_B","rad_G","rf_G","rad_R","rf_R")
saveRDS(res_lst,paste0(file_name_day[i],"_lst_parma_day.rds"))

XLConnect::writeWorksheetToFile(paste0(file_name_day[i],"_lst_day.xls"),res_lst,sheet="data")

}

for ( i in seq_along(file_name_night)) {

res_lst=data.frame(id=fabbricati_civilires_all$FABCD,
                   cbind(lst_fabbricati_b13_night_E90[[i]],
                         lst_fabbricati_b13_night_E95[[i]],
                         lst_fabbricati_b14_night_E90[[i]],
                         lst_fabbricati_b14_night_E95[[i]]))
names(res_lst)==c("id","B13_90","B13_95","B14_90","B14_95")

saveRDS(res_lst,file=paste0(file_name_night[i],"_lst_parma_night.rds"))
                           
XLConnect::writeWorksheetToFile(paste0(file_name_night[i],"_lst_night.xls"),res_lst,sheet="data")
}

###########################################################################################################################
consumo_suolo_2015=raster("parma_consumo_suolo_2015.tif")
consumo_suolo_2012=raster("parma_consumo_suolo_2012.tif")

consumo_suolo_2015_vx <- velox(consumo_suolo_2015)
consumo_suolo_2012_vx <- velox(consumo_suolo_2012)

fabbricati_civilires_centroid_200 <- gBuffer(fabbricati_civilires_centroid, width=100, byid=TRUE)
fabbricati_civilires_centroid_1ha <- gBuffer(fabbricati_civilires_centroid, width=57, byid=TRUE)

fabbricati_civilires_centroid_100 <- SpatialPolygonsDataFrame(fabbricati_civilires_centroid_200, data.frame(id=1:length(fabbricati_civilires_centroid_200)), FALSE)
fabbricati_civilires_centroid_1ha<- SpatialPolygonsDataFrame(fabbricati_civilires_centroid_1ha, data.frame(id=1:length(fabbricati_civilires_centroid_1ha)), FALSE)

cs100_fabbricati=consumo_suolo_2015_vx$extract(fabbricati_civilires_centroid_200,fun=function(x) sum(x,na.rm=T))
cs1ha_fabbricati=consumo_suolo_2015_vx$extract(fabbricati_civilires_centroid_1ha,fun=function(x) sum(x,na.rm=T))

###########################################################################################################################

banda_vnir_aster=raster("band_merge.tif")
png("aster_image_parma.png")
plotRGB(stack(banda_vnir_aster),r=3,g=2,b=1, scale=800, stretch = "Lin")
dev.off()
###########################################################################################################################

landsat_tb=raster("parma_lst_2015.tif")
landsat_parma_lst=calc_lst(landsat_tb+273.15,emis=0.93,Lwave=10.9)


###########################################################################################################################

parma_urbanizzato=parma_urbanizzato_ls[parma_urbanizzato_ls$UUR_1=="urbanizzato",]
parma_rurale=parma_urbanizzato_ls[parma_urbanizzato_ls$UUR_1=="rurale",]
parma_parco=parma_urbanizzato_ls[parma_urbanizzato_ls$UUR_1=="Parco",]
parma_urbanizzabile=parma_urbanizzato_ls[parma_urbanizzato_ls$UUR_1=="urbanizzabile",]

###########################################################################################################################

writeOGR(parma_urbanizzato,".","parma_urbanizzato",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(parma_rurale,".","parma_rurale",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(parma_parco,".","parma_parco",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(parma_urbanizzabile,".","parma_urbanizzabile",driver="ESRI Shapefile",overwrite_layer = T)

###########################################################################################################################à
# estract polygon con velox and point con raster::extract


lst_fabbricati_b13=lst_fabbricati_b13_day_E90[[1]]
lst_fabbricati_b14=lst_fabbricati_b14_day_E90[[1]]
lst_fabbricati_b13_night_E90=lst_fabbricati_b13_night_E90[[1]]
lst_fabbricati_b14=lst_fabbricati_b14_night_E90[[1]]

lst_fabbricati_landsat=extract(landsat_parma_lst,
                               fabbricati_civilires_centroid,
                               fun=function(x) mean(x,na.rm=T))


###########################################################################################################################à
# create layers

fabbricati_civilires_all$cs200=cs100_fabbricati
fabbricati_civilires_all$cs100=cs100_fabbricati
fabbricati_civilires_all$cs57=cs1ha_fabbricati

fabbricati_civilires_all$lst_b13=lst_fabbricati_b13
fabbricati_civilires_all$lst_b14=lst_fabbricati_b14

fabbricati_civilires_all$lst_landsat=lst_fabbricati_landsat


###########################################################################################################################à
# load censuary data

table_2014=readRDS("table_2014_parma_class.rds")
table_2015=readRDS("table_2015_parma_class.rds")

table_2014$perc_u_5_14=t(apply(table_2014[c(1,2)],1,function(x) sum(x,na.rm=T))/table_2014$sum2014)
table_2014$perc_o_65=apply(table_2014[c(9,10)],1,function(x) sum(x,na.rm=T))/table_2014$sum2014
table_2014$perc_o_75=apply(table_2014[c(10)],1,function(x) sum(x,na.rm=T))/table_2014$sum2014

table_2015$perc_u_5=apply(table_2015[c(1,2)],1,function(x) sum(x,na.rm=T))/table_2015$sum2015
table_2015$perc_o_65=apply(table_2015[c(9,10)],1,function(x) sum(x,na.rm=T))/table_2015$sum2015
table_2015$perc_o_75=apply(table_2015[c(10)],1,function(x) sum(x,na.rm=T))/table_2015$sum2015

saveRDS(table_2014,"table_2014_parma_class.rds")
saveRDS(table_2015,"table_2015_parma_class.rds")

table_2014=readRDS("table_2014_parma_class.rds")
table_2015=readRDS("table_2015_parma_class.rds")

XLConnect::writeWorksheetToFile("table_2014.xls",table_2014,sheet="table_2014")
XLConnect::writeWorksheetToFile("table_2015.xls",table_2015,sheet="table_2015")
###########################################################################################################################à


fabbricati_civilires_all$"p_u_5_14"=NA   
fabbricati_civilires_all$"p_o_65_14"=NA   
fabbricati_civilires_all$"p_o_75_14"=NA   
fabbricati_civilires_all$"p_u_5_15"=NA   
fabbricati_civilires_all$"p_o_65_15"=NA  
fabbricati_civilires_all$"p_o_75_15"=NA

for ( i in 1:length(table_2014$SEZCENS))
  {
  
  temp=which(fabbricati_civilires_all$SEZ_civici==table_2014$SEZCENS[i])
  fabbricati_civilires_all$"p_u_5_14"[temp]=table_2014$perc_u_5[i]   
  fabbricati_civilires_all$"p_o_65_14"[temp]=table_2014$perc_o_65[i]    
  fabbricati_civilires_all$"p_o_75_14"[temp]=table_2014$perc_o_75[i]    
  
  }
  
for ( i in 1:length(table_2015$SEZCENS)){
  
  temp=which(fabbricati_civilires_all$SEZ_civici==table_2015$SEZCENS[i])
  
  fabbricati_civilires_all$"p_u_5_15"[temp]=table_2015$perc_u_5[i]   
  fabbricati_civilires_all$"p_o_65_15"[temp]=table_2015$perc_o_65[i]    
  fabbricati_civilires_all$"p_o_75_15"[temp]=table_2015$perc_o_75[i]    
  
}  

saveRDS(fabbricati_civilires_all,"fabbricati_civilires_all.rds")

########################################################################################à
# Analisi del rischio

fabbricati_civilires_all=readRDS("fabbricati_civilires_all.rds")

fabbricati_civilires_all$cs200_scaled=as.numeric(scale(fabbricati_civilires_all$cs200,center=min(fabbricati_civilires_all$cs200,na.rm=T),scale=diff(range(fabbricati_civilires_all$cs200))))

range_lst_b13=as.numeric(c(quantile(fabbricati_civilires_all@data$lst_b13,0,na.rm=T),quantile(fabbricati_civilires_all@data$lst_b13,0.99,na.rm=T)))
fabbricati_civilires_all$lst_b13_scaled=as.numeric(scale(fabbricati_civilires_all$lst_b13,center=min(fabbricati_civilires_all$lst_b13,na.rm=T),scale=diff(range_lst_b13)))


range_lst_b14=as.numeric(c(quantile(fabbricati_civilires_all@data$lst_b14,0,na.rm=T),quantile(fabbricati_civilires_all@data$lst_b14,0.99,na.rm=T)))
fabbricati_civilires_all$lst_b14_scaled=as.numeric(scale(fabbricati_civilires_all$lst_b14,center=min(fabbricati_civilires_all$lst_b14,na.rm=T),scale=diff(range_lst_b14)))

range_lst_landsat=as.numeric(c(quantile(fabbricati_civilires_all@data$lst_landsat,0,na.rm=T),quantile(fabbricati_civilires_all@data$lst_landsat,0.99,na.rm=T)))
fabbricati_civilires_all$lst_landsat_scaled=as.numeric(scale(fabbricati_civilires_all$lst_landsat,center=min(fabbricati_civilires_all$lst_landsat,na.rm=T),scale=diff(range_lst_landsat)))

range_RES15=as.numeric(c(0,quantile(fabbricati_civilires_all@data$RES15,0.95,na.rm=T)))
fabbricati_civilires_all$RES15_scaled=as.numeric(scale(fabbricati_civilires_all$RES15,center=min(fabbricati_civilires_all$RES15,na.rm=T),scale=diff(range_RES15)))
fabbricati_civilires_all$RES15_scaled[which(fabbricati_civilires_all$RES15_scaled>1)]=1

range_RES14=as.numeric(c(0,quantile(fabbricati_civilires_all@data$RES14,0.95,na.rm=T)))
fabbricati_civilires_all$RES14_scaled=as.numeric(scale(fabbricati_civilires_all$RES14,center=min(fabbricati_civilires_all$RES14,na.rm=T),scale=diff(range_RES14)))
fabbricati_civilires_all$RES14_scaled[which(fabbricati_civilires_all$RES14_scaled>1)]=1

##############################################################################################################################################################

fabbricati_civilires_all$risk_u5_b13=(0.5*fabbricati_civilires_all$lst_b13_scaled)+((1/3)*fabbricati_civilires_all$p_u_5_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o65_b13=(0.5*fabbricati_civilires_all$lst_b13_scaled)+((1/3)*fabbricati_civilires_all$p_o_65_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o75_b13=(0.5*fabbricati_civilires_all$lst_b13_scaled)+((1/3)*fabbricati_civilires_all$p_o_75_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;

fabbricati_civilires_all$risk_u5_b14=(0.5*fabbricati_civilires_all$lst_b14_scaled)+((1/3)*fabbricati_civilires_all$p_u_5_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o65_b14=(0.5*fabbricati_civilires_all$lst_b14_scaled)+((1/3)*fabbricati_civilires_all$p_o_65_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o75_b14=(0.5*fabbricati_civilires_all$lst_b14_scaled)+((1/3)*fabbricati_civilires_all$p_o_75_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;

fabbricati_civilires_all$risk_u5_landsat=(0.5*fabbricati_civilires_all$lst_landsat_scaled)+((1/3)*fabbricati_civilires_all$p_u_5_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o65_landsat=(0.5*fabbricati_civilires_all$lst_landsat_scaled)+((1/3)*fabbricati_civilires_all$p_o_65_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;
fabbricati_civilires_all$risk_o75_landsat=(0.5*fabbricati_civilires_all$lst_landsat_scaled)+((1/3)*fabbricati_civilires_all$p_o_75_14+(1/3)*fabbricati_civilires_all$RES14_scaled+(1/3)*fabbricati_civilires_all$cs200_scaled)*0.5;

indna_risk_u5_b13=which(!is.na(fabbricati_civilires_all$risk_u5_b13))
temp_risk_u5_b13=classIntervals(fabbricati_civilires_all$risk_u5_b13[indna_risk_u5_b13], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_u5_b13_factor=cut(fabbricati_civilires_all$risk_u5_b13[indna_risk_u5_b13], breaks = temp_risk_u5_b13$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ru5_b13_N=NA
fabbricati_civilires_all@data$Ru5_b13_L=NA
fabbricati_civilires_all@data$Ru5_b13_L[indna_risk_u5_b13]=as.character(temp_class_risk_u5_b13_factor)
fabbricati_civilires_all@data$Ru5_b13_N[indna_risk_u5_b13]=as.numeric(temp_class_risk_u5_b13_factor)

indna_risk_u5_b14=which(!is.na(fabbricati_civilires_all$risk_u5_b14))
temp_risk_u5_b14=classIntervals(fabbricati_civilires_all$risk_u5_b14[indna_risk_u5_b14], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_u5_b14_factor=cut(fabbricati_civilires_all$risk_u5_b14[indna_risk_u5_b14], breaks = temp_risk_u5_b14$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ru5_b14_N=NA
fabbricati_civilires_all@data$Ru5_b14_L=NA
fabbricati_civilires_all@data$Ru5_b14_L[indna_risk_u5_b14]=as.character(temp_class_risk_u5_b14_factor)
fabbricati_civilires_all@data$Ru5_b14_N[indna_risk_u5_b14]=as.numeric(temp_class_risk_u5_b14_factor)

indna_risk_u5_landsat=which(!is.na(fabbricati_civilires_all$risk_u5_landsat))
temp_risk_u5_landsat=classIntervals(fabbricati_civilires_all$risk_u5_landsat[indna_risk_u5_landsat], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_u5_landsat_factor=cut(fabbricati_civilires_all$risk_u5_landsat[indna_risk_u5_landsat], breaks = temp_risk_u5_landsat$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ru5_landsat_N=NA
fabbricati_civilires_all@data$Ru5_landsat_L=NA
fabbricati_civilires_all@data$Ru5_landsat_L[indna_risk_u5_landsat]=as.character(temp_class_risk_u5_landsat_factor)
fabbricati_civilires_all@data$Ru5_landsat_N[indna_risk_u5_landsat]=as.numeric(temp_class_risk_u5_landsat_factor)


indna_risk_o65_b13=which(!is.na(fabbricati_civilires_all$risk_o65_b13))
temp_risk_o65_b13=classIntervals(fabbricati_civilires_all$risk_o65_b13[indna_risk_o65_b13], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o65_b13_factor=cut(fabbricati_civilires_all$risk_o65_b13[indna_risk_o65_b13], breaks = temp_risk_o65_b13$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro65_b13_N=NA
fabbricati_civilires_all@data$Ro65_b13_L=NA
fabbricati_civilires_all@data$Ro65_b13_L[indna_risk_o65_b13]=as.character(temp_class_risk_o65_b13_factor)
fabbricati_civilires_all@data$Ro65_b13_N[indna_risk_o65_b13]=as.numeric(temp_class_risk_o65_b13_factor)

indna_risk_o65_b14=which(!is.na(fabbricati_civilires_all$risk_o65_b14))
temp_risk_o65_b14=classIntervals(fabbricati_civilires_all$risk_o65_b14[indna_risk_o65_b14], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o65_b14_factor=cut(fabbricati_civilires_all$risk_o65_b14[indna_risk_o65_b14], breaks = temp_risk_o65_b14$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro65_b14_N=NA
fabbricati_civilires_all@data$Ro65_b14_L=NA
fabbricati_civilires_all@data$Ro65_b14_L[indna_risk_o65_b14]=as.character(temp_class_risk_o65_b14_factor)
fabbricati_civilires_all@data$Ro65_b14_N[indna_risk_o65_b14]=as.numeric(temp_class_risk_o65_b14_factor)

indna_risk_o65_landsat=which(!is.na(fabbricati_civilires_all$risk_o65_landsat))
temp_risk_o65_landsat=classIntervals(fabbricati_civilires_all$risk_o65_landsat[indna_risk_o65_landsat], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o65_landsat_factor=cut(fabbricati_civilires_all$risk_o65_landsat[indna_risk_o65_landsat], breaks = temp_risk_o65_landsat$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro65_landsat_N=NA
fabbricati_civilires_all@data$Ro65_landsat_L=NA
fabbricati_civilires_all@data$Ro65_landsat_L[indna_risk_o65_landsat]=as.character(temp_class_risk_o65_landsat_factor)
fabbricati_civilires_all@data$Ro65_landsat_N[indna_risk_o65_landsat]=as.numeric(temp_class_risk_o65_landsat_factor)

indna_risk_o75_b13=which(!is.na(fabbricati_civilires_all$risk_o75_b13))
temp_risk_o75_b13=classIntervals(fabbricati_civilires_all$risk_o75_b13[indna_risk_o75_b13], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o75_b13_factor=cut(fabbricati_civilires_all$risk_o75_b13[indna_risk_o75_b13], breaks = temp_risk_o75_b13$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro75_b13_N=NA
fabbricati_civilires_all@data$Ro75_b13_L=NA
fabbricati_civilires_all@data$Ro75_b13_L[indna_risk_o75_b13]=as.character(temp_class_risk_o75_b13_factor)
fabbricati_civilires_all@data$Ro75_b13_N[indna_risk_o75_b13]=as.numeric(temp_class_risk_o75_b13_factor)

indna_risk_o75_b14=which(!is.na(fabbricati_civilires_all$risk_o75_b14))
temp_risk_o75_b14=classIntervals(fabbricati_civilires_all$risk_o75_b14[indna_risk_o75_b14], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o75_b14_factor=cut(fabbricati_civilires_all$risk_o75_b14[indna_risk_o75_b14], breaks = temp_risk_o75_b14$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro75_b14_N=NA
fabbricati_civilires_all@data$Ro75_b14_L=NA
fabbricati_civilires_all@data$Ro75_b14_L[indna_risk_o75_b14]=as.character(temp_class_risk_o75_b14_factor)
fabbricati_civilires_all@data$Ro75_b14_N[indna_risk_o75_b14]=as.numeric(temp_class_risk_o75_b14_factor)

indna_risk_o75_landsat=which(!is.na(fabbricati_civilires_all$risk_o75_landsat))
temp_risk_o75_landsat=classIntervals(fabbricati_civilires_all$risk_o75_landsat[indna_risk_o75_landsat], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_o75_landsat_factor=cut(fabbricati_civilires_all$risk_o75_landsat[indna_risk_o75_landsat], breaks = temp_risk_o75_landsat$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
fabbricati_civilires_all@data$Ro75_landsat_N=NA
fabbricati_civilires_all@data$Ro75_landsat_L=NA
fabbricati_civilires_all@data$Ro75_landsat_L[indna_risk_o75_landsat]=as.character(temp_class_risk_o75_landsat_factor)
fabbricati_civilires_all@data$Ro75_landsat_N[indna_risk_o75_landsat]=as.numeric(temp_class_risk_o75_landsat_factor)


##############################################################################################################################################################
fabbricati_civilires_all$cs100=as.numeric(fabbricati_civilires_all$cs100)
fabbricati_civilires_all$cs57=as.numeric(fabbricati_civilires_all$cs57)

writeOGR(fabbricati_civilires_all,".","risk_parma_fcivres",driver="ESRI Shapefile",overwrite_layer = T)

saveRDS(fabbricati_civilires_all,"fabbricati_civilires_all.rds")

XLConnect::writeWorksheetToFile("analisi_rischio_ed.xls",data.frame(fabbricati_civilires_all@data),sheet="risk_parma")



