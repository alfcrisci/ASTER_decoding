##################################################################################################
#devtools::install_github("hunzikp/velox")
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
library(doBy)
######################################################################################################

source("aux_aster.r")

######################################################################################################

confini_comune_ls=readRDS("data/confini_comune_landsat.rds")
parma_urbanizzato_ls=readRDS("data/parma_urbanizzato_ls.rds")

######################################################################################################


fabbricati_stack=readRDS("data/vectors_list_new.rds")
ncivici=fabbricati_stack[[7]]
sez_censuarie=fabbricati_stack[[6]]
fabbcivres=fabbricati_stack[[1]]



###########################################################################################
# Elaborazioni sui fabbricati filtro e creazione centroidi

#fabbri161014_EDI2011=readRDS("data/fabbri161014_EDI2011.rds")


##########################################################################################################################################################
# check


over_edi=over(EDI2011_1grez,fabbri161014_EDI2011)
new_fabbri161014_EDI2011=fabbri161014_EDI2011[na.omit(over_edi$index),]
new_fabbri161014_EDI2011$index_building=floor(as.numeric(rownames(new_fabbri161014_EDI2011@data)))
data_building=new_fabbri161014_EDI2011@data

RES15ed=summaryBy(RES15~as.factor(data_building$index_building),data=data_building,fun=sum())
RES14ed=summaryBy(RES14~as.factor(data_building$index_building),data=data_building,fun=sum())
FAM15ed=summaryBy(RES15~as.factor(data_building$index_building),data=data_building,fun=sum())
FAM14ed=summaryBy(RES14~as.factor(data_building$index_building),data=data_building,fun=sum())


for ( i in 1:nrow(new_fabbri161014_EDI2011@data)) {
            
  
  new_fabbri161014_EDI2011@data$RES15[i]=RES15ed$RES15.0[which(RES15ed$index_building==new_fabbri161014_EDI2011@data$index_building[i])]
  new_fabbri161014_EDI2011@data$RES14[i]=RES14ed$RES14.0[which(RES15ed$index_building==new_fabbri161014_EDI2011@data$index_building[i])]
  new_fabbri161014_EDI2011@data$FAM15[i]=RES15ed$FAM15.0[which(RES15ed$index_building==new_fabbri161014_EDI2011@data$index_building[i])]
  new_fabbri161014_EDI2011@data$FAM14[i]=RES14ed$FAM14.0[which(RES15ed$index_building==new_fabbri161014_EDI2011@data$index_building[i])]
  
}

rownames(new_fabbri161014_EDI2011_data)=1:nrow(new_fabbri161014_EDI2011_data)
pols <-  new_fabbri161014_EDI2011@polygons
regions <-  SpatialPolygonsDataFrame(regions, sdf)

ids = sapply(slot(new_fabbri161014_EDI2011, "polygons"), function(i) slot(i, "ID"))

rownames(new_fabbri161014_EDI2011@data)=make.unique(ids)

saveRDS(fabbri161014_EDI2011,"data/new_fabbri161014_EDI2011.rds")

##########################################################################################################################################################

fabbri161014_EDI2011=readRDS("data/fabbri161014_EDI2011.rds")

fabbricati_all=spTransform(fabbri161014_EDI2011,proj4string(confini_comune_ls))



saveRDS(fabbricati_all,"data/fabbricati_all_ls.rds")

fabbricati_civili_all=fabbricati_all[which(fabbricati_all$DESC_TIPO  =="FABBRICATO CIVILE"),]

saveRDS(fabbricati_civili_all,"data/fabbricati_civili_all_ls.rds")

fabbricati_civilires_all=fabbricati_civili_all[which(fabbricati_civili_all$RES14 >0),]
fabbricati_civilires_all$D15_14=fabbricati_civilires_all$RES15-fabbricati_civilires_all$RES14
fabbricati_civilires_all$DFAM15_14=fabbricati_civilires_all$FAM15-fabbricati_civilires_all$FAM14

fabbricati_civilires_centroid=gCentroid(fabbricati_civilires_all,byid=T)


fabbricati_civilires_centroid_df=data.frame(fabbricati_civilires_centroid)
fabbricati_civilires_centroid_df$FABCD=fabbricati_civilires_all@data$FABCD
coordinates(fabbricati_civilires_centroid_df) =~x+y
proj4string(fabbricati_civilires_centroid_df)=proj4string(confini_comune_ls)
fabbricati_civilires_centroid=SpatialPointsDataFrame(fabbricati_civilires_centroid,data.frame(fabbricati_civilires_centroid))

saveRDS(fabbricati_civilires_centroid,"data/fabbricati_civilires_centroid.rds")


q=over(fabbricati_civilires_centroid_df,parma_urbanizzato_ls)
fabbricati_civilires_all$UUR_1=q$UUR_1

saveRDS(fabbricati_civilires_all,"data/fabbricati_civilires_all.rds")



