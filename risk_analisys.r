library(rgdal)
library(raster)
library(leaflet)
library(rgeos)
library(osmar)
library(doBy)
library(RColorBrewer)
library(mapview)
library(htmlwidgets)
library(htmltools)
library(webshot)
library(geosphere)
library(cleangeo)
library(classInt)
library(mgcv)

setwd("/home/alf/Scrivania/lav_parma_root/final")

########################################################################################################################################
# population data treatment

# parma_pop=readRDS("parma_pop.rds")
# 
# parma_pop[[1]]$class_pop=cut(parma_pop[[1]]$ETA,
#                            breaks=c(0,2,5,15,25,35,45,55,65,75,120),
#                            labels = c("00-02","03-05","06-15","16-25","26-35","36-45","46-55","56-65","66-75","76-120"),
#                            include.lowest = T)
# 
# parma_pop[[2]]$class_pop=cut(parma_pop[[2]]$ETA,
#                            breaks=c(0,2,5,15,25,35,45,55,65,75,120),
#                            labels = c("00-02","03-05","06-15","16-25","26-35","36-45","46-55","56-65","66-75","76-120"),
#                            include.lowest = T)
# 
# parma_pop[[2]]$count=1
# parma_pop[[1]]$count=1
# 
# saveRDS(parma_pop,"parma_pop.rds")
# 
# ########################################################################################################################################
# # create stats by age classes and save.
# 
# table_2014_full_sez=summaryBy(count~class_pop+SEZCENS, data=parma_pop[[1]],FUN=c(sum))
# table_2015_full_sez=summaryBy(count~class_pop+SEZCENS, data=parma_pop[[2]],FUN=c(sum))
# table_2014=as.data.frame(reshape2::acast(table_2014_full_sez,SEZCENS~class_pop,value.var="count.sum"))
# table_2015=as.data.frame(reshape2::acast(table_2015_full_sez,SEZCENS~class_pop,value.var="count.sum"))
# diff_table=table_2015[1:1158,]-table_2014
# table_2014$SEZCENS=rownames(table_2014)
# table_2015$SEZCENS=rownames(table_2015)
# diff_table$SEZCENS=rownames(table_2014)
# table_2014$sum2014=apply(table_2014[,1:10],1,function(x) sum(x,na.rm=T))
# table_2015$sum2015=apply(table_2015[,1:10],1,function(x) sum(x,na.rm=T))
# 
# write.csv(table_2014,"table_2014_parma_class.csv",row.names = F)
# write.csv(table_2015,"table_2015_parma_class.csv",row.names = F)
# write.csv(diff_table,"diff_table_2015_2014_parma_class.csv",row.names = F)
# 
# saveRDS(table_2014,"table_2014_parma_class.rds")
# saveRDS(table_2015,"table_2015_parma_class.rds")
# saveRDS(diff_table,"diff_table_2015_2014_parma_class.rds")

########################################################################################################################################
#

table_2014=readRDS("table_2014_parma_class.rds")
table_2015=readRDS("table_2015_parma_class.rds")
diff_table=readRDS("diff_table_2015_2014_parma_class.rds")

 sez_censuarie_istat=readOGR(".","R08_11_WGS84")
 sez_censuarie_parma=sez_censuarie_istat[which(sez_censuarie_istat@data$PRO_COM==34027),]
 saveRDS(sez_censuarie_parma,"sez_censuarie_parma_32N.rds")

sez_censuarie_parma=readRDS("sez_censuarie_parma_32N.rds")
sez_censuarie_geo_full_ls=sez_censuarie_parma

building_caire_geo_civ_ls=readRDS("building_caire_geo_civ_ls.rds")
sez_id_df=as.data.frame(over(building_caire_geo_civ_ls,sez_censuarie_geo_full_ls))

building_caire_geo_civ_ls@data$SEZID=sez_id_df$SEZ


sez_censuarie_geo_full_ls$"Y1400-02"=NA   
sez_censuarie_geo_full_ls$"Y1403-05"=NA   
sez_censuarie_geo_full_ls$"Y1406-15"=NA   
sez_censuarie_geo_full_ls$"Y1416-25"=NA   
sez_censuarie_geo_full_ls$"Y1426-35"=NA   
sez_censuarie_geo_full_ls$"Y1436-45"=NA   
sez_censuarie_geo_full_ls$"Y1446-55"=NA   
sez_censuarie_geo_full_ls$"Y1456-65"=NA   
sez_censuarie_geo_full_ls$"Y1466-75"=NA   
sez_censuarie_geo_full_ls$"Y1476-120"=NA  
sez_censuarie_geo_full_ls$"Y14sum2014"=NA


for ( i in 1:length(table_2014$SEZCENS)){

temp=which(sez_censuarie_geo_full_ls$SEZ==table_2014$SEZCENS[i])

sez_censuarie_geo_full_ls$"Y1400-02"[temp]=table_2014[i,1]  
sez_censuarie_geo_full_ls$"Y1403-05"[temp]=table_2014[i,2]    
sez_censuarie_geo_full_ls$"Y1406-15"[temp]=table_2014[i,3]  
sez_censuarie_geo_full_ls$"Y1416-25"[temp]=table_2014[i,4]  
sez_censuarie_geo_full_ls$"Y1426-35"[temp]=table_2014[i,5]  
sez_censuarie_geo_full_ls$"Y1436-45"[temp]=table_2014[i,6]   
sez_censuarie_geo_full_ls$"Y1446-55"[temp]=table_2014[i,7]   
sez_censuarie_geo_full_ls$"Y1456-65"[temp]=table_2014[i,8]   
sez_censuarie_geo_full_ls$"Y1466-75"[temp]=table_2014[i,9]  
sez_censuarie_geo_full_ls$"Y1476-120"[temp]=table_2014[i,10] 
sez_censuarie_geo_full_ls$"Y14sum2014"[temp]=table_2014[i,12]



}

sez_censuarie_geo_full_ls$"Y1500-02"=NA   
sez_censuarie_geo_full_ls$"Y1503-05"=NA   
sez_censuarie_geo_full_ls$"Y1506-15"=NA   
sez_censuarie_geo_full_ls$"Y1516-25"=NA   
sez_censuarie_geo_full_ls$"Y1526-35"=NA   
sez_censuarie_geo_full_ls$"Y1536-45"=NA   
sez_censuarie_geo_full_ls$"Y1546-55"=NA   
sez_censuarie_geo_full_ls$"Y1556-65"=NA   
sez_censuarie_geo_full_ls$"Y1566-75"=NA   
sez_censuarie_geo_full_ls$"Y1576-120"=NA  
sez_censuarie_geo_full_ls$"Y15sum2015"=NA  

for ( i in 1:length(table_2015$SEZCENS)){
  
  temp=which(sez_censuarie_geo_full_ls$SEZ==table_2015$SEZCENS[i])
  
  sez_censuarie_geo_full_ls$"Y1500-02"[temp]=table_2015[i,1]  
  sez_censuarie_geo_full_ls$"Y1503-05"[temp]=table_2015[i,2]    
  sez_censuarie_geo_full_ls$"Y1506-15"[temp]=table_2015[i,3]  
  sez_censuarie_geo_full_ls$"Y1516-25"[temp]=table_2015[i,4]  
  sez_censuarie_geo_full_ls$"Y1526-35"[temp]=table_2015[i,5]  
  sez_censuarie_geo_full_ls$"Y1536-45"[temp]=table_2015[i,6]   
  sez_censuarie_geo_full_ls$"Y1546-55"[temp]=table_2015[i,7]   
  sez_censuarie_geo_full_ls$"Y1556-65"[temp]=table_2015[i,8]   
  sez_censuarie_geo_full_ls$"Y1566-75"[temp]=table_2015[i,9]  
  sez_censuarie_geo_full_ls$"Y1576-120"[temp]=table_2015[i,10] 
  sez_censuarie_geo_full_ls$"Y15sum2015"[temp]=table_2015[i,12]
  
  
  
}


sez_censuarie_geo_full_ls$"YD00-02"=NA   
sez_censuarie_geo_full_ls$"YD03-05"=NA   
sez_censuarie_geo_full_ls$"YD06-15"=NA   
sez_censuarie_geo_full_ls$"YD16-25"=NA   
sez_censuarie_geo_full_ls$"YD26-35"=NA   
sez_censuarie_geo_full_ls$"YD36-45"=NA   
sez_censuarie_geo_full_ls$"YD46-55"=NA   
sez_censuarie_geo_full_ls$"YD56-65"=NA   
sez_censuarie_geo_full_ls$"YD66-75"=NA   
sez_censuarie_geo_full_ls$"YD76-120"=NA  

for ( i in 1:length(diff_table$SEZCENS)){
  
  temp=which(sez_censuarie_geo_full_ls$SEZ==diff_table$SEZCENS[i])
  sez_censuarie_geo_full_ls$"YD00-02"[temp]=diff_table[i,1]  
  sez_censuarie_geo_full_ls$"YD03-05"[temp]=diff_table[i,2]    
  sez_censuarie_geo_full_ls$"YD06-15"[temp]=diff_table[i,3]  
  sez_censuarie_geo_full_ls$"YD16-25"[temp]=diff_table[i,4]  
  sez_censuarie_geo_full_ls$"YD26-35"[temp]=diff_table[i,5]  
  sez_censuarie_geo_full_ls$"YD36-45"[temp]=diff_table[i,6]   
  sez_censuarie_geo_full_ls$"YD46-55"[temp]=diff_table[i,7]   
  sez_censuarie_geo_full_ls$"YD56-65"[temp]=diff_table[i,8]   
  sez_censuarie_geo_full_ls$"YD66-75"[temp]=diff_table[i,9]  
  sez_censuarie_geo_full_ls$"YD76-120"[temp]=diff_table[i,10] 
  
  
  
}

sez_censuarie_geo_full_ls@data$over65_15=sez_censuarie_geo_full_ls@data$`Y1566-75`+sez_censuarie_geo_full_ls@data$`Y1576-120`
sez_censuarie_geo_full_ls@data$over65_D=sez_censuarie_geo_full_ls@data$`YD66-75`+sez_censuarie_geo_full_ls@data$`YD76-120`

saveRDS(sez_censuarie_geo_full_ls,"sez_censuarie_data_ls.rds")

# mapview(sez_censuarie_geo_full_ls,zcol="YD76-120",color=heat.colors(15),alpha.regions = 0.4,legend=T)
###########################################################################################################à


building_caire_geo_civ_ls$"Y1400-02"=NA   
building_caire_geo_civ_ls$"Y1403-05"=NA   
building_caire_geo_civ_ls$"Y1406-15"=NA   
building_caire_geo_civ_ls$"Y1416-25"=NA   
building_caire_geo_civ_ls$"Y1426-35"=NA   
building_caire_geo_civ_ls$"Y1436-45"=NA   
building_caire_geo_civ_ls$"Y1446-55"=NA   
building_caire_geo_civ_ls$"Y1456-65"=NA   
building_caire_geo_civ_ls$"Y1466-75"=NA   
building_caire_geo_civ_ls$"Y1476-120"=NA  
building_caire_geo_civ_ls$"Y14sum2014"=NA


for ( i in 1:length(table_2014$SEZCENS)){
  
  temp=which(building_caire_geo_civ_ls$SEZID==table_2014$SEZCENS[i])
  
  building_caire_geo_civ_ls$"Y1400-02"[temp]=table_2014[i,1]  
  building_caire_geo_civ_ls$"Y1403-05"[temp]=table_2014[i,2]    
  building_caire_geo_civ_ls$"Y1406-15"[temp]=table_2014[i,3]  
  building_caire_geo_civ_ls$"Y1416-25"[temp]=table_2014[i,4]  
  building_caire_geo_civ_ls$"Y1426-35"[temp]=table_2014[i,5]  
  building_caire_geo_civ_ls$"Y1436-45"[temp]=table_2014[i,6]   
  building_caire_geo_civ_ls$"Y1446-55"[temp]=table_2014[i,7]   
  building_caire_geo_civ_ls$"Y1456-65"[temp]=table_2014[i,8]   
  building_caire_geo_civ_ls$"Y1466-75"[temp]=table_2014[i,9]  
  building_caire_geo_civ_ls$"Y1476-120"[temp]=table_2014[i,10] 
  building_caire_geo_civ_ls$"Y14sum2014"[temp]=table_2014[i,12]
  
  
  
}

building_caire_geo_civ_ls$"Y1500-02"=NA   
building_caire_geo_civ_ls$"Y1503-05"=NA   
building_caire_geo_civ_ls$"Y1506-15"=NA   
building_caire_geo_civ_ls$"Y1516-25"=NA   
building_caire_geo_civ_ls$"Y1526-35"=NA   
building_caire_geo_civ_ls$"Y1536-45"=NA   
building_caire_geo_civ_ls$"Y1546-55"=NA   
building_caire_geo_civ_ls$"Y1556-65"=NA   
building_caire_geo_civ_ls$"Y1566-75"=NA   
building_caire_geo_civ_ls$"Y1576-120"=NA  
building_caire_geo_civ_ls$"Y15sum2015"=NA  

for ( i in 1:length(table_2015$SEZCENS)){
  
  temp=which(building_caire_geo_civ_ls$SEZID==table_2015$SEZCENS[i])
  
  building_caire_geo_civ_ls$"Y1500-02"[temp]=table_2015[i,1]  
  building_caire_geo_civ_ls$"Y1503-05"[temp]=table_2015[i,2]    
  building_caire_geo_civ_ls$"Y1506-15"[temp]=table_2015[i,3]  
  building_caire_geo_civ_ls$"Y1516-25"[temp]=table_2015[i,4]  
  building_caire_geo_civ_ls$"Y1526-35"[temp]=table_2015[i,5]  
  building_caire_geo_civ_ls$"Y1536-45"[temp]=table_2015[i,6]   
  building_caire_geo_civ_ls$"Y1546-55"[temp]=table_2015[i,7]   
  building_caire_geo_civ_ls$"Y1556-65"[temp]=table_2015[i,8]   
  building_caire_geo_civ_ls$"Y1566-75"[temp]=table_2015[i,9]  
  building_caire_geo_civ_ls$"Y1576-120"[temp]=table_2015[i,10] 
  building_caire_geo_civ_ls$"Y15sum2015"[temp]=table_2015[i,12]
  
  
  
}


building_caire_geo_civ_ls$"YD00-02"=NA   
building_caire_geo_civ_ls$"YD03-05"=NA   
building_caire_geo_civ_ls$"YD06-15"=NA   
building_caire_geo_civ_ls$"YD16-25"=NA   
building_caire_geo_civ_ls$"YD26-35"=NA   
building_caire_geo_civ_ls$"YD36-45"=NA   
building_caire_geo_civ_ls$"YD46-55"=NA   
building_caire_geo_civ_ls$"YD56-65"=NA   
building_caire_geo_civ_ls$"YD66-75"=NA   
building_caire_geo_civ_ls$"YD76-120"=NA  

for ( i in 1:length(diff_table$SEZCENS)){
  
  temp=which(building_caire_geo_civ_ls$SEZID==diff_table$SEZCENS[i])
  
  building_caire_geo_civ_ls$"YD00-02"[temp]=diff_table[i,1]  
  building_caire_geo_civ_ls$"YD03-05"[temp]=diff_table[i,2]    
  building_caire_geo_civ_ls$"YD06-15"[temp]=diff_table[i,3]  
  building_caire_geo_civ_ls$"YD16-25"[temp]=diff_table[i,4]  
  building_caire_geo_civ_ls$"YD26-35"[temp]=diff_table[i,5]  
  building_caire_geo_civ_ls$"YD36-45"[temp]=diff_table[i,6]   
  building_caire_geo_civ_ls$"YD46-55"[temp]=diff_table[i,7]   
  building_caire_geo_civ_ls$"YD56-65"[temp]=diff_table[i,8]   
  building_caire_geo_civ_ls$"YD66-75"[temp]=diff_table[i,9]  
  building_caire_geo_civ_ls$"YD76-120"[temp]=diff_table[i,10] 
  
  
  
}




building_caire_geo_civ_ls@data$over65_15=building_caire_geo_civ_ls@data$"Y1566-75"+building_caire_geo_civ_ls@data$"Y1576-120" 
building_caire_geo_civ_ls@data$E_D_tot=(building_caire_geo_civ_ls@data$Y15sum2015/building_caire_geo_civ_ls@data$S_EDCens)*building_caire_geo_civ_ls@data$SHAPE_AREA
building_caire_geo_civ_ls@data$E_Dover65=(building_caire_geo_civ_ls@data$"over65_15"/building_caire_geo_civ_ls@data$S_EDCens)*building_caire_geo_civ_ls@data$SHAPE_AREA
building_caire_geo_civ_ls@data$E_Dover75=(building_caire_geo_civ_ls@data$"Y1576-120"/building_caire_geo_civ_ls@data$S_EDCens)*building_caire_geo_civ_ls@data$SHAPE_AREA
building_caire_geo_civ_ls@data$E_Dunder2=(building_caire_geo_civ_ls@data$"Y1500-02"/building_caire_geo_civ_ls@data$S_EDCens)*building_caire_geo_civ_ls@data$SHAPE_AREA

saveRDS(building_caire_geo_civ_ls,"building_caire_civ_data_ls.rds")


building_caire_geo_civ_ls$cs200_sum_scaled=as.numeric(scale(building_caire_geo_civ_ls$cs200_sum,center=min(building_caire_geo_civ_ls$cs200_sum,na.rm=T),scale=diff(range(building_caire_geo_civ_ls$cs200_sum))))

range_lstmean=as.numeric(c(quantile(building_caire_geo_civ_ls@data$lstmean,0,na.rm=T),quantile(building_caire_geo_civ_ls@data$lstmean,0.99,na.rm=T)))
building_caire_geo_civ_ls$lstmean_scaled=as.numeric(scale(building_caire_geo_civ_ls$lstmean,center=min(building_caire_geo_civ_ls$lstmean,na.rm=T),scale=diff(range_lstmean)))
building_caire_geo_civ_ls$lstmean_scaled[which(building_caire_geo_civ_ls$lstmean_scaled>1)]=1


range_E_D_tot=as.numeric(c(0,quantile(building_caire_geo_civ_ls@data$E_D_tot,0.95,na.rm=T)))
range_E_Dover65=as.numeric(c(0,quantile(building_caire_geo_civ_ls@data$E_Dover65,0.95,na.rm=T)))
range_E_Dover75=as.numeric(c(0,quantile(building_caire_geo_civ_ls@data$E_Dover75,0.95,na.rm=T)))
range_E_Dunder2=as.numeric(c(0,quantile(building_caire_geo_civ_ls@data$E_Dunder2,0.95,na.rm=T)))

building_caire_geo_civ_ls$dpoptot_scaled=as.numeric(scale(building_caire_geo_civ_ls$E_D_tot,center=min(building_caire_geo_civ_ls$E_D_tot,na.rm=T),scale=diff(range_E_D_tot)))
building_caire_geo_civ_ls$dpoptot_scaled[which(building_caire_geo_civ_ls$dpoptot_scaled>1)]=1

building_caire_geo_civ_ls$dpop65_scaled=as.numeric(scale(building_caire_geo_civ_ls$E_Dover65,center=min(building_caire_geo_civ_ls$E_Dover65,na.rm=T),scale=diff(range_E_Dover65)))
building_caire_geo_civ_ls$dpop65_scaled[which(building_caire_geo_civ_ls$dpop65_scaled>1)]=1

building_caire_geo_civ_ls$dpop75_scaled=as.numeric(scale(building_caire_geo_civ_ls$E_Dover75,center=min(building_caire_geo_civ_ls$E_Dover75,na.rm=T),scale=diff(range_E_Dover75)))
building_caire_geo_civ_ls$dpop75_scaled[which(building_caire_geo_civ_ls$dpop75_scaled>1)]=1

building_caire_geo_civ_ls$dpop2_scaled=as.numeric(scale(building_caire_geo_civ_ls$E_Dunder2,center=min(building_caire_geo_civ_ls$E_Dunder2,na.rm=T),scale=diff(range_E_Dunder2)))
building_caire_geo_civ_ls$dpop2_scaled[which(building_caire_geo_civ_ls$dpop2_scaled>1)]=1

saveRDS(building_caire_geo_civ_ls,"building_caire_civ_data_ls.rds")



##############################################################################################################################################################
building_caire_geo_civ_ls=readRDS("building_caire_civ_data_ls.rds")

building_caire_geo_civ_ls$riskvalue_75=(0.5*building_caire_geo_civ_ls$lstmean_scaled)+((1/3)*building_caire_geo_civ_ls$dpop75_scaled+(1/3)*building_caire_geo_civ_ls$dpoptot_scaled+(1/3)*building_caire_geo_civ_ls$cs200_sum_scaled)*0.5;
building_caire_geo_civ_ls$riskvalue_65=(0.5*building_caire_geo_civ_ls$lstmean_scaled)+((1/3)*building_caire_geo_civ_ls$dpop65_scaled+(1/3)*building_caire_geo_civ_ls$dpoptot_scaled+(1/3)*building_caire_geo_civ_ls$cs200_sum_scaled)*0.5;
building_caire_geo_civ_ls$riskvalue_2=(0.5*building_caire_geo_civ_ls$lstmean_scaled)+((1/3)*building_caire_geo_civ_ls$dpop2_scaled+(1/3)*building_caire_geo_civ_ls$dpoptot_scaled+(1/3)*building_caire_geo_civ_ls$cs200_sum_scaled)*0.5;


indna_75=which(!is.na(building_caire_geo_civ_ls$riskvalue_75))
temp_class_risk_75=classIntervals(building_caire_geo_civ_ls$riskvalue_75[indna_75], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_75_factor=cut(building_caire_geo_civ_ls$riskvalue_75[indna_75], breaks = temp_class_risk_75$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
building_caire_geo_civ_ls@data$RC_75_N=NA
building_caire_geo_civ_ls@data$RC_75_L=NA
building_caire_geo_civ_ls@data$RC_75_L[indna_75]=as.character(temp_class_risk_75_factor)
building_caire_geo_civ_ls@data$RC_75_N[indna_75]=as.numeric(temp_class_risk_75_factor)

indna_65=which(!is.na(building_caire_geo_civ_ls$riskvalue_65))
temp_class_risk_65=classIntervals(building_caire_geo_civ_ls$riskvalue_65[indna_65], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_65_factor=cut(building_caire_geo_civ_ls$riskvalue_65[indna_65], breaks = temp_class_risk_65$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
building_caire_geo_civ_ls@data$RC_65_N=NA
building_caire_geo_civ_ls@data$RC_65_L=NA
building_caire_geo_civ_ls@data$RC_65_L[indna_65]=as.character(temp_class_risk_65_factor)
building_caire_geo_civ_ls@data$RC_65_N[indna_65]=as.numeric(temp_class_risk_65_factor)

indna_2=which(!is.na(building_caire_geo_civ_ls$riskvalue_2))
temp_class_risk_2=classIntervals(building_caire_geo_civ_ls$riskvalue_2[indna_2], n =5,style="fixed",fixedBreaks=c(0,0.20,0.40,0.60,0.80,1))
temp_class_risk_2_factor=cut(building_caire_geo_civ_ls$riskvalue_2[indna_2], breaks = temp_class_risk_2$brks, labels=c("Very Low", "Low","Moderate","High","Very High"))
building_caire_geo_civ_ls@data$RC_2_N=NA
building_caire_geo_civ_ls@data$RC_2_L=NA
building_caire_geo_civ_ls@data$RC_2_L[indna_2]=as.character(temp_class_risk_2_factor)
building_caire_geo_civ_ls@data$RC_2_N[indna_2]=as.numeric(temp_class_risk_2_factor)

saveRDS(building_caire_geo_civ_ls,"building_caire_civ_data_ls.rds")


########################################################################################################################################
# reprojection

griglia_geo=readRDS("griglia_geo.rds")

confini_comune_geo=readRDS("confini_comune_geo.rds")
confini_comune_ls=readRDS("confini_comune_landsat.rds")
building_caire_geo_civ_ls=readRDS("building_caire_civ_data_ls.rds")
sez_censuarie_geo_full_ls=readRDS("sez_censuarie_data_ls.rds")


building_caire_geo_civ_ls@data$perc_csum=(building_caire_geo_civ_ls@data$cs200_sum/1241)*100
aa=mgcv::gam(lstmean~s(perc_csum),data=subset(building_caire_geo_civ_ls@data,cs200_sum>600))
plot(aa,ylab="RR - LST mean",main="Relation between relative surface temperature increase \nand degree of soil sealing (% < 200m)")
abline(h=0)


proj_caire_CRS=readRDS("proj_caire.rds")

sez_censuarie_EM_data=spTransform(sez_censuarie_geo_full_ls,CRS(proj_caire_CRS))
buildcivparma_EM_data=spTransform(building_caire_geo_civ_ls,CRS(proj_caire_CRS))
sez_censuarie_geo_data=spTransform(sez_censuarie_geo_full_ls,proj4string(confini_comune_geo))
buildcivparma_geo_data=spTransform(building_caire_geo_civ_ls,proj4string(confini_comune_geo))

saveRDS(sez_censuarie_EM_data,"sez_censuarie_EM_data.rds")
saveRDS(buildcivparma_EM_data,"buildcivparma_EM_data.rds")
saveRDS(sez_censuarie_geo_data,"sez_censuarie_geo_data.rds")
saveRDS(buildcivparma_geo_data,"buildcivparma_geo_data.rds")

writeOGR(sez_censuarie_EM_data,".","sez_censuarie_EM_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(buildcivparma_EM_data,".","buildcivparma_EM_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(buildcivparma_geo_data,".","buildcivparma_geo_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(sez_censuarie_geo_data,".","sez_censuarie_geo_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(building_caire_geo_civ_ls,".","buildcivparma_32N_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(sez_censuarie_geo_full_ls,".","sez_censuarie_32N_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(sez_censuarie_parma,".","sez_censuarie_32N_data_bis",driver="ESRI Shapefile",overwrite_layer = T)


writeOGR(griglia_geo,".","griglia_geo",driver="ESRI Shapefile",overwrite_layer = T)
griglia_32N=spTransform(griglia_geo,proj4string(sez_censuarie_geo_full_ls))
writeOGR(griglia_32N,".","griglia_32N",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(confini_comune_ls,".","confini_comune_32N",driver="ESRI Shapefile",overwrite_layer = T)

lottizzazioni=readRDS("Lottizzazioni_ls.rds")
PEEP=readRDS("PEEP_ls.rds")
writeOGR(lottizzazioni,".","lottizzazioni_32N_data",driver="ESRI Shapefile",overwrite_layer = T)
writeOGR(PEEP,".","PEEP_32N_data",driver="ESRI Shapefile",overwrite_layer = T)


########################################################################################################################################
# read vulnerability and hazard layer tmerc32 WGS84 

parma_consumo_suolo12=raster("parma_consumo_suolo_2012.tif")
parma_consumo_suolo15=raster("parma_consumo_suolo_2015.tif")

parma_lst_2015=readRDS("parma_city_lst_15_07_2015.rds")

parma_lst_2015_crop <- crop(parma_lst_2015, extent(confini_comune_ls))
parma_consumo_suolo15_crop  <- crop(parma_consumo_suolo15, extent(confini_comune_ls))
parma_consumo_suolo12_crop  <- crop(parma_consumo_suolo12, extent(confini_comune_ls))

parma_lst_2015_mask <- mask(parma_lst_2015_crop , confini_comune_ls)
parma_consumo_suolo12_mask <- mask(parma_consumo_suolo12_crop, confini_comune_ls)
parma_consumo_suolo15_mask <- mask(parma_consumo_suolo15_crop, confini_comune_ls)

writeRaster(parma_lst_2015_mask, filename="parma_lst_2015_mask.tif", overwrite=TRUE)
writeRaster(parma_consumo_suolo12_mask, filename="parma_consumo_suolo2012_mask.tif", overwrite=TRUE)
writeRaster(parma_consumo_suolo15_mask, filename="parma_consumo_suolo2015_mask.tif", overwrite=TRUE)

parma_consumo_suolo15_mask[is.na(parma_consumo_suolo15_mask)] <- 0
parma_consumo_suolo12_mask[is.na(parma_consumo_suolo12_mask)] <- 0
diff_parma_consumo_suolo=parma_consumo_suolo15_mask-parma_consumo_suolo12_mask

writeRaster(diff_parma_consumo_suolo, filename="diff15-12_parma_consumo_suolo_mask.tif", overwrite=TRUE)



########################################################################################################################################


potcom2=readOGR(".","potcom2")
proj4string(potcom2)=CRS("+proj=tmerc +lat_0=0 +lon_0=9 +k=0.9996 +x_0=500000 +y_0=-4000000 +ellps=intl +towgs84=-87,-98,-121,0,0,0,0 +units=m +no_defs")
potcom2_geo=spTransform(potcom2,proj4string(building_caire_geo_civ))
mapview(potcom2_geo)
saveRDS(potcom2_geo,"potcom2_geo.rds")


##############################################################################################################################################################àà
# palette

morab_pal_VGARV=c("#00CC00","#FFFF00","#FF7519","#FF0000","#FF00FF")
morabcols <- colorRampPalette( c(morab_pal_VGARV))(5)


# Normalize layers

