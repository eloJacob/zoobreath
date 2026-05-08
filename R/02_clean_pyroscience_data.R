#' Clean Pyroscience Data
#'
#' @param data A data.frame imported from a Pyroscience .csv file
#' @param unit Output unit for oxygen. One of: "umol.l", "%airsat", or "hPa"
#' @param salinity Salinity (PSU). Default = 38
#' @param dphi Logical. If TRUE, keeps dphi (phase) columns. Default = FALSE
#'
#' @return A cleaned tibble with standardized column names and oxygen in the chosen unit
#' @export
#'
Clean.pyroscience.data=function(data,unit,salinity=38, dphi=F){
  require(marelac)
  #press en pression lors de la mesure dbar 
  #unit = "%airsat" or "umol.l"  
  #Pour ancien logiciel
  if(str_detect(names(data)[1],"X")==T){
    #Merge firsts rows with colnames to have the correct colnames
    names(data)=str_c(names(data),data[1,],sep="_")
    #Keep only date, time, oxygen temperature and dphi id it's asked
    channel_used=data[-1,which(str_detect(colnames(data),"Oxygen")==T)]
    if(dphi==T){channel_used=cbind(channel_used,data[-1,c(18,19,20,21)])}
    data=cbind(data[-1,1:3],channel_used,data[-1,c("Temp..Probe_('C)","Pressure_(mbar)")])
    #Supprimate cols empty
    if(length(which(is.na(data[1,])==T))>0){
      data=data[,-which(is.na(data[1,])==T)]}
    #Vecteur temps
    day_start=as.numeric(str_split(data[1,"X_Date"],"/")[[1]][1])
    month_start=as.numeric(str_split(data[1,"X_Date"],"/")[[1]][2])
    year_start=as.numeric(str_split(data[1,"X_Date"],"/")[[1]][3])
    data$Time.h=time_to_hours(data[,"X_Date"],data[,"X.1_Time (HH:MM:SS)"],day_start,month_start,year_start)
    #Formate date
    data[,"X_Date"]=as.POSIXct(str_c(data[,"X_Date"],data[,"X.1_Time (HH:MM:SS)"],sep=" "),format = "%d/%m/%Y %H:%M:%S", tz="GMT")
    #Format numérique
    for(i in which(str_detect(names(data),"Ch")==T)[1]:length(names(data))){
      data[,i]<- as.numeric(as.character(data[,i]))}
    
    #rename_col
    names(data)[which(str_detect(colnames(data),"Date"))]="Date"
    names(data)[which(str_detect(colnames(data),"(HH:MM:SS)"))]="Time[H:M:S]"
    names(data)[which(str_detect(colnames(data),"X.2"))]="Time[s]"
    names(data)[which(str_detect(colnames(data),"Temp"))]="Temperature[C]"
    names(data)[which(str_detect(colnames(data),"Pressure"))]="Pressure[mbar]"
  }
  #Pour nouveau logiciel
  if(str_detect(names(data)[1],"Main")==T){
    data=data[-1,]
    #remote empty colomne
    if(length(which(is.na(data[1,])==T))>0){
      data=data[,-which(is.na(data[1,])==T)]}
    ##Compute information if ther are redondante 
    #if the time betwwen each measurements is the same => date, time , dt are the same 
    if(all(data[,which(str_detect(colnames(data),"dt"))]==data[,which(str_detect(colnames(data),"dt"))],na.rm=T)){
      data=data[,-which(str_detect(colnames(data),"dt"))[-1]]
      names(data)[which(str_detect(colnames(data),"dt"))]="Time[s]"
      #date
      data=data[,-which(str_detect(colnames(data),"Date"))[-1]]
      names(data)[which(str_detect(colnames(data),"Date"))]="Date"
      #time(h:m:s)
      data=data[,-which(str_detect(colnames(data),"Time"))[-1]]
      names(data)[which(str_detect(colnames(data),"Time"))]="Time[H:M:S]"
    }
    #if the temperature for all channel are measure with the same probe, there is no needs to keep a temperature per chanel
    if(all(data[,which(str_detect(colnames(data),"Temp"))]==data[,which(str_detect(colnames(data),"Temp"))],na.rm=T)){
      data=data[,-which(str_detect(colnames(data),"Temp"))[-1]]
      names(data)[which(str_detect(colnames(data),"Temp"))]="Temperature[C]"
    }
    #if the air pressure  were the same for all channel, there is no needs to keep a temperature per chanel
    if(all(data[,which(str_detect(colnames(data),"Pressure"))]==data[,which(str_detect(colnames(data),"Pressure"))],na.rm=T)){
      data=data[,-which(str_detect(colnames(data),"Pressure"))[-1]]
      names(data)[which(str_detect(colnames(data),"Pressure"))]="Pressure[mbar]"
    }
    #keep only the oxygen concentration measurement 
    if(dphi==T){data=data[,-which(str_detect(colnames(data),"Signal|Light|Status"))]}else{
      data=data[,-which(str_detect(colnames(data),"dphi|Signal|Light|Status"))]}
    ##Vecteur temps (h) 
    day_start=as.numeric(str_split(data[1,"Date"],"-")[[1]][1])
    month_start=as.numeric(str_split(data[1,"Date"],"-")[[1]][2])
    year_start=as.numeric(str_split(data[1,"Date"],"-")[[1]][3])
    data$Time.h=time_to_hours(data[,"Date"],data[,"Time[H:M:S]"],day_start,month_start,year_start, date_format = "d-m-y")
    ##Formate date
    data[,"X_Date"]=as.POSIXct(str_c(data[,"Date"],data[,"Time[H:M:S]"],sep=" "),format = "%d-%m-%Y h:m:%OS", tz="GMT")  
  }
  
  #Convertion µmol/l <--> air_sat
  if(unit=="umol.l"){
    id=0
    for (i in which(str_detect(colnames(data),"Oxygen"))){
      id=id+1
      name=paste0("airsat.Ch",id)
      # Calcule la saturation d'air et crée la nouvelle colonne avec le nom dynamique
      data[[name]] = data[, i] / marelac::gas_satconc(
        S = salinity,
        t = data[1,"Temperature[C]"],
        P = data[2,"Pressure[mbar]"]/1000,
        species = 'O2') * 100}}
  if(unit=="%airsat"){
    id=0
    for (i in which(str_detect(colnames(data),"Oxygen"))){
      id=id+1
      name=paste0("umol.l.Ch",id)
      # Calcule la saturation d'air et crée la nouvelle colonne avec le nom dynamique
      data[[name]] = (marelac::gas_satconc(S=salinity,
                                           t=data[1,"Temperature[C]"],
                                           P=data[2,"Pressure[mbar]"]/1000,
                                           species='O2')*data[, i] )/100}}
  if(unit=="hPa"){
    id=0
    for (i in which(str_detect(colnames(data),"Oxygen"))){
      id=id+1
      name=paste0("umol.l.Ch",id)
      # Calcule la saturation d'air et crée la nouvelle colonne avec le nom dynamique
      conv=presens::o2_unit_conv(data[, i] ,"hPa", "umol_per_l", 
                                 salinity, data[1,"Temperature[C]"], 
                                 data[2,"Pressure[mbar]"] / 1000)
      data[[name]] =conv$umol_per_l
    }}
  return(data)
}
