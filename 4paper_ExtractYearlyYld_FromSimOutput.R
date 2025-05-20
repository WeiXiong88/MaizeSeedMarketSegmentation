#######These scripts are for identifying area suitability for maize and sorghum in Africa#########
########################This work is for a high impact papaers####################################
##Scripts was created and maintained by Wei Xiong (w.xiong@cgiar.org) on 12/2/2024###############
######################################KEY STEPS###################################################
#Novelty  1: A unique approach to convert simulated yields to suitability 0~1 by considering mean yield and CV
          #  a. covert yield by a logistic function, L/(1+exp(-10*cv*(x-x0)))
          #   L=1, maximum suitability, cv: yield cv, x: mean yield, x0:median of mean yield
          #   function: y=x/5 when cv=0
          #             y=1/(1+exp(-10*cv*(x-x0)))
          #  b. suitability is the mean of SI with yields that sowing in the growing season month, averaged across crop maturity 
          #  c. suitability is calculated seperately for major growing and minor growing seasons
          #  d. computed suitability is compared to other gridded products, such as SPAM, Monfrda, Mirca to see relationship between
          #     yield/area with SI
    ####Scripts: 1. compute SI from simulated yield, 2. compare SI with gridded products.
    ####Plots: Comparison with gridded products, correlation between reported yield/area share with suitability
#Novelty  2: Suitability area for the current climate and its changes over time
          #  a. how the suitability has changed over time, where has pronounced changes
          #  b. what has attributed to this change, changes in growing season - onset, duration, changes in climate, 
          #  a extracted table to identify the contribution for SI variations - SI,yld,cv,onset,duration,p,t,drought,heat
    ####Scripts: 1. identify suitability map, visualize the suitability change, and the reason for such change
    ####Plots: suitability change over time, and key contributors
#Novelty  3: Maize cultivars and maize vs. Sorghum - 
          #  a. transitional change for cultivar maturity groups
          #  a. transitional change for crop replacement
          #  b. identify the suitability areas or possibility to rotate crops
    ####Scripts: 1. identify the space and time for the transitional change - maturity group, and crops, identify rotation hotspots
    ####Plots: space and time with transitional change

#function to compute smooth mean of simulated outcomes for each block txt file
library(pracma) #library for function findpeaks()
library(zoo)
library(dplyr)
Find_season_mon0<-function(yld){
  out<-matrix(NA,nrow(yld),6)                                       #out for return
  series<-as.numeric(colMeans(yld))                                 #identify seasons by mean yield
  series[is.na(series)] <- 0                                        #replace na yield as 0
  series<-c(series[12],rep(series,3),series[1])  #create a three year series
  series<-rollmean(series,3)                                        #smooth by 3 months
  peaks<-findpeaks(series,minpeakdistance=4,sortstr=TRUE) #minpeakheight=1,
  sea1<-peaks[(peaks[,2]>12)&(peaks[,2]<25),2][1]
  #sea1<-peaks[(peaks[,2]>12)&(peaks[,2]<25),2][1]#%%12                 #select the peak in one year
  #peaks[peaks==0]<-12  
  #for case the peak at the last year
  if(length(sea1)>0){
    for(i in 1:nrow(yld)){  #change peaks matrix start and end to yield greater than 1t/ha
      series<-as.numeric(yld[i,]) #three year series
      series[is.na(series)] <- 0
      series<-c(series[12],rep(series,3),series[1])
      series<-rollmean(series,3)
      peaks<-findpeaks(series,minpeakdistance=4,sortstr=TRUE) #find peak for each row, minpeakheight=0.9,
      if(length(peaks)>3){
        peaks<-cbind(peaks,abs(peaks[,2]-sea1))  #add the peak distance to season 1 peak
        peaks<-peaks[order(peaks[,5]),]    #sort by distance from nearest to farest
        peaks<-peaks[(peaks[,2]>12)&(peaks[,2]<25),2:4]  #slice peaks in the middle year
        for(r in 1:min(2,length(peaks)%/%3)){ #in case there is only one peak
          if(length(peaks)>3){     #if there are more than one peaks
            row<-peaks[r,]
          }else{row<-peaks}
          #while((series[row[2]]<1)&(row[2]!=row[1])){row[2]<-row[2]+1} #move the concession left if yield less tyldhan 1
          #while((series[row[3]]<1)&(row[2]!=row[1])){row[3]<-row[3]-1} #move the concession right if yield less than 1
          #if(min(abs(row[2]-row[3]),12-abs(row[2]-row[3]))>=1)
          out[i,((r-1)*3+1):((r-1)*3+3)]<-row
      }}}
  out<-out%%12
  out[out==0]<-12}
  colnames(out)<-c("1s_m","1s_ms","1s_me","2s_m","2s_ms","2s_me")
  return(as.data.frame(out)) #peaks[,2:4]
}

Find_season_mon<-function(series){
  result<-rep(NA,2)                                   #out for return                               #identify seasons by mean yield
  series[is.na(series)] <- 0                                        #replace na yield as 0
  series<-rep(series,3)                                             #create a three year series
  peaks<-findpeaks(series,minpeakdistance=3,minpeakheight=0.1,sortstr=TRUE) #minpeakheight=1,
  if(length(peaks)>0){
    peaks<-as.data.frame(peaks)%>%filter((V2>12)&(V2<25)) #peak month
    peaks<-peaks$V2%%12
    peaks<-ifelse(peaks==0,12,peaks)
    result<-peaks}
  return(result) #peaks[,2:4]
}

####################THIS SECTION IS FOR EXTRACTING YEARLY YIELD FOR THREE MATURITY GROUP AND SEASON WITH DIFFERENT SOWING SCENARIOS########
#Type: Reported - reported sowing month: reported sowing month sea1-ISIMIP, sea2-SAGA same and fixed for all maturity group
#      Fixed - Optimum and fixed sowing month: sowing month is determined by mean yield (3month rolling) of 50 years, yield is mean of three maturities
#      Optimumfixed - Optimum and fixed sowing month: sowing month is determined by mean yield (3month rolling) of 50 years, differ in maturity
#      OptimumDecade - Decaded optimum sowing month: sowing month varies every 10 years depending on mean yield/CV, differ in maturity in decade
#      OptimumYear   - Optimum yearly sowing month: sowing month varies every year given simulated annual yield, differ in maturity and year

GetyearlyYld<-function(yld,mon,type){
  out <- data.frame(matrix(ncol =2+nrow(yld), nrow = 0))  #create a empty dataframe
  #############################Reported calendar##########################################################################
  #reported sowing month is same for the three matuirty groups
  if((length(mon)>0)&(type=="Reported")){
    mon<-mon[!is.na(mon)]
    if(length(mon)>0){for(ma in 1:3){for(se in 1:length(mon)){out[nrow(out)+1,]<-c(ma,se,yld[,(ma-1)*12+mon[se]])}}}}
  ##############################Fixed sowing month, groups have the same month############################################
  if(type=="Fixed"){
    df<-as.numeric(colMeans(yld))
    df<-as.numeric(colMeans(matrix(df,3,12,byrow=T))) #average across the three maturity groups
    mon<-Find_season_mon(df)
    mon<-mon[!is.na(mon)]
    if(length(mon)>0){for(ma in 1:3){for(se in 1:length(mon)){out[nrow(out)+1,]<-c(ma,se,yld[,(ma-1)*12+mon[se]])}}}}
  #########optimum and fixed sowing month for the three maturity group, each group has different month####################
  if(type=="OptimumFixed"){  
    for(ma in 1:3){
      df<-as.numeric(colMeans(yld[,((ma-1)*12+1):((ma-1)*12+12)]))
      mon<-Find_season_mon(df)
      mon<-mon[!is.na(mon)]
      if(length(mon)>0){for(se in 1:length(mon)){out[nrow(out)+1,]<-c(ma,se,yld[,(ma-1)*12+mon[se]])}}}}
  ###############decade optimum sowing month for the three maturity groups, each group has different months##############
  if(type=="OptimumDecade"){ 
      nyear<-10 #how many years you want to change the sowing performance
      n_times<-round(nrow(df)/nyear)
      for(ma in 1:3){
        df<-yld[,((ma-1)*12+1):((ma-1)*12+12)]
        mean<-aggregate(df, list(rep(1:(nrow(df) %/% nyear + 1), each = nyear, len = nrow(df))), mean)[-1][1:n_times,]
        mon<-Find_season_mon(mean)[,c("1s_m","2s_m")]
        for(se in 1:2){
          mo<-mon[,se]
          line<-c(ma,se,rep(NaN,nrow(df)))
          for(n in 1:n_times){
            if(!is.na(mo[n])){line[(3+(n-1)*nyear):(n*nyear+2)]<-df[(1+(n-1)*nyear):(n*nyear),mo[n]]}} #replace line value if sow in season
          if((nrow(df)>n*nyear)&(!is.na(mo[n]))){line[(3+n*nyear):length(line)]<-df[(1+(n_times*nyear)):nrow(df),mo[n]]} #fill the last row with previous sowing month
          out[nrow(out)+1,]<-line}}}
  ###############yearly optimum sowing month for the three maturity groups, each group has different months##############
  if(type=="OptimumYear"){
    for(ma in 1:3){
      df<-yld[,((ma-1)*12+1):((ma-1)*12+12)]
      sea1<-c(ma,1)
      sea2<-c(ma,2)
      for(i in 1:nrow(df)){
        peaks<-findpeaks(rep(as.numeric(df[i,]),3),peakpat = "[+]{1,}[0]*[-]{1,}",minpeakdistance=4,sortstr=TRUE)
        if(!is.null(peaks)){
          peaks<-as.data.frame(peaks)%>%filter((V2>12)&(V2<25))%>%select(V1)
          if(nrow(peaks)>0){sea1<-c(sea1,peaks[1,1])}else{sea1<-c(sea1,0)}
          if(nrow(peaks)>1){sea2<-c(sea2,peaks[2,1])}else{sea2<-c(sea2,0)}
        }else{
          sea1<-c(sea1,0)
          sea2<-c(sea2,0)
        }}
      out[nrow(out)+1,]<-sea1
      out[nrow(out)+1,]<-sea2}}
  return(out)
}

ExtracYld_csv<-function(runingdir,crop,ssp,fer,startyear,endyear){
  setwd(runningdir)
  grid<-read.csv("./input/Africa_SimGrid_Confirmed_5min_4calibration.csv",header=T) #SIMUNIT 
  grid['se1m']<-grid['ReportedSow_se1m']  #the sowing month of first season
  grid['se2m']<-grid['ReportedSow_se2m']  #the sowing month of second season
  type_n<-c("Reported","Fixed","OptimumFixed","OptimumYear")  #"OptimumDecade"
  for(t in 1:length(type_n)){
    outfile<-paste("./output/",crop,"_yield_",type_n[t],"Sow_",fer,".txt",sep="")
    title<-c("SIMUNIT","matu","sea",paste("yld",seq(startyear,endyear),sep="_")) #totally 61 colnames
    cat(paste(title,collapse=" "),sep="\n",file=outfile,append=FALSE)} #print the colnames
  #loop through all grids to extract yield over years
  df<-read.csv(paste("./simout/mz_result_20crv3_obs_",fer,".csv",sep=""),header=T)
  for(rowi in 1:nrow(grid)){
    yld<-df[df$SIMUNIT==grid$SIMUNIT[rowi],3:38]
    mon<-c(grid$se1m[rowi],grid$se2m[rowi])
    if(nrow(yld)>1){
      for(n in 1:length(type_n)){
        if(((n==1)&(length(mon[!is.na(mon)]>0)))|(n>1)){  #Some rows have no reported sowing month
          out<-GetyearlyYld(yld,mon,type=type_n[n])
          if(nrow(out)>0){write.table(cbind(grid$SIMUNIT[rowi],out),
                      file=paste("./output/",crop,"_yield_",type_n[n],"Sow_",fer,".txt",sep=""),
                      sep=" ", row.names=F, col.name=F,append=T)}}}
    }}}
######################################################################################################################################
crop<-"mz"
runningdir<-"D:/works/AfricaMzSg/"
ssp<-"20crv3_obs"
fer<-"wofer_wa38hi40"
#ExtracYld(runingdir,crop,ssp,startyear=1971,endyear=2021)      #yield is multiple txt file
ExtracYld_csv(runingdir,crop,ssp,fer,startyear=1971,endyear=2021)  #yield is a one csv table
######################################################################################################################################