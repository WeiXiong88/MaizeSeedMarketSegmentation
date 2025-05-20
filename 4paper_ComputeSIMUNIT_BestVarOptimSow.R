
#######These scripts are for identifying best varietie and its dominant sowing for each simunit#########
########################This work is for a high impact papaers####################################
##Scripts was created and maintained by Wei Xiong (w.xiong@cgiar.org) on 12/2/2024###############
######################################KEY STEPS###################################################

library(dplyr)
library(pracma) #library for function findpeaks()
library(zoo)
library(DescTools) #Mode function
library("cluster")

normalize<-function(x,min_val,max_val){
  #min_val<-1
  #max_val<-8
  if(is.na(x)){return(NA)}else{
    if(x<min_val){
      return(0)
    }else if(x>max_val){
      return(1)
    }else{
      return((x-min_val)/(max_val-min_val))
    }}}


#Function to identify optimum cultivar with dominant sowing month by maximizing Ym/CV namely Ym^2/(100*SD)
optimCultivar0<-function(df){ #input: dataframe or matrix
  result<-rep(NA,4) #season 1 cultivar group, dominant month, season 2 cultivar, dominant month
  #identify season peak month using mean yield
  meanyld<-colMeans(matrix(colMeans(df),nrow=3,byrow=TRUE))  #monthly mean yield of three groups over years
  meanyld<-rollmean(c(meanyld[12],rep(meanyld,3),meanyld[1]),3)  #create a three year series
  peaks<-findpeaks(meanyld,minpeakdistance=4,sortstr=TRUE) #minpeakheight=1,
  if(length(peaks)>0){
    peaks<-as.data.frame(peaks)%>%filter((V2>12)&(V2<25))%>%select(V2) #peak month
    for(i in 1:min(2,nrow(peaks))){  #from season 1 to season 2
      mg<-matrix(NA,3,2) #maturity group Ym2/100SD and dominant sowing month
      pm<-peaks$V2[i]    #peak sowing month
      for(j in 1:3){ #loop through 3 maturity group
        temp<-df[,(1+(j-1)*12):(12+(j-1)*12)]
        temp=cbind(temp,temp,temp)[,(pm-3):(pm+3)]
        maxyld<-apply(temp,1,max)
        dommon<-apply(temp,1,function(row) index(temp)[which.max(row)])+(pm-3)
        dommon<-if (Mode(dommon)[1]%%12==0) 12 else Mode(dommon)[1]%%12
        mg[j,]<-c(mean(maxyld)^2/sd(maxyld),dommon)
      }
      result[(1+(i-1)*2):(2+(i-1)*2)]<-c(which.max(mg[,1])-2,mg[which.max(mg[,1]),2])
    }}
  return(result)
}
#Function to identify optimum cultivar with dominant sowing month by maximizing Ym/CV namely Ym^2/(100*SD) don't limite to mean peaks
optimCultivar1<-function(df){ #input: dataframe or matrix
  result<-rep(NA,4) #season 1 cultivar group, dominant month, season 2 cultivar, dominant month
  #identify season peak month using mean yield
  #meanyld<-colMeans(matrix(colMeans(df),nrow=3,byrow=TRUE))  #monthly mean yield of three groups over years
  #meanyld<-rollmean(c(meanyld[12],rep(meanyld,3),meanyld[1]),3)  #create a three year series
  #peaks<-findpeaks(meanyld,minpeakdistance=4,sortstr=TRUE) #minpeakheight=1,
  mg<-matrix(NA,3,4) #index1,mon1,index2,mon2
  for(i in 1:3){ #loop through 3 maturity group
    yld<-df[,(1+(i-1)*12):(12+(i-1)*12)]
    meanyld<-colMeans(yld) #yld of the maturity group
    meanyld<-rep(meanyld,3)
    peaks<-findpeaks(meanyld,minpeakdistance=2+i,minpeakheight=0.1,sortstr=TRUE)  #early maturing group has 3 month period, 
    if(length(peaks)>0){
      peaks<-as.data.frame(peaks)%>%filter((V2>12)&(V2<25)) #peak month
      peaks<-peaks$V2%%12
      peaks<-ifelse(peaks==0,12,peaks)
      for(j in 1:min(2,length(peaks))){
        mg[i,(j-1)*2+1]<-mean(yld[,peaks[j]])^2/sd(yld[,peaks[j]]) #index
        mg[i,(j-1)*2+2]<-peaks[j]}}}
    #deal with the minor season, if minor season overlap with the primary season, don't consider it
    if(length(which.max(mg[,1]))>0){
      #mg[,3]<-ifelse((pmin(abs(mg[,4]-mg[which.max(mg[,1]),2]),(mg[,4]-mg[which.max(mg[,1]),2])%%12)<4),NA,mg[,3])}  #making sure the month diff is greater than 4
      sea1<-(seq(mg[which.max(mg[,1]),2],mg[which.max(mg[,1]),2]+which.max(mg[,1])+1))%%12.1 #season 1 growth month
      for(j in 1:nrow(mg)){
        if(!is.na(mg[j,4])){
          sea2<-ceiling(seq(mg[j,4],mg[j,4]+j+1)%%12.1) #season 2 growth month
          if(length(intersect(sea1,sea2))>0){mg[j,3]<NA} #if there is common month between season 1 and 2, give the index NA
        }}}
      #mg[,3]<-ifelse((pmin(abs(mg[,4]-mg[which.max(mg[,1]),2]),12-(mg[,4]-mg[which.max(mg[,1]),2]))<=4),NA,mg[,3])}  #making sure the month diff is greater than 4
  result<-c(which.max(mg[,1])-2,mg[which.max(mg[,1]),2],which.max(mg[,3])-2,mg[which.max(mg[,3]),4])
  return(result)
}

#Function to identify optimum cultivar with dominant sowing month by maximizing Ym/CV namely Ym^2/(100*SD) don't limite to mean peaks
optimCultivar<-function(df){ #input: dataframe or matrix
  result<-rep(NA,8) #season 1 cultivar group, dominant month, season 2 cultivar, dominant month
  #find season 1
  mg<-matrix(NA,3,4)
  for(i in 1:3){
    yld=df[,(1+(i-1)*12):(12+(i-1)*12)]
    meanyld<-as.numeric(colMeans(yld))
    se1m<-which.max(meanyld)
    #if(meanyld[se1m]>=1){                        #make sure the primary season yield is greater than 1
      mg[i,]<-c(mean(yld[,se1m])^2/(100*sd(yld[,se1m])),se1m,mean(yld[,se1m]),sum(yld[,se1m]>=0.2)) #1season SI
  #}
  }
   if(length(which.max(mg[,1]))>0) {result[1:4]<-c(which.max(mg[,1])-2,mg[which.max(mg[,1]),2],mg[which.max(mg[,1]),3],mg[which.max(mg[,1]),4])}
  #find season 2
  if(!is.na(result[1]) & !is.na(result[2])){ #there is primary season
    se1p<-seq(result[2],result[2]+result[1]+3)%%12.1  #season 1 period
    mg<-matrix(NA,3,4)
    for(i in 1:3){
      yld<-df[,(1+(i-1)*12):(12+(i-1)*12)]
      meanyld<-as.numeric(colMeans(yld))
      se2p<-seq(se1p[1],se1p[1]-(i+1))%%12
      se2p<-replace(se2p,se2p==0,12)
      se2p<-c(se1p,se2p)
      se2p<-se2p[!duplicated(se2p)]
      meanyld[se2p]<-NA
      se2m<-which.max(meanyld)                 #2 season month
      #if((meanyld[se2m]>=1)){                             #make sure the minor season yield is greater than 1
        mg[i,]<-c(mean(yld[,se2m])^2/sd(yld[,se2m]),se2m,mean(yld[,se2m]),sum(yld[,se2m]>=0.2))
      #}
    }
    if(length(which.max(mg[,1]))>0) {result[5:8]<-c(which.max(mg[,1])-2,mg[which.max(mg[,1]),c(2,3,4)])}
  }
  return(result)
}
#For Figure 1 - cultivar maturity distribution based on historical climate (1971-2021). 
#Optimum cultivar is selected by maximizing (y^2/sd) y is mean yield and sd is standard deviation
#We compute values for first and second season. As sowing might change over years, yield at the sowing month in the 
#   season exhibiting highest yield is selected.
#Return One data
segmentation<-function(infile,outfile,crop,startyear,endyear){
  df<-read.csv(infile,header=T)%>%filter((year>=startyear)&(year<=endyear)) #all grid file
  cat(paste(c("SIMUNIT","sea1_matugroup","sea1_dominantmon","sea1_meanyld","sea1_nGT1",
              "sea2_matugroup","sea2_dominantmon","sea2_meanyld","sea2_nGT1"),collapse=","),sep="\n",file=outfile,append=FALSE)
  for(unit in unique(df$SIMUNIT)){
    df0<-df[df['SIMUNIT']==unit,3:38]
    write(paste(c(unit,optimCultivar(df0)),collapse=","),sep="\n",file=outfile, append=TRUE)
  }
}

crop<-"mz"
start<-c(1971,1981,1991,2001,2011)
end<-c(1980,1990,2000,2010,2021)
fer<-"wofer_wa38hi40"
G1<-"_GT0.2"
infile<-paste("D:/works/AfricaMzSg/simout/mz_result_20crv3_obs_",fer,".csv",sep="")
for (i in 1:5){
  outfile<-paste("D:/works/AfricaMzSg/results/",crop,"/",fer,"/",crop,"_culseg_",fer,"_",as.character(start[i]),
                 "-",as.character(end[i]),G1,".csv",sep="")  #
  segmentation(infile,outfile,"mz",start[i],end[i])
}
outfile<-paste("D:/works/AfricaMzSg/results/",crop,"/",fer,"/",crop,"_culseg_",fer,"_1971-2021",G1,".csv",sep="")  #
segmentation(infile,outfile,crop,1971,2021)
