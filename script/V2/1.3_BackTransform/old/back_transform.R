#------------------------------------------------------------
# define path
#------------------------------------------------------------
is.it.on.cluster=FALSE
if(is.it.on.cluster){
  #  setwd("/..")
  #  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("Net","Groups","BGI","work_1","2016_GapFilling")
  origin2=file.path("work_1","2016_GapFilling","_2021","data","_runs")}
if(!is.it.on.cluster){
  setwd("/..")
  origin = "Volumes/bgi/work_1/2016_GapFilling"
}
Version_now="V1"
  Version_now="V2"
  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3



  td=1
  TDno=1
  RepNum=3
  rn=1
  p=2
  t_choice <- t_choices[td]
  Percent = gappercents[p] 
  ObsOrTD <- TDnos[TDno]
for(RepNum in 1:repnums){
  for(td in 1:2){
    t_choice <- t_choices[td]
    for(p in 2:length(gappercents)){
    Percent = gappercents[p] 
    for(TDno in 1){
      ObsOrTD <- TDnos[TDno]
      
      print(paste(RepNum,ObsOrTD,t_choice,Percent))
      list.files(file.path(origin,"_2021","data","_runs"))
      traitInfo_pred_zlog <- as.matrix(read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),
                                                                 t_choice,paste0("p_",Percent),ObsOrTD,"data/mean.csv"),sep="\t",header=TRUE))
      traitInfo_pred_zlog <- traitInfo_pred_zlog[,-1]
      traitInfo_obs_zlog <- as.matrix(read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),
                                                                 t_choice,paste0("p_",Percent),ObsOrTD,"data/traitInfo.csv"),header=TRUE))
      traitInfo_obs_zlog <- traitInfo_obs_zlog[,-1]
      head(traitInfo_pred)
      head(traitInfo_obs_zlog)
      t=2
      traitInfo_pred <- exp((traitInfo_pred_zlog*apply(traitInfo_obs_zlog[,2:ncol(traitInfo_obs_zlog)],2,sd,na.rm=TRUE))+ 
                        apply(traitInfo_obs_zlog[,2:ncol(traitInfo_obs_zlog)],2,mean,na.rm=TRUE))
    }
}
}
}

