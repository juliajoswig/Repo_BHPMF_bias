
# get per observed point value: 
# - error (pred-obs)
# - average distance to all other points observed 
# - Silhouette index for this group (or calculate somewhere else...?)  
# - ORIGINAL Silhouette index for this group (or calculate somewhere else...?)  
# - DEVIATION Silhouette index for this group (or calculate somewhere else...?)  
# - number of values available
# - original number of values available

#------------------------------------------------------------
# define path
#------------------------------------------------------------
is.it.on.cluster=FALSE
if(is.it.on.cluster){
  setwd("/..")
  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("work_1","2016_GapFilling")}
if(!is.it.on.cluster){
  setwd("/..")
  origin = "Volumes/bgi/work_1/2016_GapFilling"
}
Version_now="V2"
list.files(file.path(origin,"_2021","script",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
  t_choices <- out$tsubs
  TDnos = out$TD_choices
  repnums = out$repnums
  gappercents = out$gappercents
  whichDataSet = out$whichDataSet
  ObsSpec = out$ObsSpec
  obsspec = ObsSpec
  preparation = out$preparation
  trait_guido = out$trait_guido
  trait_rainfor = out$trait_rainfor
  colz1 = out$colz1
  colz2 = out$colz2
  new.mean.fun = out$new.mean.fun
  new.sd.fun = out$new.sd.fun
  add_col_to_res <- out$add_col_to_res
  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  


  repnums=3
  TDno=2
  RepNum=1
  repnums=3
  t_choice="data"
  ObsOrTD="Obs_obs"
  ObsOrTD="Obs_obs_TD"
  TDnos=c("Obs_obs_TD","Obs_obs")
  ObsOrTD <- TDnos[TDno]

        

Percent=80
gappercents=c(1,5,10,20,30,40,50,60,70,80)


print(paste(RepNum,ObsOrTD,t_choice,Percent))
gap_m <- as.data.frame(matrix(NA,ncol=2,nrow=length(gappercents)*3))
gap_m[,1] <- rep(gappercents,3)
gap_m[,2] <- rep(1:3,length(gappercents))
colnames(gap_m) <- c("Gaps","Reps")
for(RepNum in 1:3){
  for(p in 1:length(gappercents)){
    Percent <- gappercents[p]
    print("------ % gaps now ----------")
    print(Percent)
    
    traitInfo <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                         ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
    head(traitInfo)
    t2=1
    t=1
for(t in 1:ncol(traitInfo)){
  for(t2 in 1:ncol(traitInfo)){
    if(t!=t2){
      nm=paste0(colnames(traitInfo)[c(t)],"_",colnames(traitInfo)[c(t2)])
      gap_m <- add_col_to_res(input = gap_m,new.col.names = nm)
      gap_m[gap_m[,1]==Percent&gap_m[,2]==RepNum,which(colnames(gap_m)%in%nm)]<- sum(complete.cases(traitInfo[,c(t,t2)]))
    }
  }
  }
  }
}
gap_m
gap_m[order(gap_m[,2]),]

gap_md <- aggregate(x = gap_m[,3:ncol(gap_m)],by = list(gap_m[,1]),FUN = median)
gap_sd <-aggregate(x = gap_m[,3:ncol(gap_m)],by = list(gap_m[,1]),FUN = sd)

require(xtable)
xtable(gap_md,include.rownames=FALSE)
print(xtable(tab_compare), include.rownames=FALSE)

apply(gap_md,1,median)






print(paste(RepNum,ObsOrTD,t_choice,Percent))
gap_m <- matrix(NA,ncol=8,nrow=length(gappercents)*3)
gap_m[,1] <- rep(gappercents,3)
gap_m[,2] <- rep(1:3,length(gappercents))
for(RepNum in 1:3){
for(p in 1:length(gappercents)){
  Percent <- gappercents[p]
  print("------ % gaps now ----------")
  print(Percent)
  
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
#    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
#    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
    # observed 
    traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                          ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
    head(traitInfo_obs_sparse)
    
    gap_m[gap_m[,1]==Percent&gap_m[,2]==RepNum,3:ncol(gap_m)]<- colSums(is.na(traitInfo_obs_sparse))
    
  }
}


gap_m[order(gap_m[,1]),]
