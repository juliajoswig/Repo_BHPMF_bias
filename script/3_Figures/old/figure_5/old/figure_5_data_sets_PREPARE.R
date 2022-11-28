

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
Version_now="V3"
list.files(file.path(origin,"script","analysis",Version_now))
#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions()

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
TD_choices =  out$TD_choices
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido = out$trait_guido
trait_rainfor = out$trait_rainfor

# for 0% gaps and for 70% gaps:

# load the RMSE 
# calculate per cluster an average
# load the Silhouette data per cluster
# load the mean(sd) data per cluster within silhouettes I think
# load the coefficience of variance data per cluster
# load the number of observations somehow
GapPercent1="org"
RepNum=1
TD_choice="Obs_obs_TD"
trait_sub="guido"

colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")
#COVA
  # load the coefficience of variance data per cluster GUIDO org
cova_now <- rep(NA,19)
repnums=1:10
GapPercent=80
RepNum=1
TD_choice="Obs_obs_TD"
for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
for(GapPercent in gappercents){
  if(GapPercent!=-1){GapPercent1=GapPercent}
  if(GapPercent==-1){GapPercent1="org"}
  for(RepNum in repnums){
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now2)){
      print(GapPercent1)
      cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
      colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
      
      cova_now <- rbind(cova_now,
                        cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust)),
                        rep(TD_choice,nrow(cova_0_rawclust))))
    }
  }
}
}
  print(dim(cova_now))
  head(cova_now)
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  covaG <- cova_now
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,21)
  repnums=1:10
  GapPercent=80
  TD_choice="Obs_obs_TD"
  for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
    for(GapPercent in gappercents){
    if(GapPercent!=-1){GapPercent1=GapPercent}
    if(GapPercent==-1){GapPercent1="org"}
    for(RepNum in repnums){
      path_now1 <- file.path(origin,"data_output","CoVa",
                             "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
      path_now2=file.path(origin,"data_output","CoVa",
                          "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
      if(file.exists(path_now2)){
        print(GapPercent1)
        cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
        colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
        
        cova_now <- rbind(cova_now,
                          cbind(cova_0_rawclust,
                                rep(GapPercent1,nrow(cova_0_rawclust)),
                                rep(TD_choice,nrow(cova_0_rawclust))))
      }
    }
    }
  }
  print(dim(cova_now))
  head(cova_now)
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  
  covaR <- cova_now

  #-------------------------------------
  new_mean_fun <- function(input){
      out=mean(as.numeric(input),na.rm = TRUE)
    return(out)
  }
  
  save(covaR,file=file.path(origin, "data_output","CoVa","rainfor","cova.RData"))
  write.csv(covaR,file=file.path(origin, "data_output","CoVa","rainfor","cova.csv"))

  save(covaG,file=file.path(origin, "data_output","CoVa","guido","cova.RData"))
  write.csv(covaG,file=file.path(origin, "data_output","CoVa","guido","cova.csv"))


  