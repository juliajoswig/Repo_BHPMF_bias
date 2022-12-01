
# 20210824
# Reorganize the outputs harmonize.
# create harmonic output:
# Un-transformed!

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
  originData = "/Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data/"
  originScript = "/Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git/"
}

list.files(file.path(originScript,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(originScript,"script","helper_scripts","fn_load_functions.R"))
load_functions(originScript)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)
  tsubs <- out$tsubs
  ObsOrTDs = out$ObsOrTDs
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
  
  Percent=0
  RepNum=3
  t_choice="data"
  ObsOrTD="Obs_obs"
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
#  gappercents=c(1,5,10,20,30,40,50,60,70)
  
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=1:3
  TDno=2
  p=1
  
  
  
  #-----------------------------------------------------------
  # create tab to check data availability
  # ---------------------------------------------------------
  l <- list()
  l[[1]] <- gappercents
  l[[2]] <- repnums
  l[[3]] <- t_choices
  l[[4]] <- TDnos
  l[[5]] <- c("taxInfo","traitInfo","traitInfo_zlog","zlog_transform","mean",
              "traitInfo_pred","traitInfo_pred_zlog","TDzlog_transform",
              "traitInfoTD_obs_REzlog","traitInfoTD_pred_REzlog","traitInfoTD_obs",
              "traitInfoTD_obs_zlog","traitInfoTD_pred_zlog","traitInfoTD_pred")
  l[[6]] <- NA
  
  mat <- as.data.frame(expand.grid(l))
  dim(mat)
  
  for(RepNum in 1:3){
    for(td in 1:2){
      t_choice <- t_choices[td]
      
      for(p in length(gappercents):1){
        Percent = gappercents[p] 
        
        for(TDno in 1:2){
          ObsOrTD <- TDnos[TDno]
          
          print(paste(RepNum,ObsOrTD,t_choice,Percent))
          
      #-------------------------------------------------------------------
      # load trait data   
      #-------------------------------------------------------------------
      n =  mat[,1] == Percent &
            mat[,2] == RepNum &
            mat[,3] == t_choice &
            mat[,4] == ObsOrTD &
            mat[,5] == "taxInfo"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfo"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfo_zlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "zlog_transform"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","zlog_transform.RData"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "mean"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","mean.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfo_pred"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfo_pred_zlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred_zlog.csv"))

      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "TDzlog_transform"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","TDzlog_transform.RData"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_obs_REzlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_pred_REzlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_obs"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_obs_zlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_zlog.csv"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_pred_zlog"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_zlog.csv"))
      
      n =  mat[,1] == Percent &
        mat[,2] == RepNum &
        mat[,3] == t_choice &
        mat[,4] == ObsOrTD &
        mat[,5] == "traitInfoTD_pred"
      mat[n,6] <- file.exists(file.path(originData ,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
    
  
      print(paste(RepNum,ObsOrTD,t_choice,Percent," done :)"))
        }
      }
    }
  }
  
  
  mat
  mat[mat[,6]==FALSE,]
  