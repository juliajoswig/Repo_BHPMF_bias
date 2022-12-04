
#------------------------------------------------------------
# define path
#------------------------------------------------------------
setwd("/..")
origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git"
originData = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data"
list.files(file.path(origin,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","helper_scripts","fn_load_functions.R"))
load_functions(origin)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)

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
  
  Percent=1
  t_choice="data"
  ObsOrTD="Obs_obs"
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  gappercents=c(1,5,10,20,30,40,50,60,70)
  
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=2
  TDno=2
  p=7
  td=1
  RepNum=1
  
  for(RepNum in 1){
    for(td in 1:2){
      t_choice <- t_choices[td]
      
      for(p in length(gappercents):1){
        Percent = gappercents[p] 
        
        for(TDno in 2){
          ObsOrTD <- TDnos[TDno]
          
          print(paste(RepNum,ObsOrTD,t_choice,Percent))
          
      #-------------------------------------------------------------------
      # load trait data   
      #-------------------------------------------------------------------
      taxInfo <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"), sep=",")
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clade")
      head(taxInfo)
      traitInfo_obs <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"), sep=",", dec=".",header = TRUE)[,-1]
      head(traitInfo_obs)  
      traitInfo_obs_zlog <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"), sep=",", dec=".",header = TRUE)[,-1]
      head(traitInfo_obs)
      traitInfo_obs_zlog <- cbind(traitInfo_obs[,1],traitInfo_obs_zlog);colnames(traitInfo_obs_zlog)[1] <- "ObservationID"
      head(traitInfo_obs_zlog)  
      load(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","zlog_transform.RData"))
      
  # load the output for trait predictions: mean.csv
      pred_path=file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","mean.csv")
      if(Percent!=0|file.exists(pred_path)){
        traitInfo_pred_zlog <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","mean.csv"), sep="\t", dec=".",header=TRUE)
        traitInfo_pred_zlog <- cbind(traitInfo_obs[,1],traitInfo_pred_zlog);colnames(traitInfo_pred_zlog)[1] <- "ObservationID"
        head(traitInfo_pred_zlog)
        traitInfo_pred_log <- (traitInfo_pred_zlog[,-1]*zlog_trans$sdM)+zlog_trans$meanM
        traitInfo_pred_log <- cbind(traitInfo_obs[,1],traitInfo_pred_log);colnames(traitInfo_pred_log)[1] <- "ObservationID"
        traitInfo_pred <- exp(traitInfo_pred_log[,-1])
        traitInfo_pred <- cbind(traitInfo_obs[,1],traitInfo_pred);colnames(traitInfo_pred)[1] <- "ObservationID"
        head(traitInfo_pred)
        
      #plot(traitInfo_pred[,2],traitInfo_obs[,2])
      #abline(0,1)
      #plot(traitInfo_pred_zlog[,2],traitInfo_obs_zlog[,2])
      #abline(0,1)

    #-------------------------------------------------------------------
    # save pred in both transformations
    #-------------------------------------------------------------------
        print("save pred zlog")
      write.csv(traitInfo_pred,file = 
                  file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred.csv"))
      write.csv(traitInfo_pred_zlog,file = 
                  file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred_zlog.csv"))
      }
      
  #-------------------------------------------------------------------
  # cut to test data only if necessary
  #-------------------------------------------------------------------
      if(t_choice=="data_2"){trait_now=out$trait_guido;obsIDs=out$data_2}
      if(t_choice=="data"){trait_now=out$trait_rainfor;obsIDs=out$data}
  
      # save TD version of mean and sd for Envelope TD
      zlog_trans$meanM <- zlog_trans$meanM[1:length(obsIDs),which(colnames(traitInfo_obs_zlog)[2:ncol(traitInfo_obs_zlog)]%in%trait_now)]
      zlog_trans$sdM <- zlog_trans$sdM[1:length(obsIDs),which(colnames(traitInfo_obs_zlog)[2:ncol(traitInfo_obs_zlog)]%in%trait_now)]
      save(zlog_trans,file=file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","TDzlog_transform.RData"))
      
      traitInfoTD_obs      <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_obs)%in%trait_now))] 
      traitInfoTD_obs_zlog <- traitInfo_obs_zlog[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_obs_zlog)%in%trait_now))] 
  if(Percent!=0|file.exists(pred_path)){ 
      traitInfoTD_pred <- traitInfo_pred[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_pred)%in%trait_now))]  
      traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_pred_zlog)%in%trait_now))] 
  }  
      # Re zlog to TD mean and sd
      load(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),"Obs_obs_TD","data","zlog_transform.RData"))
      traitInfoTD_obs_REzlog <- cbind(traitInfoTD_obs[,1],(log(traitInfoTD_obs[,-1])-zlog_trans$meanM)/zlog_trans$sdM)
      if(Percent!=0|file.exists(pred_path)){ 
        traitInfoTD_pred_REzlog <- cbind(traitInfoTD_pred[,1],(log(traitInfoTD_pred[,-1])-zlog_trans$meanM)/zlog_trans$sdM)}
    
    #-------------------------------------------------------------------
    # save pred in both transformations
    #-------------------------------------------------------------------
    
      write.csv(traitInfoTD_obs_REzlog,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
      if(Percent!=0|file.exists(pred_path)){ 
        write.csv(traitInfoTD_pred_REzlog,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
        }
      
      write.csv(traitInfoTD_obs,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv"))
      write.csv(traitInfoTD_obs_zlog,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_zlog.csv"))
      if(Percent!=0|file.exists(pred_path)){
      write.csv(traitInfoTD_pred_zlog,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_zlog.csv"))
      write.csv(traitInfoTD_pred,file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
    }
  
      print(paste(RepNum,ObsOrTD,t_choice,Percent," done :)"))
        }
      }
    }
  }
  
  