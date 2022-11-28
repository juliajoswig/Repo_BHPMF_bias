

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
  
  Percent=80
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs"
  ObsOrTD="Obs_obs_TD"
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  gappercents=c(80)
  
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  TDno=2
  p=1
  
  for(RepNum in 1:3){
    for(td in 1){
      t_choice <- t_choices[td]
      
      for(p in length(gappercents):1){
        Percent = gappercents[p] 
        
        for(TDno in 2){
          ObsOrTD <- TDnos[TDno]
          
          print(paste(RepNum,ObsOrTD,t_choice,Percent))
          
      #-------------------------------------------------------------------
      # load trait data   
      #-------------------------------------------------------------------
      taxInfo <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"), sep=",")
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clade")
      head(taxInfo)
      # replace with one species only
      taxInfoMF <- matrix(NA,nrow=nrow(taxInfo),ncol=5)
      taxInfoMF[,1] <- 1:nrow(taxInfo)
      taxInfoMF[,2] <- rep(1:16,round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,3] <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,4] <- rep(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,5] <- rep(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfo <- taxInfoMF
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clade")
      
      traitInfo_obs <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"), sep=",", dec=".",header = TRUE)[,-1]
      head(traitInfo_obs)  
      traitInfo_obs_zlog <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"), sep=",", dec=".",header = TRUE)[,-1]
      head(traitInfo_obs)
      traitInfo_obs_zlog <- cbind(traitInfo_obs[,1],traitInfo_obs_zlog);colnames(traitInfo_obs_zlog)[1] <- "ObservationID"
      head(traitInfo_obs_zlog)  
      load(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","zlog_transform.RData"))
      
  # load the output for trait predictions: mean.csv
      pred_path=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","meanMF.csv")
      if(Percent!=0|file.exists(pred_path)){
        traitInfo_pred_zlog <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","meanMF.csv"), sep="\t", dec=".",header=TRUE)
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
                  file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_predMF.csv"))
      write.csv(traitInfo_pred_zlog,file = 
                  file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred_zlogMF.csv"))
      }
      
  #-------------------------------------------------------------------
  # cut to test data only if necessary
  #-------------------------------------------------------------------
      if(t_choice=="data_2"){trait_now=out$trait_guido;obsIDs=out$data_2}
      if(t_choice=="data"){trait_now=out$trait_rainfor;obsIDs=out$data}
  
      # save TD version of mean and sd for Envelope TD
      zlog_trans$meanM <- zlog_trans$meanM[1:length(obsIDs),which(colnames(traitInfo_obs_zlog)[2:ncol(traitInfo_obs_zlog)]%in%trait_now)]
      zlog_trans$sdM <- zlog_trans$sdM[1:length(obsIDs),which(colnames(traitInfo_obs_zlog)[2:ncol(traitInfo_obs_zlog)]%in%trait_now)]
      save(zlog_trans,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","TDzlog_transformMF.RData"))
      
      traitInfoTD_obs      <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_obs)%in%trait_now))] 
      traitInfoTD_obs_zlog <- traitInfo_obs_zlog[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_obs_zlog)%in%trait_now))] 
  if(Percent!=0|file.exists(pred_path)){ 
      traitInfoTD_pred <- traitInfo_pred[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_pred)%in%trait_now))]  
      traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(obsIDs),c(1,which(colnames(traitInfo_pred_zlog)%in%trait_now))] 
  }  
      # Re zlog to TD mean and sd
      load(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),"Obs_obs_TD","data","zlog_transform.RData"))
      traitInfoTD_obs_REzlog <- cbind(traitInfoTD_obs[,1],(log(traitInfoTD_obs[,-1])-zlog_trans$meanM)/zlog_trans$sdM)
      if(Percent!=0){ 
        traitInfoTD_pred_REzlog <- cbind(traitInfoTD_pred[,1],(log(traitInfoTD_pred[,-1])-zlog_trans$meanM)/zlog_trans$sdM)}
    
    #-------------------------------------------------------------------
    # save pred in both transformations
    #-------------------------------------------------------------------
    
      write.csv(traitInfoTD_obs_REzlog,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlogMF.csv"))
      if(Percent!=0|file.exists(pred_path)){ 
        write.csv(traitInfoTD_pred_REzlog,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlogMF.csv"))}
      
      write.csv(traitInfoTD_obs,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obsMF.csv"))
      write.csv(traitInfoTD_obs_zlog,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_zlogMF.csv"))
      if(Percent!=0|file.exists(pred_path)){
      write.csv(traitInfoTD_pred_zlog,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_zlogMF.csv"))
      write.csv(traitInfoTD_pred,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_predMF.csv"))
    }
  
        }
      }
    }
  }
  
  