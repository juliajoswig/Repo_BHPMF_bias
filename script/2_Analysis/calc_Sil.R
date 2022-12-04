
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
  

GapPercent=50
RepNum=1

gappercents=c(1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
#-------------------------------------------------------------------
# Loop through all runs and calculate the indices
#-------------------------------------------------------------------

GapPercent=50
RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
td=1
p=3

for(RepNum in 1:3){
  for(td in 1:2){
    t_choice <- t_choices[td]
    
    for(p in 1:length(gappercents)){
      Percent = gappercents[p] 
      
      for(TDno in 1:2){
        ObsOrTD <- TDnos[TDno]
        
        print(paste(RepNum,ObsOrTD,t_choice,Percent))
        
        # define dat path
        list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
        
        path_obs0_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")
        #path_obs0_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")
        path_obs0 <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs.csv")
        path_obs <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv")
        path_pred_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")
        path_pred <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv")
        

        if(file.exists(path_pred)&file.exists(path_pred_zlog)&
           file.exists(path_obs)&file.exists(path_obs0)&file.exists(path_obs0_zlog)){
          
          # load data
          traitInfo_obs_sparse <- as.matrix(read.table(file = path_obs, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          traitInfo_obs <- as.matrix(read.table(file = path_obs0, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          traitInfo_obs_zlog <- as.matrix(read.table(file = path_obs0_zlog, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          traitInfo_pred <- as.matrix(read.table(file = path_pred, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          traitInfo_pred_zlog <- as.matrix(read.table(file = path_pred_zlog, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          
        if(sum(!is.na(traitInfo_pred_zlog))>0){# remove, once all runs are correct
          
          # taxonomy 
          tmp <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                  ObsOrTD,"data","traitInfoTD_obs.csv"),header=TRUE))[,-1]
          taxInfo <- read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                        ObsOrTD,"data","taxInfo.csv"),header=FALSE)
          colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clades")
          taxInfo <- taxInfo[which(taxInfo$ObservationID%in%tmp$ObservationID),]
          
          funInfo <- read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                        ObsOrTD,"data","funInfo.csv"),header=TRUE)[,-1]
          funInfo <- funInfo[which(funInfo$ObservationID%in%tmp$ObservationID),]
          taxInfo <- cbind(taxInfo,funInfo[,-1])
          write.csv(taxInfo,file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                      ObsOrTD,"data","taxInfoTD.csv"))
          dim(taxInfo)
          rm("tmp")
          
          library(cluster) 
          sil_l <- list()
          cl=2
          for(cl in 2:ncol(taxInfo)){
          # Silhouettes
          cats <- factor(taxInfo[,cl])
          ranks <- rank(-table(cats), ties.method="first")
          DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
          taxnow <- DF[,2]
          ixG <- !is.na(taxnow)
          si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
          si <- summary(si2)
          sil_l[[(cl-1)]] <- cbind(as.vector(unique(taxInfo[ixG,cl])),si$clus.avg.widths,si$clus.sizes)
          si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
          si <- summary(si2)
          sil_l[[(cl-1)]] <- cbind(sil_l[[(cl-1)]],si$clus.avg.widths)
          colnames(sil_l[[(cl-1)]]) <- c(colnames(taxInfo)[cl],"obs","nb","pred")
          }
          
          # -------------------------------------------------
          # create folder
          # -------------------------------------------------
          
          if(!file.exists(file.path(originData,"analyses"))){
            dir.create(file.path(originData,"analyses"))}
          if(!file.exists(file.path(originData,"analyses","Silhouette"))){
            dir.create(file.path(originData,"analyses","Silhouette"))}
          if(!file.exists(file.path(originData,"analyses","Silhouette",t_choice))){
            dir.create(file.path(originData,"analyses","Silhouette",t_choice))}
          if(!file.exists(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD))){
            dir.create(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD))}
          if(!file.exists(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD,Percent))){
            dir.create(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD,Percent))}
          if(!file.exists(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD,Percent,RepNum))){
            dir.create(file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD,Percent,RepNum))}
          
          save(sil_l,file=file.path(originData,"analyses","Silhouette",t_choice,ObsOrTD,Percent,RepNum,"Silhouette.RData"))

        }else{print("redo run ! Something went wrong")}# remove, once all runs are correct
          }else{
          print("does not exist")
        }
      }
    }
  }
}



