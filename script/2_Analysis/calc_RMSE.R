
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
    ObsOrTDs =  out$ObsOrTDs
    repnums = out$repnums
    gappercents = out$gappercents
    whichDataSet = out$whichDataSet
    ObsSpec = out$ObsSpec
    obsspec = ObsSpec
    preparation = out$preparation
    trait_guido =out$trait_guido
    trait_rainfor =out$trait_rainfor
    colz1 = out$colz1
    colz2 = out$colz2
    gappercents= c("1","5","10","20","30","40","50","60")
    repnums=1:3
    rmse_fun <- function(x,y){
      sqrt(mean((x-y)^2,na.rm=TRUE))
    }
    rmse_funMD <- function(x,y){
      sqrt(median((x-y)^2,na.rm=TRUE))
    }
    
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
t_choice="data"
ObsOrTD="Obs_obs_TD"
gappercents=c(1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
RepNum=3
td=2
p=9
TDno=2
gappercents=80


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
        list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
        
        path_obs0 <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs.csv")
        path_obs <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv")
        path_pred <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv")
        
        path_obs0_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs_zlog.csv")
        path_pred_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_pred_zlog.csv")
        if(ObsOrTD=="Obs_obs"){
          path_obs0_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")
          path_pred_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")
        }

        if(file.exists(path_pred)&file.exists(path_pred_zlog)&
           file.exists(path_obs)&file.exists(path_obs0)&file.exists(path_obs0_zlog)){
          
          # load data
          data_obs <- as.matrix(read.table(file = path_obs, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          data_obs0 <- as.matrix(read.table(file = path_obs0, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          data_obs0_zlog <- as.matrix(read.table(file = path_obs0_zlog, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          data_pred <- as.matrix(read.table(file = path_pred, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          data_pred_zlog <- as.matrix(read.table(file = path_pred_zlog, sep=",", dec=".",header=TRUE))[,-c(1:2)]
          
          # gappyfy
          data_obs0NA <- data_obs0
          data_obs0_zlogNA <- data_obs0_zlog
          data_predNA <- data_pred
          data_pred_zlogNA <- data_pred_zlog
          
          data_obs0NA[!is.na(data_obs)] <- NA
          data_predNA[!is.na(data_obs)] <- NA
          data_obs0_zlogNA[!is.na(data_obs)] <- NA
          data_pred_zlogNA[!is.na(data_obs)] <- NA
          
          #create output rmse matrix with all transformations
          rmse_m <- matrix(NA,ncol=(ncol(data_obs)+1),nrow=4)
          colnames(rmse_m) <- c("Total",colnames(data_obs))
          rownames(rmse_m) <- c("zlog","no_trans","zlogNA","no_transNA")
          rmse_m[1,1] <- rmse_funMD(data_obs0_zlog,data_pred_zlog)
          rmse_m[2,1] <- rmse_funMD(data_obs0,data_pred)
          rmse_m[3,1] <- rmse_funMD(data_obs0_zlogNA,data_pred_zlogNA)
          rmse_m[4,1] <- rmse_funMD(data_obs0NA,data_predNA)
#          t=2
 #         plot(data_obs0_zlog[,t],data_pred_zlog[,t])
#          abline(0,1)
          for(t in 1:ncol(data_obs0_zlog)){
            rmse_m[1,(t+1)] <- rmse_funMD(data_obs0_zlog[,t],data_pred_zlog[,t])
            rmse_m[2,(t+1)] <- rmse_funMD(data_obs0[,t],data_pred[,t])
            rmse_m[3,(t+1)] <- rmse_funMD(data_obs0_zlogNA[,t],data_pred_zlogNA[,t])
            rmse_m[4,(t+1)] <- rmse_funMD(data_obs0NA[,t],data_predNA[,t])
          }
          colnames(rmse_m) <- c("Total",colnames(data_obs))
          rownames(rmse_m) <- c("zlog","no_trans","zlogNA","no_transNA")
          print(rmse_m)
        
          # -------------------------------------------------
          # create folder and save
          # -------------------------------------------------
          
          if(!file.exists(file.path(originData,"analyses"))){
            dir.create(file.path(originData,"analyses"))}
          if(!file.exists(file.path(originData,"analyses","RMSE"))){
            dir.create(file.path(originData,"analyses","RMSE"))}
          if(!file.exists(file.path(originData,"analyses","RMSE",t_choice))){
            dir.create(file.path(originData,"analyses","RMSE",t_choice))}
          if(!file.exists(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD))){
            dir.create(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD))}
          if(!file.exists(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD,Percent))){
            dir.create(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD,Percent))}
          if(!file.exists(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD,Percent,RepNum))){
            dir.create(file.path(originData,"analyses","RMSE",t_choice,ObsOrTD,Percent,RepNum))}
          
          write.table(rmse_m, file = file.path(originData,"analyses","RMSE",t_choice,ObsOrTD,Percent,RepNum,paste0("RMSE.csv")),
                      sep = ",",row.names = TRUE,col.names = TRUE)

        }
      }
    }
  }
}



