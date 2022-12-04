
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
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)

  GapPercent=50
  RepNum=1

  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
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
        
        path_obs_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")
        path_pred_zlog <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")
        path_obs <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv")
        path_pred <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv")
        
        if(Percent==0){file_now=path_obs;file_zlog=path_obs_zlog}
        if(Percent!=0){file_now=path_pred;file_zlog=path_pred_zlog}
        
        if(file.exists(file_now)&file.exists(file_zlog)){
          
        # load predicted data
          data_now <- as.matrix(read.table(file = file_now, sep=",", dec=".",header=TRUE))[,-1]
          data_zlog <- as.matrix(read.table(file = file_zlog, sep=",", dec=".",header=TRUE))[,-1]
          
        require(cluster)

          cor_all <- cor(data_now[,-1],method = "pearson")
          cor_zlog <- cor(data_zlog[,-1],method = "pearson")
        
        # -------------------------------------------------
        # create folder
        # -------------------------------------------------
        # create folder

              if(!file.exists(file.path(originData,"analyses"))){
                dir.create(file.path(originData,"analyses"))}
              if(!file.exists(file.path(originData,"analyses","Correl"))){
                dir.create(file.path(originData,"analyses","Correl"))}
              if(!file.exists(file.path(originData,"analyses","Correl",t_choice))){
                dir.create(file.path(originData,"analyses","Correl",t_choice))}
              if(!file.exists(file.path(originData,"analyses","Correl",t_choice,ObsOrTD))){
                dir.create(file.path(originData,"analyses","Correl",t_choice,ObsOrTD))}
              if(!file.exists(file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent))){
                dir.create(file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent))}
              if(!file.exists(file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent,RepNum))){
                dir.create(file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent,RepNum))}
        print(Percent)
        print(cor_all)
        
        write.table(cor_all, file = file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent,RepNum,paste0("Correl.csv")),
                    sep = ",",row.names = TRUE,col.names = TRUE)
        write.table(cor_zlog, file = file.path(originData,"analyses","Correl",t_choice,ObsOrTD,Percent,RepNum,paste0("Correl_zlog.csv")),
                    sep = ",",row.names = TRUE,col.names = TRUE)
        
      }
      }
      }
      }
    }


  
  