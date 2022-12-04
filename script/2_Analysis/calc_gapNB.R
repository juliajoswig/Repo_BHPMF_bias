
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
# Loop through all runs and calculate the missing entries
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
        list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
        
        path_obs <- file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs.csv")

        if(file.exists(path_obs)){
          
          # load data
          data_obs <- as.matrix(read.table(file = path_obs, sep=",", dec=".",header=TRUE))[,-c(1:2)]
           
          
          GapNb <- matrix(NA,ncol=(ncol(data_obs)+1),nrow=2)
          GapNb[1,] <- c(sum(is.na(data_obs)),colSums(is.na(data_obs)))
          GapNb[2,] <- c(sum(is.na(data_obs))/(nrow(data_obs)*ncol(data_obs))*100,
                         colSums(is.na(data_obs))/nrow(data_obs)*100)
          
          colnames(GapNb) <- c("Total",colnames(data_obs))
          rownames(GapNb) <- c("Nb_gaps","%_gaps")
          print(GapNb)
          
          # -------------------------------------------------
          # create folder
          # -------------------------------------------------
          
          if(!file.exists(file.path(originData,"analyses"))){
            dir.create(file.path(originData,"analyses"))}
          if(!file.exists(file.path(originData,"analyses","GapNb"))){
            dir.create(file.path(originData,"analyses","GapNb"))}
          if(!file.exists(file.path(originData,"analyses","GapNb",t_choice))){
            dir.create(file.path(originData,"analyses","GapNb",t_choice))}
          if(!file.exists(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD))){
            dir.create(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD))}
          if(!file.exists(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD,Percent))){
            dir.create(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD,Percent))}
          if(!file.exists(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD,Percent,RepNum))){
            dir.create(file.path(originData,"analyses","GapNb",t_choice,ObsOrTD,Percent,RepNum))}
          
          write.table(GapNb, file = file.path(originData,"analyses","GapNb",t_choice,ObsOrTD,Percent,RepNum,paste0("GapNb.csv")),
                      sep = ",",row.names = TRUE,col.names = TRUE)
          
        }
      }
    }
  }
}



