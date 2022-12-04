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
load_functions(origin = origin,Version_now = "V2")

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
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
  colz1 =out$colz1
  colz2 =out$colz2

gappercents= c("1","5","10","20","30","40","50","60","70","80")
  #-------------------------------------------------------------------
  # Loop through all runs and calculate the indices
  #-------------------------------------------------------------------
GapPercent=50
RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
gappercents=c(1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

td=1
TDno=1
RepNum=1
rn=1
p=10
t_choice <- t_choices[td]
Percent = gappercents[p] 
ObsOrTD <- TDnos[TDno]


for(RepNum in 1:3){
  for(td in 1:2){
    t_choice <- t_choices[td]
    
    for(p in 2:length(gappercents)){
      Percent = gappercents[p] 
      
      for(TDno in 1){
        ObsOrTD <- TDnos[TDno]
        
      
          print(paste(RepNum,ObsOrTD,t_choice,Percent))
        
list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
            dat_path <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")

            if(file.exists(dat_path)){

            # load predicted data
            traitInfo_pred_zlog <- as.matrix(read.table(file = dat_path, sep=",", dec=".",head=TRUE))[,-1]
            head(traitInfo_pred_zlog)
            
            #install.packages("farver")
            require(stats)
            pca_data <- prcomp(traitInfo_pred_zlog[,-1])
            pca_sry<- summary(pca_data)
            rownms <- c("Variance explained",rownames(pca_sry$rotation))
            pca_out<- matrix(NA,ncol=5,nrow=length(rownms))
            pca_coord <- matrix(NA,ncol=5,nrow=length(rownms))
            rownames(pca_out) <- rownms
            colnames(pca_out) <- c("Axis1","Axis2","Axis3","Axis4","Axis5")
            pca_out[2:nrow(pca_out),] <- pca_sry$rotation[,1:5]
            pca_out[1,] <- pca_sry$importance[2,1:5]

            #-------------------------------------------------
            # create folder
            if(!file.exists(file.path(origin,"_2021","data","analyses"))){
              dir.create(file.path(origin,"_2021","data","analyses"))}
            if(!file.exists(file.path(origin,"_2021","data","analyses","PCA"))){
                dir.create(file.path(origin,"_2021","data","analyses","PCA"))}
              if(!file.exists(file.path(origin,"_2021","data","analyses","PCA",t_choice))){
              dir.create(file.path(origin,"_2021","data","analyses","PCA",t_choice))}
            if(!file.exists(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD))){
              dir.create(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD))}
            if(!file.exists(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent))){
              dir.create(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent))}
            if(!file.exists(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent,RepNum))){
              dir.create(file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent,RepNum))}
            print(Percent)
            write.table(pca_out,file= 
                          file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent,RepNum,paste0("PCA.csv")),
                        sep = ",",row.names = TRUE,col.names = TRUE)
            
   
            
            
      }
      }
      }
      }
    }


  
  
  