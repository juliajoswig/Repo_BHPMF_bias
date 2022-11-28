
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

#----------------------------------------------
# define functions
#----------------------------------------------
fun_cluster_mean <- function(x){
  if(sum(!is.na(x))!=0){
    return(mean(x,na.rm=TRUE))
  }else{
    return(NA)
  }
}


x=c(1:10,NA,11)
fun_cluster_dist <- function(x){
  if(sum(!is.na(x))!=0&length(x[!is.na(x)])>1){
    y=NA
    for(j in 1:length(x)){
      if(is.na(x[j])){
        y <- c(y,NA)
      }else{
        y <- c(y,mean(as.matrix(dist(c(x[j],x[-j])))[,1],na.rm=TRUE))
      }
    }      
    y=y[-1]
    y[is.na(x)] = NA
    return(y[-1])
  }else{
    return(rep(NA,length(x)))
  }
}

fun_cluster_CV <- function(x){
  if(sum(!is.na(x))!=0){
    CV=sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100
    return(CV)
  }else{
    return(NA)
  }
}

fun_count <- function(x){
    return(sum(!is.na(x)))
}


repnums=3
TDno=1
RepNum=1
repnums=3
ObsOrTD="Obs_obs_TD"
ObsOrTD="Obs_obs"
TDnos=c("Obs_obs_TD","Obs_obs","Obs_obsTot")



ObsOrTD <- "Obs_obs"
ObsOrTD <- "Obs_obs_TD"
t_choice="data"
ObsOrTD <- "Obs_obsTot"
ObsOrTD <- "Obs_obs_TD"

Percent=80
gappercents=c(1,5,10,20,30,40,50,60,70,80)
gappercents=c(80)
#gappercents=c(80)
RepNum=1
print(paste(RepNum,ObsOrTD,t_choice,Percent))

resdo=TRUE
cordo=TRUE
taxdo=TRUE
distdo=TRUE
MFdo=TRUE

for(RepNum in 1:3){
  for(Percent in gappercents){
    
    
    print("------ info now ----------")
    print(paste(RepNum,ObsOrTD,t_choice,Percent))
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
    if((ObsOrTD=="Obs_obs"|ObsOrTD=="Obs_obs_TD")&MFdo==FALSE){
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
      # taxonomy 
      taxInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD,"data","taxInfoTD.csv"),header=TRUE)[,-1]
      dim(taxInfo)
      
      # observed 
      traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                           ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
      dim(traitInfo_obs_sparse)
      # completely observed  untransformed
      traitInfo_obs <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                    ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
      traitInfo_obs_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         ObsOrTD,"data","traitInfoTD_obs_zlog.csv")))[,-c(1,2)]
      traitInfo_pred <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                     ObsOrTD,"data","traitInfoTD_pred.csv")))[,-c(1,2)]
      traitInfo_pred_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
      dim(traitInfo_obs)
    }
    if((ObsOrTD=="Obs_obs"|ObsOrTD=="Obs_obs_TD")&MFdo){
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
      # taxonomy 
      taxInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD,"data","taxInfoTD.csv"),header=TRUE)[,-1]
      head(taxInfo)
      # replace with one species only
      taxInfoMF <- matrix(NA,nrow=nrow(taxInfo),ncol=5)
      taxInfoMF[,1] <- 1:nrow(taxInfo)
      taxInfoMF[,2] <- rep(1:16,round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,3] <- rep(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,4] <- rep(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfoMF[,5] <- rep(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2),round(nrow(taxInfo)/16,0)+1)[1:nrow(taxInfo)]
      taxInfo <- cbind(taxInfoMF,taxInfo[,c(6,7)])
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clade","PFT","PlantGrowthForm")
      dim(taxInfo)
      
      # observed 
      traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                           ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
      dim(traitInfo_obs_sparse)
      # completely observed  untransformed
      traitInfo_obs <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                    ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
      traitInfo_obs_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         ObsOrTD,"data","traitInfoTD_obs_zlog.csv")))[,-c(1,2)]
      traitInfo_pred <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                     ObsOrTD,"data","traitInfoTD_predMF.csv")))[,-c(1,2)]
      traitInfo_pred_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_pred_REzlogMF.csv")))[,-c(1,2)]
      dim(traitInfo_pred)
    }
    
    if(ObsOrTD=="Obs_obsTot"&MFdo==FALSE){
      ObsOrTD_now="Obs_obs"
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD_now,"data"))
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD_now,"data"))
      # taxonomy 
      taxInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD_now,"data","taxInfo.csv"),header=FALSE)
      funInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD_now,"data","funInfo.csv"),header=TRUE)[,-1]
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clades")
      head(taxInfo)
      head(funInfo)
      taxInfo <- cbind(taxInfo,funInfo[,-1])
      
      # observed 
      traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                           ObsOrTD_now,"data","traitInfo.csv")))[,-c(1,2)]
      head(traitInfo_obs_sparse)
      dim(traitInfo_obs_sparse)
      # completely observed  untransformed
      traitInfo_obs <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                    ObsOrTD_now,"data","traitInfo.csv")))[,-c(1,2)]
      traitInfo_obs_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         ObsOrTD_now,"data","traitInfo_zlog.csv")))[,-c(1)]
      dim(traitInfo_obs)
      traitInfo_pred <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                     ObsOrTD_now,"data","traitInfo_pred.csv")))[,-c(1,2)]
      head(traitInfo_pred)
      dim(traitInfo_pred)
      traitInfo_pred_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD_now,"data","traitInfo_pred_zlog.csv")))[,-c(1,2)]
      dim(traitInfo_pred_zlog)
    }
    if(ObsOrTD=="Obs_obsTot"&MFdo){
      ObsOrTD_now="Obs_obs"
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD_now,"data"))
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD_now,"data"))
      # taxonomy 
      taxInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD_now,"data","taxInfo.csv"),header=FALSE)
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
      
      funInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                    ObsOrTD_now,"data","funInfo.csv"),header=TRUE)[,-1]
      colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clades")
      head(taxInfo)
      head(funInfo)
      taxInfo <- cbind(taxInfo,funInfo[,-1])
      
      # observed 
      traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                           ObsOrTD_now,"data","traitInfo.csv")))[,-c(1,2)]
      head(traitInfo_obs_sparse)
      dim(traitInfo_obs_sparse)
      # completely observed  untransformed
      traitInfo_obs <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                    ObsOrTD_now,"data","traitInfo.csv")))[,-c(1,2)]
      traitInfo_obs_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         ObsOrTD_now,"data","traitInfo_zlog.csv")))[,-c(1)]
      dim(traitInfo_obs)
      traitInfo_pred <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                     ObsOrTD_now,"data","traitInfo_predMF.csv")))[,-c(1,2)]
      head(traitInfo_pred)
      traitInfo_pred_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD_now,"data","traitInfo_pred_zlogMF.csv")))[,-c(1,2)]
      dim(traitInfo_pred_zlog)
    }
    
    print("data loaded")
    
    whichtraits = which(colnames(traitInfo_obs_zlog)%in%unique(c(out$trait_guido,out$trait_rainfor)))
    nb_rows = nrow(traitInfo_obs) + 
      (ncol(traitInfo_obs_zlog)*nrow(traitInfo_obs)) 
    res <- as.data.frame(matrix(NA,ncol=125,nrow=nb_rows))
    colnames(res) <- c("error","Gap","missingness","value_obs","value_pred","value_obs_zlog","value_pred_zlog","nrep","DataSet",
                       "dist_spec_obs_zlog","dist_gen_obs_zlog","dist_fam_obs_zlog","dist_clad_obs_zlog","dist_GF_obs_zlog","dist_PFT_obs_zlog",
                       "dist_spec_pred_zlog","dist_gen_pred_zlog","dist_fam_pred_zlog","dist_clad_pred_zlog","dist_GF_pred_zlog","dist_PFT_pred_zlog",
                       "dist_spec_obs","dist_gen_obs","dist_fam_obs","dist_clad_obs","dist_GF_obs","dist_PFT_obs",
                       "dist_spec_pred","dist_gen_pred","dist_fam_pred","dist_clad_pred","dist_GF_pred","dist_PFT_pred",
                       "nb_spec","nb_gen","nb_fam","nb_clad","nb_GF","nb_PFT",
                       "mean_spec_zlog","mean_gen_zlog","mean_fam_zlog","mean_clad_zlog","mean_GF_zlog","mean_PFT_zlog",#input mean
                       "mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",#input mean
                       "mean_spec_obs","mean_gen_obs","mean_fam_obs","mean_clad_obs","mean_GF_obs","mean_PFT_obs",#input mean
                       "mean_spec_sparse","mean_gen_sparse","mean_fam_sparse","mean_clad_sparse","mean_GF_sparse","mean_PFT_sparse",#input mean
                       "Sil_spec_obs_zlog","Sil_gen_obs_zlog","Sil_fam_obs_zlog","Sil_clad_obs_zlog","Sil_GF_obs_zlog","Sil_PFT_obs_zlog",#input mean
                       "Sil_spec_obs","Sil_gen_obs","Sil_fam_obs","Sil_clad_obs","Sil_GF_obs","Sil_PFT_obs",#input mean
                       "Sil_spec_pred","Sil_gen_pred","Sil_fam_pred","Sil_clad_pred","Sil_GF_pred","Sil_PFT_pred",#input mean
                       "Sil_spec_pred_zlog","Sil_gen_pred_zlog","Sil_fam_pred_zlog","Sil_clad_pred_zlog","Sil_GF_pred_zlog","Sil_PFT_pred_zlog",#input mean
                       "CV_spec_obs","CV_gen_obs","CV_fam_obs","CV_clad_obs","CV_GF_obs","CV_PFT_obs",
                       "CV_spec_pred","CV_gen_pred","CV_fam_pred","CV_clad_pred","CV_GF_pred","CV_PFT_pred",
                       "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap","nb_GF_gap","nb_PFT_gap",
                       "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred","dist_AVGF_pred","dist_AVPFT_pred",
                       "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs","dist_AVGF_obs","dist_AVPFT_obs",
                       "trait","ID","Species","Genus","Family","Clade","GF","PFT")
    
    print("res created")
    #-----------------------------------------
    # correlations deviation to distance
    #-----------------------------------------
    if(cordo&ObsOrTD!="Obs_obsTot"){
      t=1
      for(t in 1:ncol(traitInfo_obs)){
        # -----------------------------------------
        # define the place in the output matrix
        ix_trait1 = 1:nrow(traitInfo_obs)
        ix_trait1 = ix_trait1  + 
          ((t-1)*nrow(traitInfo_obs)) 
        ix_trait = ix_trait1      
        #----------------------------------------------
        
        t1=2
        for(t1 in 1:ncol(traitInfo_obs)){
          for(t in t){
            if(t1!=t){
              colnm4 = paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_pred")
              res <- add_col_to_res(new.col.names = colnm4,input = res)
              # build model from observed data
              dat_now=as.data.frame(traitInfo_obs[,c(t1,t)])
              if(sum(complete.cases(dat_now))>10){
                names(dat_now) <- c("t1","t")
                lm_now <- lm(t1~t,data = dat_now) # lm produced from NONsparse data!
                t_lm <- (traitInfo_pred[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
                res[ix_trait,colnames(res)%in%colnm4] <- t_lm
              }
            }
            if(t1!=t){
              colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_pred_zlog")
              res <- add_col_to_res(new.col.names = colnm4,input = res)
              # build model from observed data
              dat_now = as.data.frame(traitInfo_obs_zlog[,c(t1,t)])
              if(sum(complete.cases(dat_now))>10){
                names(dat_now) <- c("t1","t")
                lm_now <- lm(t1~t,data = dat_now) # lm produced from NONsparse data!
                t_lm <- (traitInfo_pred_zlog[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
                res[ix_trait,colnames(res)%in%colnm4] <- t_lm
              }
            }
            
            if(t1!=t){
              colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_obs")
              res <- add_col_to_res(new.col.names = colnm4,input = res)
              #observed
              dat_now=as.data.frame(traitInfo_obs[,c(t1,t)])
              if(sum(complete.cases(dat_now))>10){
                names(dat_now) <- c("t1","t")
                lm_now <- lm(t1~t,data = dat_now,na.action = "na.omit") # lm produced from NONsparse data!
                t_lm_obs <- (traitInfo_obs[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
                res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
              }
            }
            if(t1!=t){
              colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_obs_zlog")
              res <- add_col_to_res(new.col.names = colnm4,input = res)
              #observed
              dat_now = as.data.frame(traitInfo_obs_zlog[,c(t1,t)])
              if(sum(complete.cases(dat_now))>10){
                names(dat_now) <- c("t1","t")
                lm_now <- lm(t1~t,data = dat_now,na.action = "na.omit") # lm produced from NONsparse data!
                t_lm_obs <- (traitInfo_obs_zlog[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
                res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
              }
            }
          }
        }
        
        
      }  
    }
    #-----------------------------------------
    # taxonomy
    #-----------------------------------------
    
    if(distdo){
      
      #----------------------------------------------
      # Individuals within species
      mean_pred_zlog <- list()
      mean_obs_zlog <- list()
      mean_pred <- list()
      mean_obs <- list()
      
      dist_pred_zlog <- list()
      dist_obs_zlog <- list()
      dist_pred <- list()
      dist_obs <- list()
      
      nb_obs <- list()

      CV_pred_zlog <- list()
      CV_obs_zlog <- list()
      CV_pred <- list()
      CV_obs <- list()
      
      Sil_pred_zlog <- list()
      Sil_obs_zlog <- list()
      Sil_pred <- list()
      Sil_obs <- list()
      
      i=2
      for(i in 2:7){   
        print(i)
        
        mean_pred_zlog[[i]] <- aggregate(traitInfo_pred_zlog,list(taxInfo[,i]),fun_cluster_mean)
        mean_pred[[i]] <- aggregate(traitInfo_pred_zlog,list(taxInfo[,i]),fun_cluster_mean)
        mean_obs_zlog[[i]] <- aggregate(traitInfo_obs_zlog,list(taxInfo[,i]),fun_cluster_mean)
        mean_obs[[i]] <- aggregate(traitInfo_obs,list(taxInfo[,i]),fun_cluster_mean)
        
        CV_pred[[i]] <- aggregate(traitInfo_pred,list(taxInfo[,i]),fun_cluster_CV)
        CV_obs[[i]] <- aggregate(traitInfo_obs,list(taxInfo[,i]),fun_cluster_CV)
        
        nb_obs[[i]] <- aggregate(traitInfo_obs,list(taxInfo[,i]),fun_count)
        
        t=1
        for(t in whichtraits){
          
          print("------------trait:")
          print(t)
          ix_trait1 = 1:nrow(traitInfo_obs_zlog)
          ix_trait1 = ix_trait1  + 
            ((t-1)*nrow(traitInfo_obs_zlog)) 
          ix_trait = ix_trait1   
          print("------------from to:")
          print(min(ix_trait,na.rm = TRUE))
          print(max(ix_trait,na.rm = TRUE))
          
          res$trait[ix_trait] <- colnames(traitInfo_obs_zlog)[t]
          res$ID[ix_trait] <- as.vector(taxInfo$ObservationID)
          res$Species[ix_trait] <- as.vector(taxInfo$Species)
          res$Genus[ix_trait] <- as.vector(taxInfo$Genus)
          res$Family[ix_trait] <- as.vector(taxInfo$Family)
          res$Clade[ix_trait] <- as.vector(taxInfo$Clades)
          res$GF[ix_trait] <- as.vector(taxInfo$PlantGrowthForm)
          res$PFT[ix_trait] <- as.vector(taxInfo$PFT)
          res$error[ix_trait] <- traitInfo_pred_zlog[,t] - traitInfo_obs_zlog[,t]
          res$Gap[ix_trait] <- is.na(traitInfo_obs_sparse[,t])
          res$value_obs_zlog[ix_trait] <- traitInfo_obs_zlog[,t]
          res$value_pred_zlog[ix_trait] <- traitInfo_pred_zlog[,t]
          res$value_obs[ix_trait] <- traitInfo_obs[,t]
          res$value_pred[ix_trait] <- traitInfo_pred[,t]
          res$missingness[ix_trait] <- Percent
          res$DataSet[ix_trait] <- t_choice
          res$nrep[ix_trait] <- as.vector(RepNum)
          
          print("put into res")
          length(match(x=taxInfo[,i],table=unique(taxInfo[,i])))
          if(i==2){
            res <- add_col_to_res(new.col.names = paste0("CV_spec_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_spec_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[ix_trait,colnames(res)==paste0("CV_spec_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)==paste0("CV_spec_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_spec"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
          if(i==3){
            res <- add_col_to_res(new.col.names = paste0("CV_gen_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_gen_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[,colnames(res)==paste0("CV_gen_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[,colnames(res)==paste0("CV_gen_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_gen"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
          if(i==4){
            res <- add_col_to_res(new.col.names = paste0("CV_fam_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_fam_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[ix_trait,colnames(res)==paste0("CV_fam_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)==paste0("CV_fam_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_fam"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
          if(i==4){
            res <- add_col_to_res(new.col.names = paste0("CV_clad_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_clad_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[ix_trait,colnames(res)==paste0("CV_clad_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)==paste0("CV_clad_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_clad"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
          if(i==6){
            res <- add_col_to_res(new.col.names = paste0("CV_GF_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_GF_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[ix_trait,colnames(res)==paste0("CV_GF_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)==paste0("CV_GF_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_GF"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
          if(i==7){
            res <- add_col_to_res(new.col.names = paste0("CV_PFT_obs_",colnames(traitInfo_obs)[t]),input = res)
            res <- add_col_to_res(new.col.names = paste0("CV_PFT_pred_",colnames(traitInfo_obs)[t]),input = res)
            res[ix_trait,colnames(res)==paste0("CV_PFT_obs_",colnames(traitInfo_obs)[t])] <- CV_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)==paste0("CV_PFT_pred_",colnames(traitInfo_obs)[t])]  <- CV_pred[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
            res[ix_trait,colnames(res)=="nb_PFT"]  <- nb_obs[[i]][match(x=taxInfo[,i],table=unique(taxInfo[,i])),(t+1)]
          }
        }
      }
      
      
    }
    
    
    if(!file.exists(file.path(origin,"_2021","data","analyses"))){
      dir.create(file.path(origin,"_2021","data","analyses"))}
    if(!file.exists(file.path(origin,"_2021","data","analyses","Point_wise"))){
      dir.create(file.path(origin,"_2021","data","analyses","Point_wise"))}
    if(!file.exists(file.path(origin,"_2021","data","analyses","Point_wise",RepNum))){
      dir.create(file.path(origin,"_2021","data","analyses","Point_wise",RepNum))}
    if(!file.exists(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice))){
      dir.create(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice))}
    if(!file.exists(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD))){
      dir.create(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD))}
    if(!file.exists(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent))){
      dir.create(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent))}
    res <- res[,colSums(!is.na(res))!=0]
    res <- res[rowSums(!is.na(res))!=0,]
    
    write.csv(res,file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_indMF.csv")))
    
    print("res written.")

    }#Percent
}

#TD <- res$Sil_spec_obs_zlog
#ENV <- res$Sil_spec_obs_zlog
head(res$Sil_spec_obs_zlog)
head(resENV$Sil_spec_obs_zlog)

RepNum=1
t_choice="data"
ObsOrTD="Obs_obsMF"
Percent=80
resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_indMF.csv")))
dim(resENV)
plot(resENV$value_obs_zlog,resENV$value_pred_zlog)

RepNum=1
t_choice="data"
ObsOrTD="Obs_obsTot"
Percent=80
resENV2 <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_indMF.csv")))
head(resENV$Sil_spec_obs_zlog)
plot(resENV2$value_obs_zlog,resENV2$value_pred_zlog)
plot(resENV$value_obs,resENV$value_pred)


par(mfrow=c(1,3))
boxplot(resENV$CV_spec_obs_SLA)
boxplot(res$CV_spec_obs_SLA)

par(mfrow=c(1,3))
boxplot(resENV$CV_spec_obs_SLA,ylim=c(0,120))
boxplot(resENV$CV_spec_pred_SLA,ylim=c(0,120))
boxplot(res$CV_spec_pred_SLA,ylim=c(0,120))

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021MF.csv")))
head(res$Sil_spec_obs_zlog)
#install.packages("cluster")
#    library(cluster)
#    tax_now1 <- taxInfo[,i]
#    tax_now1[is.na(tax_now1)] <- 1:sum(is.na(tax_now1))
#    cats <- factor(tax_now1)
#    ranks <- rank(-table(cats), ties.method="first")
#    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
#    taxnow <- DF[,2]
#    si2 <- silhouette(taxnow, dist(traitInfo_obs,method = "canberra"))
#    si <- summary(si2)
#    Sil_obs[[i]] <- cbind(as.vector(unique(tax_now1)),as.vector(si$clus.avg.widths))
#    si2 <- silhouette(taxnow, dist(traitInfo_pred,method = "canberra"))
#    si <- summary(si2)
#    Sil_pred[[i]] <- cbind(as.vector(unique(tax_now1)),as.vector(si$clus.avg.widths))
#    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
#    si <- summary(si2)
#    Sil_obs_zlog[[i]] <- cbind(as.vector(unique(tax_now1)),as.vector(si$clus.avg.widths))
#    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
#    si <- summary(si2)
#    Sil_pred_zlog[[i]] <- cbind(as.vector(unique(tax_now1)),as.vector(si$clus.avg.widths))

