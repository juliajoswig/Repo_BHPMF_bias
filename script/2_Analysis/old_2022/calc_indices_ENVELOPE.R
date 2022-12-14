
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
  


  repnums=3
  TDno=1
  RepNum=1
  repnums=3
  ObsOrTD="Obs_obs_TD"
  ObsOrTD="Obs_obs"
  TDnos=c("Obs_obs_TD","Obs_obs")
  
  
  
  ObsOrTD <- "Obs_obs_TD"
  t_choice="data"

Percent=80
gappercents=c(1,5,10,20,30,40,50,60,70,80)
gappercents=c(60,70,80)
#gappercents=c(80)
print(paste(RepNum,ObsOrTD,t_choice,Percent))

resdo=TRUE
cordo=TRUE
taxdo=TRUE
distdo=TRUE

for(RepNum in 2:3){
for(Percent in gappercents){
  
  print("------ % gaps now ----------")
  print(Percent)
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
  {
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

    nb_rows = nrow(traitInfo_obs) + 
      ((ncol(traitInfo_obs))*nrow(traitInfo_obs)) 
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
    
  #-----------------------------------------
  # correlations deviation to distance
  #-----------------------------------------
  if(cordo){
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
  if(taxdo){
  {
    library(cluster)
    cats <- factor(taxInfo$Species)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_spec_obs_zlog <- cbind(as.vector(unique(taxInfo$Species)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_spec_pred_zlog <- cbind(as.vector(unique(taxInfo$Species)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Genus)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_gen_obs_zlog<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_gen_pred_zlog<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Family)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_fam_obs_zlog<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_fam_pred_zlog<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Clades)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_clad_obs_zlog<- cbind(as.vector(unique(taxInfo$Clades)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_clad_pred_zlog<- cbind(as.vector(unique(taxInfo$Clades)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PlantGrowthForm)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_obs_zlog <- cbind(as.vector(unique(taxInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_pred_zlog<- cbind(as.vector(unique(taxInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PFT)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_obs_zlog<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_pred_zlog<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
  }
  {
    library(cluster)
    cats <- factor(taxInfo$Species)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs,method = "canberra"))
    si <- summary(si2)
    Sil_spec_obs<- cbind(as.vector(unique(taxInfo$Species)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred,method = "canberra"))
    si <- summary(si2)
    Sil_spec_pred <- cbind(as.vector(unique(taxInfo$Species)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Genus)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs,method = "canberra"))
    si <- summary(si2)
    Sil_gen_obs<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred,method = "canberra"))
    si <- summary(si2)
    Sil_gen_pred<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Family)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs,method = "canberra"))
    si <- summary(si2)
    Sil_fam_obs<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred,method = "canberra"))
    si <- summary(si2)
    Sil_fam_pred<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Clades)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs,method = "canberra"))
    si <- summary(si2)
    Sil_clad_obs<- cbind(as.vector(unique(taxInfo$Clades)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred,method = "canberra"))
    si <- summary(si2)
    Sil_clad_pred<- cbind(as.vector(unique(taxInfo$Clades)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PlantGrowthForm)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_obs <- cbind(as.vector(unique(taxInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_pred<- cbind(as.vector(unique(taxInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PFT)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_obs<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_pred<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
  }
  }
  
  if(distdo){
  t=1
  for(t in 1:ncol(traitInfo_obs_zlog)){
    
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
        return(sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)*100)
      }else{
        return(NA)
      }
    }
    
    
    #----------------------------------------------
        #Individuals within species
    i=2
    mean_pred_zlog <- list()
    mean_obs_zlog <- list()
    mean_pred <- list()
    mean_obs <- list()
    
    dist_pred_zlog <- list()
    dist_obs_zlog <- list()
    dist_pred <- list()
    dist_obs <- list()
    
    CV_pred_zlog <- list()
    CV_obs_zlog <- list()
    CV_pred <- list()
    CV_obs <- list()
    
    mean_pred_zlog[[i]] <- aggregate(traitInfo_pred_zlog,list(taxInfo[,i]),fun_cluster_mean)
    mean_pred[[i]] <- aggregate(traitInfo_pred_zlog,list(taxInfo[,i]),fun_cluster_mean)
    mean_obs_zlog[[i]] <- aggregate(traitInfo_obs_zlog,list(taxInfo[,i]),fun_cluster_mean)
    mean_obs[[i]] <- aggregate(traitInfo_obs,list(taxInfo[,i]),fun_cluster_mean)
    
    CV_obs_zlog[[i]] <- aggregate(traitInfo_obs_zlog,list(taxInfo[,i]),fun_cluster_CV)
    CV_pred[[i]] <- aggregate(traitInfo_pred_zlog,list(taxInfo[,i]),fun_cluster_CV)
    CV_obs_zlog[[i]] <- aggregate(traitInfo_obs_zlog,list(taxInfo[,i]),fun_cluster_CV)
    CV_obs[[i]] <- aggregate(traitInfo_obs,list(taxInfo[,i]),fun_cluster_CV)
    #----------------------------------------------------------------------
    # DO THIS HERE
    t=1
    Sil_dat <- list(Sil_spec_obs,Sil_gen_obs,Sil_fam_obs,Sil_clad_obs,Sil_GF_obs,Sil_PFT_obs)
    colnms_Sil_pred=c("ObsID","Sil_spec_pred","Sil_gen_pred","Sil_fam_pred","Sil_clad_pred","Sil_GF_pred","Sil_PFT_pred")
    
    mylist = aggregate(traitInfo_obs_zlog[,t],list(taxInfo[,i]),fun_cluster_dist)[,2]
    names(mylist) <- unique(taxInfo[,i])
    
    cl=468
    sum(taxInfo[,i]==unique(taxInfo[,i])[cl])
    
    for(cl in 1:length(unique(taxInfo$Species))){
      ix = taxInfo[,i] == unique(taxInfo[,i])[cl]
      print("---------------------------")
      print(cl)
      print(sum(ix))
      print(length(unlist(strsplit(as.character(mylist[[cl]]),split = ","))))
      print(names(mylist)[cl])
      print(unique(taxInfo[,i])[cl])
    }
    
    traitInfo_obs_zlog[taxInfo[,i]==unique(taxInfo[,i])[cl],t]
    cl=1
    length(unlist(strsplit(as.character(mylist[[cl]]),split = ",")))
    vec=NA
    for(cl in 1:length(mylist)){
      vec <- c(vec,unlist(strsplit(as.character(mylist[[cl]]),split = ",")))
      print(cl)
      print("----")
    }
    
    length(vec)
    nrow(taxInfo)
    length(unique(taxInfo$Species))
    dat_now <- cbind(taxInfo[,i],vec[-1])
    
    
    
    #-------------------------------------------------------------------------------------------------------------------
      res$dist_gen_obs_zlog[ix_trait][ix] <- aggregate(traitInfo_obs_zlog,list(1:nrow(taxInfo),Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_gen_pred[ix_trait][ix] <- aggregate(traitInfo_pred,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_gen_pred_zlog[ix_trait][ix] <- aggregate(traitInfo_pred_zlog,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_gen_obs[ix_trait][ix] <- aggregate(traitInfo_obs,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      
      res$Sil_gen_obs[ix_trait][ix] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_gen_pred[ix_trait][ix] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_gen_obs_zlog[ix_trait][ix] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_spec_pred_zlog[ix_trait][ix] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      
      #-------------------------------------------------------------------------------------------------------------------
      res$dist_fam_obs_zlog[ix_trait][ix] <- aggregate(traitInfo_obs_zlog,list(1:nrow(taxInfo),Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_fam_pred[ix_trait][ix] <- aggregate(traitInfo_pred,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_fam_pred_zlog[ix_trait][ix] <- aggregate(traitInfo_pred_zlog,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_fam_obs[ix_trait][ix] <- aggregate(traitInfo_obs,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      
      res$Sil_fam_obs[ix_trait[i]] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_fam_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_fam_obs_zlog[ix_trait[i]] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_fam_pred_zlog[ix_trait[i]] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      
      #-------------------------------------------------------------------------------------------------------------------
      res$dist_clad_obs_zlog[ix_trait[i]] <- aggregate(traitInfo_obs_zlog,list(1:nrow(taxInfo),Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_clad_pred[ix_trait[i]] <- aggregate(traitInfo_pred,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_clad_pred_zlog[ix_trait[i]] <- aggregate(traitInfo_pred_zlog,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_clad_obs[ix_trait[i]] <- aggregate(traitInfo_obs,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)

      res$Sil_clad_obs[ix_trait[i]] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_clad_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_clad_obs_zlog[ix_trait[i]] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_clad_pred_zlog[ix_trait[i]] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      
      #-------------------------------------------------------------------------------------------------------------------
      res$dist_GF_obs_zlog[ix_trait[i]] <- aggregate(traitInfo_obs_zlog,list(1:nrow(taxInfo),Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_GF_pred[ix_trait[i]] <- aggregate(traitInfo_pred,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_GF_pred_zlog[ix_trait[i]] <- aggregate(traitInfo_pred_zlog,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_GF_obs[ix_trait[i]] <- aggregate(traitInfo_obs,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      
      res$Sil_GF_obs[ix_trait[i]] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_GF_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_GF_obs_zlog[ix_trait[i]] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_GF_pred_zlog[ix_trait[i]] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      
      #-------------------------------------------------------------------------------------------------------------------
      res$dist_PFT_obs_zlog[ix_trait[i]] <- aggregate(traitInfo_obs_zlog,list(1:nrow(taxInfo),Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_PFT_pred[ix_trait[i]] <- aggregate(traitInfo_pred,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_PFT_pred_zlog[ix_trait[i]] <- aggregate(traitInfo_pred_zlog,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      res$dist_PFT_obs[ix_trait[i]] <- aggregate(traitInfo_obs,list(ObsID=taxInfo[,1],Group=taxInfo[,i]),fun_cluster_dist)
      
      res$Sil_PFT_obs[ix_trait[i]] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_PFT_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_PFT_obs_zlog[ix_trait[i]] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_PFT_pred_zlog[ix_trait[i]] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      
    }
    
  }#t
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
  
  write.csv(res,file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_ind.csv")))
  
  if(!file.exists(file.path(origin,"_2021","data","analyses","CV"))){
    dir.create(file.path(origin,"_2021","data","analyses","CV"))}
  if(!file.exists(file.path(origin,"_2021","data","analyses","CV",RepNum))){
    dir.create(file.path(origin,"_2021","data","analyses","CV",RepNum))}
  if(!file.exists(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice))){
    dir.create(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice))}
  if(!file.exists(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD))){
    dir.create(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD))}
  if(!file.exists(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent))){
    dir.create(file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent))}
  
  save(CV_pred_zlog,file=file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent,paste0("CV_pred_zlog.RData")))
  save(CV_pred,file=file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent,paste0("CV_pred.RData")))
  save(CV_obs_zlog,file=file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent,paste0("CV_obs_zlog.RData")))
  save(CV_obs,file=file.path(origin,"_2021","data","analyses","CV",RepNum,t_choice,ObsOrTD,Percent,paste0("CV_obs.RData")))
  
  }#Percent
}

  #TD <- res$Sil_spec_obs_zlog
  #ENV <- res$Sil_spec_obs_zlog
  head(res$Sil_spec_obs_zlog)
  head(resENV$Sil_spec_obs_zlog)
  
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs"
  Percent=80
  resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  head(resENV$Sil_spec_obs_zlog)
  
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  Percent=80
  res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  head(res$Sil_spec_obs_zlog)

