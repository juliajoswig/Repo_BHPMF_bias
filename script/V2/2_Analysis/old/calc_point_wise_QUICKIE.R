
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
  


RepNum=1
repnums=3
t_choice="data"
ObsOrTD="Obs_obs_TD"
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
TDno=1
ObsOrTD <- TDnos[TDno]

print(paste(RepNum,ObsOrTD,t_choice,Percent))
        


gappercents=c(1,5,10,20,30,40,50,60,70,80)
gappercents=c(80,1)

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
    tmp <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                            ObsOrTD,"data","traitInfoTD_obs.csv"),header=TRUE))[,-1]
    taxInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                  ObsOrTD,"data","taxInfo.csv"),header=FALSE)
    colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clades")
    taxInfo <- taxInfo[which(taxInfo$ObservationID%in%tmp$ObservationID),]
    dim(taxInfo)
    funInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                  ObsOrTD,"data","funInfo.csv"),header=TRUE)[,-1]
    funInfo <- funInfo[which(funInfo$ObservationID%in%tmp$ObservationID),]
    taxInfo <- cbind(taxInfo,funInfo[,-1])
    write.csv(taxInfo,file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                ObsOrTD,"data","taxInfoTD.csv"))
    head(taxInfo)
    rm("tmp")
    #list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))

    # observed 
    # sparse untransformed
    traitInfo_obs_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                          ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
    head(traitInfo_obs_sparse)
    # completely observed  untransformed
    traitInfo_obs <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                  ObsOrTD,"data","traitInfo.csv")))[,-c(1,2)]
    head(traitInfo_obs)
    # sparse zlog-transformed
    traitInfo_obs_zlog_sparse <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_obs_zlog.csv")))[,-c(1,2)]
    head(traitInfo_obs_zlog_sparse)
    # completely observed zlog-transformed
    traitInfo_log <- log(traitInfo_obs)
    load(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","zlog_transform.RData"))
    traitInfo_obs_zlog <- (traitInfo_log - zlog_trans$meanM) /zlog_trans$sdM
    head(traitInfo_obs_zlog)
    write.csv(cbind(traitInfo_obs_zlog,taxInfo[,1]),file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_obsTot_zlog.csv"))
    traitInfo_obs_zlog[traitInfo_obs_sparse] <- NA
    
    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
    # predicted 
    # complete untransformed
    traitInfo_pred <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                         ObsOrTD,"data","traitInfoTD_pred.csv")))[,-c(1,2)]
    head(traitInfo_pred)
    # sparse untransformed
    traitInfo_pred_sparse <- traitInfo_pred
    traitInfo_pred_sparse[is.na(traitInfo_obs_sparse)] <- NA
    # complete zlog-transformed
    traitInfo_pred_zlog <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_pred_zlog.csv")))[,-c(1,2)]
    head(traitInfo_pred_zlog)
    # sparse zlog-transformed
    traitInfo_pred_zlog_sparse <- traitInfo_pred_zlog
    traitInfo_pred_zlog_sparse[is.na(traitInfo_obs_sparse)] <- NA
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
            colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_pred")
            res <- add_col_to_res(new.col.names = colnm4,input = res)
            #observed
            dat_now=as.data.frame(traitInfo_pred[,c(t1,t)])
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now) # lm produced from NONsparse data!
            t_lm_obs <- (traitInfo_obs[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
          }
          if(t1!=t){
            colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_pred_zlog")
            res <- add_col_to_res(new.col.names = colnm4,input = res)
            #observed
            dat_now=as.data.frame(traitInfo_pred_zlog[,c(t1,t)])
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now) # lm produced from NONsparse data!
            t_lm_obs <- (traitInfo_obs[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
          }
          if(t1!=t){
            colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_obs_sparse")
            res <- add_col_to_res(new.col.names = colnm4,input = res)
            #observed
            dat_now=as.data.frame(traitInfo_obs_sparse[,c(t1,t)])
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now,na.action = "na.omit") # lm produced from NONsparse data!
            t_lm_obs <- (traitInfo_obs[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
          }
          
          if(t1!=t){
            colnm4=paste0(colnames(traitInfo_obs)[t],"_from_",colnames(traitInfo_obs)[t1],"_lm_obs_zlog_sparse")
            res <- add_col_to_res(new.col.names = colnm4,input = res)
            #observed
            dat_now=as.data.frame(traitInfo_obs_zlog_sparse[,c(t1,t)])
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now,na.action = "na.omit") # lm produced from NONsparse data!
            t_lm_obs <- (traitInfo_obs[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
          }
        }
    }
    
    
  }  
  
  #-----------------------------------------
  # taxonomy
  #-----------------------------------------
  {
    library(cluster)
    cats <- factor(taxInfo$Species)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_spec_obs_zlog<- cbind(as.vector(unique(taxInfo$Species)),as.vector(si$clus.avg.widths))
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
    Sil_GF_obs_zlog <- cbind(as.vector(unique(funInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_pred_zlog<- cbind(as.vector(unique(funInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PFT)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_obs_zlog<- cbind(as.vector(unique(funInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_pred_zlog<- cbind(as.vector(unique(funInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
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
    Sil_GF_obs <- cbind(as.vector(unique(funInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_pred<- cbind(as.vector(unique(funInfo$PlantGrowthForm[ixG])),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PFT)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_obs<- cbind(as.vector(unique(funInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_pred<- cbind(as.vector(unique(funInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
  }
  
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
    res$GF[ix_trait] <- as.vector(funInfo$PlantGrowthForm)
    res$PFT[ix_trait] <- as.vector(funInfo$PFT)
    res$error[ix_trait] <- traitInfo_pred_zlog[,t] - traitInfo_obs_zlog[,t]
    res$Gap[ix_trait] <- is.na(traitInfo_obs_sparse[,t])
    res$value_obs_zlog[ix_trait] <- traitInfo_obs_zlog[,t]
    res$value_pred_zlog[ix_trait] <- traitInfo_pred_zlog[,t]
    res$value_obs[ix_trait] <- traitInfo_obs[,t]
    res$value_pred[ix_trait] <- traitInfo_pred[,t]
    res$missingness[ix_trait] <- Percent
    res$DataSet[ix_trait] <- t_choice
    res$nrep[ix_trait] <- as.vector(RepNum)
    
    print("------------missingness---------------")
    print(unique(res$missingness))

    i=80
    for(i in 1: nrow(traitInfo_obs_zlog)){
      print(paste(t,"-",i))    
      {
      #----------------------------------------------
      #Individuals within species
      #----------------------------------------------
      {
        ix1 = taxInfo$Species==as.vector(taxInfo$Species[i])
      res$mean_spec_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_spec_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_spec_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_spec[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_spec[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$dist_spec_obs_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_spec_pred_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred_zlog[i,t],
                 traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_spec_obs[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs[i,t],
                 traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
        res$dist_spec_pred[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred[i,t],
                 traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
        res$CV_spec_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_spec_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      res$Sil_spec_obs[ix_trait[i]] <- as.numeric(Sil_spec_obs[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_spec_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_pred[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_spec_obs_zlog[ix_trait[i]] <- as.numeric(Sil_spec_obs_zlog[which(Sil_spec_obs_zlog[,1]==res$Species[ix_trait[i]]),2])
      res$Sil_spec_pred_zlog[ix_trait[i]] <- as.numeric(Sil_spec_pred_zlog[which(Sil_spec_pred_zlog[,1]==res$Species[ix_trait[i]]),2])
      }
       #----------------------------------------------
      #Species within genera
      #----------------------------------------------
      {
        ix1 = taxInfo$Genus==as.vector(taxInfo$Genus[i])
      res$mean_gen_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_gen_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_gen[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      res$mean_gen_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_gen[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$dist_gen_obs_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_gen_pred_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred_zlog[i,t],
                 traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_gen_obs[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs[i,t],
                 traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
        res$dist_gen_pred[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred[i,t],
                 traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
        res$CV_gen_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_gen_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      res$Sil_gen_obs[ix_trait[i]] <- as.numeric(Sil_gen_obs[which(Sil_gen_obs[,1]==res$Genus[ix_trait[i]]),2])
      res$Sil_gen_pred[ix_trait[i]] <- as.numeric(Sil_gen_pred[which(Sil_gen_pred[,1]==res$Genus[ix_trait[i]]),2])
      res$Sil_gen_obs_zlog[ix_trait[i]] <- as.numeric(Sil_gen_obs_zlog[which(Sil_gen_obs_zlog[,1]==res$Genus[ix_trait[i]]),2])
      res$Sil_gen_pred_zlog[ix_trait[i]] <- as.numeric(Sil_gen_pred_zlog[which(Sil_gen_pred_zlog[,1]==res$Genus[ix_trait[i]]),2])
      }
      #----------------------------------------------
      #Genera within families
      #----------------------------------------------
        {
        ix1 = taxInfo$Family==as.vector(taxInfo$Family[i])
      res$mean_fam_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_fam_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_fam_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_fam[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_fam[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$dist_fam_obs_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_fam_pred_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred_zlog[i,t],
                 traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_fam_obs[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs[i,t],
                 traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
        res$dist_fam_pred[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred[i,t],
                 traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
        res$CV_fam_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_fam_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      res$Sil_fam_obs[ix_trait[i]] <- as.numeric(Sil_fam_obs[which(Sil_fam_obs[,1]==res$Family[ix_trait[i]]),2])
      res$Sil_fam_pred[ix_trait[i]] <- as.numeric(Sil_fam_pred[which(Sil_fam_pred[,1]==res$Family[ix_trait[i]]),2])
      res$Sil_fam_obs_zlog[ix_trait[i]] <- as.numeric(Sil_fam_obs_zlog[which(Sil_fam_obs_zlog[,1]==res$Family[ix_trait[i]]),2])
      res$Sil_fam_pred_zlog[ix_trait[i]] <- as.numeric(Sil_fam_pred_zlog[which(Sil_fam_pred_zlog[,1]==res$Family[ix_trait[i]]),2])
        }
        #----------------------------------------------
      #Families within Clades
      #----------------------------------------------
      {
        ix1 = taxInfo$Clades==as.vector(taxInfo$Clades[i])
      res$mean_clad_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_clad_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_clad_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_clad[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_clad[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$dist_clad_obs_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_clad_pred_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred_zlog[i,t],
                 traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_clad_obs[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs[i,t],
                 traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
        res$dist_clad_pred[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred[i,t],
                 traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
        res$CV_clad_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_clad_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      res$Sil_clad_obs[ix_trait[i]] <- as.numeric(Sil_clad_obs[which(Sil_clad_obs[,1]==res$Clade[ix_trait[i]]),2])
      res$Sil_clad_pred[ix_trait[i]] <- as.numeric(Sil_clad_pred[which(Sil_clad_pred[,1]==res$Clade[ix_trait[i]]),2])
      res$Sil_clad_obs_zlog[ix_trait[i]] <- as.numeric(Sil_clad_obs_zlog[which(Sil_clad_obs_zlog[,1]==res$Clade[ix_trait[i]]),2])
      res$Sil_clad_pred_zlog[ix_trait[i]] <- as.numeric(Sil_clad_pred_zlog[which(Sil_clad_pred_zlog[,1]==res$Clade[ix_trait[i]]),2])
      }
      #---------------------------------------------------------------------
      # PlantGrowthForm
      #---------------------------------------------------------------------
        {
        ix1 = taxInfo$PlantGrowthForm==as.vector(funInfo$PlantGrowthForm[i])
      ix1[is.na(ix1)] <- FALSE
      res$mean_GF_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix1[is.na(ix1)] <- FALSE
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_GF[ix_trait[i]] <- sum(ix1,na.rm = TRUE)
      if(!is.na(as.vector(funInfo$PlantGrowthForm[i]))){
      if(sum(ix1)>1){
        res$dist_GF_obs_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_GF_pred_zlog[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred_zlog[i,t],
                 traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
        res$dist_GF_obs[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs[i,t],
                 traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
        res$dist_GF_pred[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_pred[i,t],
                 traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
        res$CV_GF_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_GF_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
        }
      if(!is.na(res$GF[ix_trait[i]])){
        try(res$Sil_GF_obs[ix_trait[i]] <- as.numeric(Sil_GF_obs[which(Sil_GF_obs[,1]==res$GF[ix_trait[i]]),2]))
        try(res$Sil_GF_pred[ix_trait[i]] <- as.numeric(Sil_GF_pred[which(Sil_GF_pred[,1]==res$GF[ix_trait[i]]),2]))
        try(res$Sil_GF_obs_zlog[ix_trait[i]] <- as.numeric(Sil_GF_obs_zlog[which(Sil_GF_obs_zlog[,1]==res$GF[ix_trait[i]]),2]))
        try(res$Sil_GF_pred_zlog[ix_trait[i]] <- as.numeric(Sil_GF_pred_zlog[which(Sil_GF_pred_zlog[,1]==res$GF[ix_trait[i]]),2]))
      }
        }
      #---------------------------------------------------------------------
      # PFT
      #---------------------------------------------------------------------
      {
        ix1 = taxInfo$PFT==as.vector(funInfo$PFT[i])
      ix1[is.na(ix1)] <- FALSE
      res$mean_PFT_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix1[is.na(ix1)] <- FALSE
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_PFT[ix_trait[i]] <- sum(ix1,na.rm = TRUE)
      if(!is.na(as.vector(funInfo$PFT[i]))){
        if(sum(ix1)>1){
          res$dist_PFT_obs_zlog[ix_trait[i]] <- mean(as.matrix(
            dist(c(traitInfo_obs_zlog[i,t],
                   traitInfo_obs_zlog[ix,t])))[,1],na.rm = TRUE)
          res$dist_PFT_pred_zlog[ix_trait[i]] <- mean(as.matrix(
            dist(c(traitInfo_pred_zlog[i,t],
                   traitInfo_pred_zlog[ix,t])))[,1],na.rm = TRUE)
          res$dist_PFT_obs[ix_trait[i]] <- mean(as.matrix(
            dist(c(traitInfo_obs[i,t],
                   traitInfo_obs[ix,t])))[,1],na.rm = TRUE)
          res$dist_PFT_pred[ix_trait[i]] <- mean(as.matrix(
            dist(c(traitInfo_pred[i,t],
                   traitInfo_pred[ix,t])))[,1],na.rm = TRUE)
          res$CV_PFT_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
          res$CV_PFT_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
        }}
      
      if(!is.na(res$PFT[ix_trait[i]])){
        res$Sil_PFT_obs[ix_trait[i]] <- as.numeric(Sil_PFT_obs[which(Sil_PFT_obs[,1]==res$PFT[ix_trait[i]]),2])
        res$Sil_PFT_pred[ix_trait[i]] <- as.numeric(Sil_PFT_pred[which(Sil_PFT_pred[,1]==res$PFT[ix_trait[i]]),2])
        res$Sil_PFT_obs_zlog[ix_trait[i]] <- as.numeric(Sil_PFT_obs_zlog[which(Sil_PFT_obs_zlog[,1]==res$PFT[ix_trait[i]]),2])
        res$Sil_PFT_pred_zlog[ix_trait[i]] <- as.numeric(Sil_PFT_pred_zlog[which(Sil_PFT_pred_zlog[,1]==res$PFT[ix_trait[i]]),2])
      }
      }
        #----------------------------------------------------------------------
      }
  
    }#i
  }#t


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
  write.csv(res,file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  
  }#Percent

res_save <- res
summary(res_save)

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

list.files(file.path(origin,"_2021","data","analyes"))
write.csv(res,file=file.path(origin,"_2021","data","analyses","Point_wise",paste0("res_2021.csv")))
res<- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",paste0("res_2021.csv")))
summary(res$missingness)


summary(res)





if(sum(colnames(res)%in%colnm)==0){
  newcol=rep(NA,nrow(res))
  res= cbind(res,newcol)
  colnames(res)[ncol(res)] <- colnm
}

