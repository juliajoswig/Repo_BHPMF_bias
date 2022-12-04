
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
  

  RepNum=1
  t_choice="data_2"
  ObsOrTD="Obs_obs_TD"
  ObsOrTD="Obs_obs"
  Percent=80
  
  repnums=3
  RepNum=2
  repnums=3
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  ObsOrTD="Obs_obs"
  TDnos=c("Obs_obs_TD","Obs_obs")
  TDno=2
  ObsOrTD <- TDnos[TDno]

        

Percent = 80
gappercents=c(1,5,10,20,30,40,50,60,70,80)
gappercents=c(80)
resdo=TRUE
cordo=TRUE
taxdo=TRUE
distdo=TRUE

print(paste(RepNum,ObsOrTD,t_choice,Percent))

for(Percent in gappercents){
  
  print("------ % gaps now ----------")
  print(Percent)
  
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
  {
    list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
    list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
    list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
    # taxonomy 
    tmp <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                            ObsOrTD,"data","traitInfoTD_obs.csv"),header=TRUE))[,-1]
    taxInfo <- read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                  ObsOrTD,"data","taxInfo.csv"),header=FALSE)
    colnames(taxInfo) <- c("ObservationID","Species","Genus","Family","Clades")
    taxInfo <- taxInfo[which(taxInfo$ObservationID%in%tmp$ObservationID),]
    dim(taxInfo)
    funInfo <- read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                  ObsOrTD,"data","funInfo.csv"),header=TRUE)[,-1]
    funInfo <- funInfo[which(funInfo$ObservationID%in%tmp$ObservationID),]
    taxInfo <- cbind(taxInfo,funInfo[,-1])
    write.csv(taxInfo,file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                ObsOrTD,"data","taxInfoTD.csv"))
    head(taxInfo)
    rm("tmp")
    
    #-----------------------------
    # observed 
    #-----------------------------
    # completely observed  untransformed
    traitInfo_obs <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                  "Obs_obs_TD","data","traitInfoTD_obs.csv")))[,-c(1,2)]
    head(traitInfo_obs)
    # sparse untransformed
    traitInfo_obs_sparse <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                          ObsOrTD,"data","traitInfoTD_obs.csv")))[,-c(1,2)]
    head(traitInfo_obs_sparse)
    dim(traitInfo_obs_sparse)
    # sparse zlog-transformed
    traitInfo_obs_zlog_sparse <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                                ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")))[,-c(1,2)]
    # completely observed zlog-transformed
    traitInfo_obs_zlog <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,"p_0",
                                                              ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")))[,-c(1,2)]
    list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
    #-----------------------------
    # imputed / predicted 
    #-----------------------------
    # complete predicted untransformed
    traitInfo_pred <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                   ObsOrTD,"data","traitInfoTD_pred.csv")))[,-c(1,2)]
    dim(traitInfo_pred)
    # sparse predicted untransformed
    traitInfo_pred_sparse <- traitInfo_pred
    traitInfo_pred_sparse[is.na(traitInfo_obs_sparse)] <- NA
    # complete predicted zlog-transformed
    traitInfo_pred_zlog <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),
                                                          ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
    head(traitInfo_pred_zlog)
    # sparse predicted zlog-transformed
    traitInfo_pred_zlog_sparse <- traitInfo_pred_zlog
    traitInfo_pred_zlog_sparse[is.na(traitInfo_obs_sparse)] <- NA
    par(mfrow=c(1,4))
    plot(traitInfo_pred[,1],traitInfo_obs[,1])
    abline(0,1)
    plot(traitInfo_pred_zlog[,1],traitInfo_obs_zlog[,1])
    abline(0,1)
    plot(traitInfo_pred_zlog_sparse[,1],traitInfo_obs_zlog_sparse[,1])
    abline(0,1)
    plot(traitInfo_pred_sparse[,1],traitInfo_obs_sparse[,1])
    abline(0,1)
    
    # we now have loaded:
    # traitInfo_obs
    # traitInfo_obs_sparse
    # traitInfo_obs_zlog
    # traitInfo_obs_zlog_sparse
    # traitInfo_pred
    # traitInfo_pred_sparse
    # traitInfo_pred_zlog
    # traitInfo_pred_zlog_sparse
  }

  
    if(resdo|!file.exists(file.path(originData,"analyses","Point_wise",
                                    RepNum,t_choice,ObsOrTD,Percent,paste0("res_2022_12_2.csv")))){
    nb_rows = nrow(traitInfo_obs) + 
      ((ncol(traitInfo_obs))*nrow(traitInfo_obs)) 
    nms <- c("error","Gap","missingness",
             "trait","ID","Species","Genus","Family","Clade","GF","PFT",
             "value_obs","value_pred","value_obs_zlog","value_pred_zlog",
             "value_obs_sparse","value_pred_sparse","value_obs_zlog_sparse","value_pred_zlog_sparse",
             
             "nrep","DataSet",
             
             "nb_spec","nb_gen","nb_fam","nb_clad","nb_GF","nb_PFT",

             "mean_spec_obs_zlog","mean_gen_obs_zlog","mean_fam_obs_zlog","mean_clad_obs_zlog","mean_GF_obs_zlog","mean_PFT_obs_zlog",#input mean
             "mean_spec_obs_zlog_sparse","mean_gen_obs_zlog_sparse","mean_fam_obs_zlog_sparse","mean_clad_obs_zlog_sparse","mean_GF_obs_zlog_sparse","mean_PFT_obs_zlog_sparse",#input mean
             "mean_spec_obs","mean_gen_obs","mean_fam_obs","mean_clad_obs","mean_GF_obs","mean_PFT_obs",#input mean
             "mean_spec_obs_sparse","mean_gen_obs_sparse","mean_fam_obs_sparse","mean_clad_obs_sparse","mean_GF_obs_sparse","mean_PFT_obs_sparse",#input mean
             
             "mean_spec_pred_zlog","mean_gen_pred_zlog","mean_fam_pred_zlog","mean_clad_pred_zlog","mean_GF_pred_zlog","mean_PFT_pred_zlog",#input mean
             "mean_spec_pred_zlog_sparse","mean_gen_pred_zlog_sparse","mean_fam_pred_zlog_sparse","mean_clad_pred_zlog_sparse","mean_GF_pred_zlog_sparse","mean_PFT_pred_zlog_sparse",#input mean
             "mean_spec_pred","mean_gen_pred","mean_fam_pred","mean_clad_pred","mean_GF_pred","mean_PFT_pred",#input mean
             "mean_spec_pred_sparse","mean_gen_pred_sparse","mean_fam_pred_sparse","mean_clad_pred_sparse","mean_GF_pred_sparse","mean_PFT_pred_sparse",#input mean
             
             "dist_spec_obs_zlog","dist_gen_obs_zlog","dist_fam_obs_zlog","dist_clad_obs_zlog","dist_GF_obs_zlog","dist_PFT_obs_zlog",#input dist
             "dist_spec_obs_zlog_sparse","dist_gen_obs_zlog_sparse","dist_fam_obs_zlog_sparse","dist_clad_obs_zlog_sparse","dist_GF_obs_zlog_sparse","dist_PFT_obs_zlog_sparse",#input dist
             "dist_spec_obs","dist_gen_obs","dist_fam_obs","dist_clad_obs","dist_GF_obs","dist_PFT_obs",#input dist
             "dist_spec_obs_sparse","dist_gen_obs_sparse","dist_fam_obs_sparse","dist_clad_obs_sparse","dist_GF_obs_sparse","dist_PFT_obs_sparse",#input dist
             
             "dist_spec_pred_zlog","dist_gen_pred_zlog","dist_fam_pred_zlog","dist_clad_pred_zlog","dist_GF_pred_zlog","dist_PFT_pred_zlog",#input dist
             "dist_spec_pred_zlog_sparse","dist_gen_pred_zlog_sparse","dist_fam_pred_zlog_sparse","dist_clad_pred_zlog_sparse","dist_GF_pred_zlog_sparse","dist_PFT_pred_zlog_sparse",#input dist
             "dist_spec_pred","dist_gen_pred","dist_fam_pred","dist_clad_pred","dist_GF_pred","dist_PFT_pred",#input dist
             "dist_spec_pred_sparse","dist_gen_pred_sparse","dist_fam_pred_sparse","dist_clad_pred_sparse","dist_GF_pred_sparse","dist_PFT_pred_sparse",#input dist
             
             "CV_spec_obs","CV_gen_obs","CV_fam_obs","CV_clad_obs","CV_GF_obs","CV_PFT_obs",
             "CV_spec_pred","CV_gen_pred","CV_fam_pred","CV_clad_pred","CV_GF_pred","CV_PFT_pred"
             )
    res <- as.data.frame(matrix(NA,ncol=length(nms),nrow=nb_rows))
    colnames(res) <- nms
    
  }else{
    res <- read.csv(file=file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2022_12_2.csv")))
  }
  
  #-------------------------------------------------------------------
  # distance to taxon mean
  #-------------------------------------------------------------------
  
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
    res$GF[ix_trait] <- as.vector(funInfo$PlantGrowthForm)
    res$PFT[ix_trait] <- as.vector(funInfo$PFT)
    res$error[ix_trait] <- traitInfo_pred_zlog[,t] - traitInfo_obs_zlog[,t]
    res$Gap[ix_trait] <- is.na(traitInfo_obs_sparse[,t])
    
    res$value_obs_zlog[ix_trait] <- traitInfo_obs_zlog[,t]
    res$value_pred_zlog[ix_trait] <- traitInfo_pred_zlog[,t]
    res$value_obs[ix_trait] <- traitInfo_obs[,t]
    res$value_pred[ix_trait] <- traitInfo_pred[,t]
    res$value_obs_zlog_sparse[ix_trait] <- traitInfo_obs_zlog_sparse[,t]
    res$value_pred_zlog_sparse[ix_trait] <- traitInfo_pred_zlog_sparse[,t]
    res$value_obs_sparse[ix_trait] <- traitInfo_obs_sparse[,t]
    res$value_pred_sparse[ix_trait] <- traitInfo_pred_sparse[,t]
    
    res$missingness[ix_trait] <- Percent
    res$DataSet[ix_trait] <- t_choice
    res$nrep[ix_trait] <- as.vector(RepNum)
  }
  
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
    print("------------missingness---------------")
    print(unique(res$missingness))

    i=80
    for(i in 1: nrow(traitInfo_obs_zlog)){
      print(paste(t,"-",i))    
      {
      #----------------------------------------------
      #Individuals within SPECIES
      #----------------------------------------------
      {
        ix1 = taxInfo$Species==as.vector(taxInfo$Species[i])
        
        res$mean_spec_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
       
        res$mean_spec_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
        res$mean_spec_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
        res$nb_spec[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$CV_spec_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_spec_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      }
       #----------------------------------------------
      #Species within GENERA
      #----------------------------------------------
      {
        ix1 = taxInfo$Genus==as.vector(taxInfo$Genus[i])
        
        res$mean_gen_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        res$mean_gen_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
        res$mean_gen_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        
        ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_gen[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$CV_gen_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_gen_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
      }
      #----------------------------------------------
      #Genera within FAMILIES
      #----------------------------------------------
      {
        ix1 = taxInfo$Family==as.vector(taxInfo$Family[i])
        res$mean_fam_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        res$mean_fam_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
        res$mean_fam_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        
        ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
        res$nb_fam[ix_trait[i]] <- sum(ix1)
        if(sum(ix1)>1){
          res$CV_fam_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
          res$CV_fam_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
        }
      }
      #----------------------------------------------
      #Families within CLADES
      #----------------------------------------------
      {
        ix1 = taxInfo$Clades==as.vector(taxInfo$Clades[i])
        res$mean_clad_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        res$mean_clad_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
        res$mean_clad_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        
        ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
        res$nb_clad[ix_trait[i]] <- sum(ix1)
        if(sum(ix1)>1){
          res$CV_clad_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
          res$CV_clad_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
        }
      }
      #---------------------------------------------------------------------
      # PlantGrowthForm
      #---------------------------------------------------------------------
      {
      ix1 = taxInfo$PlantGrowthForm==as.vector(funInfo$PlantGrowthForm[i])
      ix1[is.na(ix1)] <- FALSE

      res$mean_GF_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
      
      res$mean_GF_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
      
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_GF[ix_trait[i]] <- sum(ix1)
      if(sum(ix1)>1){
        res$CV_GF_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
        res$CV_GF_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
      }
  }
      #---------------------------------------------------------------------
      # PFT
      #---------------------------------------------------------------------
      {
        ix1 = taxInfo$PFT==as.vector(funInfo$PFT[i])
        ix1[is.na(ix1)] <- FALSE
        
        res$mean_PFT_obs_zlog[ix_trait[i]] <- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_obs_zlog_sparse[ix_trait[i]] <- mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_obs_sparse[ix_trait[i]] <- mean(traitInfo_obs_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        res$mean_PFT_pred_zlog[ix_trait[i]] <- mean(traitInfo_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_pred_zlog_sparse[ix_trait[i]] <- mean(traitInfo_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_pred[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
        res$mean_PFT_pred_sparse[ix_trait[i]] <- mean(traitInfo_pred_sparse[ix1,t],na.rm = TRUE) #input species mean
        
        ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
        res$nb_PFT[ix_trait[i]] <- sum(ix1)
        if(sum(ix1)>1){
          res$CV_PFT_obs[ix_trait[i]] <- sd(traitInfo_obs[ix,t])/mean(traitInfo_obs[ix,t])
          res$CV_PFT_pred[ix_trait[i]] <- sd(traitInfo_pred[ix,t])/mean(traitInfo_pred[ix,t])
        }
      }
        #----------------------------------------------------------------------
      }
  
    }#i
  }#t
  
  {
    # Value distance to species mean
    res$dist_spec_obs                    <- res$value_obs -             res$mean_spec_obs
    res$dist_spec_obs_sparse             <- res$value_obs_sparse -      res$mean_spec_obs_sparse
    res$dist_spec_obs_zlog               <- res$value_obs_zlog -        res$mean_spec_obs_zlog
    res$dist_spec_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_spec_obs_zlog_sparse
    res$dist_spec_pred                   <- res$value_pred -            res$mean_spec_pred
    res$dist_spec_pred_sparse            <- res$value_pred_sparse -     res$mean_spec_pred_sparse
    res$dist_spec_pred_zlog              <- res$value_pred_zlog -       res$mean_spec_pred_zlog
    res$dist_spec_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_spec_pred_zlog_sparse
    
    # Value distance to genus mean
    res$dist_gen_obs                    <- res$value_obs -             res$mean_gen_obs
    res$dist_gen_obs_sparse             <- res$value_obs_sparse -      res$mean_gen_obs_sparse
    res$dist_gen_obs_zlog               <- res$value_obs_zlog -        res$mean_gen_obs_zlog
    res$dist_gen_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_gen_obs_zlog_sparse
    res$dist_gen_pred                   <- res$value_pred -            res$mean_gen_pred
    res$dist_gen_pred_sparse            <- res$value_pred_sparse -     res$mean_gen_pred_sparse
    res$dist_gen_pred_zlog              <- res$value_pred_zlog -       res$mean_gen_pred_zlog
    res$dist_gen_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_gen_pred_zlog_sparse
    
    # Value distance to family mean
    res$dist_fam_obs                    <- res$value_obs -             res$mean_fam_obs
    res$dist_fam_obs_sparse             <- res$value_obs_sparse -      res$mean_fam_obs_sparse
    res$dist_fam_obs_zlog               <- res$value_obs_zlog -        res$mean_fam_obs_zlog
    res$dist_fam_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_fam_obs_zlog_sparse
    res$dist_fam_pred                   <- res$value_pred -            res$mean_fam_pred
    res$dist_fam_pred_sparse            <- res$value_pred_sparse -     res$mean_fam_pred_sparse
    res$dist_fam_pred_zlog              <- res$value_pred_zlog -       res$mean_fam_pred_zlog
    res$dist_fam_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_fam_pred_zlog_sparse
    
    # Value distance to clades mean
    res$dist_clad_obs                    <- res$value_obs -             res$mean_clad_obs
    res$dist_clad_obs_sparse             <- res$value_obs_sparse -      res$mean_clad_obs_sparse
    res$dist_clad_obs_zlog               <- res$value_obs_zlog -        res$mean_clad_obs_zlog
    res$dist_clad_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_clad_obs_zlog_sparse
    res$dist_clad_pred                   <- res$value_pred -            res$mean_clad_pred
    res$dist_clad_pred_sparse            <- res$value_pred_sparse -     res$mean_clad_pred_sparse
    res$dist_clad_pred_zlog              <- res$value_pred_zlog -       res$mean_clad_pred_zlog
    res$dist_clad_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_clad_pred_zlog_sparse
    
    # Value distance to GF mean
    res$dist_GF_obs                    <- res$value_obs -             res$mean_GF_obs
    res$dist_GF_obs_sparse             <- res$value_obs_sparse -      res$mean_GF_obs_sparse
    res$dist_GF_obs_zlog               <- res$value_obs_zlog -        res$mean_GF_obs_zlog
    res$dist_GF_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_GF_obs_zlog_sparse
    res$dist_GF_pred                   <- res$value_pred -            res$mean_GF_pred
    res$dist_GF_pred_sparse            <- res$value_pred_sparse -     res$mean_GF_pred_sparse
    res$dist_GF_pred_zlog              <- res$value_pred_zlog -       res$mean_GF_pred_zlog
    res$dist_GF_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_GF_pred_zlog_sparse
    
    # Value distance to PFT mean
    res$dist_PFT_obs                    <- res$value_obs -             res$mean_PFT_obs
    res$dist_PFT_obs_sparse             <- res$value_obs_sparse -      res$mean_PFT_obs_sparse
    res$dist_PFT_obs_zlog               <- res$value_obs_zlog -        res$mean_PFT_obs_zlog
    res$dist_PFT_obs_zlog_sparse        <- res$value_obs_zlog_sparse - res$mean_PFT_obs_zlog_sparse
    res$dist_PFT_pred                   <- res$value_pred -            res$mean_PFT_pred
    res$dist_PFT_pred_sparse            <- res$value_pred_sparse -     res$mean_PFT_pred_sparse
    res$dist_PFT_pred_zlog              <- res$value_pred_zlog -       res$mean_PFT_pred_zlog
    res$dist_PFT_pred_zlog_sparse       <- res$value_pred_zlog_sparse- res$mean_PFT_pred_zlog_sparse
  }
  }

  
  
  if(!file.exists(file.path(originData,"analyses"))){
    dir.create(file.path(originData,"analyses"))}
  if(!file.exists(file.path(originData,"analyses","Point_wise"))){
    dir.create(file.path(originData,"analyses","Point_wise"))}
  if(!file.exists(file.path(originData,"analyses","Point_wise",RepNum))){
    dir.create(file.path(originData,"analyses","Point_wise",RepNum))}
  if(!file.exists(file.path(originData,"analyses","Point_wise",RepNum,t_choice))){
    dir.create(file.path(originData,"analyses","Point_wise",RepNum,t_choice))}
  if(!file.exists(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD))){
    dir.create(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD))}
  if(!file.exists(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent))){
    dir.create(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent))}
  write.csv(res,file=file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2022_12_2.csv")))
  
  }#Percent



res <- read.csv(file=file.path(originData,"analyses","Point_wise",RepNum,
                               t_choice,ObsOrTD,Percent,paste0("res_2022_12_2.csv")))
head(res)





























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
}

