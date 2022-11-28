
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
Version_now="V1"
list.files(file.path(origin,"_2021","script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
  tsubs <- out$tsubs
  TD_choices = out$TD_choices
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
  gappercents=c(1,5,10,20,30,40,50,60)
  
rainfor_nrow=1136
rainfor_ncol=6
GapPercent=60
nb_rows = rainfor_nrow + 
  (rainfor_ncol*rainfor_nrow) + 
  (length(gappercents)*rainfor_nrow) 
res <- as.data.frame(matrix(NA,ncol=113,nrow=nb_rows))
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
#                   "mean_specRM","mean_genRM","mean_famRM","mean_cladRM","mean_GFRM","mean_PFTRM",#input mean
                   "Sil_spec_obs","Sil_gen_obs","Sil_fam_obs","Sil_clad_obs","Sil_GF_obs","Sil_PFT_obs",#input mean
                   "Sil_spec_pred","Sil_gen_pred","Sil_fam_pred","Sil_clad_pred","Sil_GF_pred","Sil_PFT_pred",#input mean
                   "CV_spec_obs","CV_gen_obs","CV_fam_obs","CV_clad_obs","CV_GF_obs","CV_PFT_obs",
                   "CV_spec_pred","CV_gen_pred","CV_fam_pred","CV_clad_pred","CV_GF_pred","CV_PFT_pred",
#                   "dist_specRM","dist_genRM","dist_famRM","dist_cladRM","dist_GFRM","dist_PFTRM",
#                   "nb_specRM","nb_genRM","nb_famRM","nb_cladRM","nb_GFRM","nb_PFTRM",
                   # missing: nb_GFRM_gap
                   "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap","nb_GF_gap","nb_PFT_gap",
 #                  "nb_specRM_gap","nb_genRM_gap","nb_famRM_gap","nb_cladRM_gap","nb_GFRM_gap","nb_PFTRM_gap",
                   "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred","dist_AVGF_pred","dist_AVPFT_pred",
                   "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs","dist_AVGF_obs","dist_AVPFT_obs",
#                   "dist_AVspecRM_pred","dist_AVgenRM_pred","dist_AVfamRM_pred","dist_AVcladRM_pred","dist_AVGFRM_pred","dist_AVPFTRM_pred",
#                   "dist_AVspecRM_obs","dist_AVgenRM_obs","dist_AVfamRM_obs","dist_AVcladRM_obs","dist_AVGFRM_obs","dist_AVPFTRM_obs",
                   "trait","ID","Species","Genus","Family","Clade","GF","PFT")
res_Obs_obs_TD <- res
res_Obs_obs <- res

TD_choice="Obs_obs"
TD_choice="Obs_obs_TD"

if(TD_choice=="Obs_obs_TD"){res=res_Obs_obs_TD}
if(TD_choice=="Obs_obs"){res=res_Obs_obs}
RepNum=1
trait_sub="rainfor"
GapPercent=60
for(GapPercent in gappercents){
  print(summary(res))
  print("------ % gaps now ----------")
  print(GapPercent)
  
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
  {
    # completely observed
    if(TD_choice=="Obs_obs_TD"&trait_sub=="rainfor"){traitInfo_obs=out$rainforTD_observed;taxInfo=out$rainforTD_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs_TD}
    if(TD_choice=="Obs_obs_TD"&trait_sub=="guido"){traitInfo_obs=out$guidoTD_observed;taxInfo=out$guidoTD_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs_TD}
    if(TD_choice=="Obs_obs"&trait_sub=="guido"){traitInfo_obs=out$guido_observed;taxInfo=out$guidoTD_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs}
    if(TD_choice=="Obs_obs"&trait_sub=="rainfor"){traitInfo_obs=out$rainfor_observed;taxInfoTD=out$rainforTD_tax;taxInfo=out$rainfor_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs}
    
    #sparse observed
    #  traitInfo_sparse <- load(file = file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","TestData_org.RData"))
    traitInfo_sparse <- read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","traitInfo.csv"), sep=",", dec=".")

    if(TD_choice=="Obs_obs"){ traitInfo_sparse <- traitInfo_sparse[which(as.numeric(taxInfoTD$IDs)%in%as.numeric(taxInfo$IDs)),colnames(traitInfo_sparse)%in%colnames(traitInfo_obs)]}
    
    # load the output for trait predictions: mean.csv
    traitInfo_pred_zlog <- as.matrix(read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), 
                                                          paste0("p",GapPercent,"_",trait_sub),TD_choice,"data/mean.csv"),
                                                sep="\t", dec=".",header=TRUE))
    if(TD_choice=="Obs_obs"){ traitInfo_pred_zlog <- traitInfo_pred_zlog[which(as.numeric(taxInfoTD$IDs)%in%as.numeric(taxInfo$IDs)),colnames(traitInfo_pred_zlog)%in%colnames(traitInfo_obs)]
    traitInfoTD_pred_zlog=traitInfo_pred_zlog}
    dim(traitInfoTD_pred_zlog)
    #-------------------------------------------------------------------
    # cut to test data only if necessary
    #-------------------------------------------------------------------
#    if(trait_sub=="guido"){ traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
#    if(trait_sub=="rainfor"){
#      traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]
#      traitInfo_sparse  <- traitInfo_sparse[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}
#    if(trait_sub=="guido"){ traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
#    if(trait_sub=="rainfor"){traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  
    
    #-------------------------------------------------------------------
    # back transformation necessary.
    #-------------------------------------------------------------------
    summary(traitInfo_sparse)
    summary(traitInfo_obs)
    summary(traitInfo_pred)#tbd
    summary(traitInfo_obs_zlog)#tbd
    summary(traitInfoTD_pred_zlog)
#    traitInfo_obs <- traitInfo_obs_c#NEW 2021
    t=3
    plot(traitInfoTD_pred_zlog[,t],traitInfo_obs_zlog[,t])
#    abline(0,1)
    
    traitInfo_obs_zlog <- traitInfo_obs
    for(i in 1:ncol(traitInfo_obs)){
      traitInfo_obs_zlog[,i] <- (log(traitInfo_obs[,i]) - scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),1])/
        scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),2]
    }
    
    summary(traitInfoTD_pred_zlog)
    summary(traitInfo_obs_zlog)
    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_zlog[,1])
    plot(traitInfo_pred[,1],traitInfo_obs[,1])
    abline(0,1)
    
    traitInfo_pred <- traitInfoTD_pred_zlog
    for(i in 1:ncol(traitInfoTD_pred_zlog)){
      traitInfo_pred[,i] <- exp((traitInfoTD_pred_zlog[,i]*scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),2])+ 
                                  scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),1])
    }
    
    plot(traitInfo_pred[,1],traitInfo_obs_c[,1])
#    abline(0,1)
    traitInfoTD_pred_zlog_sparse <- traitInfoTD_pred_zlog
    traitInfoTD_pred_zlog_sparse[is.na(traitInfo_sparse)] <- NA
    traitInfo_obs_zlog_sparse <- traitInfo_obs_zlog
    traitInfo_obs_zlog_sparse[is.na(traitInfo_sparse)] <- NA
    
  }

  
  #-----------------------------------------
  # correlations deviation to distance
  #-----------------------------------------
  #-----------------------------------------
  traitInfo_obs_zlogCOR <- traitInfo_obs_zlog
  traitInfo_obs_zlog_sparseCOR <- traitInfo_obs_zlog_sparse
  traitInfo_pred_zlogCOR <- traitInfo_pred_zlog
  
  t=1  
  for(t in 1:ncol(traitInfo_obs_zlog_sparse)){
    # -----------------------------------------
    # define the place in the output matrix
    ix_trait1 = 1:nrow(traitInfo_obs_zlog)
    ix_trait1 = ix_trait1  + 
      ((t-1)*nrow(traitInfo_obs_zlog)) 
    ix_trait = ix_trait1   + 
      ((which(GapPercent==gappercents)-1)*nrow(traitInfo_obs_zlog)) 
    #----------------------------------------------
    
    t1=2
    for(t1 in 1:ncol(traitInfo_obs_zlog_sparseCOR)){
        for(t in t){
          if(t1!=t){
            colnm4=paste0(colnames(traitInfo_obs_zlog_sparseCOR)[t],"_from_",colnames(traitInfo_obs_zlog_sparseCOR)[t1],"_lm")
            res <- add_col_to_res(new.col.names = colnm4,input = res)
            #observed
            dat_now=traitInfo_obs_zlogCOR[,c(t1,t)]
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now) # lm produced from NONsparse data!
            t_lm_obs <- (traitInfo_obs_zlogCOR[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
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
    cats <- factor(taxInfo$AccSpeciesName)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_spec_obs<- cbind(as.vector(unique(taxInfo$AccSpeciesName)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_spec_pred<- cbind(as.vector(unique(taxInfo$AccSpeciesName)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Genus)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_gen_obs<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_gen_pred<- cbind(as.vector(unique(taxInfo$Genus)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$Family)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_fam_obs<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_fam_pred<- cbind(as.vector(unique(taxInfo$Family)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PhylogeneticGroup)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    si2 <- silhouette(taxnow, dist(traitInfo_obs_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_clad_obs<- cbind(as.vector(unique(taxInfo$PhylogeneticGroup)),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow, dist(traitInfo_pred_zlog,method = "canberra"))
    si <- summary(si2)
    Sil_clad_pred<- cbind(as.vector(unique(taxInfo$PhylogeneticGroup)),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$GF)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_obs<- cbind(as.vector(unique(taxInfo$GF[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_GF_pred<- cbind(as.vector(unique(taxInfo$GF[ixG])),as.vector(si$clus.avg.widths))
    
    cats <- factor(taxInfo$PFT)
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow <- DF[,2]
    ixG <- !is.na(taxnow)
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_obs_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_obs<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
    si2 <- silhouette(taxnow[ixG], dist(traitInfo_pred_zlog[ixG,],method = "canberra"))
    si <- summary(si2)
    Sil_PFT_pred<- cbind(as.vector(unique(taxInfo$PFT[ixG])),as.vector(si$clus.avg.widths))
  }
  
  t=1
  for(t in 1:ncol(traitInfo_obs_zlog)){
    print("------------trait:")
    print(t)
    ix_trait1 = 1:nrow(traitInfo_obs_zlog)
    ix_trait1 = ix_trait1  + 
      ((t-1)*nrow(traitInfo_obs_zlog)) 
    ix_trait = ix_trait1   + 
      ((which(GapPercent==gappercents)-1)*nrow(traitInfo_obs_zlog)) 
    ix_trait=ix_trait1
    print("------------from to:")
    print(min(ix_trait,na.rm = TRUE))
    print(max(ix_trait,na.rm = TRUE))
    
    res$trait[ix_trait] <- colnames(traitInfo_obs_zlog)[t]
    res$ID[ix_trait] <- as.vector(taxInfo$IDs)
    res$Species[ix_trait] <- as.vector(taxInfo$AccSpeciesName)
    res$Genus[ix_trait] <- as.vector(taxInfo$Genus)
    res$Family[ix_trait] <- as.vector(taxInfo$Family)
    res$Clade[ix_trait] <- as.vector(taxInfo$PhylogeneticGroup)
    res$GF[ix_trait] <- as.vector(taxInfo$GF)
    res$PFT[ix_trait] <- as.vector(taxInfo$PFT)
    res$error[ix_trait] <- traitInfoTD_pred_zlog[,t]-traitInfo_obs_zlog[,t]
    res$Gap[ix_trait] <- is.na(traitInfo_sparse[,t])
    res$value_obs_zlog[ix_trait] <- traitInfo_obs_zlog[,t]
    res$value_pred_zlog[ix_trait] <- traitInfoTD_pred_zlog[,t]
    res$value_obs[ix_trait] <- traitInfo_obs[,t]
    res$value_pred[ix_trait] <- traitInfo_pred[,t]
    res$missingness[ix_trait] <- GapPercent
    res$DataSet[ix_trait] <- trait_sub
    res$nrep[ix_trait] <- as.vector(RepNum)
    
    print("------------missingness---------------")
    print(unique(res$missingness))

    i=80
    for(i in 1: nrow(traitInfo_obs_zlog)){
     # print(i)    
      {
      #----------------------------------------------
      #Individuals within species
      #----------------------------------------------

      ix1 = taxInfo$AccSpeciesName==as.vector(taxInfo$AccSpeciesName[i])
      res$mean_spec_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_spec_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
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
      res$Sil_spec_pred[ix_trait[i]] <- as.numeric(Sil_spec_pred[which(Sil_spec_obs[,1]==res$Species[ix_trait[i]]),2])
      #----------------------------------------------
      #Species within genera
      #----------------------------------------------
      ix1 = taxInfo$Genus==as.vector(taxInfo$Genus[i])
      res$mean_gen_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_gen_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
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
      #----------------------------------------------
      #Genera within families
      #----------------------------------------------
      ix1 = taxInfo$Family==as.vector(taxInfo$Family[i])
      res$mean_fam_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_fam_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
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
      #----------------------------------------------
      #Families within Clades
      #----------------------------------------------
      ix1 = taxInfo$PhylogeneticGroup==as.vector(taxInfo$PhylogeneticGroup[i])
      res$mean_clad_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_clad_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
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
      
      #---------------------------------------------------------------------
      # GF
      #---------------------------------------------------------------------
      ix1 = taxInfo$GF==as.vector(taxInfo$GF[i])
      ix1[is.na(ix1)] <- FALSE
      res$mean_GF_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_GF[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix1[is.na(ix1)] <- FALSE
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_GF[ix_trait[i]] <- sum(ix1,na.rm = TRUE)
      if(!is.na(as.vector(taxInfo$GF[i]))){
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
      }}
      if(!is.na(res$GF[ix_trait[i]])){
      try(res$Sil_GF_obs[ix_trait[i]] <- as.numeric(Sil_GF_obs[which(Sil_GF_obs[,1]==res$GF[ix_trait[i]]),2]))
      try(res$Sil_GF_pred[ix_trait[i]] <- as.numeric(Sil_GF_pred[which(Sil_GF_pred[,1]==res$GF[ix_trait[i]]),2]))
      }
      #---------------------------------------------------------------------
      # PFT
      #---------------------------------------------------------------------
      ix1 = taxInfo$PFT==as.vector(taxInfo$PFT[i])
      ix1[is.na(ix1)] <- FALSE
      res$mean_PFT_sparse[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT_zlog[ix_trait[i]] <- mean(traitInfoTD_pred_zlog[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT_obs[ix_trait[i]] <- mean(traitInfo_obs[ix1,t],na.rm = TRUE) #input species mean
      res$mean_PFT[ix_trait[i]] <- mean(traitInfo_pred[ix1,t],na.rm = TRUE) #input species mean
      ix1[is.na(ix1)] <- FALSE
      ix <- ix1;ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      res$nb_PFT[ix_trait[i]] <- sum(ix1,na.rm = TRUE)
      if(!is.na(as.vector(taxInfo$PFT[i]))){
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
        }
      #----------------------------------------------------------------------
      }
      # basta

    }#i
    
  }#t

  }#GapPercent

res_save <- res
summary(res_save)

list.files(file.path(origin,"_2021","data","analyes"))
write.csv(res,file=file.path(origin,"_2021","data","analyses","Point_wise",paste0("res",,".csv")))
res<- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",paste0("res",,".csv")))
summary(res$missingness)


summary(res)





if(sum(colnames(res)%in%colnm)==0){
  newcol=rep(NA,nrow(res))
  res= cbind(res,newcol)
  colnames(res)[ncol(res)] <- colnm
}

