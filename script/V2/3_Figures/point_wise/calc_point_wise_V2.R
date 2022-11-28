
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
gappercents <- c(1,60,5,50,10,40,20,30)

GapPercent=60
res <- as.data.frame(matrix(NA,ncol=86,nrow=ncol(traitInfo_obs_zlog)*nrow(traitInfo_obs_zlog)*length(gappercents)))
colnames(res) <- c("error","Gap","missingness","value_obs","value_pred","nrep",
                   "dist_spec","dist_gen","dist_fam","dist_clad","dist_GF","dist_PFT",
                   "nb_spec","nb_gen","nb_fam","nb_clad","nb_GF","nb_PFT",
                   "mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",#input mean
                   "dist_specRM","dist_genRM","dist_famRM","dist_cladRM","dist_GFRM","dist_PFTRM",
                   "nb_specRM","nb_genRM","nb_famRM","nb_cladRM","nb_GFRM","nb_PFTRM",
                   "mean_specRM","mean_genRM","mean_famRM","mean_cladRM","mean_GFRM","mean_PFTRM",#input mean
                   # missing: nb_GFRM_gap
                   "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap","nb_GF_gap","nb_PFT_gap",
                   "nb_specRM_gap","nb_genRM_gap","nb_famRM_gap","nb_cladRM_gap","nb_GFRM_gap","nb_PFTRM_gap",
                   "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred","dist_AVGF_pred","dist_AVPFT_pred",
                   "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs","dist_AVGF_obs","dist_AVPFT_obs",
                   "dist_AVspecRM_pred","dist_AVgenRM_pred","dist_AVfamRM_pred","dist_AVcladRM_pred","dist_AVGFRM_pred","dist_AVPFTRM_pred",
                   "dist_AVspecRM_obs","dist_AVgenRM_obs","dist_AVfamRM_obs","dist_AVcladRM_obs","dist_AVGFRM_obs","dist_AVPFTRM_obs",
                   "trait","ID","Species","Genus","Family","Clade","GF","PFT")

RepNum=1
trait_sub="rainfor"
TD_choice="Obs_obs_TD"
for(GapPercent in gappercents){
  print("----------------------")
  print(GapPercent)
  print("----------------------")
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
  {
    # completely observed
    if(TD_choice=="Obs_obs_TD"&trait_sub=="rainfor"){traitInfo_obs=out$rainforTD_observed;taxInfo=out$rainforTD_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs_TD}
    if(TD_choice=="Obs_obs_TD"&trait_sub=="guido"){traitInfo_obs=out$guidoTD_observed;taxInfo=out$guidoTD_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs_TD}
    if(TD_choice=="Obs_obs"&trait_sub=="rainfor"){traitInfo_obs=out$rainfor_observed;taxInfo=out$rainfor_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs}
    if(TD_choice=="Obs_obs"&trait_sub=="guido"){traitInfo_obs=out$guido_observed;taxInfo=out$guido_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs}
    #sparse observed
    #  traitInfo_sparse <- load(file = file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","TestData_org.RData"))
    traitInfo_sparse <- read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","traitInfo.csv"), sep=",", dec=".")
    # load the output for trait predictions: mean.csv
    traitInfo_pred_zlog <- as.matrix(read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), 
                                                          paste0("p",GapPercent,"_",trait_sub),TD_choice,"data/mean.csv"),
                                                sep="\t", dec=".",header=TRUE))
    
    
    #-------------------------------------------------------------------
    # cut to test data only if necessary
    #-------------------------------------------------------------------
    if(trait_sub=="guido"){ traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
    if(trait_sub=="rainfor"){traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  
    if(trait_sub=="guido"){ traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
    if(trait_sub=="rainfor"){traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  
    
    #-------------------------------------------------------------------
    # back transformation necessary.
    #-------------------------------------------------------------------
    summary(traitInfo_sparse)
    summary(traitInfo_obs)
    summary(traitInfo_obs_c)
    summary(traitInfoTD_pred_zlog)
    
#    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_c[,1])
#    abline(0,1)
    
    traitInfo_obs_zlog <- traitInfo_obs
    for(i in 1:ncol(traitInfo_obs)){
      traitInfo_obs_zlog[,i] <- (log(traitInfo_obs[,i]) - scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),1])/
        scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),2]
    }
    
#    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_zlog[,1])
#    abline(0,1)
    
    traitInfo_pred <- traitInfoTD_pred_zlog
    for(i in 1:ncol(traitInfoTD_pred_zlog)){
      traitInfo_pred[,i] <- exp((traitInfoTD_pred_zlog[,i]*scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),2])+ 
                                  scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),1])
    }
    
#    plot(traitInfo_pred[,1],traitInfo_obs_c[,1])
#    abline(0,1)
    
  }
  
  #---------------------------------------------------------------------------------------------------------------
  traitInfoTD_pred_zlog_sparse <- traitInfoTD_pred_zlog
  traitInfoTD_pred_zlog_sparse[is.na(traitInfo_sparse)] <- NA
  traitInfo_obs_zlog_sparse <- traitInfo_obs_zlog
  traitInfo_obs_zlog_sparse[is.na(traitInfo_sparse)] <- NA
  
  
  #-----------------------------------------
  # correlations
  #-----------------------------------------
t=1  
  for(t in 1:ncol(traitInfo_obs_zlog)){
    #----------------------------------------------
    # correlations deviation to distance
    #----------------------------------------------
    ix_trait1 = 1:nrow(traitInfo_obs_zlog) + 
      ((t-1)*nrow(traitInfo_obs_zlog)) 
    ix_trait = ix_trait1 + (ix_trait1*(1-which(GapPercent%in%gappercents)))
    for(t1 in 1:ncol(traitInfo_obs_zlog)){
        for(t in t){
          if(t1!=t){
            colnm1=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESobs")
            colnm2=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESpred")
            colnm3=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESdev")
            colnm4=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_lm")
            res <- add_col_to_res(input = res,new.col.names = colnm1)
            res <- add_col_to_res(input = res,new.col.names = colnm2)
            res <- add_col_to_res(input = res,new.col.names = colnm3)
            res <- add_col_to_res(input = res,new.col.names = colnm4)
          #observed
            dat_now=traitInfo_obs_zlog[,c(t1,t)]
            names(dat_now) <- c("t1","t")
            lm_now <- lm(t1~t,data = dat_now)
            #          plot(dat_now)
            #          abline(lm_now)
            t_lm_obs <- (traitInfo_obs_zlog[,1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            t_res_obs <- (t_lm_obs-traitInfo_obs_zlog[,2])
            #         points(traitInfo_obs_zlog[,1],t_lm_obs,col="orange",pch=16)
            
            #HPMF predicted
            dat_now=data.frame(t1=traitInfo_pred_zlog[,t1],
                               t=traitInfo_pred_zlog[,t])
            t_lm_pred <- (traitInfo_pred_zlog[,t1]*lm_now$coefficients[2])+lm_now$coefficients[1]
            #         points(traitInfo_pred_zlog[,1],t_lm_pred,col="hotpink",pch=16)
            t_res_pred <- (t_lm_pred-traitInfo_obs_zlog[,t])
            deviation=t_res_pred-t_res_obs
#            boxplot(deviation,ylim=c(-4,4))
            res[ix_trait,colnames(res)%in%colnm1] <- t_res_obs
            res[ix_trait,colnames(res)%in%colnm2] <- t_res_pred
            res[ix_trait,colnames(res)%in%colnm3] <- deviation
            res[ix_trait,colnames(res)%in%colnm4] <- t_lm_obs
          }
        }
      }
  }  
  
  #-----------------------------------------
  # taxonomy
  #-----------------------------------------
  # Randomize the taxonomic information :)
  taxInfoRM <- taxInfo[sample(1:nrow(taxInfo),nrow(taxInfo)),]
  t=1
  for(t in 1:ncol(traitInfo_obs_zlog)){
    print(t)
    ix_trait1 = 1:nrow(traitInfo_obs_zlog) + 
      ((t-1)*nrow(traitInfo_obs_zlog)) 
    ix_trait = ix_trait1 + (ix_trait1*(1-which(GapPercent%in%gappercents)))
    
    
    res$trait[ix_trait] <- colnames(traitInfo_obs_zlog)[t]
#    res$rnep[ix_trait] <- RepNum
    res$ID[ix_trait] <- as.vector(taxInfo$IDs)
    res$Species[ix_trait] <- as.vector(taxInfo$AccSpeciesName)
    res$Genus[ix_trait] <- as.vector(taxInfo$Genus)
    res$Family[ix_trait] <- as.vector(taxInfo$Family)
    res$Clade[ix_trait] <- as.vector(taxInfo$PhylogeneticGroup)
    res$GF[ix_trait] <- as.vector(taxInfo$GF)
    res$PFT[ix_trait] <- as.vector(taxInfo$PFT)
    res$error[ix_trait] <- traitInfoTD_pred_zlog[,t]-traitInfo_obs_zlog[,t]
    res$Gap[ix_trait] <- is.na(traitInfo_sparse[,t])
    res$value_obs[ix_trait] <- traitInfo_obs_zlog[,t]
    res$value_pred[ix_trait] <- traitInfoTD_pred_zlog[,t]
    res$missingness[ix_trait] <- GapPercent
    
    i=1
    for(i in 1: nrow(traitInfo_obs_zlog)){
      print(i)    
      {
      #----------------------------------------------
      #Individuals within species
      #----------------------------------------------
      ix1 = taxInfo$AccSpeciesName==as.vector(taxInfo$AccSpeciesName[i])
      res$mean_spec[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input species mean
      res$nb_spec_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_spec[ix_trait[i]] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_spec[ix_trait[i]] <- 0
        res$dist_AVspec_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVspec_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_spec[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVspec_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVspec_obs[ix_trait[i]]  <- traitInfo_obs_zlog[i,t]- mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #RANDOM
      ix1 = taxInfoRM$AccSpeciesName==as.vector(taxInfoRM$AccSpeciesName[i])
      res$mean_specRM[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input species mean
      res$nb_spec_gapRM[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_specRM[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_specRM[ix_trait[i]] <- 0
        res$dist_AVspecRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVspecRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_specRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVspecRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVspecRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      
      #----------------------------------------------
      #Species within genera
      #----------------------------------------------
      ix1 = taxInfo$Genus==as.vector(taxInfo$Genus[i])
      res$mean_gen[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_gen_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_gen[ix_trait[i]] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_gen[ix_trait[i]] <- 0
        res$dist_AVgen_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVgen_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_gen[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVgen_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVgen_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #RANDOM
      ix1 = taxInfoRM$Genus==as.vector(taxInfoRM$Genus[i])
      res$mean_genRM[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_genRM_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_genRM[ix_trait[i]] = sum(ix1)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_genRM[ix_trait[i]] <- 0
        res$dist_AVgenRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVgenRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_genRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVgenRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVgenRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #----------------------------------------------
      #Genera within families
      #----------------------------------------------
      ix1 = taxInfo$Family==as.vector(taxInfo$Family[i])
      res$mean_fam[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_fam_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_fam[ix_trait[i]] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_fam[ix_trait[i]] <- 0
        res$dist_AVfam_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVfam_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_fam[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVfam_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVfam_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #RANDOM
      ix1 = taxInfoRM$Family==as.vector(taxInfoRM$Family[i])
      res$mean_famRM[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_famRM_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_famRM[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_famRM[ix_trait[i]] <- 0
        res$dist_AVfamRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVfamRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_famRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVfamRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVfamRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #----------------------------------------------
      #Families within Clades
      #----------------------------------------------
      ix1 = taxInfo$PhylogeneticGroup==as.vector(taxInfo$PhylogeneticGroup[i])
      res$mean_clad[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_clad_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_clad[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_clad[ix_trait[i]] <- 0
        res$dist_AVclad_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVclad_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_clad[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVclad_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVclad_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #RANDOM
      ix1 = taxInfoRM$PhylogeneticGroup==as.vector(taxInfoRM$PhylogeneticGroup[i])
      res$mean_cladRM[ix_trait[i]]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_cladRM_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_cladRM[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_cladRM[ix_trait[i]] <- 0
        res$dist_AVcladRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVcladRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_cladRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVcladRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVcladRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      #---------------------------------------------------------------------
      # GF
      #---------------------------------------------------------------------
      ix1 = taxInfo$GF==as.vector(taxInfo$GF[i])
      ix1[is.na(ix1)] <- FALSE
      if(!is.na(taxInfo$GF[i])){
      res$mean_GF[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input GFies mean
      res$nb_GF_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_GF[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_GF[ix_trait[i]] <- 0
        res$dist_AVGF_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVGF_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_GF[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVGF_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVGF_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }}
      #RANDOM
      ix1 = taxInfoRM$GF==as.vector(taxInfoRM$GF[i])
      ix1[is.na(ix1)] <- FALSE
      if(!is.na(taxInfo$GF[i])){
      res$mean_GFRM[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input GFies mean
      res$nb_GFRM_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_GFRM[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_GFRM[ix_trait[i]] <- 0
        res$dist_AVGFRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVGFRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_GFRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVGFRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVGFRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }}
      #---------------------------------------------------------------------
      # GF
      #---------------------------------------------------------------------
      ix1 = taxInfo$PFT==as.vector(taxInfo$PFT[i])
      ix1[is.na(ix1)] <- FALSE
      if(!is.na(taxInfo$PFT[i])){
        res$mean_PFT[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input PFTies mean
      res$nb_PFT_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_PFT[ix_trait[i]] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_PFT[ix_trait[i]] <- 0
        res$dist_AVPFT_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVPFT_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_PFT[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVPFT_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVPFT_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }}
      #RANDOM
      ix1 = taxInfoRM$PFT==as.vector(taxInfoRM$PFT[i])
      ix1[is.na(ix1)] <- FALSE
      if(!is.na(taxInfo$PFT[i])){
        res$mean_PFTRM[ix_trait[i]] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input PFTies mean
      res$nb_PFTRM_gap[ix_trait[i]] = sum(!is.na(traitInfo_sparse[ix1,t]),na.rm = TRUE)
      res$nb_PFTRM[ix_trait[i]] = sum(ix1,na.rm = TRUE)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_PFTRM[ix_trait[i]] <- 0
        res$dist_AVPFTRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVPFTRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }else{
        res$dist_PFTRM[ix_trait[i]] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1],na.rm = TRUE)
        res$dist_AVPFTRM_pred[ix_trait[i]] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog_sparse[ix1,t],na.rm = TRUE)
        res$dist_AVPFTRM_obs[ix_trait[i]] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t],na.rm = TRUE)
      }
      }
      
      
      }
      #----------------------------------------------------------------------
      # basta

    }
  }

  }

res_save <- res
head(res_save)

list.files(file.path(origin,"_2021","data","analyes"))
write.csv(res,file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))
res<- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))
head(res)








if(sum(colnames(res)%in%colnm)==0){
  newcol=rep(NA,nrow(res))
  res= cbind(res,newcol)
  colnames(res)[ncol(res)] <- colnm
}

