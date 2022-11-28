
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
gappercents <- c(1,5,10,20,30 , 40 , 50 , 60)

GapPercent=60
res <- as.data.frame(matrix(NA,ncol=49,nrow=ncol(traitInfo_obs_zlog)*nrow(traitInfo_obs_zlog)*length(gappercents)))
colnames(res) <- c("error","Gap","missingness","value_obs","value_pred",
                   "dist_spec","dist_gen","dist_fam","dist_clad",
                   "dist_specRM","dist_genRM","dist_famRM","dist_cladRM",
                   "nb_spec","nb_gen","nb_fam","nb_clad",
                   "mean_spec","mean_gen","mean_fam","mean_clad",#input mean
                   "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap",
                   "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred",
                   "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs",
                   "dist_AVspecRM_pred","dist_AVgenRM_pred","dist_AVfamRM_pred","dist_AVcladRM_pred",
                   "dist_AVspecRM_obs","dist_AVgenRM_obs","dist_AVfamRM_obs","dist_AVcladRM_obs",
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
    
    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_c[,1])
    abline(0,1)
    
    traitInfo_obs_zlog <- traitInfo_obs
    for(i in 1:ncol(traitInfo_obs)){
      traitInfo_obs_zlog[,i] <- (log(traitInfo_obs[,i]) - scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),1])/
        scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),2]
    }
    
    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_zlog[,1])
    abline(0,1)
    
    traitInfo_pred <- traitInfoTD_pred_zlog
    for(i in 1:ncol(traitInfoTD_pred_zlog)){
      traitInfo_pred[,i] <- exp((traitInfoTD_pred_zlog[,i]*scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),2])+ 
                                  scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),1])
    }
    
    plot(traitInfo_pred[,1],traitInfo_obs_c[,1])
    abline(0,1)
    
  }
  
  #---------------------------------------------------------------------------------------------------------------
  traitInfoTD_pred_zlog_sparse <- traitInfoTD_pred_zlog
  traitInfoTD_pred_zlog_sparse[is.na(traitInfo_sparse)] <- NA
  
  
  #-----------------------------------------
  # Randomize the taxonomic information :)
  #-----------------------------------------
  taxInfoRM <- taxInfo[sample(1:nrow(taxInfo),nrow(taxInfo)),]
  
  for(t in 1:ncol(traitInfo_obs_zlog)){
    head(res)
    #----------------------------------------------
    # correlations deviation to distance
    #----------------------------------------------
    ix_trait = 1:nrow(traitInfo_obs_zlog) + ((t-1)*nrow(traitInfo_obs_zlog))
    for(t1 in 1:ncol(traitInfo_obs_zlog)){
        for(t in t){
          if(t1!=t){
            colnm1=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESobs")
            colnm2=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESpred")
            colnm3=paste0(colnames(traitInfo_obs_zlog)[t],"_from_",colnames(traitInfo_obs_zlog)[t1],"_RESdev")
            res <- add_col_to_res(res = res,new.col.names = colnm1)
            res <- add_col_to_res(res = res,new.col.names = colnm2)
            res <- add_col_to_res(res = res,new.col.names = colnm3)
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
            boxplot(deviation,ylim=c(-4,4))
            res[ix_trait,colnames(res)%in%colnm1] <- t_res_obs
            res[ix_trait,colnames(res)%in%colnm2] <- t_res_pred
            res[ix_trait,colnames(res)%in%colnm3] <- deviation
          }
        }
      }
  }  
  
  head(res)
  n=1
  for(t in 1:ncol(traitInfo_obs_zlog)){
    
    res[ix_trait,colnames(res)%in%"trait"] <- colnames(traitInfo_obs_zlog)[t]
    res[ix_trait,colnames(res)%in%"ID"] <- as.vector(taxInfo$IDs)
    res[ix_trait,colnames(res)%in%"Species"] <- as.vector(taxInfo$AccSpeciesName)
    res[ix_trait,colnames(res)%in%"Genus"] <- as.vector(taxInfo$Genus)
    res[ix_trait,colnames(res)%in%"Family"] <- as.vector(taxInfo$Family)
    res[ix_trait,colnames(res)%in%"Clade"] <- as.vector(taxInfo$PhylogeneticGroup)
    res[ix_trait,colnames(res)%in%"GF"] <- as.vector(taxInfo$GF)
    res[ix_trait,colnames(res)%in%"PFT"] <- as.vector(taxInfo$PFT)
    res[ix_trait,colnames(res)%in%"error"] <- traitInfoTD_pred_zlog[,t]-traitInfo_obs_zlog[,t]
    res[ix_trait,colnames(res)%in%"Gap"] <- is.na(traitInfo_sparse[,t])
    res[ix_trait,colnames(res)%in%"value_obs"] <- traitInfo_obs_zlog[,t]
    res[ix_trait,colnames(res)%in%"value_pred"] <- traitInfoTD_pred_zlog[,t]
    res[ix_trait] <- GapPercent
    i=1
    for(i in 1: nrow(traitInfo_obs_zlog)){
      print(n)    
      {
      #----------------------------------------------
      #Individuals within species
      #----------------------------------------------
      ix1 = taxInfo[,colnames(taxInfo)%in%"AccSpeciesName"]==as.vector(taxInfo[,colnames(taxInfo)%in%"AccSpeciesName"][i])
      res[n,colnames(res)%in%"mean_spec"] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input species mean
      res[n,colnames(res)%in%"nb_spec_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res[n,colnames(res)%in%"nb_spec"] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_spec"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_spec"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVspec_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVspec_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #RANDOM
      ix1 = taxInfoRM[,colnames(taxInfoRM)%in%"AccSpeciesName"]==as.vector(taxInfoRM[,colnames(taxInfoRM)%in%"AccSpeciesName"][i])
      #    res$mean_spec[n] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input species mean
      #    res$nb_spec_gap[n] = sum(!is.na(traitInfo_sparse[ix1,t]))
      #    res$nb_spec[n] = sum(ix1)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_specRM"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_specRM"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVspecRM_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVspecRM_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      
      #----------------------------------------------
      #Species within genera
      #----------------------------------------------
      ix1 = taxInfo[,colnames(taxInfo)%in%"Genus"]==as.vector(taxInfo[,colnames(taxInfo)%in%"Genus"][i])
      res[n,colnames(res)%in%"mean_gen"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res[n,colnames(res)%in%"nb_gen_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res[n,colnames(res)%in%"nb_gen"] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_gen"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_gen"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVgen_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVgen_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #RANDOM
      ix1 = taxInfoRM[,colnames(taxInfoRM)%in%"Genus"]==as.vector(taxInfoRM[,colnames(taxInfoRM)%in%"Genus"][i])
#      res[n,colnames(res)%in%"mean_genRM"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
#      res[n,colnames(res)%in%"nb_genRM_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
#      res[n,colnames(res)%in%"nb_genRM"] = sum(ix1)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_genRM"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_genRM"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVgenRM_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVgenRM_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #----------------------------------------------
      #Genera within families
      #----------------------------------------------
      ix1 = taxInfo[,colnames(taxInfo)%in%"Family"]==as.vector(taxInfo[,colnames(taxInfo)%in%"Family"][i])
      res[n,colnames(res)%in%"mean_fam"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res[n,colnames(res)%in%"nb_fam_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res[n,colnames(res)%in%"nb_fam"] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_fam"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_fam"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVfam_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVfam_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #RANDOM
      ix1 = taxInfoRM[,colnames(taxInfoRM)%in%"Family"]==as.vector(taxInfoRM[,colnames(taxInfoRM)%in%"Family"][i])
      #      res[n,colnames(res)%in%"mean_famRM"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      #      res[n,colnames(res)%in%"nb_famRM_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      #      res[n,colnames(res)%in%"nb_famRM"] = sum(ix1)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_famRM"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_famRM"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVfamRM_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVfamRM_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #----------------------------------------------
      #Families within Clades
      #----------------------------------------------
      ix1 = taxInfo[,colnames(taxInfo)%in%"PhylogeneticGroup"]==as.vector(taxInfo[,colnames(taxInfo)%in%"PhylogeneticGroup"][i])
      res[n,colnames(res)%in%"mean_clad"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res[n,colnames(res)%in%"nb_clad_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res[n,colnames(res)%in%"nb_clad"] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_clad"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_clad"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVclad_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVclad_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      #RANDOM
      ix1 = taxInfoRM[,colnames(taxInfoRM)%in%"PhylogeneticGroup"]==as.vector(taxInfoRM[,colnames(taxInfoRM)%in%"PhylogeneticGroup"][i])
      #      res[n,colnames(res)%in%"mean_cladRM"]  <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      #      res[n,colnames(res)%in%"nb_cladRM_gap"] = sum(!is.na(traitInfo_sparse[ix1,t]))
      #      res[n,colnames(res)%in%"nb_cladRM"] = sum(ix1)
      ix <- ix1
      ix[taxInfoRM$ID==as.vector(taxInfoRM$ID[i])] <- FALSE
      if(sum(ix)==0){
        res[n,colnames(res)%in%"dist_cladRM"] <- 0
      }else{
        res[n,colnames(res)%in%"dist_cladRM"] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res[n,colnames(res)%in%"dist_AVcladRM_pred"] <- traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
        res[n,colnames(res)%in%"dist_AVcladRM_obs"] <- traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t])
      }
      }
      #----------------------------------------------------------------------
      # basta
      n=n+1
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

