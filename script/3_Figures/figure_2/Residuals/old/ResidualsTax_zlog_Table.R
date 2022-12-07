
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
  GapPercent=80
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  mean.fun <- function(input){
    if(sum(!is.na(input))==0){return(NA)}else{
      return(mean(input[!is.na(input)]))}
  }
  count.fun <- function(input){
    return(sum(!is.na(input)))
  }
  range.fun <- function(input){
    if(sum(!is.na(input))==0){return(NA)}else{
      return(max(input[!is.na(input)])-min(input[!is.na(input)]))
    }
  }
  


#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
t_choice="data"
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}
if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}
if(t_choice=="data"){trait_names=trait_rainfor}
if(t_choice=="data_2"){trait_names=trait_guido}

txs=c("indiv","spec","gen","fam","clad")
tx=1 #(0=indiv,1=spec,2=genus,3=family,4=clades)
#for(tx in 0:4){
  tx=tx+1
  par(mfrow=c(4,4))
    
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
    file.path(originData,"analyes","Point_wise","res.csv")
    colz=colz1
    
    #load data
    {
      # load Envelope data
      # load TDext
      ObsOrTD="Obs_obs"
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
      
      # load TD data
      # total trait data 
      TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-c(1,2)]
      TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                    "Obs_obs_TD","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-c(1,2)]
      TDtd <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                               "Obs_obs_TD","data","traitInfoTD_pred_REzlog.csv"),header=TRUE))[,-c(1,2)]
      TD_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                     "Obs_obs_TD","data","taxInfo.csv"), sep=",")
      colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
      head(TD_tax)
      colnames(TD_tax)
      dim(TD_tax)
      # load Envelope data
      # total trait data 
      ExTD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                              "Obs_obs","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-1]

      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
      # predicted 
      TDext <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs","data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
      head(TDext)
      head(TDtd)
      head(TDext)
      # plot(TD[,1],TDtd[,1])
      #  abline(0,1)
      #  plot(TD[,1],TDext[,1])
      #  abline(0,1)
    }
    
    trait_names=colnames(TD)
    
    #_________________________________________________
    # aggregate TD to taxon mean
    #-------------------------------------------------

    TD_tx <- aggregate(TD,by=list(TD_tax[,tx]),FUN = mean.fun)
    TD_sparse_tx <- aggregate(TD_sparse,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDtd_tx <- aggregate(TDtd,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDext_tx <- aggregate(TDext,by=list(TD_tax[,tx]),FUN = mean.fun)
    
    TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
    TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
    TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
    TDext_tx <- TDext_tx[,-which(colnames(TDext_tx)=="Group.1")]
    

    #_________________________________________________
    # taxon mean: residuals (same aggregation)
    #-------------------------------------------------
      dist_spec_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_gen_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_fam_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_clad_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      
      # calc distance to taxon mean
      t=1
      for(t in 1:ncol(TD)){
        tax <- TD_tax[,2]
        for(i in 1:nrow(TD)){
          s = TD[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_spec_mean[i,t] <- TD[i,t]-sm
        }
        tax <- TD_tax[,3]
        for(i in 1:nrow(TD)){
          s = TD[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_gen_mean[i,t] <- TD[i,t]-sm
        }
        tax <- TD_tax[,4]
        for(i in 1:nrow(TD)){
          s = TD[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_fam_mean[i,t] <- TD[i,t]-sm
        }
        tax <- TD_tax[,5]
        for(i in 1:nrow(TD)){
          s = TD[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_clad_mean[i,t] <- TD[i,t]-sm
        }
      }
    
    par(mfrow=c(1,1))  
    
    #_________________________________________________
    # aggregate OBS or TD to taxon mean
    #-------------------------------------------------
    dist_spec_mean_tx <- aggregate(abs(dist_spec_mean),by=list(TD_tax[,2]),FUN = mean.fun)
      dist_spec_mean_tx <- dist_spec_mean_tx[,-which(colnames(dist_spec_mean_tx)=="Group.1")]

    #_________________________________________________
    # calculate residuals for single values
    #-------------------------------------------------
    # for OBS and IMP_obs
    resid <- abs(cbind(TD-TDtd))
    resid_tx <- abs(cbind(TD_tx-TDtd_tx))
    dim(resid)
    # for OBS and IMP_obsExt
    residE <- abs(cbind(TD-TDext))
    residE_tx <- abs(cbind(TD_tx-TDext_tx))
    dim(residE)
    t=1
    
    #_________________________________________________
    # calculate average residuals for whole taxa
    #-------------------------------------------------
    # average distance to species mean
    distTDtd_TD <- aggregate(abs(TDtd-TD),by=list(TD_tax[,2]),FUN = mean.fun)
    distTDtd_TD <- distTDtd_TD[,-which(colnames(distTDtd_TD)=="Group.1")]
    colMeans(distTDtd_TD)
    
    # average distance to species mean IMP_obsExt
    distTDext_TD <- aggregate(abs(TDext-TD),by=list(TD_tax[,2]),FUN = mean.fun)
    distTDext_TD <- distTDtd_TD[,-which(colnames(distTDtd_TD)=="Group.1")]
    colMeans(distTDtd_TD)
    
    plot(unlist(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t])),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t]))
    print(cor(unlist(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t])),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t])))
    abline(0,1)
    
    for(t in 1:ncol(TD)){
#      print(cor(abs(dist_spec_mean[dist_spec_mean[,t]!=0,t]),resid[dist_spec_mean[,t]!=0,t]))
#      plot(abs(dist_spec_mean[dist_spec_mean[,t]!=0,t]),resid[dist_spec_mean[,t]!=0,t])
      print(cor(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t]),resid_tx[dist_spec_mean_tx[,t]!=0,t]))
      plot(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t]),resid_tx[dist_spec_mean_tx[,t]!=0,t])
      abline(0,1)
    }
    
}

