# atm S7: Relative error to trait range for test data 2. Error (absolute) calculated
# as imputed - observed for test data imputed with and without the envelope (back-transformed data 

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
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
GapPercent=50
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
for(tx in 0:4){
  tx=tx+1
  par(mfrow=c(4,4))
    
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
    file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
    colz=colz1
    
    #load data
    {
      # load Envelope data
      # load TDenvelope
      ObsOrTD="Obs_obs"
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
      
      # load TD data
      # total trait data 
      TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
      TD_sparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                    "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
      TDtd <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                               "Obs_obs_TD","data","traitInfoTD_pred.csv"),header=TRUE))[,-c(1,2)]
      TD_tax <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                     "Obs_obs_TD","data","taxInfo.csv"), sep=",")
      colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
      head(TD_tax)
      colnames(TD_tax)
      dim(TD_tax)
      # load Envelope data
      # total trait data 
      Env <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                              "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
      summary(Env)
      Env <- Env[,colnames(Env)%in%colnames(TD)]
      Env_tax <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                      "Obs_obs","data","taxInfo.csv"), sep=",")
      colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
      length(Env_tax[,tx]%in%TD_tax[,tx])
      sum(Env_tax[,tx]%in%TD_tax[,tx])
      dim(Env)
      Env <- Env[Env_tax[,tx]%in%TD_tax[,tx],]
      Env_tax_c <- Env_tax[Env_tax[,tx]%in%TD_tax[,tx],]
      length(unique(Env_tax_c[,tx]))
      length(unique(TD_tax[,tx]))
      
      list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
      # predicted 
      TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
      TDenv
      head(TDtd)
      head(TDenv)
      # plot(TD[,1],TDtd[,1])
      #  abline(0,1)
      #  plot(TD[,1],TDenv[,1])
      #  abline(0,1)
    }
    
    trait_names=colnames(TD)
    

    #_________________________________________________
    # zlog
    #-------------------------------------------------
    t=1
    zlog=FALSE
    if(zlog){
      res <- list()
      res$TD <- TD
      res$TDtd <- TDtd
      res$TDenv <- TDenv
      res$TD_sparse <- TD_sparse
      res$Env_sparse <- Env
      t=1
      for(t in 1:length(trait_names)){
        whichcol=which(colnames(TD)==trait_names[t])
        
        TDmn=mean(log(TD[,whichcol]),na.rm = TRUE)
        TDsd=mean(log(TD[,whichcol]),na.rm = TRUE)
        
        # zlog transform on the basis of TD        
        res$TDtd[,whichcol] <- (log(TDtd[,whichcol])-TDmn)/TDsd #zlog TDtd
        res$TD[,whichcol] <- (log(TD[,whichcol])-TDmn)/TDsd #zlog TD
        res$TDenv[,whichcol] <- (log(TDenv[,whichcol])-TDmn)/TDsd #zlog TDenv
        res$TD_sparse[,whichcol] <- (log(TD_sparse[,whichcol])-TDmn)/TDsd #zlog TD_sparse
        res$Env_sparse[,whichcol] <- (log(Env[,whichcol])-TDmn)/TDsd #zlog TD_sparse
      }
      TD_zlog <- res$TD
      TDtd_zlog <- res$TDtd
      TDtax_zlog <- res$TD
      TDenv_zlog <- res$TDenv
      TD_sparse_zlog <- res$TD_sparse
      Env_sparse_zlog <- res$Env_sparse
    }
    
    #_________________________________________________
    # aggregate TD to taxon mean
    #-------------------------------------------------
    if(zlog){
    TD_tx <- aggregate(TD_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TD_sparse_tx <- aggregate(TD_sparse_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDtd_tx <- aggregate(TDtd_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDenv_tx <- aggregate(TDenv_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    Env_sparse_tx <- aggregate(Env_sparse_zlog,by=list(Env_tax_c[,tx]),FUN = mean.fun)
   
    TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
    TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
    TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
    TDenv_tx <- TDenv_tx[,-which(colnames(TDenv_tx)=="Group.1")]
    Env_sparse_tx <- Env_sparse_tx[,-which(colnames(Env_sparse_tx)=="Group.1")]

    }else{
  TD_tx <- aggregate(TD,by=list(TD_tax[,tx]),FUN = mean.fun)
  TD_sparse_tx <- aggregate(TD_sparse,by=list(TD_tax[,tx]),FUN = mean.fun)
  TDtd_tx <- aggregate(TDtd,by=list(TD_tax[,tx]),FUN = mean.fun)
  TDenv_tx <- aggregate(TDenv,by=list(TD_tax[,tx]),FUN = mean.fun)
  Env_sparse_tx <- aggregate(Env,by=list(Env_tax_c[,tx]),FUN = mean.fun)
  
  TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
  TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
  TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
  TDenv_tx <- TDenv_tx[,-which(colnames(TDenv_tx)=="Group.1")]
  Env_sparse_tx <- Env_sparse_tx[,-which(colnames(Env_sparse_tx)=="Group.1")]
}
    
    #_________________________________________________
    # taxon mean: residuals (same aggregation)
    #-------------------------------------------------
    if(zlog){
    { 
      dist_spec_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_gen_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_fam_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      dist_clad_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
      # calc distance to species mean
      t=1
      for(t in 1:ncol(TD)){
        tax <- TD_tax[,2]
        for(i in 1:nrow(TD)){
          s = TD_zlog[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_spec_mean[i,t] <- TD_zlog[i,t]-sm
        }
        tax <- TD_tax[,3]
        for(i in 1:nrow(TD)){
          s = TD_zlog[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_gen_mean[i,t] <- TD_zlog[i,t]-sm
        }
        tax <- TD_tax[,4]
        for(i in 1:nrow(TD)){
          s = TD_zlog[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_fam_mean[i,t] <- TD_zlog[i,t]-sm
        }
        tax <- TD_tax[,5]
        for(i in 1:nrow(TD)){
          s = TD_zlog[tax==tax[i],t]
          sm = mean(s,na.rm = TRUE)
          dist_clad_mean[i,t] <- TD_zlog[i,t]-sm
        }
      }
    }
    }else{
      {
        { 
          dist_spec_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_gen_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_fam_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_clad_mean <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          # calc distance to species mean
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
        }
      }
      {
        res_TDtd_TD <- TDtd - TD
        { 
          dist_spec_meantd <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_gen_meantd <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_fam_meantd <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          dist_clad_meantd <- matrix(NA,ncol=ncol(TD), nrow=nrow(TD))
          # calc distance to species mean
          t=1
          for(t in 1:ncol(TD)){
            tax <- TD_tax[,2]
            for(i in 1:nrow(TD)){
              s = TD[tax==tax[i],t]
              sm = mean(s,na.rm = TRUE)
              dist_spec_meantd[i,t] <- TDtd[i,t]-sm
            }
            tax <- TD_tax[,3]
            for(i in 1:nrow(TD)){
              s = TD[tax==tax[i],t]
              sm = mean(s,na.rm = TRUE)
              dist_gen_meantd[i,t] <- TDtd[i,t]-sm
            }
            tax <- TD_tax[,4]
            for(i in 1:nrow(TD)){
              s = TD[tax==tax[i],t]
              sm = mean(s,na.rm = TRUE)
              dist_fam_meantd[i,t] <- TDtd[i,t]-sm
            }
            tax <- TD_tax[,5]
            for(i in 1:nrow(TD)){
              s = TD[tax==tax[i],t]
              sm = mean(s,na.rm = TRUE)
              dist_clad_meantd[i,t] <- TDtd[i,t]-sm
            }
          }
        }
      }
    }
    par(mfrow=c(1,1))  
    
    #_________________________________________________
    # aggregate TD to taxon mean
    #-------------------------------------------------
    dist_spec_mean_tx <- aggregate(abs(dist_spec_mean),by=list(TD_tax[,2]),FUN = mean.fun)
    dist_spec_mean_tx <- dist_spec_mean_tx[,-which(colnames(dist_spec_mean_tx)=="Group.1")]
    
    dist_spec_meantd_tx <- aggregate(abs(dist_spec_meantd),by=list(TD_tax[,2]),FUN = mean.fun)
    dist_spec_meantd_tx <- dist_spec_meantd_tx[,-which(colnames(dist_spec_meantd_tx)=="Group.1")]
    
    resid <- abs(cbind(TD_zlog-TDtd_zlog))
    resid_tx <- abs(cbind(TD_tx-TDtd_tx))
    dim(resid)
    t=1
    
    
    (TDtd-TD)
    distTDtd_TD <- aggregate(abs(TDtd-TD),by=list(TD_tax[,2]),FUN = mean.fun)
    distTDtd_TD <- distTDtd_TD[,-which(colnames(distTDtd_TD)=="Group.1")]
    #plot(distTDtd_TD[,t],resid[,t])
    
    dim(dist_spec_mean)
    dim(dist_spec_meantd)
    
    plot(abs(dist_spec_meantd_tx[,t]-dist_spec_mean_tx[,t]),resid_tx[,t])
    
    print(cor(unlist(abs()),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t])))
    plot(unlist(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t])),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t]))
    print(cor(unlist(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t])),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t])))
    plot(unlist(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t])),unlist(resid_tx[dist_spec_mean_tx[,t]!=0,t]))
    abline(0,1)
    
    for(t in 1:ncol(TD)){
#      print(cor(abs(dist_spec_mean[dist_spec_mean[,t]!=0,t]),resid[dist_spec_mean[,t]!=0,t]))
#      plot(abs(dist_spec_mean[dist_spec_mean[,t]!=0,t]),resid[dist_spec_mean[,t]!=0,t])
      print(cor(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t]),resid_tx[dist_spec_mean_tx[,t]!=0,t]))
      plot(abs(dist_spec_mean_tx[dist_spec_mean_tx[,t]!=0,t]),resid_tx[dist_spec_mean_tx[,t]!=0,t])
      abline(0,1)
    }
    
}

