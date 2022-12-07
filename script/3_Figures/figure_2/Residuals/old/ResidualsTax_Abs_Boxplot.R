# atm S7: Relative error to trait range for test data 2. Error (absolute) calculated
# as imputed - observed for test data imputed with and without the envelope (back-transformed data 


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
  GapPercent=50
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3


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
tx=0 #(0=indiv,1=spec,2=genus,3=family,4=clades)
#for(tx in 0:4){
  tx=tx+1
  par(mfrow=c(4,4))
  {
    
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
    file.path(originData,"analyes","Point_wise","res.csv")
    colz=colz1
    
    #load data

    
    {
      # load Envelope data
      # load TDenvelope
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
      # load TD data
      # total trait data 
      TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
      TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                    "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
      TD_tax <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         "Obs_obs_TD","data","taxInfo.csv"),
                                        sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
      ID_TD <- TD_tax[,1]
      head(TD_tax)
      dim(TD_tax)
      dim(TD)
      dim(TD_sparse)
      # load Envelope data
      # total trait data 
      EnvTot <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                                 "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
      Envsparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                                    "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
      taxEnv <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                                 "Obs_obs","data","taxInfo.csv"),header=TRUE))
      taxEnv <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                          "Obs_obs","data","taxInfo.csv"),
                                         sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
      ID_env <- taxEnv[,1]
      dim(taxEnv)
      #cut down to columns
      Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
      Env_sp <- Envsparse[,colnames(Envsparse)%in%colnames(TD)]
      #cut down to taxa necessary
      Env_tax_c <- taxEnv[taxEnv[,tx]%in%TD_tax[,tx],]
      Env <- Env[taxEnv[,tx]%in%TD_tax[,tx],]
      Env_sp <- Env_sp[taxEnv[,tx]%in%TD_tax[,tx],]
      
      
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      # predicted 
      TDtd <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                           "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
      TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
      head(TDtd)
      head(TDenv)
#      plot(TD[,1],TDtd[,1])
#      abline(0,1)
#      plot(TD_sparse[,1],Env[match(TD_tax$ID,taxEnv_c$ObservationID),1])
#      abline(0,1)
#      plot(TD[,1],Env_sp[match(TD_tax$ID,taxEnv_c$ObservationID),1])
#      plot(TD[,1],TDenv[,1])
#      abline(0,1)
#      plot(TD_sparse[,1],TD[,1])
#      abline(0,1)
    }
    
    trait_names=colnames(TD)
    
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
    
    #_________________________________________________
    # aggregate TD to taxon mean
    #-------------------------------------------------
    TD_tx <- aggregate(TD,by=list(TD_tax[,tx]),FUN = mean.fun)
    TD_sparse_tx <- aggregate(TD_sparse,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDtd_tx <- aggregate(TDtd,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDenv_tx <- aggregate(TDenv,by=list(TD_tax[,tx]),FUN = mean.fun)
    Env_sparse_tx <- aggregate(Env_sp,by=list(Env_tax_c[,tx]),FUN = mean.fun)
    
    TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
    TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
    TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
    TDenv_tx <- TDenv_tx[,-which(colnames(TDenv_tx)=="Group.1")]
    Env_sparse_tx <- Env_sparse_tx[,-which(colnames(Env_sparse_tx)=="Group.1")]
    
    
    #_________________________________________________
    # zlog based on TD
    #-------------------------------------------------
    {
      res <- list()
      res$TDtd_tx <- TDtd_tx
      res$TD_tx <- TD_tx
      res$TDenv_tx <- TDenv_tx
      res$TD_sparse_tx <- TD_sparse_tx
      res$Env_sparse_tx <- Env_sparse_tx
      for(t in 1:length(trait_names)){
        whichcol=which(colnames(TD_tx)==trait_names[t])
        
        TDmn=mean(log(TD[,whichcol]),na.rm = TRUE)
        TDsd=mean(log(TD[,whichcol]),na.rm = TRUE)
        
        # zlog transform on the basis of TD        
        res$TD_tx[,whichcol] <-         (log(TD_tx[,whichcol])        -TDmn)/TDsd #zlog TD
        res$TDtd_tx[,whichcol] <-       (log(TDtd_tx[,whichcol])      -TDmn)/TDsd #zlog TDtd
        res$TDenv_tx[,whichcol] <-      (log(TDenv_tx[,whichcol])     -TDmn)/TDsd #zlog TDenv
        res$TD_sparse_tx[,whichcol] <-  (log(TD_sparse_tx[,whichcol]) -TDmn)/TDsd #zlog TD_sparse
        res$Env_sparse_tx[,whichcol] <- (log(Env_sparse_tx[,whichcol])-TDmn)/TDsd #zlog TD_sparse
      }
}
      TDtd_now <- res$TDtd_tx
      TD_now <- res$TD_tx
      TDenv_now <- res$TDenv_tx
      TD_sparse_now <- res$TD_sparse_tx
      Env_sparse_now <- res$Env_sparse_tx

#      plot(TD_sparse_tx[,1],Env_sparse_tx[match(TD_tax$ID,Env_tax_c$ID),1])
      abline(0,1)
      
      # residuals of TDtd to TD: 
      bxplTDtd_TD= TDtd_now - TD_now
      #        plot(TDtd_now,TD_now);abline(0,1)
      # residuals of TDtd to TDsparse: 
      bxplTDtd_TDsparse=TDtd_now- TD_sparse_now
      #        plot(TDtd_now,TD_sparse_now);abline(0,1)
      # residuals of TDtd to Envsparse: 
      bxplTDtd_Envsparse=TDtd_now- Env_sparse_now[TD_tax$ID%in%Env_tax_c$ID,]
      #        plot(TDtd_now,Env_sparse_now);abline(0,1)
      # residuals of TDenv to Envsparse: 
      bxplTDenv_Envsparse=TDenv_now- Env_sparse_now[TD_tax$ID%in%Env_tax_c$ID,]
      #        plot(TDenv_now,Env_sparse_now);abline(0,1)
      # residuals of TDenv to TD: 
      bxplTDenv_TD = TDenv_now - TD_now
      #        plot(TD_now,TDenv_now);abline(0,1)
    

    

   # par(mfrow=c(1,1),mar=c(12,4,1,1))
    pdf(file=file.path(origin,"figures","Figure_2",paste0(txs[tx],"_boxplots_",t_choice,".pdf")),width=8,height=2.3,pointsize = 9)
    par(mfrow=c(1,6),mar=c(12,4,1,1))
    t=1

    for(t in 1:ncol(TD)){
      head(bxplTDenv_TD)
      
      # (bxplTDenv_Envsparse),  (bxplTDenv_TD))[seq(from = t,to = ncol(TD)*5,by = ncol(TD))]
      dat_plot=cbind((bxplTDtd_TD),(bxplTDtd_TDsparse),(bxplTDtd_Envsparse),(bxplTDenv_Envsparse),(bxplTDenv_TD))[seq(from = t,to = ncol(TD)*5,by = ncol(TD))]
      colnames(dat_plot) =  c("IMP_obs - OBS","IMP_obs - OBS_sparse", 
                             "IMP_obs - TRY17_sparse", "IMP_obsExt - TRY17_sparse",  "IMP_obsExt - OBS")
      dat_plot <- abs(dat_plot)
      
      boxplot(dat_plot,main=trait_names[t],ylim=c(0,
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.95,na.rm=TRUE))),col=colz,las=2,ylab=paste("abs(zlog) distance"))
      rect(xleft = 0,ybottom = min(apply(dat_plot,MARGIN = 2,quantile,prob=0,na.rm=TRUE))-1,xright = 6,ytop = median(dat_plot[,1],na.rm = TRUE),col=colz2[1])
      rect(xleft = 0,ybottom = median(dat_plot[,1],na.rm = TRUE),xright = 6,ytop = max(apply(dat_plot,MARGIN = 2,quantile,prob=1,na.rm=TRUE))+1,col=colz2[3])
      abline(h = median(dat_plot[,1],na.rm = TRUE),col=colz2[8],lwd=1)
      abline(v = 3.5,col="white",lwd=1.5)
      boxplot(dat_plot,main=trait_names[t],ylim=c(min(apply(dat_plot,MARGIN = 2,quantile,prob=.01,na.rm=TRUE)),
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.99,na.rm=TRUE))),
              col=colz,las=2,ylab=paste("zlog(",units[t],")"),
              add=TRUE)
    }   
    dev.off()
    
  }
  
}


dev.off()



{
  # load Envelope data
  # load TDenvelope
  ObsOrTD="Obs_obs"
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
  TDtd <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                           "Obs_obs_TD","data","traitInfoTD_pred.csv"),header=TRUE))[,-c(1,2)]
  TD_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                 "Obs_obs_TD","data","taxInfo.csv"), sep=",")
  colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  head(TD_tax)
  colnames(TD_tax)
  dim(TD_tax)
  # load Envelope data
  # total trait data 
  Env <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                          "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  Env <- Env[,colnames(Env)%in%colnames(TD)]
  sum(colnames(Env)==colnames(TD))==ncol(Env)
  
  Env_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                  "Obs_obs","data","taxInfo.csv"), sep=",")
  colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  dim(Env_tax)==dim(Env)
  Env <- Env[Env_tax[,tx]%in%TD_tax[,tx],]
  Env_tax_c <- Env_tax[Env_tax[,tx]%in%TD_tax[,tx],]
  
  mtch <- match(TD_tax[,tx],Env_tax_c[,tx])
  print(sum(TD_tax[,tx]==Env_tax_c[mtch,tx])==length(TD_tax[,tx]))
  Env_tax_c<- Env_tax_c[mtch,]
  Env<- Env[mtch,]
  
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  # predicted 
  TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv
  head(TDtd)
  head(TDenv)
  head(TD_sparse)
  head(Env)
  par(mfrow=c(1,1))
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],Env[,1])      
  abline(0,1)
}