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
tx=tx+1
print(tx)

  par(mfrow=c(4,4))
    
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
    file.path(originData,"analyes","Point_wise","res.csv")
    colz=colz1
    
    # load Envelope data
      # load TDenvelope
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
      
      # load TD data directly zlog according to OBS
      # total trait data 
      TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfoTD_obs_zlog.csv"),header=TRUE))[,-c(1,2)]
      TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                    "Obs_obs_TD","data","traitInfoTD_obs_zlog.csv"),header=TRUE))[,-c(1,2)]
      TD_tax <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         "Obs_obs_TD","data","taxInfo.csv"),
                                        sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
      ID_TD <- TD_tax[,1]
      head(TD_tax)
      dim(TD_tax)
      dim(TD)
      dim(TD_sparse)
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),"Obs_obs","data"))
      EnvTot <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                                 "Obs_obs","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-1]
      Envsparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                                    "Obs_obs","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-1]
      taxEnv <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                                 "Obs_obs","data","taxInfoTD.csv"),header=TRUE))
      ID_env <- taxEnv[,1]


      # predicted 
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      TDtd <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                           "Obs_obs","data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
      TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs_TD","data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
      dim(TDtd)
      dim(TDenv)

    trait_names=colnames(TD)
    
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
    
    
    TDtd_now <- TDtd_tx
    TD_now <- TD_tx
    TDenv_now <- TDenv_tx
    TD_sparse_now <- TD_sparse_tx
    Env_sparse_now <- Env_sparse_tx

    try(dev.off())
    plot(TD_sparse_tx[,1],Env_sparse_tx[,1])
    abline(0,1)
      
    # residuals of TDtd to TD: 
    bxplTDtd_TD= TDtd_now - TD_now
    #        plot(TDtd_now,TD_now);abline(0,1)
    # residuals of TDtd to TDsparse: 
    bxplTDtd_TDsparse=TDtd_now- TD_sparse_now
      #        plot(TDtd_now,TD_sparse_now);abline(0,1)
      # residuals of TDtd to Envsparse: 
      bxplTDtd_Envsparse=TDtd_now- Env_sparse_now
      #        plot(TDtd_now,Env_sparse_now);abline(0,1)
      # residuals of TDenv to Envsparse: 
      bxplTDenv_Envsparse=TDenv_now- Env_sparse_now
      #        plot(TDenv_now,Env_sparse_now);abline(0,1)
      # residuals of TDenv to TD: 
      bxplTDenv_TD = TDenv_now - TD_now
      #        plot(TD_now,TDenv_now);abline(0,1)
    dev.off()

    

   # par(mfrow=c(1,1),mar=c(12,4,1,1))
    pdf(file=file.path(origin,"figures","Residuals",paste0(txs[tx],"_boxplots_",t_choice,".pdf")),width=8,height=2.3,pointsize = 9)
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
  

