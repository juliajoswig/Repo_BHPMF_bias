
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
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  

  #------------------------------------------------------------
  # chose data set
  #------------------------------------------------------------
  repnums=3
  t_choice="data"
  GapPercent=80
  if(t_choice=="data"){RepNum=1}
  if(t_choice=="data_2"){RepNum=2}
  
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
  if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
  if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","mm2")}
  
colz=colz1

{
  # load Envelope data
  # load TDenvelope
  ObsOrTD="Obs_obs_TD"
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))

  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  head(TD)
  # load Envelope data
  # total trait data 
  EnvTot <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                          "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  summary(EnvTot)
  Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
  head(Env)
  
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
}

{
  # load Envelope data
  # load TDenvelope
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  taxTD <- as.data.frame(read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                     "Obs_obs_TD","data","taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_TD <- taxTD[,1]
  head(taxTD)
  dim(taxTD)
  dim(TD)
  dim(TD_sparse)
  # load Envelope data
  # total trait data 
  EnvTot <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                             "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  Envsparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                                "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  taxEnv <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs","data","taxInfo.csv"),header=TRUE))
  taxEnv <- as.data.frame(read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                      "Obs_obs","data","taxInfo.csv"),
                                     sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_env <- taxEnv[,1]
  dim(taxEnv)
  summary(EnvTot)
  Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
  Env_sp <- Envsparse[,colnames(Envsparse)%in%colnames(TD)]
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],Env[match(taxTD$ID,ID_env),1])
  abline(0,1)
  plot(TD[,1],Env_sp[match(taxTD$ID,ID_env),1])
  plot(TD[,1],TDenv[,1])
  abline(0,1)
  plot(TD_sparse[,1],TD[,1])
  abline(0,1)
}

colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
trait_names=colnames(TD)

defense_prep=TRUE
par(mfrow=c(1,1),mar=c(5,6,1,1))
if(t_choice=="data"){pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_Density.pdf"),width=18,height=4)}
if(t_choice=="data_2"){pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_S2_Density.pdf"),width=18,height=4)}
if(t_choice=="data"&defense_prep==TRUE){pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_Density_only2.pdf"),width=18,height=4)}
if(t_choice=="data_2"&defense_prep==TRUE){pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_S2_Density_only2.pdf"),width=18,height=4)}

{
  if(t_choice=="data"){
  limmin <- c(0,-1,-.01,0,0,0)
  limmax <- c(60,60,1,40,5,7)
  li <- c(0,0,0,0,0,0)
  lx <- c(.12,.16,3,.1,1.1,.8)
  lx <- c(.14,.53,3.1,.1,1.2,.85)
  par(mfrow=c(1,6),mar=c(5,6,3,1))
  col_below=colz[5]#"gray"
  }
  if(t_choice=="data_2"){
    limmin <- c(0,0,0,0,-100)
    limmax <- c(60,3,8,1,5000)
    li <- c(0,0,0,0,0,0)
    lx <- c(.14,1.5,40,5,.0015)
    
    par(mfrow=c(1,1),mar=c(5,6,3,1))
    par(mfrow=c(1,5),mar=c(5,6,3,1))
    col_below=colz[5]#"gray"
  }
  if(t_choice=="data"&defense_prep){
    limmin <- c(0,-1,-.01,0,0,0)
    limmax <- c(60,60,1,40,5,7)
    li <- c(0,0,0,0,0,0)
    lx <- c(.12,.16,3,.1,1.1,.8)
    lx <- c(.14,.25,3.1,.1,1.2,.85)
    par(mfrow=c(1,6),mar=c(5,6,3,1))
    col_below=colz[5]#"gray"
  }
  t=1
  for(t in 1:length(trait_names)){
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1.5,bw = "SJ")
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1,bw = "SJ")
    plot(de,col=col_below,lwd=5,cex.main=3,cex.lab=2.5,
         cex.axis=2,main = trait_names[t],xlab=units[t],ylim=c(li[t],lx[t]),
         xlim=c(limmin[t],limmax[t]),lty=1)
    if((t!=3&t_choice=="data_2")|t_choice=="data"){polygon(de, col=col_below, border=colz[6],lwd=4)}
    if(!defense_prep){lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col="white",lwd=4)}
    if(!defense_prep){lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col=colz[1],lwd=4)}
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col="white",lwd=4)
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col=colz[2],lwd=4)
    if(!defense_prep){lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col="white",lwd=4)}
    if(!defense_prep){lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col=colz[7],lwd=4)}
    #mean(TDenv[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    #median(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    #abline(v=mean(TDenv[,colnames(TDenv)==trait_names[t]]),col=colz[1],lty=1,lwd=4)  
    #abline(v=mean(TDtd[,colnames(TDtd)==trait_names[t]]),col=colz[2],lty=1,lwd=4)  
    #abline(v=mean(TD[,colnames(TD)==trait_names[t]]),col=colz[6],lty=1,lwd=4)  
    #abline(v=mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE),col=colz[7],lty=1,lwd=4)  
    
    if(!defense_prep){
      if(t==1){
      legend(11,.15,legend = c("TD","ExTD","TD_td","TD_ext"),col = colz[c(6,7,2,1)],
             lwd=4,lty=1,cex=1.9,bty = "n")
    }
    }
    if(defense_prep){
        legend("topright",legend = c("TD","TDtd"),col = colz[c(6,2)],
               lwd=4,lty=1,cex=1.9,bty = "n")
      }
    }
    
  }
  

dev.off()



