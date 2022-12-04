
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
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  
  
  #-------------------------------------------------------------------
  # chose trait data   
  #-------------------------------------------------------------------
  repnums=3
  t_choice="data"
  GapPercent=80
  if(t_choice=="data"){RepNum=1}
  if(t_choice=="data_2"){RepNum=2}
  
  
  
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
units <- c("mm2 mg-1","m","","","mm2")
colz=colz1

#load data

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
#colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
#colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))
trait_names=colnames(TD)


par(mfrow=c(1,1),mar=c(5,6,1,1))
pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_S2_Density.pdf"),width=18,height=4)

{
  limmin <- c(0,0,0,0,-100)
  limmax <- c(60,3,8,1,5000)
  li <- c(0,0,0,0,0,0)
  lx <- c(.14,1.5,40,5,.0015)

  par(mfrow=c(1,1),mar=c(5,6,3,1))
  par(mfrow=c(1,5),mar=c(5,6,3,1))
  col_below=colz[5]#"gray"
  t=3
  for(t in 1:length(trait_names)){
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1.5,bw = "SJ")
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1,bw = "SJ")
   # plot(de)
    mean(TD[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    mean(TDtd[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    mean(TDenv[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
    
    plot(de,col=col_below,lwd=5,cex.main=3,cex.lab=2.5,
         cex.axis=2,main = trait_names[t],xlab=units[t],ylim=c(li[t],lx[t]),
         xlim=c(limmin[t],limmax[t]),lty=1)
    if(t!=3){polygon(de, col=col_below, border=colz[6],lwd=4)}
    lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col="white",lwd=4)
    lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col=colz[1],lwd=4)
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col="white",lwd=4)
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col=colz[2],lwd=4)
    lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col="white",lwd=4)
    lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col=colz[7],lwd=4)
    #summary(Env[,colnames(Env)==trait_names[t]])
    #de <- density(c(Env[,colnames(Env)==trait_names[t]],a[[t]]),na.rm = T,adjust = 1,bw = "SJ")
    #plot(de)
    #summary(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw) 
    #abline(v=mean(TDenv[,colnames(TDenv)==trait_names[t]]),col=colz[1],lty=1,lwd=4)  
    #abline(v=mean(TDtd[,colnames(TDtd)==trait_names[t]]),col=colz[2],lty=1,lwd=4)  
    #abline(v=mean(TD[,colnames(TD)==trait_names[t]]),col=colz[6],lty=1,lwd=4)  
    #abline(v=mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE),col=colz[7],lty=1,lwd=4)  
    
    if(t==1){
      legend(11,.15,legend = c("TD","Env","TD_td","TD_env"),col = colz[c(6,7,2,1)],
             lwd=4,lty=1,cex=1.9,bty = "n")
    }
    
    
  }
  
}
dev.off()





par(mfrow=c(1,1),mar=c(5,6,1,1))
pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_S2_Density.pdf"),width=18,height=4)

{
  limmin <- c(0,-1,-.01,0,0,0)
  limmax <- c(60,60,1,40,5,7)
  li <- c(0,0,0,0,0,0)
  lx <- c(.12,.16,3,.1,1.1,.8)
  lx <- c(.14,.53,3.1,.1,1.2,.85)
  
  par(mfrow=c(1,6),mar=c(5,6,3,1))
  col_below=colz[5]#"gray"
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl <- rep(NA,sum(ix_trait))
  t=2
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    sum(ix_trait,na.rm = T)
    ix_trait[is.na(ix_trait)] <- FALSE
    
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1.5,bw = "SJ")
    de <- density(TD[,colnames(TD)==trait_names[t]],na.rm = T,adjust = 1,bw = "SJ")
    plot(de,col=col_below,lwd=5,cex.main=3,cex.lab=2.5,
         cex.axis=2,main = trait_names[t],xlab=units[t],ylim=c(li[t],lx[t]),
         xlim=c(limmin[t],limmax[t]),lty=1)
    polygon(de, col=col_below, border=col_below) 
    lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col="white",lwd=4)
    lines(density(TDenv[,colnames(TD)==trait_names[t]],bw = de$bw),col=colz[1],lwd=4)
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col="white",lwd=4)
    lines(density(TDtd[,colnames(TDtd)==trait_names[t]],bw = de$bw),col=colz[2],lwd=4)
    lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col="white",lwd=4)
    lines(density(Env[,colnames(Env)==trait_names[t]],na.rm=TRUE,bw = de$bw),col=colz[7],lwd=4)
    abline(v=mean(TDenv[,colnames(TDenv)==trait_names[t]]),col=colz[1],lty=2,lwd=2)  
    abline(v=mean(TDtd[,colnames(TDtd)==trait_names[t]]),col=colz_solid[2],lty=1,lwd=2)  
    abline(v=mean(TD[,colnames(TD)==trait_names[t]]),col="white",lty=1,lwd=3)  
    abline(v=mean(TD[,colnames(TD)==trait_names[t]]),col=col_below,lty=1,lwd=2)  
    abline(v=mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE),col=colz[7],lty=2,lwd=2)  
  }
}
dev.off()
