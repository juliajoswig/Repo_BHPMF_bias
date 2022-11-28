

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
t_choices <- out$t_choices
TDnos = out$TDnos
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
gappercents= c("1","5","10","20","30","40","50","60")
repnums=1:3

GapPercent=50
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load processed data   
#-------------------------------------------------------------------
  colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
  colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255),"#b2182b",rgb(178/255,24/255,43/255,alpha = .8),"#b97880","#b2182b")
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  
  
  #-------------------------------------------------------------------
  # load data
  #-------------------------------------------------------------------
  RepNum=1
  t_choice="data"
  # TD
  ObsOrTD="Obs_obs_TD"
  path_TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs.csv")
  path_TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfoTD_obs_REzlog.csv")
  path_TDtax <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","taxInfo.csv")
  path_TDfun <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","funInfo.csv")
  TD <- as.matrix(read.table(file = path_TD, sep=",", dec=".",header=TRUE))[,-c(1:2)]
  TDtax <- as.matrix(read.table(file = path_TDtax, sep=",", dec=".",col.names = c("ObservationID","Species","Genus","Family","Clade")))
  TDfun <- as.matrix(read.table(file = path_TDfun, sep=",", dec=".",header = TRUE))
  TDtax <- cbind(TDtax,TDfun[,c(3,4)])
  dim(TDtax)
  head(TDtax)
  # TDtd
  path_TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_80",ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")
  TDtd <- as.matrix(read.table(file = path_TD, sep=",", dec=".",header=TRUE))[,-c(1:2)]
  dev.off()
  head(TDtd)
  plot(TD[,1],TDtd[,1])
  # TDenv
  ObsOrTD="Obs_obs"
  path_TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_80",ObsOrTD,"data","traitInfoTD_pred_REzlog.csv")
  TDenv <- as.matrix(read.table(file = path_TD, sep=",", dec=".",header=TRUE))[,-c(1:2)]
  head(TDtd)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
  # Env
  # ObsOrTD="Obs_obs"
  # path_Env <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfo.csv")
  # path_Envtax <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_80",ObsOrTD,"data","taxInfo.csv")
  # path_PFT <- file.path(origin,"_2021","data","helper_data","functional_information","PFT_Obs.csv")
  # path_GF <- file.path(origin,"_2021","data","helper_data","functional_information","GF_Obs.csv")
  # Env <- as.matrix(read.table(file = path_Env, sep=",", dec=".",header=TRUE))[,-c(1:2)]
  # Envtax <- as.data.frame(read.table(file = path_Envtax, sep=",", dec=".",col.names = c("ObservationID","Species","Genus","Family","Clade")))
  # head(PFT)
  # head(Env)
  # GF <- as.data.frame(read.table(file = path_GF, sep=",", dec=".",header = TRUE))
  # PFT <- as.data.frame(read.table(file = path_PFT, sep=",", dec=".",header = TRUE))
  # EnvFun1 <- cbind(PFT$ObservationID,PFT$PFT,GF$GF)
  # colnames(EnvFun1) <- c("ObservationID","PFT","GF")
  # Envfun <- merge(x = Envtax,y=EnvFun1, by="ObservationID",all.x = TRUE)
  # head(EnvFun1)
  # dim(Envfun)
  # dim(Envtax)
  # dim(Env)
  # Env <- Env[,which(colnames(Env)%in%colnames(TD))]
#  Env_sp <- Env[which(Envtax[,2]%in%TDtax[,2]),]
#  Env_gen <- Env[which(Envtax[,3]%in%TDtax[,3]),]
#  Env_fam <- Env[which(Envtax[,4]%in%TDtax[,4]),]
#  Env_clad <- Env[which(Envtax[,5]%in%TDtax[,5]),]
#  Env_GF <- Env[which(Envfun[,6]%in%TDtax[,6]),]
#  Env_PFT <- Env[which(Envfun[,7]%in%TDtax[,7]),]
  
  
  #-----------------------------------------------------------------------------
  # Silhouette index
  #-----------------------------------------------------------------------------
  library(cluster) 
  sil_TD <- list()
  sil_TDtd <- list()
  sil_TDenv <- list()
  sil_Env <- list()
  
  cl=3
 for(cl in 2:7){
   # TD TDtd and TDenv
  cats <- factor(TDtax[,cl])
  ranks <- rank(-table(cats), ties.method="first")
  DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
  taxnow <- DF[,2]
  ixG <- !is.na(taxnow)
  si2 <- silhouette(taxnow[ixG], dist(TD[ixG,],method = "canberra"))
  si=summary(si2)
  sil_TD[[(cl-1)]] <- cbind(as.vector(unique(TDtax[ixG,cl])),si$clus.avg.widths,si$clus.sizes)
  si2 <- silhouette(taxnow[ixG], dist(TDtd[ixG,],method = "canberra"))
  si=summary(si2)
  sil_TDtd[[(cl-1)]] <- cbind(as.vector(unique(TDtax[ixG,cl])),si$clus.avg.widths,si$clus.sizes)
  si2 <- silhouette(taxnow[ixG], dist(TDenv[ixG,],method = "canberra"))
  si=summary(si2)
  sil_TDenv[[(cl-1)]] <- cbind(as.vector(unique(TDtax[ixG,cl])),si$clus.avg.widths,si$clus.sizes)
  
  # Env
  # if(cl==2){cats <- factor(Envtax[which(Envtax[,2]%in%TDtax[,2]),cl]);tax_Env=Envtax[which(Envtax[,cl]%in%TDtax[,cl]),cl]}
  # if(cl==3){cats <- factor(Envtax[which(Envtax[,3]%in%TDtax[,3]),cl]);tax_Env=Envtax[which(Envtax[,cl]%in%TDtax[,cl]),cl]}
  # if(cl==4){cats <- factor(Envtax[which(Envtax[,4]%in%TDtax[,4]),cl]);tax_Env=Envtax[which(Envtax[,cl]%in%TDtax[,cl]),cl]}
  # if(cl==5){cats <- factor(Envtax[which(Envtax[,5]%in%TDtax[,5]),cl]);tax_Env=Envtax[which(Envtax[,cl]%in%TDtax[,cl]),cl]}
  # if(cl==6){cats <- factor(Envfun[which(Envfun[,6]%in%TDtax[,6]),cl]);tax_Env=Envfun[which(Envfun[,cl]%in%TDtax[,cl]),cl]}
  # if(cl==7){cats <- factor(Envfun[which(Envfun[,7]%in%TDtax[,7]),cl]);tax_Env=Envfun[which(Envfun[,cl]%in%TDtax[,cl]),cl]}
  # ranks <- rank(-table(cats), ties.method="first")
  # DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
  # taxnow <- DF[,2]
  # ixG <- !is.na(taxnow)
  # if(cl==2){dat_now <- Env[which(Envtax[,cl]%in%TDtax[,cl]),]}
  # if(cl==3){dat_now <- Env[which(Envtax[,cl]%in%TDtax[,cl]),]}
  # if(cl==4){dat_now <- Env[which(Envtax[,cl]%in%TDtax[,cl]),]}
  # if(cl==5){dat_now <- Env[which(Envtax[,cl]%in%TDtax[,cl]),]}
  # if(cl==6){dat_now <- Env[which(Envfun[,cl]%in%TDtax[,cl]),]}
  # if(cl==7){dat_now <- Env[which(Envfun[,cl]%in%TDtax[,cl]),]}
  # 
  # try(rm("si2"))
  # print(dim(dat_now))
  # 
  # try(si2 <- silhouette(taxnow[ixG], dist(dat_now[ixG,],method = "canberra")))
  # if(exists("si2")){
  # si=summary(si2)
  # sil_Env[[(cl-1)]] <- cbind(as.vector(unique(tax_Env[ixG])),si$clus.avg.widths,si$clus.sizes)
  # }else{sil_Env[[(cl-1)]] <- matrix(NA,ncol=3,nrow=nrow(tax_Env))}
}
  str(sil_TD)
    
  #-----------------------------------------------------------------------------
  # plot it now
  #-----------------------------------------------------------------------------
  
  pdf(file=file.path(origin,"_2021","figures","Figure_5",paste0("Figure_5_Silouette.pdf")),width=20,height=10)
  
  par(mar=c(17,0,4,0),mfrow=c(1,7))
  {
    ymax=1
    ymin=-1
    plot(1:10,col="white",ylim=c(0,ymax),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    axis(side = 2,line = -20,cex.axis=2)
    axis(side = 2,at = .5,line = -14,labels = "Silhouette index",tick = FALSE,cex.axis=7,cex.lab=5)
    
    names_now=c("Species","Genus","Family","Clades","GF","PFT")
    t=1
    for(t in 1:6){
      bxplTD=as.numeric(as.vector(sil_TD[[t]][,2]))
      bxplTDtd=as.numeric(as.vector(sil_TDtd[[t]][,2]))
      bxplTDenv=as.numeric(as.vector(sil_TDenv[[t]][,2]))
      
      dat_now=cbind(bxplTD,bxplTDtd,bxplTDenv)
      dat_now[dat_now==0]=NA
      boxplot(dat_now,ylim=c(ymin,ymax),ylab="Silhouette index",col=colz_solid,main=names_now[t],cex.main=5,
              las=2,cex.axis=4.5,cex.lab=4.5,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
      axis(1,at = 1:ncol(dat_now),labels = c("TD","TD_td","TD_ext"),las=2,cex.axis=4.5,tick=FALSE)
      yBot_now=median(dat_now[,1],na.rm = TRUE)
      i=1
      rect(xleft = .5,ybottom = ymin,xright = 4.5,ytop = yBot_now[i],col=colz[6],border = NA)
      rect(xleft = .5,ybottom = yBot_now,xright = 4.5,ytop = ymax,col=colz[10],border = NA) 
      boxplot(dat_now,ylim=c(ymin,ymax),ylab="Silhouette index",col=colz_solid,main=names_now[t],cex.main=5,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n",add=TRUE)
      
    }
    
  }
  
  dev.off()
  
  
  
  
  
  
  
  
  
  
  
  
  
  