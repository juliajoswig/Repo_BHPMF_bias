
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

colz2=c("#d53e4f","#f46d43","#fdae61",
        "#fee08b","#ffffbf","#e6f598","#abdda4",
        "#66c2a5","#3288bd")
#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
RepNum=2
t_choice="data_2"


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
  ObsOrTD="Obs_obs"
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfoTD_obs_REzlog.csv"),header=TRUE))[,-c(1,2)]
  TDtd <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                         "Obs_obs_TD","data","traitInfoTD_pred_REzlog.csv"),header=TRUE))[,-c(1,2)]
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
  Env_tax <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                 "Obs_obs","data","taxInfo.csv"), sep=",")
  colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  head(Env_tax)
  
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  # predicted 
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs","data","traitInfoTD_pred_REzlog.csv")))[,-c(1,2)]
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
}

trait_names=colnames(TD)

count.fun <- function(input){
  return(sum(!is.na(input)))
}
# aggregate TD to species mean
TD_sp <- aggregate(TD,by=list(TD_tax$Species),FUN = mean)

TDtd_sp <- aggregate(TDtd,by=list(TD_tax$Species),FUN = mean)
TDenv_sp <- aggregate(TDenv,by=list(TD_tax$Species),FUN = mean)
TDenv_sp <- TDenv_sp[TDenv_sp$Group.1%in%TD_sp$Group.1,]

#taxonomic information
TD_spNB <- aggregate(TD,by=list(TD_tax$Species),FUN = count.fun)
TDenv_spNB <- aggregate(Env,by=list(Env_tax$Species),FUN = count.fun)
TDenv_spNB <- TDenv_spNB[TDenv_spNB$Group.1%in%TD_sp$Group.1,]
dim(TD_sp)
dim(TDenv_sp)
dim(TDtd_sp)

#create color according to sample number
x <- 1:40
dat <- data.frame(x = x,y = x^2 + 1)
#rbPal <- colorRampPalette(c("#ef8a62","#67a9cf"))
rbPal <- colorRampPalette(c(rgb(146,197,222,maxColorValue = 255,alpha=255/2),
                            rgb(247,247,247,maxColorValue = 255,alpha=255/2),
                            rgb(202,0,32,maxColorValue = 255,alpha=255/2)))
dat$Col <- rbPal((max(x)-1))[as.numeric(cut(dat$y,breaks = (max(x)-1)))]
plot(dat$x,dat$y,pch = 20,col = dat$Col)
colpal=dat$Col



  
par(mfrow=c(1,1),mar=c(5,6,1,1))
pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_xy_species.pdf"),width=10,height=5)

{
  par(mfrow=c(1,2),mar=c(5,6,3,1))
  # TDtd
  nbs=TD_spNB$SLA
  colz=rep("black",length(nbs))
  for(i in 1:length(colpal)){
    colz[nbs==i] <- colpal[i]
  }
  plot_data=c(NA,NA)
  t=2
  plot_data <- rbind(plot_data,cbind(TD_sp[,t],TDtd_sp[,t]))
  plot(TD_sp[,t],TDtd_sp[,t],col=colz,pch=16,xlim=c(-3.5,2.5),ylim=c(-3,5),
       xlab="TD, species mean",cex=.3,
       ylab="TDtd, species mean")
  head(plot_data)
  for(t in 3:ncol(TD_sp)){
    plot_data <- rbind(plot_data,cbind(TD_sp[,t],TDtd_sp[,t]))
    points(TD_sp[,t],TDtd_sp[,t],col=colz,pch=16,
         xlab="TD, species mean",cex=.6,
         ylab="TDtd, species mean")
  }
  colnames(plot_data) = c("obs","pred")
  dim(plot_data)
  plot_data <- as.data.frame(plot_data)
  text(x = -1,y = 4, labels = paste0("Pearson cor=",round(cor(plot_data[complete.cases(plot_data),])[2,1],digits=2)))
  abline(0,1)

  
  
  # TDenv
  plot_data=c(NA,NA)

  t=2
  nbs=TDenv_spNB[,t]
  colz=rep("black",length(nbs))
  for(i in 1:length(colpal)){
    colz[nbs==i] <- colpal[i]
  }
  plot_data <- rbind(plot_data,cbind(TD_sp[,t],TDenv_sp[,t]))
  plot(TD_sp[,t],TDenv_sp[,t],col=colz,pch=16,xlim=c(-3.5,2.5),ylim=c(-3,5),
       xlab="TD, species mean",cex=.3,
       ylab="TDenv, species mean")
  for(t in 3:ncol(TD_sp)){
    
    nbs=TDenv_spNB[,t]
    colz=rep("black",length(nbs))
    for(i in 1:length(colpal)){
      colz[nbs==i] <- colpal[i]
    }
    plot_data <- rbind(plot_data,cbind(TD_sp[,t],TDenv_sp[,t]))
    points(TD_sp[,t],TDenv_sp[,t],col=colz,pch=16,
           xlab="TD, species mean",cex=.6,
           ylab="TDenv, species mean")
  }
  colnames(plot_data) = c("obs","pred")
  dim(plot_data)
  plot_data <- as.data.frame(plot_data)
  text(x = -1,y = 4, labels = paste0("Pearson cor=",round(cor(plot_data[complete.cases(plot_data),])[2,1],digits=2)))
  abline(0,1)
}
dev.off()


