
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
  

GapPercent=50
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

units=c("mm2 mg-1","m","mg mm-3","mg g-1","mg g-1","g m-2")

RepNum=1
t_choice="data_2"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

RepNum=1
t_choice="data_2"
ObsOrTD="Obs_obs"
Percent=80
resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))
res <- res[,colSums(!is.na(res))!=0]
trait_names=as.vector(unique(res$trait))
trait_names <- trait_names[!is.na(trait_names)]
missingness = unique(as.vector(res$missingness))
missingness <- missingness[!is.na(missingness)]
m=1
t=2
w=1



par(mfrow=c(1,1),mar=c(5,6,1,1))
pdf(file=file.path(origin,"_2021","figures","figure_1","Fig_S1_Density.pdf"),width=18,height=4)

{
  limmin <- c(0,-1,-.1,0,0,0)
  limmax <- c(60,40,50,1,110000,7)
  li <- c(0,0,0,0,0,0)
  lx <- c(.1,.8,.1,4,.00022,.85)
  
  par(mfrow=c(1,5),mar=c(5,6,3,1))
  col_below="gray"
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl <- rep(NA,sum(ix_trait))
  t=5
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    sum(ix_trait,na.rm = T)
    ix_trait[is.na(ix_trait)] <- FALSE
    
    de <- density(res$value_obs[ix_trait],na.rm = T,adjust = 1,bw = "SJ")
    if(t==2){de <- density(res$value_obs[ix_trait],na.rm = T,adjust=3)}
    if(t==3){de <- density(res$value_obs[ix_trait],na.rm = T,adjust=12)}
    if(t==5){de <- density(res$value_obs[ix_trait],na.rm = T,adjust=2)}
    plot(de,col=col_below,lwd=5,cex.main=3,cex.lab=2.5,
         cex.axis=2,main = trait_names[t],xlab=units[t],ylim=c(li[t],lx[t]),
         xlim=c(limmin[t],limmax[t]),lty=1)
    polygon(de, col=col_below, border=col_below) 
    lines(density(resENV$value_pred[ix_trait],bw = de$bw),col="white",lwd=4)
    lines(density(resENV$value_pred[ix_trait],bw = de$bw),col=colz[1],lwd=4)
    lines(density(res$value_pred[ix_trait],bw = de$bw),col="white",lwd=4)
    lines(density(res$value_pred[ix_trait],bw = de$bw),col=colz_solid[2],lwd=4)
    abline(v=median(resENV$value_pred[ix_trait]),col=colz[1],lty=3,lwd=2)  
    abline(v=median(res$value_pred[ix_trait]),col=colz_solid[2],lty=2,lwd=4)  
    abline(v=median(res$value_obs[ix_trait]),col="black",lty=1,lwd=2)  
  }
}
dev.off()

plot(res$value_pred_zlog[ix_trait],resENV$value_pred_zlog[ix_trait])
abline(0,1)

par(mfrow=c(1,1),mar=c(5,6,1,1))
pdf(file=file.path(origin,"_2021","figures","figure_1","Fig_S1_Density.pdf"),width=18,height=4)
par(mfrow=c(1,5),mar=c(5,6,1,1))

{
  limmin <- c(0,-1,-.1,0,0,0)
  limmax <- c(60,55,50,1,110000,7)
  li <- c(0,0,0,0,0,0)
  lx <- c(.1,1.5,.6,3.5,.0005,.85)
  
  col_below="gray"
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl <- rep(NA,sum(ix_trait))
  t=5
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    sum(ix_trait,na.rm = T)
    ix_trait[is.na(ix_trait)] <- FALSE
    max(res$value_obs[ix_trait])
    de <- density(res$value_obs[ix_trait],na.rm = T)
    
    
    plot(de,col=col_below,lwd=5,cex.main=2,cex.lab=3,cex.axis=2,xlab = trait_names[t],main="",ylim=c(li[t],lx[t]),xlim=c(limmin[t],limmax[t]),lty=1)
    #polygon(de, col=col_below, border=col_below) 
      lines(density(res$value_pred[ix_trait],bw = de$bw),col="white",lwd=4)
      lines(density(res$value_pred[ix_trait],bw = de$bw),col=colz_solid[2],lwd=4)
    abline(v=median(res$value_obs[ix_trait]),col="black",lty=1,lwd=2)  
    abline(v=median(res$value_pred[ix_trait]),col=colz_solid[2],lty=2,lwd=2)  
  }
}
dev.off()
