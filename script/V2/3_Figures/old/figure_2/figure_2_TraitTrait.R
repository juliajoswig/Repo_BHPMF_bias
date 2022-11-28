
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
units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")
colz=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#a65628")
colz=c("#b2182b","#ef8a62")

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
#res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
#TDtd <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
TDtd <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
TDtd <- TDtd[,-c(1,2)]
TDtd_cor <- cor(TDtd,use = "pairwise.complete.obs")
trait_names=colnames(TDtd)

RepNum=1 
t_choice="data"
ObsOrTD="Obs_obs"
Percent=80
#resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
#TDenv <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
TDenv <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
ID_TDenv <- TDenv[,2]
TDenv <- TDenv[,-c(1,2)]
TDenv_cor <- cor(TDenv,use = "pairwise.complete.obs")
head(TDenv)

plot(TDtd[,1],TDenv[,1])

RepNum=1 
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=0
#TD <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
TD <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
ID_TD <- TD[,2]
head(TD)
TD <- TD[,-c(1,2)]
TD_cor <- cor(TD,use = "pairwise.complete.obs")

RepNum=1 
t_choice="data"
ObsOrTD="Obs_obs"
Percent=80
Env <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
head(Env)
dim(Env)
ID_env <- Env[,2]
Env <- Env[,-c(1,2)]
Env_cor <- cor(Env[,match(colnames(TD),colnames(Env))],use = "pairwise.complete.obs")


par(mfrow=c(1,2))
heatmap(TDtd_cor-TD_cor,scale="none",Colv=NA,Rowv=NA)
heatmap((TDtd_cor-TD_cor)-(TDenv_cor-TD_cor),scale="none",Colv=NA,Rowv=NA)

dev.off()

require(FactoMineR)
dat_plot=cbind(c(TDtd_cor),c(TD_cor),c(Env_cor),c(TDenv_cor))
colnames(dat_plot)=c("TDtd_cor","TD","Env","TDenv")
rownms <- TDtd_cor
for(t in 1:nrow(rownms)){ rownms[t,] <- rownames(rownms)}
for(t in 1:nrow(rownms)){rownms[t,] <- paste0(rownms[t,],"_",rownms[t:nrow(rownms),t])}
rownms
row.names(dat_plot) <- c(rownms)

dat_plot=dat_plot[!(dat_plot[,1]==dat_plot[,2]&dat_plot[,1]==dat_plot[,3]&dat_plot[,1]==dat_plot[,4]),]
#PCA(dat_plot)

#dat_plot <- rbind(dat_plot,rep(0,ncol(dat_plot)),rep(1,ncol(dat_plot)))
#row.names(dat_plot)[c(nrow(dat_plot)-1,nrow(dat_plot))] <- c("0","1")
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_TraitTraitheat.pdf"),width=13,height=11)
par(mar=c(5,1,12,1))
heatmap(dat_plot)
dev.off()


par(mar=c(5,1,12,1))
heatmap(dat_plot[,c(1,2,4,3)])



pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_TraitTrait.pdf"),width=15,height=20)
par(mar=c(10,1,1,1))
layout(matrix(c(1, 2, 2, 2, 2, 2,3,
                1, 2, 2, 2, 2, 2,3,
                1, 2, 2, 2, 2, 2,3,
                1, 2, 2, 2, 2, 2,3,
                1, 2, 2, 2, 2, 2,3
                ), nrow=5, byrow=TRUE))
#layout.show(n=8)
plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
text(4,5,labels=c("Distance to observed envelope data (Pearson)"),srt=90,cex=7)

  dat_plot1 = c(TDtd_cor-TD_cor)
  dat_plot2 = c(TDenv_cor-TD_cor)
  dat_plot3 = c(TDtd_cor-Env_cor)
  dat_plot4 = c(TDenv_cor-Env_cor)
  ylims=c(max(c(abs(dat_plot1)),c(abs(dat_plot2)),c(abs(dat_plot3)),c(abs(dat_plot4))))*1.2
  ylims=c(-ylims,ylims)
  barplot(c(abs(dat_plot1),abs(dat_plot2)),
          col = c(rep(colz[1],length(dat_plot1)),rep(colz[2],length(dat_plot2))),
                  horiz = TRUE,xlim=ylims,xaxt="n")
  legend(.2,89,legend = c("test data", "envelope data"),
         title = "TD; filled with",pch=15,col=colz[2:1],
         cex=4,bty = "n",pt.cex =6)
  
  abline(v=mean(abs(dat_plot1)),lwd=8,lty=1,col="black")
  abline(v=mean(abs(dat_plot2)),lwd=8,lty=1,col="black")
  abline(v=mean(abs(dat_plot1)),lwd=6.5,lty=1,col=colz[2])
  abline(v=mean(abs(dat_plot2)),lwd=6.5,lty=1,col=colz[1])
  
  barplot(-c(abs(dat_plot3),abs(dat_plot4)),add=TRUE,
          col = c(rep(colz[1],length(dat_plot1)),rep(colz[2],length(dat_plot2))),
          horiz = TRUE,xlim=ylims,xaxt="n")
  abline(v=-mean(abs(dat_plot3)),lwd=8,lty=1,col="black")
  abline(v=-mean(abs(dat_plot4)),lwd=8,lty=1,col="black")
  abline(v=-mean(abs(dat_plot3)),lwd=6.5,lty=1,col=colz[1])
  abline(v=-mean(abs(dat_plot4)),lwd=6.5,lty=1,col=colz[2])
  

  axis(1,cex.axis=1,lwd=6,lwd.ticks=0,col.axis="white")
  axis(1,  lwd=0, lwd.ticks=6,col.axis="white")
  axis(1,line=1.5,cex.axis=3,lwd.ticks=0,tick=FALSE)
  axis(1,line=6,at = 0,"Pearson correlation coefficient",cex.axis=6,tick = FALSE)
  
  
  barplot(c(abs(dat_plot1),abs(dat_plot2)),add=TRUE,
          col = c(rep(colz[1],length(dat_plot1)),rep(colz[2],length(dat_plot2))),
          horiz = TRUE,xlim=ylims,xaxt="n")
  barplot(-c(abs(dat_plot3),abs(dat_plot4)),add=TRUE,
          col = c(rep(colz[1],length(dat_plot1)),rep(colz[2],length(dat_plot2))),
          horiz = TRUE,xlim=ylims,xaxt="n")
  
  abline(v=0,col="white",lwd=6)

  
  plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
  text(5,5,labels=c("Distance to observed test data (Pearson)"),srt=-90,cex=7)

  dev.off()
  





