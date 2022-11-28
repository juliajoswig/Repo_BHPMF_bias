
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

RepNum=1 
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
#TD <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
TD80 <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
head(TD80)
TD80 <- TD80[,-c(1,2)]

RepNum=1 
t_choice="data"
ObsOrTD="Obs_obs"
Percent=80
Env <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
head(Env)
dim(Env)
ID_env <- Env[,2]
Env <- Env[,-c(1,2)]

dev.off()
plot(TDtd[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDtd[,1])
abline(0,1)
plot(Env[match(ID_TD,ID_env),1],TDtd[,1])
plot(Env[match(ID_TD,ID_env),1],TDenv[,1])
abline(0,1)
print("Yaaay loaded :)")


# table
m1 <- matrix(NA,ncol=4,nrow=length(trait_names))
m2 <- matrix(NA,ncol=4,nrow=length(trait_names))
m3 <- matrix(NA,ncol=4,nrow=length(trait_names))
m4 <- matrix(NA,ncol=4,nrow=length(trait_names))
colnames(m1) <- c("TD","TD_td","Env","TD_env")
colnames(m2) <- c("TD","TD_td","Env","TD_env")
colnames(m3) <- c("TD","TD_td","Env","TD_env")
colnames(m4) <- c("TD","TD_td","Env","TD_env")
rownames(m1) <- trait_names
rownames(m2) <- trait_names
rownames(m3) <- trait_names
rownames(m4) <- trait_names

t=2
for(t in 1:length(trait_names)){
  m1[t,1]<- quantile(TD[,colnames(TD)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,2]<- quantile(TDtd[,colnames(TDtd)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .5,na.rm = TRUE)
  
  m2[t,1]<- quantile(x = TD[,colnames(TD)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,2]<- quantile(x = TDtd[,colnames(TDtd)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .9,na.rm = TRUE)
  
  m3[t,1]<- quantile(x = TD[,colnames(TD)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,2]<- quantile(x = TDtd[,colnames(TDtd)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .1,na.rm = TRUE)
  
  m4[t,1]<- mean(x = TD[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m4[t,2]<- mean(x = TDtd[,colnames(TDtd)==trait_names[t]],na.rm = TRUE)
  m4[t,3] <- mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
  m4[t,4] <- mean(TDenv[,colnames(TDenv)==trait_names[t]],na.rm = TRUE)
  
}
   
mode(m1) <- "numeric"
t=1
for(t in 1:6){
  dat_plot=rbind(m1[t,],m2[t,],m3[t,])
  rownames(dat_plot) <- c("10th quantile","median","90th quantile")
  heatmap(dat_plot,Rowv = NA)
}
m_tot=rep(NA,4)
for(t in 1:6){
  m_now <- rbind(m1[t,],m2[t,],m3[t,])
  rownames(m_now) <- paste0(trait_names[t],c(" 10qnt"," md"," 90qnt"))
  m_tot <- rbind(m_tot,m_now)
}
heatmap(m_tot,Rowv = NA)
m_tot <- m_tot[complete.cases(m_tot),]
heatmap(m_tot)

heatmap(m1)
heatmap(m2)
heatmap(m3)
heatmap(m4)

units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")
dev.off()

m_now=m4
pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_Mean.pdf"),width=15,height=20)
  t=1
  par(mar=c(7,7,1,1))
  layout(matrix(c(1, 2, 2, 2, 2, 2,8,
                  1, 3,3,3,3, 3,8,
                  1, 4, 4, 4, 4,4, 8,
                  1, 5, 5, 5, 5,5, 8,
                  1, 6,6,6,6,6,8,
                  1,7,7,7,7,7,8), nrow=6, byrow=TRUE))
  #layout.show(n=8)
  plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
  text(4,5,labels=c("Distance to observed envelope data (mean)"),srt=90,cex=7)

  for(t in 1:6){
    dat_plot=c(abs(m_now[t,4]-m_now[t,3]),
               abs(m_now[t,2]-m_now[t,3]),
               abs(m_now[t,4]-m_now[t,1]),
               abs(m_now[t,2]-m_now[t,1]))
    barplot(abs(dat_plot[c(3:4)]),
            ylab=trait_names[t],cex.lab=4,xaxt="n",
            horiz = TRUE,
            col=colz[c(1,2,1,2)],
            xlim=c(-max(abs(dat_plot))*1.25,max(abs(dat_plot))*1.25),
            xlab="")
    barplot(-dat_plot[c(1:2)],add=TRUE,
            horiz = TRUE,xaxt="n",
            col=colz[c(1,2,1,2)],
            xlim=c(-max(abs(dat_plot))*1.25,max(abs(dat_plot))*1.25))
    abline(v=0,lwd=5,col="white")
    if(t==1){
      legend(.8,2.4,legend = c("test data", "envelope data"),title = "TD; filled with",pch=15,col=colz[2:1],cex=4,bty = "n",pt.cex = 4)
    }
    
    axis(1,cex.axis=1,lwd=3,lwd.ticks=0,col.axis="white")
    axis(1,  lwd=0, lwd.ticks=3,col.axis="white")
    axis(1,line=1.5,cex.axis=3,lwd.ticks=0,tick=FALSE)
    axis(1,line=4,at = 0,units[t],cex.axis=3,tick = FALSE)
    
  }
  plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
  text(5,5,labels=c("Distance to observed test data (mean)"),srt=-90,cex=7)

dev.off()
