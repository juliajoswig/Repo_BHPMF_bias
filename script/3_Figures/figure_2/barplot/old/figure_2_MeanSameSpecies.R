
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
t_choice="data"
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")
colz=colz1


{
  # load Envelope data
  # load TDenvelope
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD",t_choice,"traitInfo.csv"),header=TRUE))[,-c(1,2)]
  taxTD <- as.data.frame(read.table(file = file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD",t_choice,"taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_TD <- taxTD[,1]
  head(taxTD)
  dim(taxTD)
  dim(TD)
  # load Envelope data
  # total trait data 
  EnvTot <- as.data.frame(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs",t_choice,"traitInfo.csv"),header=TRUE))[,-1]
  taxEnv <- as.data.frame(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs",t_choice,"taxInfo.csv"),header=TRUE))
  taxEnv <- as.data.frame(read.table(file = file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                     "Obs_obs",t_choice,"taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_env <- taxEnv[,1]
  dim(taxEnv)
  summary(EnvTot)
  Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
  head(Env)
  # chose only those Env with same species information
  Env_sp <- Env[taxEnv[,2]%in%taxTD[,2],]
  dim(Env_sp)
  Env_gen <- Env[taxEnv[,3]%in%taxTD[,3],]
  dim(Env_gen)
  Env_fam <- Env[taxEnv[,4]%in%taxTD[,4],]
  dim(Env_fam)
  Env_cl <- Env[taxEnv[,5]%in%taxTD[,5],]
  dim(Env_cl)
  
  list.files(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD",t_choice))
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs",t_choice,"traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD",t_choice,"traitInfoTD_pred.csv")))[,-c(1,2)]
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
}

#-------------------------------------------------------------------
# chose 
#-------------------------------------------------------------------



dev.off()
plot(TDtd[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDtd[,1])
abline(0,1)
length(ID_env)
length(ID_TD)

plot(Env[match(ID_TD,ID_env),2],TDtd[,2])
plot(Env[match(ID_TD,ID_env),1],TDenv[,1])
abline(0,1)
print("Yaaay loaded :)")

if(t_choice=="data"){trait_names=trait_rainfor}
if(t_choice=="data_2"){trait_names=trait_guido}

# table
m1 <- matrix(NA,ncol=5,nrow=length(trait_names))
m2 <- matrix(NA,ncol=5,nrow=length(trait_names))
m3 <- matrix(NA,ncol=5,nrow=length(trait_names))
m4 <- matrix(NA,ncol=5,nrow=length(trait_names))
colnames(m1) <- c("TD","TD_td","Env","TD_env","Env_sp")
colnames(m2) <-  c("TD","TD_td","Env","TD_env","Env_sp")
colnames(m3) <-  c("TD","TD_td","Env","TD_env","Env_sp")
colnames(m4) <-  c("TD","TD_td","Env","TD_env","Env_sp")
rownames(m1) <- trait_names
rownames(m2) <- trait_names
rownames(m3) <- trait_names
rownames(m4) <- trait_names

t=1
for(t in 1:length(trait_names)){
  m1[t,1]<- quantile(TD[,colnames(TD)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,2]<- quantile(TDtd[,colnames(TDtd)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,5] <- quantile(Env_sp[,colnames(Env_sp)==trait_names[t]],probs = .5,na.rm = TRUE)
  
  m2[t,1]<- quantile(x = TD[,colnames(TD)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,2]<- quantile(x = TDtd[,colnames(TDtd)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .9,na.rm = TRUE)
  m2[t,5] <- quantile(Env_sp[,colnames(Env_sp)==trait_names[t]],probs = .9,na.rm = TRUE)
  
  m3[t,1]<- quantile(x = TD[,colnames(TD)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,2]<- quantile(x = TDtd[,colnames(TDtd)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,3] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,4] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .1,na.rm = TRUE)
  m3[t,5] <- quantile(Env_sp[,colnames(Env_sp)==trait_names[t]],probs = .1,na.rm = TRUE)
  
  m4[t,1]<- mean(x = TD[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m4[t,2]<- mean(x = TDtd[,colnames(TDtd)==trait_names[t]],na.rm = TRUE)
  m4[t,3] <- mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
  m4[t,4] <- mean(TDenv[,colnames(TDenv)==trait_names[t]],na.rm = TRUE)
  m4[t,5] <- mean(Env_sp[,colnames(Env_sp)==trait_names[t]],na.rm = TRUE)
  
}
   
mode(m1) <- "numeric"
t=1
for(t in 1:6){
  dat_plot=rbind(m1[t,],m2[t,],m3[t,])
  rownames(dat_plot) <- c("10th quantile","median","90th quantile")
  heatmap(dat_plot,Rowv = NA,main=trait_names[t])
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

require(FactoMineR)
PCA(rbind(m1,m2,m4))
dev.off()


units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")
m_now=m4[,c(1,2,3,5)]

pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_Mean.pdf"),width=15,height=20)
par(mar=c(7,7,1,1),mfrow=c(1,7))

colz=c("#d1e5f0",#TD
       "#ef8a62",#TDtd
       "#2166ac",#Env
       #"#67a9cf",#Env_sp
       "#b2182b"#TDenv
       )
t=1
for(t in 1:length(trait_names)){
barplot(m_now[t,],col=colz,border = "white",ylab=units[t],cex.lab=2,main=trait_names[t])
#abline(h=seq(0,max(m_now[t,]),length.out=6),col="gray",lty=2)
#barplot(m_now[t,],col=colz,border = "white",add=TRUE)
}
dev.off()









#plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
#text(4,5,labels=c("Distance to observed envelope data (mean)"),srt=90,cex=7)


#########################################################
### C) Customizing and plotting the heat map
#########################################################
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

# creates a own color palette from red to green
my_palette <- colorRampPalette(c(colz[1], colz[2], colz[3], colz[4]))(n = 399)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,1,length=100),  # for red
               seq(1.01,3,length=100),           # for yellow
               seq(3.1,9,length=100),           # for yellow
               seq(9.1,23,length=100))             # for green

# creates a 5 x 5 inch image
png("../images/heatmaps_in_r.png",    # create PNG for the heat map        
    width = 5*400,        # 5 x 300 pixels
    height = 5*400,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap.2(t(m_now),
          cellnote = round(t(m_now),digits = 2),  # same data set for cell labels
          main = "mean (total)", # heat map title
          notecol="black",      # change font color of cell labels to black
          density.info="none",  # turns off density plot inside color legend
          trace="none",         # turns off trace lines inside the heat map
          margins =c(12,9),     # widens margins around plot
          col=my_palette,       # use on color palette defined earlier
          breaks=col_breaks,    # enable color transition at specified limits
          dendrogram="row",     # only draw a row dendrogram
          Colv="NA")            # turn off column clustering































t=1
for(t in 1:6){
  dat_plot=c(abs(m_now[t,4]-m_now[t,3]),
             abs(m_now[t,2]-m_now[t,3]),
             abs(m_now[t,4]-m_now[t,1]),
             abs(m_now[t,2]-m_now[t,1]))
  dat_plot <- t(matrix(dat_plot))
  colnames(dat_plot) <- c("TDenv - Env","TDtd - Env",
                          "TDenv - TD","TDtd - TD")
  barplot(abs(dat_plot),col=colz[c(1:4)],
          main=trait_names[t],cex.lab=4,las=2,
          xlab="")
  
}
plot(1:10,1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",xlab="",ylab="")
text(5,5,labels=c("Distance to observed test data (mean)"),srt=-90,cex=7)

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
