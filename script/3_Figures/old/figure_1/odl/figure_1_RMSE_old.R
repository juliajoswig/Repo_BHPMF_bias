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
Version_now="V1"
list.files(file.path(origin,"script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin = origin,Version_now = "V1")

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
ObsOrTDs =  out$ObsOrTDs
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido =out$trait_guido
trait_rainfor =out$trait_rainfor
colz1 = out$colz1
colz2 = out$colz2
gappercents= c("1","5","10","20","30","40","50","60")
repnums=1:3

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
colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
aph=.5
colz_add_alpha=c(rgb(178, 24, 43,maxColorValue = 255,alpha = aph),
                 rgb(214, 96, 77,alpha = aph,maxColorValue = 255),
                 rgb(244, 165, 130,alpha = aph,maxColorValue = 255),
                 rgb(253, 219, 199,alpha = aph,maxColorValue = 255),
                 rgb(146, 197, 222,alpha = aph,maxColorValue = 255),
                 rgb(67, 147, 195,alpha = aph,maxColorValue = 255),
                 rgb(33, 102, 172,alpha = aph,maxColorValue = 255))

colz_add_alpha <- rep(rgb(red = 1,green = 1,blue = 1,alpha = .1),6)


GapPercent=50
RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
td=1
p=3
res_matrix_name="res_20210303"
res <- read.csv(file=file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")))

res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
m=1
t=2
w=1

colnames(res)
for(m in 1:length(missingness)){}
  par(mfrow=c(1,1),mar=c(10,6,2,2))
  t=1
  ix_trait=res$trait==trait_names[t]&res$Obs_or_Spec=="Obs_obs_TD"&res$TraitChoice=="data"
  ix_trait[is.na(ix_trait)] <- FALSE
  sum(ix_trait)
  
  bxpl=NA
  bxpl_gap=NA
  bxpl_median=NA
  bxpl_gap_median=NA
  t=1
  pdf(file=file.path(origin,"_2021","Figures","Figure_Sx","Fig_Sx_Error.pdf"),width=12,height=9)
  par(mfrow=c(2,3),mar=c(4,9,4,2))
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait=res$trait==trait_names[t]&res$Obs_or_Spec=="Obs_obs_TD"&res$TraitChoice=="data"
    ix_trait[is.na(ix_trait)] <- FALSE
    print(sum(ix_trait))
    ix_trait2=ix_trait&res$Gap
    ix_trait2[is.na(ix_trait2)] <- FALSE
    ix_trait3=ix_trait!res$Gap
    ix_trait3[is.na(ix_trait2)] <- FALSE
    print(sum(ix_trait2))
    
    plot(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
         main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=2,
         ylab="Predicted (zlog)",xlab="Observed (zlog)")
#    rect(xleft = -4,ybottom = -4,xright = 4,ytop = 4,col = "black")
    abline(0,1)
    points(res$value_obs_zlog[ix_trait2],res$value_pred_zlog[ix_trait2],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
         main=trait_names[t],col="white",pch=16,cex=.5,
         ylab="Predicted zlog value",xlab="Observed zlog value")
    points(res$value_obs_zlog[ix_trait3],res$value_pred_zlog[ix_trait3],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col="black",pch=16,cex=.5,
           ylab="Predicted zlog value",xlab="Observed zlog value")
    abline(0,1)
    
  }
  
  dev.off()
  t=1
  
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    print(sum(ix_trait))
    ix_trait2=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
    ix_trait2[is.na(ix_trait2)] <- FALSE
    print(sum(ix_trait2))
    
    rmse_now(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait])
    mean_now=res$mean_gen_obs[ix_trait]
    
    bxpl=cbind(bxpl,rmse_now(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait]))
    bxpl_gap=cbind(bxpl_gap,rmse_now(res$value_obs_zlog[ix_trait2],res$value_pred_zlog[ix_trait2]))
    bxpl_median=cbind(bxpl_median,rmse_now_median(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait]))
    bxpl_gap_median=cbind(bxpl_gap_median,rmse_now_median(res$value_obs_zlog[ix_trait2],res$value_pred_zlog[ix_trait2]))
    
  }
  
  bxpl <- t(as.matrix(bxpl[!is.na(bxpl)]))
  bxpl_gap <- bxpl_gap[!is.na(bxpl_gap)]
  colnames(bxpl) <- trait_names
  bxpl_median <- t(as.matrix(bxpl_median[!is.na(bxpl_median)]))
  bxpl_gap_median <- bxpl_gap_median[!is.na(bxpl_gap_median)]
  colnames(bxpl_median) <- trait_names
  
  
  pdf(file=file.path(origin,"_2021","Figures","Figure_1","Fig_1_RMSE.pdf"),width=10,height=8)
  par(mfrow=c(1,1),mar=c(7,7,2,2))
  i=2

  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="Missingness [% gaps]",ylab="RMSE",xlim=c(0,100),ylim=c(0,1.2))
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,cex=.1,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(0,1))
  lines(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],lwd=5)
  
  for(i in 3:ncol(dat_plot2)){
    points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,cex=.1,
           main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(0,1))
    lines(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],lwd=5)
  }  
  
  
  lines(dat_plot2[dat_plot2[,1]!=(-1),1],
  rowMeans(dat_plot2[dat_plot2[,1]!=(-1),2:ncol(dat_plot2)]),col="black",lwd=6.5,lty=3)

  rect(xleft = 70,ybottom = 0,xright = 105,ytop = 1.25,col = "white",border = FALSE)
  level_names=c("Average","SLA","PlantHeight", "SSD","LeafN","LeafP","LeafNArea")
  legend(70, 1, level_names, col = c("black",colz[2:length(colz)]),bty = "n",
         text.col = "black", lty = c(3,rep(1,6)), lwd=5,cex=1.8,
         merge = TRUE, bg = "white")  
dev.off()


pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllTaxa.pdf"),width=12,height=6)
par(mfrow=c(2,4),mar=c(4,4,2,2))
i=2

for(i in 2:5){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),c(1,i)],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:5){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  rect(xleft=-.1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()

pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllGroups.pdf"),width=12,height=6)
par(mfrow=c(2,6),mar=c(4,4,2,2))
i=2

for(i in 2:7){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),c(1,i)],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:7){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  rect(xleft=-.1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()
}



#--------------------------------------
# Silhouette for maxima
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="Sil_"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,grep(colnames(dat_rmse),pattern = "_x")]
colnames(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot <- dat_plot[complete.cases(dat_plot),]

#chose rainfor
ix2=res[,colnames(res)=="TraitChoice"]=="rainfor"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
dat_rmse2 <- dat_rmse2[,grep(colnames(dat_rmse2),pattern = "_x")]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]



pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllTaxa_max.pdf"),width=12,height=6)
par(mfrow=c(2,6),mar=c(4,4,2,2))
i=2

for(i in 2:7){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(0.1,1.1))
  rect(xleft=-1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),i],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:7){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(0.1,1.1))
  rect(xleft=-1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),i],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()


#--------------------------------------
# Silhouette for minima
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="Sil_"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,grep(colnames(dat_rmse),pattern = "_i")]
colnames(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot <- dat_plot[complete.cases(dat_plot),]

#chose rainfor
ix2=res[,colnames(res)=="TraitChoice"]=="rainfor"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
dat_rmse2 <- dat_rmse2[,grep(colnames(dat_rmse2),pattern = "_i")]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]



pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllGroups_min.pdf"),width=12,height=6)
par(mfrow=c(2,6),mar=c(4,4,2,2))
i=2

for(i in 2:7){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),i],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:7){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),i],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()



#--------------------------------------
# Silhouette for 25qnt
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="Sil_"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,grep(colnames(dat_rmse),pattern = "_25")]
colnames(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot <- dat_plot[complete.cases(dat_plot),]

#chose rainfor
ix2=res[,colnames(res)=="TraitChoice"]=="rainfor"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
dat_rmse2 <- dat_rmse2[,grep(colnames(dat_rmse2),pattern = "_25")]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]



pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllGroups_25.pdf"),width=12,height=6)
par(mfrow=c(2,6),mar=c(4,4,2,2))
i=2

for(i in 2:7){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),i],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:7){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),i],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()



#--------------------------------------
# Silhouette for 25qnt
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="Sil_"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,grep(colnames(dat_rmse),pattern = "_75")]
colnames(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot <- dat_plot[complete.cases(dat_plot),]

#chose rainfor
ix2=res[,colnames(res)=="TraitChoice"]=="rainfor"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
dat_rmse2 <- dat_rmse2[,grep(colnames(dat_rmse2),pattern = "_75")]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]



pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllGroups_75.pdf"),width=12,height=6)
par(mfrow=c(2,6),mar=c(4,4,2,2))
i=2

for(i in 2:7){
  plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
       main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot[dat_plot[,1]==(-1),i],col=colz1[9-i],lwd=2)
  points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz1[9-i],pch=16,
         main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
  lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz1[9-i],lwd=3)
}


for(i in 2:7){
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
       main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-1,-0.2))
  rect(xleft=-1, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=75, ytop=3, density = NULL, angle = 45,
       border = "white",col="#fde0ef")
  rect(xleft=-1, ybottom=-2, xright=75, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col="#e6f5d0")
  abline(h=seq(-1.2,2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),i],col=colz2[9-i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[9-i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[9-i],lwd=3)
}
dev.off()
