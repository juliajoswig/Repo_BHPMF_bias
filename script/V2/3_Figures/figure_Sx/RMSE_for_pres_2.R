setwd("/..")
origin = "Volumes/bgi/work_1/2016_GapFilling"
#origin="2016_GapFilling"
res_matrix_name="res_20201020"
#res_o=res

res <- read.table(file.path(origin,"runs","META",paste0(res_matrix_name,".csv")),sep=",",dec=".")
summary(res)
list.files(file.path(origin,"runs","META"))
#colz=rainbow(7)#c("red","blue","green","magenta","turquoise","orange","yellow")#
#colz=c("#d73027","#f46d43","#fdae61","#fee08b","#d9ef8b","#a6d96a","#66bd63","#1a9850")
#colz=c("#d7191c","#2c7bb6","#fdae61","#abd9e9")#"#ffffbf",
#colz1=c("#f0f9e8","#bae4bc","#7bccc4","#43a2ca","#0868ac")#"#ffffbf",
#colz2=c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")#"#ffffbf",
output_term="2020"

# RMSE increases with gap-size
# RMSE different for traits
# RMSE different for data set

#--------------------------------------
# Silhouette increases with gap-size
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="Sil_"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,grep(colnames(dat_rmse),pattern = "_md")]
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
dat_rmse2 <- dat_rmse2[,grep(colnames(dat_rmse2),pattern = "_md")]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]


{
#--------------------------------------
# Sil along gap size
#--------------------------------------
pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size.pdf"),width=12,height=6)
par(mfrow=c(2,4),mar=c(4,4,2,2))
i=2

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
