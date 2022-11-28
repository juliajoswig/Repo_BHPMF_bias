setwd("/..")
origin = "Volumes/bgi/work_1/2016_GapFilling"
#origin="2016_GapFilling"
res_matrix_name="res_20201020"
res_matrix_name="res_20201126"#"res_20201112"
#res_o=res
colz2 <- out$colz2

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

colup="#f6e8c3"
coldown="#f5f5f5"
colup="#bdbdbd"
coldown="#f0f0f0"
colup="white"
coldown="#f0f0f0"
cxlab=1.7

#--------------------------------------
# Sil along gap size
#--------------------------------------
pdf(file=file.path(origin,"_2021","figures","Figure_2","figure_2_SilSpecRF.pdf"),width=6,height=8)
par(mfrow=c(1,1),mar=c(5.5,5.5,1,1))
i=2
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[5],pch=16,
       main="",xlab="Missingness [% gaps]",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2),
       cex.lab=2,cex.axis=2)
  rect(xleft=-10, ybottom=dat_plot2[dat_plot2[,1]==(-1),i], xright=80, ytop=.5, density = NULL, angle = 45,
       border = "white",col=colup)
  rect(xleft=-10, ybottom=-.6, xright=80, ytop=dat_plot2[dat_plot2[,1]==(-1),i], density = NULL, angle = 45,
       border = "white",col=coldown)
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col="black",lwd=3.5)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col=colz2[5],lwd=3)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz2[5],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=3.5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz2[5],lwd=3)

dev.off()
