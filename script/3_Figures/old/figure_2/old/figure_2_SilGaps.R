setwd("/..")
origin = "Volumes/bgi/work_1/2016_GapFilling"
#origin="2016_GapFilling"
res_matrix_name="res_20201020"
res_matrix_name="res_20201126"#"res_20201112"
res_matrix_name="res_20201126"#"res_20201112"
res <- read.table(file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
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
    colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
    colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
#--------------------------------------
# Sil along gap size
#--------------------------------------

dat_plot2 <- dat_plot2[1:10,]
pdf(file=file.path(origin,"_2021","figures","Figure_2","Figure_2_Sil_Missingness.pdf"),width=10,height=8)
  par(mfrow=c(1,1),mar=c(7,7,1,1))
  i=2
  plot(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,cex.axis=2,
       xlab="Missingness [% gaps]",ylab="Silhouette index",xlim=c(0,100),ylim=c(-.5,.2),cex.lab=3)
  #rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
  #     border = "white",col=colup)
  #rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
  #     border = "white",col=coldown)
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col=colz[i],lwd=2)
  points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=5)
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz[i],lwd=4)
  
  for(i in 3:ncol(dat_plot2)){
    points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
         main=colnames(dat_plot2)[i],ylab="Silhouette index",xlim=c(0,60),ylim=c(-.5,.2))
  #  rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
  #       border = "white",col=colup)
  #  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
  #       border = "white",col=coldown)
    abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
    abline(h=dat_plot2[dat_plot2[,1]==(-1),c(1,i)],col=colz[i],lwd=2)
    points(dat_plot2[dat_plot2[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
           main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
    lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col="black",lwd=5)
    lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz[i],lwd=4)
  }
  
  level_names=c("Species","Genus","Family","Clade","PFT","GF")
  rect(xleft=65, ybottom=-.5, xright=103.5, ytop=.1, density = NULL, angle = 45,
            border = "white",col="white")
       legend(70, 0, level_names, col = colz[2:length(colz)],
              text.col = "black", lty = 1, lwd=4,cex=2,
              merge = TRUE, bg = "white")
dev.off()


pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllTaxa_Guido.pdf"),width=12,height=6)
par(mfrow=c(1,1),mar=c(4,4,2,2))
i=2
    
    plot(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
         xlab="Missingness [% gaps]",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
    #rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
    #     border = "white",col=colup)
    #rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
    #     border = "white",col=coldown)
    abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
    abline(h=dat_plot[dat_plot[,1]==(-1),c(1,i)],col=colz[i],lwd=2)
    points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
           main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
    lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
    lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz[i],lwd=3)
    
    for(i in 3:5){
      points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
             main=colnames(dat_plot)[i],ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
      #  rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==(-1),i], xright=75, ytop=.5, density = NULL, angle = 45,
      #       border = "white",col=colup)
      #  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==(-1),i], density = NULL, angle = 45,
      #       border = "white",col=coldown)
      abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
      abline(h=dat_plot[dat_plot[,1]==(-1),c(1,i)],col=colz[i],lwd=2)
      points(dat_plot[dat_plot[,1]!=(-1),c(1,i)],col=colz[i],pch=16,
             main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
      lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
      lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz[i],lwd=3)
    }

dev.off()
