setwd("/..")
origin = "Volumes/bgi/work_1/2016_GapFilling"
#origin="2016_GapFilling"
res_matrix_name="res_20201126"#"res_20201112"
res_matrix_name="res_20210211"
res <- read.table(file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")

#--------------------------------------
# PCA increases with gap-size
#--------------------------------------
    sel_ObSpe1="Obs_obs_TD"
    pattern_now="PCA"
    #chose guido
    ix1=res[,colnames(res)=="TraitChoice"]=="data"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
    Percent1 <- as.vector(res[ix1,colnames(res)=="GapPercent"])
    mode(Percent1) <- "numeric"
    
    dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
    colnames(dat_rmse)
    
    dat_rmse <- as.matrix(dat_rmse)
    mode(dat_rmse)="numeric"
    colnames(dat_rmse)
    
    dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
    dat_plot <- dat_plot[complete.cases(dat_plot),]
    
    #chose rainfor
    ix2=as.vector(res[,colnames(res)=="TraitChoice"]=="data_2"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1)
    mode(Percent1) <- "numeric"
    
    dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
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
# Correl along gap size
#--------------------------------------
dat_plot2 <- dat_plot2[1:8,]
pdf(file=file.path(origin,"_2021","figures","Figure_3","Figure_3_PCA.pdf"),width=10,height=8)
  par(mfrow=c(1,1),mar=c(7,7,1,1))
  i=2
  plot(dat_plot2[dat_plot2[,1]!=0,c(1,i)],col=colz[i],pch=16,cex.axis=2,
       xlab="Missingness [% gaps]",ylab="Variance explained",xlim=c(0,60),ylim=c(0,60),cex.lab=3)
  abline(h=dat_plot2[dat_plot2[,1]==0,c(1,i)],col=colz[i],lwd=2)
  abline(v=seq(1,100,by=10),col="gray")
  points(dat_plot2[dat_plot2[,1]!=0,c(1,i)],col=colz[i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz[i],lwd=8)
  i=3
  points(dat_plot2[dat_plot2[,1]!=0,c(1,i)],col=colz[i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz[i],lwd=8)
  i=4
  points(dat_plot2[dat_plot2[,1]!=0,c(1,i)],col=colz[i],pch=16,
         main=colnames(dat_plot2)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
  lines(dat_plot2[,c(1,i)][order(dat_plot2[,1]),],col=colz[i],lwd=8)
  
  level_names=c("PCA 1","PCA 2","PCA 3")
  legend(x = 40,y = 47, level_names, col = colz[1:3],
              text.col = "black", lty = 1, lwd=4,cex=2,
              merge = TRUE, bg = "white")
dev.off()


pdf(file=file.path(origin,"plots","Silhouettes","Sil_Gap_size_AllTaxa_Guido.pdf"),width=12,height=6)
par(mfrow=c(1,1),mar=c(4,4,2,2))
i=2
    
    plot(dat_plot[dat_plot[,1]!=0,c(1,i)],col=colz[i],pch=16,
         xlab="Missingness [% gaps]",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
    #rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==0,i], xright=75, ytop=.5, density = NULL, angle = 45,
    #     border = "white",col=colup)
    #rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==0,i], density = NULL, angle = 45,
    #     border = "white",col=coldown)
    abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
    abline(h=dat_plot[dat_plot[,1]==0,c(1,i)],col=colz[i],lwd=2)
    points(dat_plot[dat_plot[,1]!=0,c(1,i)],col=colz[i],pch=16,
           main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
    lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
    lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz[i],lwd=3)
    
    for(i in 3:5){
      points(dat_plot[dat_plot[,1]!=0,c(1,i)],col=colz[i],pch=16,
             main=colnames(dat_plot)[i],ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
      #  rect(xleft=-.1, ybottom=dat_plot[dat_plot[,1]==0,i], xright=75, ytop=.5, density = NULL, angle = 45,
      #       border = "white",col=colup)
      #  rect(xleft=-.1, ybottom=-.6, xright=75, ytop=dat_plot[dat_plot[,1]==0,i], density = NULL, angle = 45,
      #       border = "white",col=coldown)
      abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
      abline(h=dat_plot[dat_plot[,1]==0,c(1,i)],col=colz[i],lwd=2)
      points(dat_plot[dat_plot[,1]!=0,c(1,i)],col=colz[i],pch=16,
             main=colnames(dat_plot)[i],xlab="% Gaps",ylab="Silhouette index",xlim=c(0,72),ylim=c(-.5,.2))
      lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col="black",lwd=3.5)
      lines(dat_plot[,c(1,i)][order(dat_plot[,1]),],col=colz[i],lwd=3)
    }

dev.off()
