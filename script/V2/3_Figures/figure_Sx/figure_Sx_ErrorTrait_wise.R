# get per observed point value: 
# - error (pred-obs)
# - average distance to all other points observed 
# - Silhouette index for this group (or calculate somewhere else...?)  
# - ORIGINAL Silhouette index for this group (or calculate somewhere else...?)  
# - DEVIATION Silhouette index for this group (or calculate somewhere else...?)  
# - number of values available
# - original number of values available

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
list.files(file.path(origin,"_2021","script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
  tsubs <- out$tsubs
  TD_choices = out$TD_choices
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
  rmse_now <- function(a,b){
    return(sqrt(mean(b-a,na.rm=TRUE)^2))
  }
  rmse_now_median <- function(a,b){
    return(sqrt(median(b-a,na.rm=TRUE)^2))
  }
  
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

res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise","res2.csv"))

res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
m=1
t=2
w=1


for(m in 1:length(missingness)){}
  par(mfrow=c(1,1),mar=c(10,6,2,2))
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  sum(ix_trait)
  
  bxpl=NA
  bxpl_gap=NA
  bxpl_median=NA
  bxpl_gap_median=NA
  
  
  
  pdf(file=file.path(origin,"_2021","Figures","Figure_Sx","Fig_Sx_Error.pdf"),width=12,height=9)
  par(mfrow=c(2,3),mar=c(4,9,4,2))
  t=1
  limmin <- c(0,0,0,0,0,0)
  limmax <- c(65,65,1.1,60,8,7)
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    #print(sum(ix_trait))
    ix_trait2=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
    ix_trait2[is.na(ix_trait2)] <- FALSE
    ix_trait3=res$trait==trait_names[t]&res$missingness==missingness[m]&!res$Gap
    ix_trait3[is.na(ix_trait2)] <- FALSE
    
    plot(res$value_obs[ix_trait],res$value_pred[ix_trait],xlim=c(limmin[t],limmax[t]),ylim=c(limmin[t],limmax[t]),
         main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=2,
         ylab="Predicted (untransformed)",xlab="Observed (untransformed)")
    #    rect(xleft = -4,ybottom = -4,xright = 4,ytop = 4,col = "black")
    rect(xleft = limmin[t]-10,ybottom = limmin[t]-10,xright = limmax[t]+10,ytop = limmax[t]+10,col = "gray")
    points(res$value_obs[ix_trait],res$value_pred[ix_trait],xlim=c(limmin[t],limmax[t]),ylim=c(limmin[t],limmax[t]),
         main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=2,
         ylab="Predicted (untransformed)",xlab="Observed (untransformed)")
    #    rect(xleft = -4,ybottom = -4,xright = 4,ytop = 4,col = "black")
    abline(0,1)
    points(res$value_obs[ix_trait2],res$value_pred[ix_trait2],#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col="white",pch=15,cex=.5,
           ylab="Predicted zlog value",xlab="Observed zlog value")
    points(res$value_obs[ix_trait3],res$value_pred[ix_trait3],#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col="black",pch=16,cex=.5,
           ylab="Predicted zlog value",xlab="Observed zlog value")
    abline(0,1)
    text(limmax[t]*.7,limmax[t]*.9,round(rmse_now(res$value_obs[ix_trait2],res$value_pred[ix_trait2]),digits = 3),cex=3,col="black")
    text(limmax[t]*.25,limmax[t]*.9,"RMSE=",cex=3,col="black")
    
  }
  dev.off()
  
  pdf(file=file.path(origin,"_2021","Figures","Figure_Sx","Fig_Sx_Errorzlog.pdf"),width=12,height=9)
  par(mfrow=c(2,3),mar=c(4,9,4,2))
  t=1
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    #print(sum(ix_trait))
    ix_trait2=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
    ix_trait2[is.na(ix_trait2)] <- FALSE
    ix_trait3=res$trait==trait_names[t]&res$missingness==missingness[m]&!res$Gap
    ix_trait3[is.na(ix_trait2)] <- FALSE
    #print(sum(ix_trait2))
    print(min(res$value_obs_zlog[ix_trait]))
    print(max(res$value_obs_zlog[ix_trait]))
    
    plot(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
         main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=2,
         ylab="Predicted (zlog)",xlab="Observed (zlog)")
    rect(xleft = -4,ybottom = -4,xright = 4,ytop = 4,col = "gray")
    points(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=2,
           ylab="Predicted (zlog)",xlab="Observed (zlog)")
    abline(0,1)
    points(res$value_obs_zlog[ix_trait2],res$value_pred_zlog[ix_trait2],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col="white",pch=16,cex=.5,
           ylab="Predicted zlog value",xlab="Observed zlog value")
    points(res$value_obs_zlog[ix_trait3],res$value_pred_zlog[ix_trait3],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col="black",pch=16,cex=.5,
           ylab="Predicted zlog value",xlab="Observed zlog value")
    abline(0,1)
    text(1,2.7,round(rmse_now(res$value_obs_zlog[ix_trait2],res$value_pred_zlog[ix_trait2]),digits = 3),cex=3,col="black")
    text(-1.7,2.7,"RMSE=",cex=3,col="black")
    
  }
  dev.off()
  
  pdf(file=file.path(origin,"_2021","Figures","Figure_Sx","Fig_Sx_ErrorzlogTOT.pdf"),width=6,height=8)
  par(mfrow=c(1,1),mar=c(6,6,1,1))
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  #print(sum(ix_trait))
  ix_trait2=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait2[is.na(ix_trait2)] <- FALSE
  ix_trait3=res$trait==trait_names[t]&res$missingness==missingness[m]&!res$Gap
  ix_trait3[is.na(ix_trait2)] <- FALSE
  #print(sum(ix_trait2))
  print(min(res$value_obs_zlog[ix_trait]))
  print(max(res$value_obs_zlog[ix_trait]))
  
  plot(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-4,4),ylim=c(-4,4),#xlim=c(-4,4),ylim=c(-4,4),
       col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=.5,
       ylab="Predicted (zlog)",xlab="Observed (zlog)")
  rect(xleft = -4,ybottom = -4,xright = 4,ytop = 4,col = "gray")
  points(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
         main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=.6,
         ylab="Predicted (zlog)",xlab="Observed (zlog)")
t=2
  for(t in 2:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    #print(sum(ix_trait))
    ix_trait2=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
    ix_trait2[is.na(ix_trait2)] <- FALSE
    ix_trait3=res$trait==trait_names[t]&res$missingness==missingness[m]&!res$Gap
    ix_trait3[is.na(ix_trait2)] <- FALSE
    #print(sum(ix_trait2))
    print(min(res$value_obs_zlog[ix_trait]))
    print(max(res$value_obs_zlog[ix_trait]))
    
    points(res$value_obs_zlog[ix_trait],res$value_pred_zlog[ix_trait],xlim=c(-3,3),ylim=c(-3,3),#xlim=c(-4,4),ylim=c(-4,4),
           main=trait_names[t],col=colz[t],pch=16,cex.main=4,cex.axis=2,cex.lab=2.2,cex=.6,
           ylab="Predicted (zlog)",xlab="Observed (zlog)")
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

