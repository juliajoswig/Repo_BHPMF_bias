

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
  TDnos = out$TDnos
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
  colz = rainbow(ncol(Correls))
  colz=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
         "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
         "#fbb4ae", "#b3cde3", "#ccebc5","#decbe4", "#fed9a6", "#ffffcc","#e5d8bd","#fddaec")

  t_choice="data_2"
  t_choice="data"
  
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                       "Obs_obs","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                       "Obs_obs_TD","data"))
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1)]
  TDID=TD[,1]
  TD <- TD[,-1]
  head(TD)
  TDtd <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","80"),
                                           "Obs_obs_TD","data","traitInfo_pred.csv"),header=TRUE))[,-c(1,2)]
  
  head(TDtd)
  TDenv <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","80"),
                                           "Obs_obs","data","traitInfoTD_pred.csv"),header=TRUE))[,-c(1,2)]
  head(TDenv)
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                       "Obs_obs","data"))
  Env <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                                          "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-c(1)]
  head(Env)
  EnvID=Env[,1]
  Env <- Env[,-1]
  
  summary(cbind(TD[,2],TDtd[,2],TDenv[,2]))
  summary(cbind(TD[,2],TDtd[,2],TDenv[,2],Env[which(EnvID%in%TDID),1]))
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
  plot(TD[,1],Env[which(EnvID%in%TDID),1])
  abline(0,1)
  
  # --- 
  #do the PCAs
  # --- 
  require(FactoMineR)
  trait_names=trait_rainfor
  pca_TDenv <- PCA(log(TDenv))
  pca_TDtd <- PCA(log(TDtd))
  pca_TD <- PCA(log(TD))
  
  #-------------------------------------------------------------------------------------------------
  # plot the PCAs
  #-------------------------------------------------------------------------------------------------
  
  #install.packages("ks")
  # packages
  library(vegan)
  library(ks) 
  library(stats) 
  library(calibrate)
  
  #----------------------------------------------------------------------
  # calculate the window size 
  #----------------------------------------------------------------------
  sc_TD = pca_TD$ind$coord[,1:2]
  sc_TDtd = pca_TDtd$ind$coord[,1:2]
  sc_TDenv = pca_TDenv$ind$coord[,1:2]
  
  if(file.exists(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TD.RData"))){
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TD.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TD.RData"))
  }else{
    sc <- sc_TD #scores(env.pca, choices=c(1,2), display=c("species"))  
    H <- Hpi(x=sc) # optimal bandwidth estimation
    save(H,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TD.RData"))
    est <- kde(x=sc, H=H, compute.cont=TRUE) # kernel density estimation
    save(est,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TD.RData"))
    print("kernel estimation done!")
  }
  
  if(file.exists(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDtd.RData"))){
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDtd.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDtd.RData"))
  }else{
    sc <- sc_TDtd #scores(env.pca, choices=c(1,2), display=c("species"))  
    H <- Hpi(x=sc) # optimal bandwidth estimation
    save(H,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDtd.RData"))
    est <- kde(x=sc, H=H, compute.cont=TRUE) # kernel density estimation
    save(est,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDtd.RData"))
    print("kernel estimation done!")
  }
  
  if(file.exists(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDenv.RData"))){
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDenv.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDenv.RData"))
  }else{
    sc <- sc_TDenv #scores(env.pca, choices=c(1,2), display=c("species"))  
    H <- Hpi(x=sc) # optimal bandwidth estimation
    save(H,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDenv.RData"))
    est <- kde(x=sc, H=H, compute.cont=TRUE) # kernel density estimation
    save(est,file=file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDenv.RData"))
    print("kernel estimation done!")
  }
  
  
  
  #----------------------------------------------------------------------
  # chose one approach
  #----------------------------------------------------------------------
  aprc="TD"
  aprcs=c("TD","TDtd","TDenv")
  
  for(a in 1:3){
    aprc=aprcs[a]
    
  if(aprc=="TD"){
    dat_now=TD
    pca_now=pca_TD
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TD.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TD.RData"))
  }
  
  if(aprc=="TDtd"){
    dat_now=TDenv
    pca_now=pca_TDenv
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDtd.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDtd.RData"))
  }
  if(aprc=="TDenv"){
    dat_now=TDenv
    pca_now=pca_TDenv
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"H_TDenv.RData"))
    load(file.path(origin,"_2021","data","helper_data","Figure_4_PCA",t_choice,"est_TDenv.RData"))
  }
  
  #----------------------------------------------------------------------
  # some adjustments for the plot
  #----------------------------------------------------------------------
  # set contour probabilities for drawing contour levels
  #  cl <- contourLevels(est, prob=c(0.5, 0.05, 0.001), approx=TRUE)
  cl <- contourLevels(est, prob=c(0.8, 0.5, 0.05, 0.01), approx=TRUE)
  
  de_fit.m <- data.frame(matrix(data=NA,nrow=ncol(dat_now),ncol=ncol(dat_now)))
  de_fit.m[,1] <- pca_now$var$coord[,1]
  de_fit.m[,2] <- pca_now$var$coord[,2]
  
  
  rownames(de_fit.m) <- colnames(dat_now)
  colnames(de_fit.m) <- c("PC1","PC2","r2","Pr(>1)")
  
  # arrows
  act.place.x <- vector(mode="numeric",length=nrow(de_fit.m))
  act.place.y <- vector(mode="numeric",length=nrow(de_fit.m))
  de_fit.m <- cbind(de_fit.m,act.place.x,act.place.y)
  
  # set contour probabilities for drawing contour levels
  cl <- contourLevels(est, prob=c(0.8, 0.5, 0.05, 0.01), approx=TRUE)
  
  de_fit.m <- data.frame(matrix(data=NA,nrow=ncol(dat_now),ncol=2))
  de_fit.m[,1] <- pca_now$var$coord[,1]
  de_fit.m[,2] <- pca_now$var$coord[,2]
  pca_now$eig[3,]
  
  rownames(de_fit.m) <- rownames(pca_now$var$coord)
  colnames(de_fit.m) <- c("PC1","PC2")
  
  # arrows
  act.place.x <- vector(mode="numeric",length=nrow(de_fit.m))
  act.place.y <- vector(mode="numeric",length=nrow(de_fit.m))
  de_fit.m <- cbind(de_fit.m,act.place.x,act.place.y)
  
  
  #----------------------------------------------------------------------
  # plot PCA
  #----------------------------------------------------------------------
  
  if(t_choice=="data_2"){    pdf(file=file.path(origin,"_2021", "figures","Figure_4",paste0("figure_S4_PCA",aprc,t_choice,".pdf")),height=14,width=16)}
  if(t_choice=="data"){    pdf(file=file.path(origin,"_2021", "figures","Figure_4",paste0("figure_4_PCA",aprc,t_choice,".pdf")),height=14,width=16)}
  {
   
  xlims=c(-10,10)
#  if(aprc=="TD"){xlims=c(10,-10)}

    par(mfrow=c(1,1),mar=c(9,10,2,2),bg="white")
    plot(est, cont=seq(1,100,by=2), display="filled.contour2", add=FALSE, cex.axis=3,cex.lab=3,frame=F,
         ylab="",xaxt="n",yaxt="n",
         xlab="",
         xlim=xlims,
         ylim=c(-10,10)
    )
    ## add the tick
    axis(1, at = seq(from=-10,to=10,by = 5), label = rep("",5) , tck = -0.02,cex.lab=3)
    ## add the labels
    axis(1, at = seq(from=-10,to=10,by = 5), line = 2, lwd = 0, cex.axis = 0.9,cex.axis=3)
    ## add the label
    axis(side = 1,at=-.2,
         labels = expression('PC'[1]),tick = FALSE,
         cex.axis=5,line = 7)
    axis(side = 1,at=3.1,
         labels = paste0("(",round(pca_now$eig[1,2],digits = 0)," %)"),tick = FALSE,
         cex.axis=5,line = 6)
    
    ## add the tick
    axis(2, at = seq(from=-10,to=10,by = 5), label = rep("",5) , tck = -0.02,cex.lab=3)
    ## add the labels
    axis(2, at = seq(from=-10,to=10,by = 5), line = 2, lwd = 0, cex.axis = 0.9,cex.axis=3)
    ## add the label
    axis(side = 2,at=-.8,
         labels = expression('PC'[2]),tick = FALSE,
         cex.axis=5,line = 4.7)
    axis(side = 2,at=3,
         labels = paste0("(",round(pca_now$eig[2,2],digits = 0)," %)"),tick = FALSE,
         cex.axis=5,line = 5.2)
    
    plot(est,abs.cont=cl[1], labels="0.2",labcex=1, add=TRUE, lwd=1)
    plot(est,abs.cont=cl[2], labels="0.5",labcex=1, add=TRUE, lwd=1)
    plot(est,abs.cont=cl[3], labels="0.95",labcex=1, add=TRUE, lwd=1)
    plot(est,abs.cont=cl[4], labels="0.99",labcex=1, add=TRUE, lwd=1)
    
    abline(h=0, lty=3, lwd=0.5)
    abline(v=0, lty=3, lwd=0.5)
    
    ##########################################################################################
    # add arrows and trait names
    
    i=1
    for(i in 1:nrow(de_fit.m)){
      x <- c(0,de_fit.m[i,1]*8)
      y <- c(0,de_fit.m[i,2]*8)
      ## draw arrows from point to point :
      s <- seq(length(x)-1)  # one shorter than data
      arrows(x[s], y[s], x[s+1], y[s+1],angle = 25, col = colz[i],lwd=2.2)
      adjx=1.1;adjy=1.1
      
        text(x[2]*adjx, y[2]*adjy, labels = trait_names[i],
             cex=3,col=colz[i]) 
    }
    }
  
  dev.off()
}  
  
  
  
  
  
  
  
  
  
  
  #-----------------------------------------------------------------------------------------
  # Procrustes
  #-----------------------------------------------------------------------------------------
  suppressPackageStartupMessages(library(vegan))
  suppressPackageStartupMessages(library(ade4))


  #-----------------------------------------------------------------------------------------
  # create input table
  #-----------------------------------------------------------------------------------------
  mat=matrix(NA, nrow=12,ncol=3+ncol(TD))
  colnames(mat) <- c("Sums of squares","Correlation in a symmetric Procrustes rotation","Significance",colnames(TD))
  rownames(mat) <- c("Traits (TD,TDtd)","x","y","Values (TD,TDtd)",
                     "Traits (TD,TDenv)","x","y","Values (TD,TDenv)",
                     "Traits (TDenv,TDtd)","x","y","Values (TDenv,TDtd)")
  pca1=pca_TD
  pca2=pca_TDtd
  
  pro <- procrustes(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)],permutations = 999)
  mat[1,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))
  mat[2,4:ncol(mat)] <-  round(abs(p$Yrot[,1]-p$X[,1]),digits = 2)
  mat[3,4:ncol(mat)] <-  round(abs(p$Yrot[,2]-p$X[,2]),digits = 2)
  #barplot(abs(p$Yrot[,1]-p$X[,1]),las=2,ylim=c(0,.1))
  barplot(cbind(abs(p$Yrot[,1]-p$X[,1]),abs(p$Yrot[,2]-p$X[,2]))[c(1,7,2,8,3,9,4,10,5,11,6,12)],las=2,ylim=c(0,.1),col=colz[c(1,1,2,2,3,3,4,4,5,5,6,6)])
  axis(1,labels = trait_names,at=seq(1.25,14,2.35),las=2,tick=FALSE)
  
  pca1=pca_TD
  pca2=pca_TDenv
  
  pro <- procrustes(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)],permutations = 999)
  mat[5,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))
  mat[6,4:ncol(mat)] <-  round(abs(p$Yrot[,1]-p$X[,1]),digits = 2)
  mat[7,4:ncol(mat)] <-  round(abs(p$Yrot[,2]-p$X[,2]),digits = 2)
  #barplot(abs(p$Yrot[,1]-p$X[,1]),las=2,ylim=c(0,.1))
  barplot(cbind(abs(p$Yrot[,1]-p$X[,1]),abs(p$Yrot[,2]-p$X[,2]))[c(1,7,2,8,3,9,4,10,5,11,6,12)],las=2,ylim=c(0,.1),col=colz[c(1,1,2,2,3,3,4,4,5,5,6,6)])
  axis(1,labels = trait_names,at=seq(1.25,14,2.35),las=2,tick=FALSE)
  
  
  pca1=pca_TDtd
  pca2=pca_TDenv
  
  pro <- procrustes(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$var$coord[,c(1,2)], Y = pca2$var$coord[,c(1,2)],permutations = 999)
  mat[9,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))
  mat[10,4:ncol(mat)] <-  round(abs(p$Yrot[,1]-p$X[,1]),digits = 2)
  mat[11,4:ncol(mat)] <-  round(abs(p$Yrot[,2]-p$X[,2]),digits = 2)
  #barplot(abs(p$Yrot[,1]-p$X[,1]),las=2,ylim=c(0,.1))
  barplot(cbind(abs(p$Yrot[,1]-p$X[,1]),abs(p$Yrot[,2]-p$X[,2]))[c(1,7,2,8,3,9,4,10,5,11,6,12)],las=2,ylim=c(0,.1),col=colz[c(1,1,2,2,3,3,4,4,5,5,6,6)])
  axis(1,labels = trait_names,at=seq(1.25,14,2.35),las=2,tick=FALSE)
  
  #-----------------------------------------------------------------------------------------
  # add single points
  #-----------------------------------------------------------------------------------------
  pca1=pca_TD
  pca2=pca_TDtd
  
  pro <- procrustes(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)],permutations = 999)
  mat[4,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))

  pca1=pca_TD
  pca2=pca_TDenv
  
  pro <- procrustes(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)],permutations = 999)
  mat[8,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))
  
  
  pca1=pca_TDtd
  pca2=pca_TDenv
  
  pro <- procrustes(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)], symmetric = FALSE)
  
  plot(pro, kind = 1, type = "text")
  plot(pro, kind = 2, type = "text")
  p  <- protest(X = pca1$ind$coord[,c(1,2)], Y = pca2$ind$coord[,c(1,2)],permutations = 999)
  mat[12,1:3] <- c(round(p$ss,digits = 3),round(p$scale,digits = 3),round(p$signif,digits = 4))
  
  require(xtable)
print(xtable(mat))  
    