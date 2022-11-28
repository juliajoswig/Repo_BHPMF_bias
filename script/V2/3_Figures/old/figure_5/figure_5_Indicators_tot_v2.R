

#------------------------------------------------------------
# define path
#------------------------------------------------------------
is.it.on.cluster=TRUE
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
  t_choices <- out$t_choices
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
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
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
  colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .5),rgb(103/255,169/255,207/255,alpha = .5))
  colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))
  colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
  colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
  
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  Percent=80
  res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  res <- res[,colSums(!is.na(res))!=0]
  dim(res)
  
  # load the for the total (envelope + test data 1 and 2)
  # takes quite long, as big
  PercentTot=1
  ObsOrTD="Obs_obsTot"
  res2 <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,PercentTot,paste0("res_ind.csv")))
  dim(res2)
  
  try(dev.off(),silent = TRUE)
  t=1
  plot(cbind(res_now2[colnames(res2)==paste0("CV_spec_pred_",trait_names[t])],
             res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])]))

  plot(cbind(res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])],
    res_now2[colnames(res2)==paste0("CV_spec_pred_",trait_names[t])]-
             res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])]))
  
  
#-------------------------------------------------------------------
# load sd
#-------------------------------------------------------------------
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
  cnfd <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","std.csv"), sep="\t",header=TRUE)

#-------------------------------------------------------------------
# colors
#-------------------------------------------------------------------
  xs <- round(res$nb_spec,digits = 3)
  xg <- round(res$nb_gen,digits = 3)
  xf <- round(res$nb_fam,digits = 3)
  xc <- round(res$nb_clad,digits = 3)
  xc2 <- round(res$nb_fam+res$nb_clad,digits = 3)
  
  xs <- res$nb_spec==1
  xg <- res$nb_gen==1
  xf <- res$nb_fam==1
  xc <- res$nb_clad==1
  xs[is.na(xs)] <- FALSE
  xg[is.na(xg)] <- FALSE
  xf[is.na(xf)] <- FALSE
  xc[is.na(xc)] <- FALSE
  perc=.60
  c_NB_sp=rep(0,nrow(res))
  c_NB_sp[xs] <- c_NB_sp[xs]+1
  
  c_NB_gen=rep(0,nrow(res))
  c_NB_gen[xg] <- c_NB_gen[xg]+1
  
  c_NB_fam=rep(0,nrow(res))
  c_NB_fam[xf] <- c_NB_fam[xf]+1
  
  c_NB_clad=rep(0,nrow(res))
  c_NB_clad[xc] <- c_NB_clad[xc]+1


  #c_NB_fc=rep(0,nrow(res))
  #ix=xc2>quantile(xc2,probs = perc,na.rm = TRUE);ix[is.na(ix)] <- FALSE
  #c_NB_fc[ix] <- c_NB_fc[ix]+1
  
  rgb_input <- cbind(c_NB_sp,c_NB_gen,c_NB_fam)
  rgb_input[rgb_input==0] <- .1
  #rgb_input[is.na(rgb_input)] <- 0
  colz_rgbinv=rgb((1-rgb_input),alpha = .9)
  colz_rgb=rgb(rgb_input,alpha = .9)
  colz_rgb
  
  try(dev.off(),silent = TRUE)
  
  plot(res$value_obs,res$value_pred,col=colz_rgbinv,pch=16)
  
  plot(res2$value_obs,res2$value_pred,col=colz_rgbinv,pch=16)# the total data set
  
  plot(res$value_obs,res$value_pred,col=colz_rgb,pch=16)
  plot(cbind(res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])],
             res_now2[colnames(res2)==paste0("CV_spec_pred_",trait_names[t])]-
               res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])]),col=colz_rgbinv,pch=16)
  
# ------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  
  print(paste(RepNum,ObsOrTD,t_choice,Percent))
  
  t=1
  dat_cor <- list()
  dat_tot=matrix(NA,ncol=32,nrow=0)
  colnames(dat_tot) <- c("Error","Confid","Dev spec","Dev gen","Dev fam","Dev clad","Dev GF","Dev PFT","DevCV spec","DevCV gen","DevCV fam","DevCV clad","DevCV GF","DevCV PFT",
    "Dist_spec","Dist_gen","Dist_fam","Dist_clad","Dist_GF","Dist_PFT",
    "Nb_spec","Nb_gen","Nb_fam","Nb_clad","Nb_GF","Nb_PFT","CV_spec","CV_gen","CV_fam","CV_clad","CV_GF","CV_PFT")
  colnames(dat_tot) <- c("Error","Confidence","Dev species","Dev genus","Dev family","Dev clade","Dev GF","Dev PFT",
                         "DevCV species","DevCV genus","DevCV family","DevCV clade","DevCV GF","DevCV PFT",
                         "Distance species","Distance genus","Distance family","Distance clade","Distance GF","Distance PFT",
                         "Number of species","Number of genera","Number of families","Number of clades","Number of GFs","Number of PFTs",
                         "CV_spec","CV_gen","CV_fam","CV_clad","CV_GF","CV_PFT")
  t=1
  for(t in 1:length(trait_names)){
    
    res_now=res[res$trait==trait_names[t],]
    res_now <- res_now[!is.na(res_now$trait),]

    res_now2=res2[res2$trait==trait_names[t],]
    res_now2 <- res_now2[!is.na(res_now2$trait),]

    print(dim(res_now))
    print(length(cnfd[,t]))
    dat_plot <- as.matrix(cbind(#dav
      
      abs(res_now$value_pred_zlog-res_now$value_obs_zlog),#Error
      cnfd[,t],#SD
      abs(res_now$dist_spec_pred_zlog-res_now$dist_spec_obs_zlog),#dev tax
      abs(res_now$dist_gen_pred_zlog-res_now$dist_gen_obs_zlog),#dev tax
      abs(res_now$dist_fam_pred_zlog-res_now$dist_fam_obs_zlog),#dev tax
      abs(res_now$dist_clad_pred_zlog-res_now$dist_clad_obs_zlog),#dev tax
      abs(res_now$dist_GF_pred_zlog-res_now$dist_GF_obs_zlog),#dev tax
      abs(res_now$dist_PFT_pred_zlog-res_now$dist_PFT_obs_zlog),#dev tax

      # CV
      abs(res_now2[colnames(res2)==paste0("CV_spec_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])]),#dev tax
      abs(res_now2[colnames(res2)==paste0("CV_gen_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_gen_obs_",trait_names[t])]),#dev tax
      abs(res_now2[colnames(res2)==paste0("CV_fam_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_fam_obs_",trait_names[t])]),#dev tax
      abs(res_now2[colnames(res2)==paste0("CV_clad_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_clad_obs_",trait_names[t])]),#dev tax
      abs(res_now2[colnames(res2)==paste0("CV_GF_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_GF_obs_",trait_names[t])]),#dev tax
      abs(res_now2[colnames(res2)==paste0("CV_PFT_pred_",trait_names[t])]-res_now2[colnames(res2)==paste0("CV_PFT_obs_",trait_names[t])]),#dev tax
      
      #dist
      res_now$dist_spec_obs_zlog,
      res_now$dist_gen_obs_zlog,
      res_now$dist_fam_obs_zlog,
      res_now$dist_clad_obs_zlog,
      res_now$dist_GF_obs_zlog,
      res_now$dist_PFT_obs_zlog,
      
      # nb
      res_now$nb_spec,
      res_now$nb_gen,
      res_now$nb_fam,
      res_now$nb_clad,
      res_now$nb_GF,
      res_now$nb_PFT,
      
      # CV
      res_now2[colnames(res2)==paste0("CV_spec_obs_",trait_names[t])],#dev tax
      res_now2[colnames(res2)==paste0("CV_gen_obs_",trait_names[t])],#dev tax
      res_now2[colnames(res2)==paste0("CV_fam_obs_",trait_names[t])],#dev tax
      res_now2[colnames(res2)==paste0("CV_clad_obs_",trait_names[t])],#dev tax
      res_now2[colnames(res2)==paste0("CV_GF_obs_",trait_names[t])],#dev tax
      res_now2[colnames(res2)==paste0("CV_PFT_obs_",trait_names[t])]#dev tax
))
  #    grep(colnames(res_now),pattern = "_lm_obs_zlog"))
    dim(dat_plot)
    head(dat_plot)
    
       colnames(dat_plot) <- c("Error","Confidence",
                              "Deviation species","Deviation genus","Deviation family","Deviation clade","Deviation GF","Deviation PFT",
                              "DevCV species","DevCV genus","DevCV family","DevCV clade","DevCV GF","DevCV PFT",
                              "Distance species","Distance genus","Distance family","Distance clade","Distance GF","Distance PFT",
                              "Number of species","Number of genera","Number of families","Number of clades","Number of GFs","Number of PFTs",
                              "CV_spec","CV_gen","CV_fam","CV_clad","CV_GF","CV_PFT")
      dat_plot <- as.data.frame(dat_plot)
      plot(cbind(dat_plot$`DevCV species`, dat_plot$CV_spec),col=colz_rgbinv,pch=16)
       
      dat_cor[[t]] <- as.matrix(cor(dat_plot,use = "pairwise.complete.obs"))
      
      dat_tot <- rbind(dat_tot,dat_plot)
  }
#dat_tot <- dat_tot[,colSums(is.na(dat_tot))!=0]

  #dat_cor <- as.matrix(cor(dat_tot,use = "pairwise.complete.obs"))
  par(mar=c(10,8,5,0))
  boxplot(as.matrix(dat_cor[[1]]),las=2)
  abline(h=0)
  
  pdf(file=file.path(origin,"_2021","figures","Figure_5","Figure_5_Indicators_V2.pdf"),width=13,height=12)
  par(mar=c(7,8,5,0))
  layout(mat = matrix(c(1,3,6,7,
                      1,4,8,10,
                      2,5,9,11),nrow=4))
{ 
  t=1
  i=1
  dat_now=rbind(dat_cor[[1]][i,],dat_cor[[2]][i,],dat_cor[[3]][i,],dat_cor[[4]][i,],dat_cor[[5]][i,],dat_cor[[6]][i,])
  
  boxplot(abs(dat_now),las=2,col=colz[c(3,5)],ylab="Pearson correlation coeff",cex.main=3,
          cex.lab=2,cex.axis=1,ylim=c(0,1),frame=FALSE,main=colnames(dat_cor[[2]])[dat_cor[[2]][i,]==1])
  abline(v=seq(.5,50,1),col=colz[c(3)],lty=1,lwd=1)
  abline(v=seq(.5,50,2),col="gray",lty=1,lwd=1)
  abline(h=1,col="gray",lty=1,lwd=2)
  abline(h=0,col="gray",lty=1,lwd=2)
  abline(h=seq(-1,1,.5),col="gray",lty=2,lwd=1)
  rect(xleft = -1,ybottom = 1,xright = 20,ytop = 2,col="white",border = NA)
  boxplot(dat_now,las=2,col=colz[c(3,5)],cex.main=3,
          cex.lab=2,cex.axis=1,ylim=c(-.5,1),frame=FALSE,add=TRUE)
  i=11
  for(i in 2:11){
    dat_now=rbind(dat_cor[[1]][i,],dat_cor[[2]][i,],dat_cor[[3]][i,],dat_cor[[4]][i,],dat_cor[[5]][i,],dat_cor[[6]][i,])
    boxplot(abs(dat_now[,13:ncol(dat_now)]),las=2,col=colz[c(3,5)],ylab="Pearson correlation coeff",cex.main=3,
            cex.lab=2,cex.axis=1,ylim=c(0,1),frame=FALSE,main=colnames(dat_cor[[2]])[dat_cor[[2]][i,]==1][1])
    abline(v=seq(.5,50,1),col=colz[c(3)],lty=1,lwd=1)
    abline(v=seq(.5,50,2),col="gray",lty=1,lwd=1)
    abline(h=1,col="gray",lty=1,lwd=2)
    abline(h=0,col="gray",lty=1,lwd=2)
    abline(h=seq(-1,1,.5),col="gray",lty=2,lwd=1)
    rect(xleft = -1,ybottom = 1,xright = 20,ytop = 2,col="white",border = NA)
    boxplot(abs(dat_now[,13:ncol(dat_now)]),las=2,col=colz[c(3,5)],ylab="Pearson correlation coeff",cex.main=3,
            cex.lab=2,cex.axis=1,ylim=c(-.5,1),frame=FALSE,main=colnames(dat_cor[[2]])[dat_cor[[2]][i,]==1],add=TRUE)
  }
  
  }
  dev.off()
  
  dat_tab=rep(NA,ncol(dat_now))
  for(i in 1:8){
    dat_now=rbind(dat_cor[[1]][i,],dat_cor[[2]][i,],dat_cor[[3]][i,],dat_cor[[4]][i,],dat_cor[[5]][i,])
    dat_tab <- rbind(dat_tab,apply(dat_now,2,median))
  }
  
  dat_tab <- dat_tab[-1,]
  rownames(dat_tab) <- colnames(dat_cor[[1]])[1:8]
  #install.packages("xtable")  
  library(xtable)
  
  
  tab <- as.data.frame(t(dat_tab))
  xtable(tab)
  #xtable(tab, type = "latex", file = file.path(origin,"_2021","Figures","Tables","Table_Indicators.tex"))
  
  