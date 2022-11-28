

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
#-------------------------------------------------------------------
# colors
#-------------------------------------------------------------------

xs <- round(res$dist_spec_obs_zlog/max(res$dist_spec_obs_zlog,na.rm = TRUE),digits = 3)
xg <- round(res$dist_gen_obs_zlog/max(res$dist_gen_obs_zlog,na.rm = TRUE),digits = 3)
xf <- round(res$dist_fam_obs_zlog/max(res$dist_fam_obs_zlog,na.rm = TRUE),digits = 3)
xc <- round(res$dist_clad_obs_zlog/max(res$dist_clad_obs_zlog,na.rm = TRUE),digits = 3)
xc2 <- round((res$dist_fam_obs_zlog+res$dist_clad_obs_zlog)/max((res$dist_fam_obs_zlog+res$dist_clad_obs_zlog)),digits = 3)

# kick out those values which have any tax info
c_DIST=rep(0,nrow(res))
ix=xs>quantile(xs,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST[ix] <- c_DIST[ix]+.25
ix=xg>quantile(xg,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST[ix] <- c_DIST[ix]+.25
ix=xg>quantile(xg,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST[ix] <- c_DIST[ix]+.25
ix=xf>quantile(xf,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST[ix] <- c_DIST[ix]+.25
ix=xc>quantile(xc,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST[ix] <- c_DIST[ix]+.25
#-----
perc=.80
c_DIST_sp=rep(0,nrow(res))
ix=xs>quantile(xs,probs = perc,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST_sp[ix] <- c_DIST_sp[ix]+1

c_DIST_gen=rep(0,nrow(res))
ix=xg>quantile(xg,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST_gen[ix] <- c_DIST_gen[ix]+1

c_DIST_fam=rep(0,nrow(res))
ix=xf>quantile(xf,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST_fam[ix] <- c_DIST_fam[ix]+1

c_DIST_clad=rep(0,nrow(res))
ix=xc>quantile(xc,probs = .90,na.rm = T);ix[is.na(ix)] <- FALSE
c_DIST_clad[ix] <- c_DIST_clad[ix]+1

c_DIST_fc=rep(0,nrow(res))
ix=xc2>quantile(xc2,probs = perc,na.rm = TRUE);ix[is.na(ix)] <- FALSE
c_DIST_fc[ix] <- c_DIST_fc[ix]+1
#-----
xs <- round(res$nb_spec/max(res$nb_spec,na.rm = TRUE),digits = 3)
xg <- round(res$nb_gen/max(res$nb_gen,na.rm = TRUE),digits = 3)
xf <- round(res$nb_fam/max(res$nb_fam,na.rm = TRUE),digits = 3)
xc <- round(res$nb_clad/max(res$nb_clad,na.rm = TRUE),digits = 3)
xc2 <- round((res$nb_fam+res$nb_clad)/max((res$nb_fam+res$nb_clad)),digits = 3)

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

rgb_input <- cbind(xs,xg,xf)
rgb_input <- cbind(c_DIST_sp,c_DIST_gen,c_DIST_fc)
rgb_input <- cbind(c_NB_sp,c_NB_gen,c_NB_fc)
rgb_input[rgb_input==0] <- .1
#rgb_input[is.na(rgb_input)] <- 0
colz_rgbinv=rgb((1-rgb_input),alpha = .9)
colz_rgb=rgb(rgb_input,alpha = .9)
colz_rgb

dev.off()
plot(res$value_obs,res$value_pred,col=colz_rgbinv,pch=16)
plot(res$value_obs,res$value_pred,col=colz_rgb,pch=16)

# ------------------------------------------------------------------------
# 
# ------------------------------------------------------------------------
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  m=1
  t=2
  w=1

  print(paste(RepNum,ObsOrTD,t_choice,Percent))
  
  pdf(file=file.path(origin,"_2021","figures","Figure_5","Figure_5_Dist_Spec.pdf"),width=20,height=5)
  par(mfrow=c(1,5),mar=c(7,4,2,2))
  colz_now=colz_rgbinv
  {
  #-------------------------------------------------------------------
  #  RMSE
  #-------------------------------------------------------------------

  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$value_pred_zlog-res_now$value_obs_zlog),res_now$dist_spec_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_now,xlab="RMSE",ylab="Distance to co-species plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(2,2,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(2,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  #-------------------------------------------------------------------
  #  Species
  #-------------------------------------------------------------------
  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_spec_pred_zlog-res_now$dist_spec_obs_zlog),res_now$dist_spec_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_now,xlab="Deviation pred-obs (sp)",ylab="Distance to co-species plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(.5,2,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(.5,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  #-------------------------------------------------------------------
  #  Genus
  #-------------------------------------------------------------------
  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_gen_pred_zlog-res_now$dist_gen_obs_zlog),res_now$dist_spec_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_now,xlab="Deviation pred-obs (gen)",ylab="Distance to co-species plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(.8,2,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(.8,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  #-------------------------------------------------------------------
  #  Family
  #-------------------------------------------------------------------
  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_fam_pred_zlog-res_now$dist_fam_obs_zlog),res_now$dist_spec_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_now,xlab="Deviation pred-obs (fam)",ylab="Distance to co-species plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(.8,2.0,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(.8,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  #-------------------------------------------------------------------
  #  Clades
  #-------------------------------------------------------------------
  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_clad_pred_zlog-res_now$dist_clad_obs_zlog),res_now$dist_spec_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_now,xlab="Deviation pred-obs (clad)",ylab="Distance to co-species plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(.8,2.0,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(.8,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  }
  dev.off()
  