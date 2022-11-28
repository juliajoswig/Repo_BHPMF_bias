

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
# load sd
#-------------------------------------------------------------------
list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
cnfd <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","std.csv"), sep="\t",header=TRUE)

  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  m=1
  t=2
  w=1

  print(paste(RepNum,ObsOrTD,t_choice,Percent))
  
  pdf(file=file.path(origin,"_2021","figures","Figure_5","Figure_5_Ind_clad.pdf"),width=10,height = 8)
  par(mfrow=c(2,3),mar=c(4,4,2,2))
  #-------------------------------------------------------------------
  #  RMSE
  #-------------------------------------------------------------------
  t=2
  dat_plot=c(NA,NA)
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$value_pred_zlog-res_now$value_obs_zlog),res_now$dist_clad_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_alpha[1],xlab="Error (absolute)",ylab="Distance to co-clade plants",cex.lab=2)
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
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_spec_pred_zlog-res_now$dist_spec_obs_zlog),res_now$dist_clad_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_alpha[1],xlab="Deviation pred-obs (sp)",ylab="Distance to co-clade plants",cex.lab=2)
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
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_gen_pred_zlog-res_now$dist_gen_obs_zlog),res_now$dist_clad_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_alpha[1],xlab="Deviation pred-obs (gen)",ylab="Distance to co-clade plants",cex.lab=2)
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
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_fam_pred_zlog-res_now$dist_fam_obs_zlog),res_now$dist_clad_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_alpha[1],xlab="Deviation pred-obs (fam)",ylab="Distance to co-clade plants",cex.lab=2)
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
    dat_plot <- rbind(dat_plot,cbind(abs(res_now$dist_clad_pred_zlog-res_now$dist_clad_obs_zlog),res_now$dist_clad_obs_zlog))
  }
  dat_lm<- data.frame(dat_plot)
  lm_now <- lm(X2~X1,dat_lm)
  plot(dat_plot,pch=16,col=colz_alpha[1],xlab="Deviation pred-obs (clad)",ylab="Distance to co-clade plants",cex.lab=2)
  abline(lm_now$coefficients,lwd=2)
  lm_smr <- summary(lm_now)
  text(.8,2.0,labels = paste0("Pearson: ",round(cor(dat_plot,use = "pairwise.complete.obs")[1,2],digits = 2)),cex = 2)
  text(.8,1.5,labels = paste0("R2 (adj): ", round(lm_smr$adj.r.squared,digits = 2)),cex = 2)
  
  
  dev.off()
  