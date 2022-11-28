

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
  #res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  #TDtd <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
  TDtd <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
  TDtd <- TDtd[,-c(1,2)]
  trait_names=colnames(TDtd)
  
  
  RepNum=1 
  t_choice="data"
  ObsOrTD="Obs_obs"
  Percent=80
  #resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  #TDenv <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred_REzlog.csv"))
  TDenv <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_pred.csv"))
  ID_TDenv <- TDenv[,2]
  TDenv <- TDenv[,-c(1,2)]
  head(TDenv)
  
  plot(TDtd[,1],TDenv[,1])
  
  RepNum=1 
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  Percent=0
  #TD <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
  TD <- as.data.frame(read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv")))
  TDtax <- read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"),sep = ",")
  colnames(TDtax) <- c("ID","Species","Genus","Family","Clades")
  dim(TDtax)
  ID_TD <- TD[,2]
  TD <- TD[,-c(1,2)]
  head(TD)
  
  RepNum=1 
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  Percent=80
  #TD <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfoTD_obs_REzlog.csv"))
  TD80 <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
  TD80 <- TD80[,-c(1,2)]
  
  RepNum=1 
  t_choice="data"
  ObsOrTD="Obs_obs"
  Percent=0
  Env <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
  Env <- Env[,-c(1,2)]
  
  RepNum=1 
  t_choice="data"
  ObsOrTD="Obs_obs"
  Percent=80
  Env80 <- read.csv(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
  Envtax <- read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"),sep = ",")
  colnames(Envtax) <- c("ID","Species","Genus","Family","Clades")
  dim(Envtax)
  ID_env <- Env80[,2]
  Env80 <- Env80[,-c(1,2)]
  
  try(dev.off(),silent = TRUE)
  colz=c("#f4a582","#4393c3")
  
  # get the CV
  new.sd.fun <- function(input){
    if(sum(input[!is.na(input)])<=1){return(NA)}
    return(sd(input,na.rm = TRUE))}
  new.count.fun<- function(input){
    return(sum(!is.na(input)))}
  new.sd.fun(c(NA,NA,1,2))  
  
  {
    l <- list()
    l_now <- list()
    l$Env <- list()
    l$TD <- list()
    l$TDtd <- list()
    l$TDenv <- list()

    res <- list()
    mean_now <-aggregate(TD,by = list(unlist(TDtax$Species)),mean,na.rm=TRUE)
    sd_now <- aggregate(TD,by = list(unlist(TDtax$Species)),FUN = new.sd.fun)
    l_now$mean_now <- mean_now
    l_now$sd_now <- sd_now
    #res$TD <- cbind(l_now$sd_now[,1],l_now$sd_now[,2:ncol(TD)]/l_now$mean_now[,2:ncol(TD)])
    
    #TD    
    TD_count <- aggregate(TD,by = list(unlist(TDtax$Species)),FUN = new.count.fun)
    mean_now <-aggregate(TD,by = list(unlist(TDtax$Species)),mean,na.rm=TRUE)
    sd_now <- aggregate(TD,by = list(unlist(TDtax$Species)),FUN = new.sd.fun)
    l_now$mean_now <- mean_now
    l_now$sd_now <- sd_now
    TD_CV <- cbind(l_now$sd_now[,1],l_now$sd_now[,2:ncol(l_now$sd_now)]/l_now$mean_now[,2:ncol(l_now$sd_now)],
                   l_now$sd_now[,2:ncol(l_now$sd_now)],
                   l_now$mean_now[,2:ncol(l_now$sd_now)])
    names(TD_CV) <- c("Species",
                      paste0(trait_names,"_CV"),
                      paste0(trait_names,"_sd"),
                      paste0(trait_names,"_mean"))
    #TDtd    
    mean_now <-aggregate(TDtd,by = list(unlist(TDtax$Species)),mean,na.rm=TRUE)
    sd_now <- aggregate(TDtd,by = list(unlist(TDtax$Species)),FUN = new.sd.fun)
    l_now$mean_now <- mean_now
    l_now$sd_now <- sd_now
    TDtd_CV <- cbind(l_now$sd_now[,1],l_now$sd_now[,2:ncol(l_now$sd_now)]/l_now$mean_now[,2:ncol(l_now$sd_now)],
                   l_now$sd_now[,2:ncol(l_now$sd_now)],
                   l_now$mean_now[,2:ncol(l_now$sd_now)])
    names(TDtd_CV) <- c("Species",
                      paste0(trait_names,"_CV"),
                      paste0(trait_names,"_sd"),
                      paste0(trait_names,"_mean"))
    #TDenv   
    mean_now <-aggregate(TDenv,by = list(unlist(TDtax$Species)),mean,na.rm=TRUE)
    sd_now <- aggregate(TDenv,by = list(unlist(TDtax$Species)),FUN = new.sd.fun)
    l_now$mean_now <- mean_now
    l_now$sd_now <- sd_now
    TDenv_CV <- cbind(l_now$sd_now[,1],l_now$sd_now[,2:ncol(l_now$sd_now)]/l_now$mean_now[,2:ncol(l_now$sd_now)],
                     l_now$sd_now[,2:ncol(l_now$sd_now)],
                     l_now$mean_now[,2:ncol(l_now$sd_now)])
    names(TDenv_CV) <- c("Species",
                        paste0(trait_names,"_CV"),
                        paste0(trait_names,"_sd"),
                        paste0(trait_names,"_mean"))
    #Env    
    Env_count <- aggregate(Env,by = list(unlist(Envtax$Species)),FUN = new.count.fun)
    mean_now <-aggregate(Env,by = list(unlist(Envtax$Species)),mean,na.rm=TRUE)
    sd_now <- aggregate(Env,by = list(unlist(Envtax$Species)),FUN = new.sd.fun)
    l_now$mean_now <- mean_now
    l_now$sd_now <- sd_now
    Env_CV <- cbind(l_now$sd_now[,1],l_now$sd_now[,2:ncol(l_now$sd_now)]/l_now$mean_now[,2:ncol(l_now$sd_now)],
                      l_now$sd_now[,2:ncol(l_now$sd_now)],
                      l_now$mean_now[,2:ncol(l_now$sd_now)])
    names(Env_CV) <- c("Species",
                         paste0(trait_names,"_CV"),
                         paste0(trait_names,"_sd"),
                         paste0(trait_names,"_mean"))
  }
  
  
  
  #------------------------------------------------------------------------------
  # define colours
  colz_list <- c("#fddbc7","#f4a582","#d6604d","#b2182b","#d1e5f0","#92c5de","#4393c3","#2166ac")[7:1]
  colz_list <-c(  "#e0f3db",  "#ccebc5",  "#a8ddb5",  "#7bccc4", "#4eb3d3",  "#2b8cbe",  "#08589e")
  barplot(rep(1,length(colz_list)),col=colz_list)

  colz_l <- list()
  colz_l[[1]] <- rep("black",nrow(TD_count))
  colz_l[[1]][TD_count[,2]==1] <- colz_list[1] 
  colz_l[[1]][TD_count[,2]==2] <- colz_list[2]
  colz_l[[1]][TD_count[,2]>=3&TD_count[,2]<5] <- colz_list[3]
  colz_l[[1]][TD_count[,2]>=5&TD_count[,2]<10] <- colz_list[4]
  colz_l[[1]][TD_count[,2]>=10&TD_count[,2]<15] <- colz_list[5]
  colz_l[[1]][TD_count[,2]>=15&TD_count[,2]<20] <- colz_list[6]
  colz_l[[1]][TD_count[,2]>=20&TD_count[,2]<25] <- colz_list[7]
  colz_l[[1]][TD_count[,2]>=25] <- colz_list[7]
  
  
   #------------------------------------------------------------------------------
  pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_CV_spec_singleTraits.pdf"),width=10,height=5)
  t=2
  for(t in 1:length(trait_names)){
  t=t+1
    par(mfrow=c(1,2))
  
  dat_plot=cbind(TD_CV[,t],TDtd_CV[,t])
  plot(dat_plot,pch=16,col=colz_l[[1]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
       xlab = "Observed test data",
       ylab = "Imputed with test data",
       main=paste0("",colnames(TD_CV)[t]))
  cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
  text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
  abline(0,1)
  
  ix=match(TD_CV$Species,Env_CV[,1])
  colz_l2 <- list()
  colz_l2[[t]] <- rep("gray",nrow(Env_count))
  tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_CV",replacement = ""))
  colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
  colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
  colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
  colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
  colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
  colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
  colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
  colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
  colz_l2[[t]] <- colz_l2[[t]][ix]
  
  dat_plot=cbind(TD_CV[,t],TDenv_CV[,t])
  plot(dat_plot,pch=16,col=colz_l2[[t]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
       xlab = "Observed test data",
       ylab = "Imputed with envelope data",
       main=paste0("",colnames(TD_CV)[t]))
  
  cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
  text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
  abline(0,1)
  }
  dev.off()
  
  
  pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_sd_spec_singleTraits.pdf"),width=10,height=5)
  t=2
  for(t in 1:length(trait_names)){
    t=t+7
    par(mfrow=c(1,2))
    
    dat_plot=cbind(TD_CV[,t],TDtd_CV[,t])
    plot(dat_plot,pch=16,col=colz_l[[1]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with test data",
         main=paste0("",colnames(TD_CV)[t]))
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    ix=match(TD_CV$Species,Env_CV[,1])
    colz_l2 <- list()
    colz_l2[[t]] <- rep("gray",nrow(Env_count))
    tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_sd",replacement = ""))
    colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
    colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
    colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
    colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
    colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
    colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
    colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
    colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
    colz_l2[[t]] <- colz_l2[[t]][ix]
    
    dat_plot=cbind(TD_CV[,t],TDenv_CV[,t])
    plot(dat_plot,pch=16,col=colz_l2[[t]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with envelope data",
         main=paste0("",colnames(TD_CV)[t]))
    
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
  }
  dev.off()
  
  
  
  pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_mean_spec_singleTraits.pdf"),width=10,height=5)
  t=2
  for(t in 1:length(trait_names)){
    t=t+13
    par(mfrow=c(1,2))
    
    dat_plot=cbind(TD_CV[,t],TDtd_CV[,t])
    plot(dat_plot,pch=16,col=colz_l[[1]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with test data",
         main=paste0("",colnames(TD_CV)[t]))
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    ix=match(TD_CV$Species,Env_CV[,1])
    colz_l2 <- list()
    colz_l2[[t]] <- rep("gray",nrow(Env_count))
    tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_mean",replacement = ""))
    colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
    colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
    colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
    colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
    colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
    colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
    colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
    colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
    colz_l2[[t]] <- colz_l2[[t]][ix]
    
    dat_plot=cbind(TD_CV[,t],TDenv_CV[,t])
    plot(dat_plot,pch=16,col=colz_l2[[t]],ylim=c(0,max(TD_CV[,t],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with envelope data",
         main=paste0("",colnames(TD_CV)[t]))
    
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(max(TD_CV[,t],na.rm = TRUE)*.2,max(TD_CV[,t],na.rm = TRUE)*.9,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
  }
  dev.off()
  
  

  #------------------------------------------------------------------------------
  pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_CV_spec.pdf"),width=10,height=5)
  {  
  par(mfrow=c(1,2))
  
  dat_plot=cbind(c(unlist(TD_CV[,2:7])),c(unlist(TDtd_CV[,2:7])))
  rownames(dat_plot) <- 1:nrow(dat_plot)
  plot(dat_plot,pch=16,col=colz_l[[1]],
       xlim=c(0,max(dat_plot,na.rm = TRUE)),
       ylim=c(0,max(dat_plot,na.rm = TRUE)),
       xlab = "Observed test data",
       ylab = "Imputed with test data",
       main=paste0("CV"))
  cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
  text(.2,1.2,paste0("Pearson=",round(cor_now,digits = 2)))
  abline(0,1)
  
  ix=match(TD_CV$Species,Env_CV[,1])
  count_l <- list()
  colz_l2 <- list()
  for(t in 2:7){
    colz_l2[[t]] <- rep("gray",nrow(Env_count))
    tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_CV",replacement = ""))
    colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
    colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
    colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
    colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
    colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
    colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
    colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
    colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
    colz_l2[[t]] <- colz_l2[[t]][ix]
    
    count_l[[t]] <- rep(0,nrow(Env_count))
    count_l[[t]][Env_count[,tenv]==1] <- 1
    count_l[[t]][Env_count[,tenv]==2] <- 2
    count_l[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- "3to5"
    count_l[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- "5to10"
    count_l[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- "10to15"
    count_l[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- "15to20"
    count_l[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- "25"
    count_l[[t]][Env_count[,tenv]>=25] <- colz_list[7]
    count_l[[t]] <- count_l[[t]][ix]
  }
  
  colz_now <- c(colz_l2[[2]])
  for(t in 3:7){colz_now <- c(colz_now,colz_l2[[t]])}
  count_now <- c(count_l[[2]])
  for(t in 3:7){count_now <- c(count_now,count_l[[t]])}
  
  dat_plot=cbind(c(unlist(TD_CV[,2:7])),c(unlist(TDenv_CV[,2:7])))
  plot(dat_plot,pch=16,col=colz_now,
       xlim=c(0,max(TD_CV[,2:7],na.rm = TRUE)),
       ylim=c(0,max(TD_CV[,2:7],na.rm = TRUE)),
       xlab = "Observed test data",
       ylab = "Imputed with envelope data",
       main=paste0("CV"),cex=.5)
  
  cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
  text(.2,1.2,paste0("Pearson=",round(cor_now,digits = 2)))
  abline(0,1)
  }
  dev.off()
  
  
  #------------------------------------------------------------------------------
  pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_CV_spec_barplot.pdf"),width=10,height=5)
  
  par(mfrow=c(1,2))
  
  dat_plot=cbind(c(unlist(TD_CV[,2:7])),c(unlist(TDtd_CV[,2:7])))
  rownames(dat_plot) <- 1:nrow(dat_plot)
  dat_plot <- dat_plot[,2]-dat_plot[,1]
  boxplot(dat_plot[TD_count==2],
          dat_plot[TD_count>3&TD_count<=5],
          dat_plot[TD_count>5&TD_count<=10],
          dat_plot[TD_count>10&TD_count<=15],
          dat_plot[TD_count>10&TD_count<=20],
          dat_plot[TD_count>20&TD_count<=25],
          dat_plot[TD_count>25],col=colz_list)
  abline(h=0)
  
  ix=match(TD_CV$Species,Env_CV[,1])
  dat_plot=cbind(c(unlist(TD_CV[,2:7])),c(unlist(TDenv_CV[,2:7])))
  rownames(dat_plot) <- 1:nrow(dat_plot)
  dat_plot <- dat_plot[,2]-dat_plot[,1]
  boxplot(dat_plot[Env_count[ix]==2],
          dat_plot[TD_count>3&TD_count<=5],
          dat_plot[TD_count>5&TD_count<=10],
          dat_plot[TD_count>10&TD_count<=15],
          dat_plot[TD_count>10&TD_count<=20],
          dat_plot[TD_count>20&TD_count<=25],
          dat_plot[TD_count>25],col=colz_list)
  abline(h=0)
  
  cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
  text(.2,1.2,paste0("Pearson=",round(cor_now,digits = 2)))
  abline(0,1)
  
  dev.off()
  
  
  #------------------------------------------------------------------------------
    pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_sd_spec.pdf"),width=10,height=5)
    
    t=2
    par(mfrow=c(1,2))
    
    dat_plot=cbind(c(unlist(TD_CV[,8:13])),c(unlist(TDtd_CV[,8:13])))
    rownames(dat_plot) <- 1:nrow(dat_plot)
    plot(dat_plot,pch=16,col=colz_l[[1]],
         xlim=c(0,max(dat_plot,na.rm = TRUE)),
         ylim=c(0,max(dat_plot,na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with test data",
         main=paste0("sd"))
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(5,20,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    ix=match(TD_CV$Species,Env_CV[,1])
    
    colz_l2 <- list()
    for(t in 8:13){
      colz_l2[[t]] <- rep("gray",nrow(res))
      tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_sd",replacement = ""))
      colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
      colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
      colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
      colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
      colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
      colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
      colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
      colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
      colz_l2[[t]] <- colz_l2[[t]][ix]
    }
    
    colz_now <- c(colz_l2[[8]])
    for(t in 9:13){
      colz_now <- c(colz_now,colz_l2[[t]])
    }
    dat_plot=cbind(c(unlist(TD_CV[,8:13])),c(unlist(TDenv_CV[,8:13])))
    rownames(dat_plot) <- 1:nrow(dat_plot)
    plot(dat_plot,pch=16,col=colz_now,
         xlim=c(0,max(TD_CV[,8:13],na.rm = TRUE)),
         ylim=c(0,max(TD_CV[,8:13],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with envelope data",
         main=paste0("sd"))
    
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(5,20,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    dev.off()
    
    #------------------------------------------------------------------------------
    pdf(file = file.path(origin,"_2021","figures","Figure_5","fig_5_mean_spec.pdf"),width=10,height=5)
    
    t=2
    par(mfrow=c(1,2))
    
    dat_plot1=cbind(c(unlist(TD_CV[,14:ncol(TD_CV)])),c(unlist(TDtd_CV[,14:ncol(TD_CV)])))
    dat_plot=cbind(c(unlist(TD_CV[,14:ncol(TD_CV)])),c(unlist(TDtd_CV[,14:ncol(TD_CV)])))
    rownames(dat_plot) <- 1:nrow(dat_plot)
    plot(dat_plot,pch=16,col=colz_l[[1]],
         xlim=c(0,max(TD_CV[,14:ncol(TD_CV)],na.rm = TRUE)),
         ylim=c(0,max(TD_CV[,14:ncol(TD_CV)],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with test data",
         main=paste0("mean"))
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(10,50,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    ix=match(TD_CV$Species,Env_CV[,1])
    
    colz_l2 <- list()
    for(t in 14:ncol(TD_CV)){
      colz_l2[[t]] <- rep("gray",nrow(res))
      tenv=which(colnames(Env_count)%in%gsub(colnames(TD_CV)[t],pattern = "_mean",replacement = ""))
      colz_l2[[t]][Env_count[,tenv]==1] <- colz_list[1]
      colz_l2[[t]][Env_count[,tenv]==2] <- colz_list[2]
      colz_l2[[t]][Env_count[,tenv]>=3&Env_count[,tenv]<5] <- colz_list[3]
      colz_l2[[t]][Env_count[,tenv]>=5&Env_count[,tenv]<10] <- colz_list[4]
      colz_l2[[t]][Env_count[,tenv]>=10&Env_count[,tenv]<15] <- colz_list[5]
      colz_l2[[t]][Env_count[,tenv]>=15&Env_count[,tenv]<20] <- colz_list[6]
      colz_l2[[t]][Env_count[,tenv]>=20&Env_count[,tenv]<25] <- colz_list[7]
      colz_l2[[t]][Env_count[,tenv]>=25] <- colz_list[7]
      colz_l2[[t]] <- colz_l2[[t]][ix]
    }
    
    colz_now <- c(colz_l2[[14]])
    for(t in 15:ncol(TD_CV)){
      colz_now <- c(colz_now,colz_l2[[t]])
    }
    dat_plot=cbind(c(unlist(TD_CV[,14:ncol(TD_CV)])),c(unlist(TDenv_CV[,14:ncol(TD_CV)])))
    rownames(dat_plot) <- 1:nrow(dat_plot)
    plot(dat_plot,pch=16,col=colz_now,
         xlim=c(0,max(TD_CV[,14:ncol(TD_CV)],na.rm = TRUE)),
         ylim=c(0,max(TD_CV[,14:ncol(TD_CV)],na.rm = TRUE)),
         xlab = "Observed test data",
         ylab = "Imputed with envelope data",
         main=paste0("mean"))
    
    cor_now <- cor(dat_plot,use = "pairwise.complete.obs")[1,2]
    text(10,50,paste0("Pearson=",round(cor_now,digits = 2)))
    abline(0,1)
    
    dev.off()
    
    
    