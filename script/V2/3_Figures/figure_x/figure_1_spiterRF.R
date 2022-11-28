
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

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
TD_choice="Obs_obs_TD"
trait_sub="rainfor"
RepNum=1
GapPercent=60
{
  # completely observed
  if(TD_choice=="Obs_obs_TD"&trait_sub=="rainfor"){traitInfo_obs=out$rainforTD_observed;taxInfo=out$rainforTD_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs_TD}
  if(TD_choice=="Obs_obs_TD"&trait_sub=="guido"){traitInfo_obs=out$guidoTD_observed;taxInfo=out$guidoTD_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs_TD}
  if(TD_choice=="Obs_obs"&trait_sub=="rainfor"){traitInfo_obs=out$rainfor_observed;taxInfo=out$rainfor_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs}
  if(TD_choice=="Obs_obs"&trait_sub=="guido"){traitInfo_obs=out$guido_observed;taxInfo=out$guido_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs}
  #sparse observed
  #  traitInfo_sparse <- load(file = file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","TestData_org.RData"))
  traitInfo_sparse <- read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","traitInfo.csv"), sep=",", dec=".")
  # load the output for trait predictions: mean.csv
  traitInfo_pred_zlog <- as.matrix(read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), 
                                                        paste0("p",GapPercent,"_",trait_sub),TD_choice,"data/mean.csv"),
                                              sep="\t", dec=".",header=TRUE))
  
  
  #-------------------------------------------------------------------
  # cut to test data only if necessary
  #-------------------------------------------------------------------
  if(trait_sub=="guido"){ traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
  if(trait_sub=="rainfor"){
    traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]
    traitInfo_sparse  <- traitInfo_sparse[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}
  if(trait_sub=="guido"){ traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
  if(trait_sub=="rainfor"){traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  
  
  #-------------------------------------------------------------------
  # back transformation necessary.
  #-------------------------------------------------------------------
  summary(traitInfo_sparse)
  summary(traitInfo_obs)
  summary(traitInfo_obs_c)
  summary(traitInfoTD_pred_zlog)
  traitInfo_obs <- traitInfo_obs_c#NEW 2021
  #    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_c[,1])
  #    abline(0,1)
  
  traitInfo_obs_zlog <- traitInfo_obs
  for(i in 1:ncol(traitInfo_obs)){
    traitInfo_obs_zlog[,i] <- (log(traitInfo_obs[,i]) - scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),1])/
      scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),2]
  }
  
  #    plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_zlog[,1])
  #    abline(0,1)
  
  traitInfo_pred <- traitInfoTD_pred_zlog
  for(i in 1:ncol(traitInfoTD_pred_zlog)){
    traitInfo_pred[,i] <- exp((traitInfoTD_pred_zlog[,i]*scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),2])+ 
                                scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),1])
  }
  
  #    plot(traitInfo_pred[,1],traitInfo_obs_c[,1])
  #    abline(0,1)
  traitInfoTD_pred_zlog_sparse <- traitInfoTD_pred_zlog
  traitInfoTD_pred_zlog_sparse[is.na(traitInfo_sparse)] <- NA
  traitInfo_obs_zlog_sparse <- traitInfo_obs_zlog
  traitInfo_obs_zlog_sparse[is.na(traitInfo_sparse)] <- NA
  
}

#res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise","res.csv"))

# colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
# res <- res[,colSums(!is.na(res))!=0]
# trait_names = as.vector(unique(res$trait))
# trait_names <- trait_names[!is.na(trait_names)]
# missingness = unique(as.vector(res$missingness))
# missingness <- missingness[!is.na(missingness)]
# summary(res_now$value_obs)
m=1
t=1

#install.packages("randomForest")
library(randomForest)
pairs(traitInfo_obs_zlog)
par(mfrow=c(5,2))
plot(traitInfo_obs_zlog[,c(1,3)],xlim=c(-2.5,3.5),ylim=c(-4,2))# bit change
plot(traitInfo_pred_zlog[,c(1,3)],xlim=c(-2.5,3.5),ylim=c(-4,2))

plot(traitInfo_obs_zlog[,c(1,2)],xlim=c(-2.5,3.5),ylim=c(-3,1))# no change?
plot(traitInfo_pred_zlog[,c(1,2)],xlim=c(-2.5,3.5),ylim=c(-3,1))

plot(traitInfo_obs_zlog[,c(1,4)],xlim=c(-2.5,3),ylim=c(-3,3))#NO change
plot(traitInfo_pred_zlog[,c(1,4)],xlim=c(-2.5,3),ylim=c(-3,3))

plot(traitInfo_obs_zlog[,c(1,5)],xlim=c(-2.5,3.5),ylim=c(-3,3))# bit change?
plot(traitInfo_pred_zlog[,c(1,5)],xlim=c(-2.5,3.5),ylim=c(-3,3))

plot(traitInfo_obs_zlog[,c(1,6)],xlim=c(-2.5,3.5),ylim=c(-3,3))# bit change?
plot(traitInfo_pred_zlog[,c(1,6)],xlim=c(-2.5,3.5),ylim=c(-3,3))

plot(traitInfo_obs_zlog[,c(2,3)],xlim=c(-2.5,1),ylim=c(-4,2))#change
plot(traitInfo_pred_zlog[,c(2,3)],xlim=c(-2.5,1),ylim=c(-4,2))


  taxnow<- matrix(NA,ncol=ncol(taxInfo),nrow=nrow(taxInfo))
  colnames(taxnow) <- colnames(taxInfo)
  for(i in 1:ncol(taxInfo)){
  cats <- factor(taxInfo[,i])
  ranks <- rank(-table(cats), ties.method="first")
  DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
  taxnow[,i] <- DF[,2]
  }
  t=1
  importances <- list()
  for(t in 1:6){
    dat_obs <- cbind(taxnow[,-1],traitInfo_obs_zlog)
    dat_pred <- cbind(taxnow[,-1],traitInfo_pred_zlog)
    dat_obs <- dat_obs[complete.cases(dat_obs),]
    dat_pred <- dat_pred[complete.cases(dat_pred),]
    #    dat_obs <- cbind(taxnow[,2:5],traitInfo_obs_zlog)
#    dat_pred <- cbind(taxnow[,2:5],traitInfo_pred_zlog)
    colnames(dat_obs)[(6+t)] <- "trait"
    colnames(dat_pred)[(6+t)] <- "trait"
    # Perform training:
    importances[[t]] <- rep(NA,(ncol(dat_obs)-1))
    for(i in 1:10){
    rf_obs = randomForest(trait ~ ., data=dat_obs[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
    rf_pred = randomForest(trait ~ ., data=dat_pred[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
    importances[[t]] <- cbind(importances[[t]],rf_obs$importance[,1],rf_pred$importance[,1])
    }
    }
  
  library(fmsb)
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,seq(3,to=ncol(importances[[t]]),by = 2)] - importances[[t]][,seq(2,to=ncol(importances[[t]]),by = 2)]
    data_plot[data_plot<0] <- -.01
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.1,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  
  
  taxnow<- matrix(NA,ncol=ncol(taxInfo),nrow=nrow(taxInfo))
  colnames(taxnow) <- colnames(taxInfo)
  for(i in 1:ncol(taxInfo)){
    cats <- factor(taxInfo[,i])
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow[,i] <- DF[,2]
  }
  
  t=1
  importances <- list()
  for(t in 1:6){
    dat_obs <- cbind(taxnow[,c(2:5)],traitInfo_obs_zlog)
    dat_pred <- cbind(taxnow[,c(2:5)],traitInfo_pred_zlog)
    dat_obs <- dat_obs[complete.cases(dat_obs),]
    dat_pred <- dat_pred[complete.cases(dat_pred),]
    #    dat_obs <- cbind(taxnow[,2:5],traitInfo_obs_zlog)
    #    dat_pred <- cbind(taxnow[,2:5],traitInfo_pred_zlog)
    colnames(dat_obs)[(4+t)] <- "trait"
    colnames(dat_pred)[(4+t)] <- "trait"
    importances[[t]] <- rep(NA,(ncol(dat_obs)-1))
    # Perform training:
    for(i in 1:10){
      rf_obs = randomForest(trait ~ ., data=dat_obs[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
      rf_pred = randomForest(trait ~ ., data=dat_pred[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
      importances[[t]] <- cbind(importances[[t]],rf_obs$importance[,1],rf_pred$importance[,1])
    }
  }
  
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF_Nofunction.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,seq(3,to=ncol(importances[[t]]),by = 2)] - importances[[t]][,seq(2,to=ncol(importances[[t]]),by = 2)]
    data_plot[data_plot<0] <- -.01
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.1,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  

  
  taxnow<- matrix(NA,ncol=ncol(taxInfo),nrow=nrow(taxInfo))
  colnames(taxnow) <- colnames(taxInfo)
  for(i in 1:ncol(taxInfo)){
    cats <- factor(taxInfo[,i])
    ranks <- rank(-table(cats), ties.method="first")
    DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
    taxnow[,i] <- DF[,2]
  }
  
  unique(taxnow[!is.na(taxnow[,6]),6])
  taxnow[is.na(taxnow[,6]),6] <- 9
#  taxnow[is.na(taxnow[,6]),6] <- c(9:(8+(length(taxnow[is.na(taxnow[,6]),6]))))
  unique(taxnow[!is.na(taxnow[,7]),7])
  taxnow[is.na(taxnow[,7]),7] <- 14
#  taxnow[is.na(taxnow[,7]),7] <- c(14:(13+(length(taxnow[is.na(taxnow[,7]),7]))))
  importances <- list()
  for(t in 1:6){
    dat_obs <- cbind(taxnow[,-1],traitInfo_obs_zlog)
    dat_pred <- cbind(taxnow[,-1],traitInfo_pred_zlog)
    dat_obs <- dat_obs[complete.cases(dat_obs),]
    dat_pred <- dat_pred[complete.cases(dat_pred),]
    #    dat_obs <- cbind(taxnow[,2:5],traitInfo_obs_zlog)
    #    dat_pred <- cbind(taxnow[,2:5],traitInfo_pred_zlog)
    colnames(dat_obs)[(6+t)] <- "trait"
    colnames(dat_pred)[(6+t)] <- "trait"
    importances[[t]] <- rep(NA,(ncol(dat_obs)-1))
    # Perform training:
    for(i in 1:10){
      rf_obs = randomForest(trait ~ ., data=dat_obs[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
      rf_pred = randomForest(trait ~ ., data=dat_pred[sample(1:nrow(dat_obs),round(nrow(dat_obs)*.9)),], ntree=100, mtry=2, importance=TRUE)
      importances[[t]] <- cbind(importances[[t]],rf_obs$importance[,1],rf_pred$importance[,1])
    }
  }
  
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF_Filledfunction.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  t=2
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,seq(3,to=ncol(importances[[t]]),by = 2)] - importances[[t]][,seq(2,to=ncol(importances[[t]]),by = 2)]
    data_plot[data_plot<0] <- -.01
#    data_plot[data_plot>.15] <- .15
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.1,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  

  
  library(fmsb)
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,3] - importances[[t]][,1]
    data_plot[data_plot<0] <- -.01
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.1,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  
  
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,3] - importances[[t]][,1]
    data_plot[data_plot<0] <- -.01
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.1,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  
  
  
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_RF_op.pdf"),width = 6,height=5)
  #mar=c(10,4,2,2)
  par(mar=c(2,2,2,2),mfrow=c(1,2))
  t=3
  for(t in 1:6){
    #  barplot(c(importances[[t]][,1],NA,importances[[t]][,3]),las=2,main=trait_names[t],ylim=c(0,1))
    
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,1]
    data_plot[data_plot<0] <- 0
    data_plot <- as.data.frame(t(data_plot))
    data_now <- rbind(rep(.3,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    
    # for data 1 = predicted - observed
    data_plot <- importances[[t]][,3]
    data_plot[data_plot<0] <- 0
    data_plot <- as.data.frame(t(data_plot))
    data_now <- rbind(rep(.3,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    
  }
  
  