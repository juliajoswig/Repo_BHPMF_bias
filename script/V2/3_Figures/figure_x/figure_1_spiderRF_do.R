
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
colz=c("#92c5de","#0571b0","#ca0020","#f4a582")
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
TD_choice="Obs_obs_TD"
TD_choice="Obs_obs"
trait_sub="rainfor"
RepNum=1
GapPercent=60

load(file=file.path(origin,"_2021","data","analyses","RF",trait_sub,"Obs_obs","importances.RData"))
imp <- importances
load(file=file.path(origin,"_2021","data","analyses","RF",trait_sub,"Obs_obs_TD","importances.RData"))
imp_TD <- importances

  library(fmsb)
pdf(file=file.path(origin,"_2021","figures","Figure_4",paste0("Figure_4_RF_singlecompare.pdf")),width = 6,height=5)
par(mar=c(2,2,2,2),mfrow=c(1,1))
for(t in 1:6){
  data_plot_obsTD <- imp_TD[[t]][,seq(2,to=ncol(imp_TD[[t]]),by = 2)]
  data_plot_predTD <-imp_TD[[t]][,seq(3,to=ncol(imp_TD[[t]]),by = 2)]
  data_plot_obs <-   imp[[t]][,seq(2,to=ncol(imp[[t]]),by = 2)]
  data_plot_pred <-  imp[[t]][,seq(3,to=ncol(imp[[t]]),by = 2)]
  
  #OBS TD
  data_plot=data_plot_obsTD
  data_plot <- as.data.frame(t(data_plot))
  colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
  data_plot_obsTD2=apply(data_plot,MARGIN = 2,FUN = mean,na.rm=TRUE)
  #PRED TD
  data_plot=data_plot_predTD
  data_plot <- as.data.frame(t(data_plot))
  colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
  data_plot_predTD2=apply(data_plot,MARGIN = 2,FUN = mean,na.rm=TRUE)
  #PRED TD
  data_plot=data_plot_pred
  data_plot <- as.data.frame(t(data_plot))
  colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
  data_plot_pred2=apply(data_plot,MARGIN = 2,FUN = mean,na.rm=TRUE)
  #Obs TD
  data_plot=data_plot_obs
  data_plot <- as.data.frame(t(data_plot))
  colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
  data_plot_obs2=apply(data_plot,MARGIN = 2,FUN = mean,na.rm=TRUE)

  data_plot=rbind(data_plot_obsTD2,data_plot_predTD2,data_plot_pred2)
  
  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
  data <- rbind(rep(.8,11) , rep(0,11) , as.data.frame(data_plot))
  radarchart( data  , axistype=1 , 
              
              #custom polygon
              pcol=colz , pfcol=rgb(0.5,0.5,0.5,0.1) , plwd=4 , 
              
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,.8,.2), cglwd=0.8,
              
              #custom labels
              vlcex=1.5 
  )
#  radarchart(data_now,pcol = colz,axistype = 2,plty = 1,plwd = 3,add=TRUE)
  text(0,0,trait_names[t])
}
dev.off()  



pdf(file=file.path(origin,"_2021","figures","Figure_4",paste0("Figure_4_RF_compare.pdf")),width = 18,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,3))
  for(t in 1:6){
    data_plot_obsTD <- imp_TD[[t]][,seq(2,to=ncol(imp_TD[[t]]),by = 2)]
    data_plot_predTD <-imp_TD[[t]][,seq(3,to=ncol(imp_TD[[t]]),by = 2)]
    data_plot_obs <-   imp[[t]][,seq(2,to=ncol(imp[[t]]),by = 2)]
    data_plot_pred <-  imp[[t]][,seq(3,to=ncol(imp[[t]]),by = 2)]
    
    #OBS TD
    data_plot=data_plot_obsTD
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.8,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
    #PRED TD
    data_plot=data_plot_predTD
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.8,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
    #PRED TD
    data_plot=data_plot_pred
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now <- rbind(rep(.8,10) , rep(0,10) , data_plot)
    radarchart(data_now,pcol = colz,axistype = 2,plty = 1)
    text(0,0,trait_names[t])
  }
  dev.off()  
  
  
  pdf(file=file.path(origin,"_2021","figures","Figure_4",paste0("Figure_4_RF.pdf")),width = 18,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    data_plot_obsTD <- imp_TD[[t]][,seq(2,to=ncol(imp_TD[[t]]),by = 2)]
    data_plot_predTD <-imp_TD[[t]][,seq(3,to=ncol(imp_TD[[t]]),by = 2)]
    data_plot_obs <-   imp[[t]][,seq(2,to=ncol(imp[[t]]),by = 2)]
    data_plot_pred <-  imp[[t]][,seq(3,to=ncol(imp[[t]]),by = 2)]
    
    #OBS TD
    data_plot=data_plot_predTD-data_plot_obsTD
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now1 <- apply(data_plot,2,mean,na.rm=TRUE)
    #PRED TD
    data_plot=data_plot_pred-data_plot_obs
    data_plot <- as.data.frame(t(data_plot))
    colnames(data_plot)[1:4] <- c("Species","Genus","Family","Clade")
    data_now2 <- apply(data_plot,2,mean,na.rm=TRUE)
    data_now <- rbind(rep(.8,10) , rep(0,10) , as.data.frame(rbind(data_now1,data_now2)))
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
  
  