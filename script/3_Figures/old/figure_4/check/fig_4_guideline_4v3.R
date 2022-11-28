

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
Version_now="V3"
list.files(file.path(origin,"script","analysis",Version_now))
#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions()

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
TD_choices =  out$TD_choices
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido = out$trait_guido
trait_rainfor = out$trait_rainfor

# for 0% gaps and for 70% gaps:

# load the RMSE 
# calculate per cluster an average
# load the Silhouette data per cluster
# load the mean(sd) data per cluster within silhouettes I think
# load the coefficience of variance data per cluster
# load the number of observations somehow
GapPercent1="org"
RepNum=1
TD_choice="Obs_obs_TD"
trait_sub="guido"

colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")
#COVA
  # load the coefficience of variance data per cluster GUIDO org
cova_now <- rep(NA,13)
repnums=1:9
GapPercent=80
TD_choice="Obs_obs_TD"
for(GapPercent in gappercents){
  if(GapPercent!=-1){GapPercent1=GapPercent}
  if(GapPercent==-1){GapPercent1="org"}
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "guido",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now2)){
      print(GapPercent1)
      cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust))))
    }
  }
}
  print(dim(cova_now))
  
  colnames(cova_now)[ncol(cova_now)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  covaG <- cova_now
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,15)
  repnums=1:9
  GapPercent=80
  TD_choice="Obs_obs_TD"
  for(GapPercent in gappercents){
    if(GapPercent!=-1){GapPercent1=GapPercent}
    if(GapPercent==-1){GapPercent1="org"}
    for(RepNum in repnums){
      path_now1 <- file.path(origin,"data_output","CoVa",
                             "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
      path_now2=file.path(origin,"data_output","CoVa",
                          "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
      if(file.exists(path_now2)){
        print(GapPercent1)
        cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
        print(dim(cova_0_rawclust))
        cova_now <- rbind(cova_now,cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust))))
      }
    }
  }
  print(dim(cova_now))
  
  colnames(cova_now)[ncol(cova_now)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  
  covaR <- cova_now
  
  #-------------------------------------
  new.mean_fun <- function(input){
      out=mean(as.numeric(input),na.rm = TRUE)
    return(out)
  }
  
  R_now=covaR
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R <- aggregate(x=ag_now,
                     by=list(Group=R_now$Group,Cluster=R_now$Cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=covaG
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  G <- aggregate(x=ag_now,
                 by=list(Group=R_now$Group,Cluster=R_now$Cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)


  
  pdf(file=file.path(origin,"plots","figure_4","Attractor_OR_Example_V2.pdf"),width=12,height=4)
  par(mfrow=c(2,5))
  ix=4:(ncol(R)/2)
  for(t in 1:length(ix)){
    
    group_nm="Species"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    gc=6
    dat_dist=cbind(R[ix_1,ix[t]]-R[ix_1,ix[t]+gc],
                   R[ix_5,ix[t]]-R[ix_5,ix[t]+gc],
                   R[ix_10,ix[t]]-R[ix_10,ix[t]+gc],
                   R[ix_20,ix[t]]-R[ix_20,ix[t]+gc],
                   R[ix_40,ix[t]]-R[ix_40,ix[t]+gc],
                   R[ix_50,ix[t]]-R[ix_50,ix[t]+gc],
                   R[ix_70,ix[t]]-R[ix_70,ix[t]+gc])
    
    dat_obs=cbind(R[ix_org,ix[t]+gc],
                  R[ix_0,ix[t]+gc],
                  R[ix_1,ix[t]+gc],
                  R[ix_5,ix[t]+gc],
                  R[ix_10,ix[t]+gc],
                  R[ix_20,ix[t]+gc],
                  R[ix_30,ix[t]+gc],
                  R[ix_40,ix[t]+gc],
                  R[ix_50,ix[t]+gc],
                  R[ix_60,ix[t]+gc],
                  R[ix_70,ix[t]+gc])
    
    dat_pred=cbind(R[ix_org,ix[t]],
                   R[ix_0,ix[t]],
                   R[ix_1,ix[t]],
                   R[ix_5,ix[t]],
                   R[ix_10,ix[t]],
                   R[ix_20,ix[t]],
                   R[ix_30,ix[t]],
                   R[ix_40,ix[t]],
                   R[ix_50,ix[t]],
                   R[ix_60,ix[t]],
                   R[ix_70,ix[t]])
    
    plot(dat_obs[,1],dat_pred[,1],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main=paste0("Observed ",colnames(R)[ix[t]]))     
    #plot(dat_obs[,2],dat_pred[,2],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="0%")
    plot(dat_obs[,3],dat_pred[,3],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="1%")
    plot(dat_obs[,4],dat_pred[,4],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="5%")
    plot(dat_obs[,5],dat_pred[,5],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="10%")
    plot(dat_obs[,6],dat_pred[,6],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="20%")
    plot(dat_obs[,7],dat_pred[,7],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="30%")
    plot(dat_obs[,8],dat_pred[,8],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="40%")
    plot(dat_obs[,9],dat_pred[,9],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="50%")
    plot(dat_obs[,10],dat_pred[,10],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="60%")
    plot(dat_obs[,11],dat_pred[,11],xlim=c(-5,5),ylim=c(-5,5),pch=16,col=colz1[8],xlab="Observed cova",ylab="Predicted cova",main="70%")
    
  }
  dev.off()
  
  
  
  
  pdf(file=file.path(origin,"plots","figure_4","Attractor_OR_Example.pdf"),width=8,height=4)
  par(mfrow=c(1,2))
  t=1
  
  {
    colnames(G)
    ix=4:(ncol(G)/2)
    for(t in 1:length(ix)){

    group_nm="Species"
    ix_org=G$Group==group_nm&G$Gaps=="org"
    ix_0=G$Group==group_nm&G$Gaps=="0"
    ix_1=G$Group==group_nm&G$Gaps=="1"
    ix_5=G$Group==group_nm&G$Gaps=="5"
    ix_10=G$Group==group_nm&G$Gaps=="10"
    ix_20=G$Group==group_nm&G$Gaps=="20"
    ix_30=G$Group==group_nm&G$Gaps=="30"
    ix_40=G$Group==group_nm&G$Gaps=="40"
    ix_50=G$Group==group_nm&G$Gaps=="50"
    ix_60=G$Group==group_nm&G$Gaps=="60"
    ix_70=G$Group==group_nm&G$Gaps=="70"
    colnames(G)
    gc=5
    dat_plot=cbind(G[ix_1,ix[t]]-G[ix_1,ix[t]+gc],
                   G[ix_5,ix[t]]-G[ix_5,ix[t]+gc],
                   G[ix_10,ix[t]]-G[ix_10,ix[t]+gc],
                   G[ix_20,ix[t]]-G[ix_20,ix[t]+gc],
                   G[ix_40,ix[t]]-G[ix_40,ix[t]+gc],
                   G[ix_50,ix[t]]-G[ix_50,ix[t]+gc],
                   G[ix_70,ix[t]]-G[ix_70,ix[t]+gc])
    
   # boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Distance")   
  #  abline(h=0,col="red")
    
    dat_plot=cbind(G[ix_org,ix[t]+gc],
                   G[ix_1,ix[t]+gc],
                   G[ix_5,ix[t]+gc],
                   G[ix_10,ix[t]+gc],
                   G[ix_40,ix[t]+gc],
                   G[ix_50,ix[t]+gc],
                   G[ix_70,ix[t]+gc])
    boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Observed")   
    abline(h=0,col="red")    
    dat_plot=cbind(G[ix_org,ix[t]],
                   G[ix_1,ix[t]],
                   G[ix_5,ix[t]],
                   G[ix_10,ix[t]],
                   G[ix_40,ix[t]],
                   G[ix_50,ix[t]],
                   G[ix_70,ix[t]])
    
    #colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Predicted")   
    abline(h=0,col="red")
    }
    
    ix=4:(ncol(R)/2)
    for(t in 1:length(ix)){

    group_nm="Species"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    gc=6
    dat_plot=cbind(R[ix_1,ix[t]]-R[ix_1,ix[t]+gc],
                   R[ix_5,ix[t]]-R[ix_5,ix[t]+gc],
                   R[ix_10,ix[t]]-R[ix_10,ix[t]+gc],
                   R[ix_20,ix[t]]-R[ix_20,ix[t]+gc],
                   R[ix_40,ix[t]]-R[ix_40,ix[t]+gc],
                   R[ix_50,ix[t]]-R[ix_50,ix[t]+gc],
                   R[ix_70,ix[t]]-R[ix_70,ix[t]+gc])
    
 #   boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Distance")   
  #  abline(h=0,col="red")
    
    dat_plot=cbind(R[ix_org,ix[t]+gc],
                   R[ix_1,ix[t]+gc],
                   R[ix_5,ix[t]+gc],
                   R[ix_10,ix[t]+gc],
                   R[ix_40,ix[t]+gc],
                   R[ix_50,ix[t]+gc],
                   R[ix_70,ix[t]+gc])
    boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Observed")   
    abline(h=0,col="red")    
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],
                   R[ix_5,ix[t]],
                   R[ix_10,ix[t]],
                   R[ix_40,ix[t]],
                   R[ix_50,ix[t]],
                   R[ix_70,ix[t]])
    
    #colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2,ylab="Predicted")   
    abline(h=0,col="red")
    
    }
    dev.off()
  }
  
    group_nm="Genus"
    ix_org=G$Group==group_nm&G$Gaps=="org"
    ix_0=G$Group==group_nm&G$Gaps=="0"
    ix_1=G$Group==group_nm&G$Gaps=="1"
    ix_5=G$Group==group_nm&G$Gaps=="5"
    ix_10=G$Group==group_nm&G$Gaps=="10"
    ix_20=G$Group==group_nm&G$Gaps=="20"
    ix_30=G$Group==group_nm&G$Gaps=="30"
    ix_40=G$Group==group_nm&G$Gaps=="40"
    ix_50=G$Group==group_nm&G$Gaps=="50"
    ix_60=G$Group==group_nm&G$Gaps=="60"
    ix_70=G$Group==group_nm&G$Gaps=="70"
    dat_plot=cbind(G[ix_org,ix[t]],
                   G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                   G[ix_20,ix[t]],G[ix_40,ix[t]],
                   G[ix_50,ix[t]],G[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
    abline(h=0,col="red")
    


    group_nm="Family"
    ix_org=G$Group==group_nm&G$Gaps=="org"
    ix_0=G$Group==group_nm&G$Gaps=="0"
    ix_1=G$Group==group_nm&G$Gaps=="1"
    ix_5=G$Group==group_nm&G$Gaps=="5"
    ix_10=G$Group==group_nm&G$Gaps=="10"
    ix_20=G$Group==group_nm&G$Gaps=="20"
    ix_30=G$Group==group_nm&G$Gaps=="30"
    ix_40=G$Group==group_nm&G$Gaps=="40"
    ix_50=G$Group==group_nm&G$Gaps=="50"
    ix_60=G$Group==group_nm&G$Gaps=="60"
    ix_70=G$Group==group_nm&G$Gaps=="70"
    dat_plot=cbind(G[ix_org,ix[t]],
                   G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                   G[ix_20,ix[t]],G[ix_40,ix[t]],
                   G[ix_50,ix[t]],G[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    
    group_nm="PG"
    ix_org=G$Group==group_nm&G$Gaps=="org"
    ix_0=G$Group==group_nm&G$Gaps=="0"
    ix_1=G$Group==group_nm&G$Gaps=="1"
    ix_5=G$Group==group_nm&G$Gaps=="5"
    ix_10=G$Group==group_nm&G$Gaps=="10"
    ix_20=G$Group==group_nm&G$Gaps=="20"
    ix_30=G$Group==group_nm&G$Gaps=="30"
    ix_40=G$Group==group_nm&G$Gaps=="40"
    ix_50=G$Group==group_nm&G$Gaps=="50"
    ix_60=G$Group==group_nm&G$Gaps=="60"
    ix_70=G$Group==group_nm&G$Gaps=="70"
    dat_plot=cbind(G[ix_org,ix[t]],
                   G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                   G[ix_20,ix[t]],G[ix_40,ix[t]],
                   G[ix_50,ix[t]],G[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    
    
    group_nm="GF"
    ix_org=G$Group==group_nm&G$Gaps=="org"
    ix_0=G$Group==group_nm&G$Gaps=="0"
    ix_1=G$Group==group_nm&G$Gaps=="1"
    ix_5=G$Group==group_nm&G$Gaps=="5"
    ix_10=G$Group==group_nm&G$Gaps=="10"
    ix_20=G$Group==group_nm&G$Gaps=="20"
    ix_30=G$Group==group_nm&G$Gaps=="30"
    ix_40=G$Group==group_nm&G$Gaps=="40"
    ix_50=G$Group==group_nm&G$Gaps=="50"
    ix_60=G$Group==group_nm&G$Gaps=="60"
    ix_70=G$Group==group_nm&G$Gaps=="70"
    dat_plot=cbind(G[ix_org,ix[t]],
                   G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                   G[ix_20,ix[t]],G[ix_40,ix[t]],
                   G[ix_50,ix[t]],G[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    }
    
    #----------------------------
    
    ix=4:ncol(R)
    for(t in 1:length(ix)){
    group_nm="Species"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    group_nm="Genus"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    group_nm="Family"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    
    
    group_nm="PG"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
    
    group_nm="GF"
    ix_org=R$Group==group_nm&R$Gaps=="org"
    ix_0=R$Group==group_nm&R$Gaps=="0"
    ix_1=R$Group==group_nm&R$Gaps=="1"
    ix_5=R$Group==group_nm&R$Gaps=="5"
    ix_10=R$Group==group_nm&R$Gaps=="10"
    ix_20=R$Group==group_nm&R$Gaps=="20"
    ix_30=R$Group==group_nm&R$Gaps=="30"
    ix_40=R$Group==group_nm&R$Gaps=="40"
    ix_50=R$Group==group_nm&R$Gaps=="50"
    ix_60=R$Group==group_nm&R$Gaps=="60"
    ix_70=R$Group==group_nm&R$Gaps=="70"
    
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    }
    dev.off()
  }
  
    
  pdf(file=file.path(origin,"plots","figure_4","Attractor_SLA.pdf"),width=10,height=10)
  par(mfrow=c(3,4))
  t=1
  {
  group_nm="Species"
  ix_org=G$Group==group_nm&G$Gaps=="org"
  ix_0=G$Group==group_nm&G$Gaps=="0"
  ix_1=G$Group==group_nm&G$Gaps=="1"
  ix_5=G$Group==group_nm&G$Gaps=="5"
  ix_10=G$Group==group_nm&G$Gaps=="10"
  ix_20=G$Group==group_nm&G$Gaps=="20"
  ix_30=G$Group==group_nm&G$Gaps=="30"
  ix_40=G$Group==group_nm&G$Gaps=="40"
  ix_50=G$Group==group_nm&G$Gaps=="50"
  ix_60=G$Group==group_nm&G$Gaps=="60"
  ix_70=G$Group==group_nm&G$Gaps=="70"
  dat_plot=cbind(G[ix_org,ix[t]],
                 G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                 G[ix_20,ix[t]],G[ix_40,ix[t]],
                 G[ix_50,ix[t]],G[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  group_nm="Species"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  dat_plot=cbind(R[ix_org,ix[t]],
                 R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                 R[ix_20,ix[t]],R[ix_40,ix[t]],
                 R[ix_50,ix[t]],R[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  
  
  
  group_nm="Genus"
  ix_org=G$Group==group_nm&G$Gaps=="org"
  ix_0=G$Group==group_nm&G$Gaps=="0"
  ix_1=G$Group==group_nm&G$Gaps=="1"
  ix_5=G$Group==group_nm&G$Gaps=="5"
  ix_10=G$Group==group_nm&G$Gaps=="10"
  ix_20=G$Group==group_nm&G$Gaps=="20"
  ix_30=G$Group==group_nm&G$Gaps=="30"
  ix_40=G$Group==group_nm&G$Gaps=="40"
  ix_50=G$Group==group_nm&G$Gaps=="50"
  ix_60=G$Group==group_nm&G$Gaps=="60"
  ix_70=G$Group==group_nm&G$Gaps=="70"
  dat_plot=cbind(G[ix_org,ix[t]],
                 G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                 G[ix_20,ix[t]],G[ix_40,ix[t]],
                 G[ix_50,ix[t]],G[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz1,ylim=c(-3,3),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  group_nm="Genus"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  dat_plot=cbind(R[ix_org,ix[t]],
                 R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                 R[ix_20,ix[t]],R[ix_40,ix[t]],
                 R[ix_50,ix[t]],R[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz,ylim=c(-3,3),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  
  group_nm="Family"
  ix_org=G$Group==group_nm&G$Gaps=="org"
  ix_0=G$Group==group_nm&G$Gaps=="0"
  ix_1=G$Group==group_nm&G$Gaps=="1"
  ix_5=G$Group==group_nm&G$Gaps=="5"
  ix_10=G$Group==group_nm&G$Gaps=="10"
  ix_20=G$Group==group_nm&G$Gaps=="20"
  ix_30=G$Group==group_nm&G$Gaps=="30"
  ix_40=G$Group==group_nm&G$Gaps=="40"
  ix_50=G$Group==group_nm&G$Gaps=="50"
  ix_60=G$Group==group_nm&G$Gaps=="60"
  ix_70=G$Group==group_nm&G$Gaps=="70"
  dat_plot=cbind(G[ix_org,ix[t]],
                 G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                 G[ix_20,ix[t]],G[ix_40,ix[t]],
                 G[ix_50,ix[t]],G[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
  abline(h=0,col="red")

  group_nm="Family"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  dat_plot=cbind(R[ix_org,ix[t]],
                 R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                 R[ix_20,ix[t]],R[ix_40,ix[t]],
                 R[ix_50,ix[t]],R[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
  abline(h=0,col="red")
  

  group_nm="PG"
  ix_org=G$Group==group_nm&G$Gaps=="org"
  ix_0=G$Group==group_nm&G$Gaps=="0"
  ix_1=G$Group==group_nm&G$Gaps=="1"
  ix_5=G$Group==group_nm&G$Gaps=="5"
  ix_10=G$Group==group_nm&G$Gaps=="10"
  ix_20=G$Group==group_nm&G$Gaps=="20"
  ix_30=G$Group==group_nm&G$Gaps=="30"
  ix_40=G$Group==group_nm&G$Gaps=="40"
  ix_50=G$Group==group_nm&G$Gaps=="50"
  ix_60=G$Group==group_nm&G$Gaps=="60"
  ix_70=G$Group==group_nm&G$Gaps=="70"
  dat_plot=cbind(G[ix_org,ix[t]],
                 G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                 G[ix_20,ix[t]],G[ix_40,ix[t]],
                 G[ix_50,ix[t]],G[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  group_nm="PG"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  
  dat_plot=cbind(R[ix_org,ix[t]],
                 R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                 R[ix_20,ix[t]],R[ix_40,ix[t]],
                 R[ix_50,ix[t]],R[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  
  group_nm="GF"
  ix_org=G$Group==group_nm&G$Gaps=="org"
  ix_0=G$Group==group_nm&G$Gaps=="0"
  ix_1=G$Group==group_nm&G$Gaps=="1"
  ix_5=G$Group==group_nm&G$Gaps=="5"
  ix_10=G$Group==group_nm&G$Gaps=="10"
  ix_20=G$Group==group_nm&G$Gaps=="20"
  ix_30=G$Group==group_nm&G$Gaps=="30"
  ix_40=G$Group==group_nm&G$Gaps=="40"
  ix_50=G$Group==group_nm&G$Gaps=="50"
  ix_60=G$Group==group_nm&G$Gaps=="60"
  ix_70=G$Group==group_nm&G$Gaps=="70"
  dat_plot=cbind(G[ix_org,ix[t]],
                 G[ix_1,ix[t]],G[ix_5,ix[t]],G[ix_10,ix[t]],
                 G[ix_20,ix[t]],G[ix_40,ix[t]],
                 G[ix_50,ix[t]],G[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz1,ylim=c(-10,10),main=paste0(group_nm," data set 1 ",colnames(G)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  
  group_nm="GF"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  
  dat_plot=cbind(R[ix_org,ix[t]],
                 R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                 R[ix_20,ix[t]],R[ix_40,ix[t]],
                 R[ix_50,ix[t]],R[ix_70,ix[t]])
  colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
  boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
  abline(h=0,col="red")
  
  dev.off()
  }

dev.off()

  group_nm="PFT"
  ix_org=R$Group==group_nm&R$Gaps=="org"
  ix_0=R$Group==group_nm&R$Gaps=="0"
  ix_1=R$Group==group_nm&R$Gaps=="1"
  ix_5=R$Group==group_nm&R$Gaps=="5"
  ix_10=R$Group==group_nm&R$Gaps=="10"
  ix_20=R$Group==group_nm&R$Gaps=="20"
  ix_30=R$Group==group_nm&R$Gaps=="30"
  ix_40=R$Group==group_nm&R$Gaps=="40"
  ix_50=R$Group==group_nm&R$Gaps=="50"
  ix_60=R$Group==group_nm&R$Gaps=="60"
  ix_70=R$Group==group_nm&R$Gaps=="70"
  
  boxplot(cbind(R[ix_org,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_40,ix[t]],
                R[ix_50,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=,xaxt="n")   
  axis(1,las=2,at=c(1:10),labels = c("observed","1","5","10","20","40","50","70"))
  abline(h=0,col="red")  
  
  

  
  
  
  
  
  
  
  
  
  pdf(file=file.path(origin,"plots","figure_4","Cova_Towards0"),width=25,height=8)
    boxplot(cbind(Rorg_specC[,t],R0_specC[,t],R1_specC[,t],R5_specC[,t],
               R10_specC[,t],R20_specC[,t],R30_specC[,t],R40_specC[,t],
               R50_specC[,t],R70_specC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    boxplot(cbind(Rorg_genC[,t],R0_genC[,t],R1_genC[,t],R5_genC[,t],
                  R10_genC[,t],R20_genC[,t],R30_genC[,t],R40_genC[,t],
                  R50_genC[,t],R70_genC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    boxplot(cbind(Rorg_famC[,t],R0_famC[,t],R1_famC[,t],R5_famC[,t],
                  R10_famC[,t],R20_famC[,t],R30_famC[,t],R40_famC[,t],
                  R50_famC[,t],R70_famC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    boxplot(cbind(Rorg_pgC[,t],R0_pgC[,t],R1_pgC[,t],R5_pgC[,t],
                  R10_pgC[,t],R20_pgC[,t],R30_pgC[,t],R40_pgC[,t],
                  R50_pgC[,t],R70_pgC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    boxplot(cbind(Rorg_gfC[,t],R0_gfC[,t],R1_gfC[,t],R5_gfC[,t],
                  R10_gfC[,t],R20_gfC[,t],R30_gfC[,t],R40_gfC[,t],
                  R50_gfC[,t],R70_gfC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    boxplot(cbind(Rorg_pftC[,t],R0_pftC[,t],R1_pftC[,t],R5_pftC[,t],
                  R10_pftC[,t],R20_pftC[,t],R30_pftC[,t],R40_pftC[,t],
                  R50_pftC[,t],R70_pftC[,t]),ylimc=c(-3,3),col=colz)
    abline(h=0)
    dev.off()
  
  pdf(file=file.path(origin,"plots","test_arrow.pdf"),width=25,height=8)
  plot(1:nrow(Rorg_spec),col="white",xaxt="n",ylab="",ylim=c(-50,50))
  i=1
  for(i in 1:nrow(Rorg_spec)){
    dat_plot=c(Rorg_specC[i,t],R0_specC[i,t],R1_specC[i,t],R1_specC[i,t],
               R1_specC[i,t],R1_specC[i,t],R1_specC[i,t],R1_specC[i,t])
    if(sum(!is.na(dat_plot))>3){
      points(rep(i,length(dat_plot)),dat_plot,col=colz[2:length(dat_plot)])
      lines(x=rep(i,length(dat_plot)),y=dat_plot,col="black")
      arrows(i, dat_plot[1], x1 = i, y1 = dat_plot[length(dat_plot)],length=.08)
    }
  }
  abline(h=0)
  dev.off()
  