

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
cova_now <- rep(NA,19)
repnums=1:2
GapPercent=80
RepNum=1
TD_choice="Obs_obs_TD"
for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
for(GapPercent in gappercents){
  if(GapPercent!=-1){GapPercent1=GapPercent}
  if(GapPercent==-1){GapPercent1="org"}
  for(RepNum in repnums){
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now2)){
      print(GapPercent1)
      cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
      colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
      
      cova_now <- rbind(cova_now,
                        cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust)),
                        rep(TD_choice,nrow(cova_0_rawclust))))
    }
  }
}
}
  print(dim(cova_now))
  head(cova_now)
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  covaG <- cova_now
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,21)
  repnums=1:2
  GapPercent=80
  TD_choice="Obs_obs_TD"
  for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
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
        colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
        
        cova_now <- rbind(cova_now,
                          cbind(cova_0_rawclust,
                                rep(GapPercent1,nrow(cova_0_rawclust)),
                                rep(TD_choice,nrow(cova_0_rawclust))))
      }
    }
    }
  }
  print(dim(cova_now))
  head(cova_now)
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  
  covaR <- cova_now

  #-------------------------------------
  new_mean_fun <- function(input){
      out=mean(as.numeric(input),na.rm = TRUE)
    return(out)
  }
  
  R_now=covaR
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent","TD_choice"))]);  mode(ag_now) <- "numeric"
  R <- aggregate(x=ag_now,
                     by=list(Group=R_now$Group,Cluster=R_now$Cluster,TD_choice=R_now$TD_choice,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  colnames(R)
  colnames(G)
  
  R_now=covaG
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent","TD_choice"))]);  mode(ag_now) <- "numeric"
  G <- aggregate(x=ag_now,
                 by=list(Group=R_now$Group,Cluster=R_now$Cluster,TD_choice=R_now$TD_choice,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  

  
  ix=4:ncol(R)
  pdf(file=file.path(origin,"figures","figure_4","figure_4_c.pdf"),width=12,height=8)
  par(mfrow=c(2,6))
  t=2
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
    
    head(G)
    
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
    

    
    group_nm="PFT"
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
    
    dat_plot=cbind(R[ix_org,ix[t]],
                   R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                   R[ix_20,ix[t]],R[ix_40,ix[t]],
                   R[ix_50,ix[t]],R[ix_70,ix[t]])
    colnames(dat_plot)=c("observed","1","5","10","20","40","50","70")
    boxplot(dat_plot,col=colz,ylim=c(-10,10),main=paste0(group_nm," data set 2 ",colnames(R)[ix[t]]),las=2)   
    abline(h=0,col="red")
    
  }
  dev.off()
  
  