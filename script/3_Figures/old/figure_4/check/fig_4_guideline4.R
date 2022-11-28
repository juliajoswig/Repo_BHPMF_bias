

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

  #COVA
  # load the coefficience of variance data per cluster GUIDO org
  cova_now <- rep(NA,8)
  repnums=1:5
  GapPercent1="org"
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
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust))))
    }
  }
  }

  names(cova_now)[ncol(cova_now)] <-"GapPercent" 

  covaG_spec <- cova_now[cova_now$level=="Species",2:ncol(cova_now)]
  covaG_gen <- cova_now[cova_now$level=="Genus",2:ncol(cova_now)]
  covaG_fam <- cova_now[cova_now$level=="Family",2:ncol(cova_now)]
  covaG_pg <- cova_now[cova_now$level=="PG",2:ncol(cova_now)]
  covaG_gf <- cova_now[cova_now$level=="GF",2:ncol(cova_now)]
  covaG_pft <- cova_now[cova_now$level=="PFT",2:ncol(cova_now)]
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,8)
  repnums=1:5
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
  
  covaR_spec <- cova_now[cova_now$level=="Species",2:ncol(cova_now)]
  covaR_gen <- cova_now[cova_now$level=="Genus",2:ncol(cova_now)]
  covaR_fam <- cova_now[cova_now$level=="Family",2:ncol(cova_now)]
  covaR_pg <- cova_now[cova_now$level=="PG",2:ncol(cova_now)]
  covaR_gf <- cova_now[cova_now$level=="GF",2:ncol(cova_now)]
  covaR_pft <- cova_now[cova_now$level=="PFT",2:ncol(cova_now)]
  
  covaR <- cova_now
  
  #All
  Rorg <- covaR[covaR$GapPercent=="org",]
  R0 <- covaR[covaR$GapPercent=="0",]
  R1 <- covaR[covaR$GapPercent=="1",]
  R5 <- covaR[covaR$GapPercent=="5",]
  R10 <- covaR[covaR$GapPercent=="10",]
  R20 <- covaR[covaR$GapPercent=="20",]
  R30 <- covaR[covaR$GapPercent=="30",]
  R40 <- covaR[covaR$GapPercent=="40",]
  R50 <- covaR[covaR$GapPercent=="50",]
  R60 <- covaR[covaR$GapPercent=="60",]
  R70 <- covaR[covaR$GapPercent=="70",]
 
  new.mean_fun <- function(input){
      out=mean(as.numeric(input),na.rm = TRUE)
    return(out)
  }
  
  R_now=covaR
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R <- aggregate(x=ag_now,
                     by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=Rorg
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  RorgA <- aggregate(x=ag_now,
                     by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R0
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R0A <- aggregate(x=ag_now,
                     by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R1
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R1A <- aggregate(x=ag_now,
                     by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R5
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R5A <- aggregate(x=ag_now,
                   by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  R_now=R10
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R10A <- aggregate(x=ag_now,
                   by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  R_now=R20
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R20A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R30
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R30A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R40
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R40A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R_now=R50
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R50A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  
  R_now=R60
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R60A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  
  R_now=R70
  ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("level","cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
  R70A <- aggregate(x=ag_now,
                    by=list(Group=R_now$level,cluster=R_now$cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)
  
  R0A$cluster==RorgA$cluster
  colnames(RorgA)
  ix=4:ncol(R0A)
  R0D <- R0A[,ix]-RorgA[,ix]
  R1D <- R1A[,ix]-RorgA[,ix]
  R5D <- R5A[,ix]-RorgA[,ix]
  R10D <- R10A[,ix]-RorgA[,ix]
  R20D <- R20A[,ix]-RorgA[,ix]
  R30D <- R30A[,ix]-RorgA[,ix]
  R40D <- R40A[,ix]-RorgA[,ix]
  R50D <- R50A[,ix]-RorgA[,ix]
  R60D <- R60A[,ix]-RorgA[,ix]
  R70D <- R70A[,ix]-RorgA[,ix]
  
  
  t=1
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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
  R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
  R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
  R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")
  
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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
                R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")  
  
  
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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
                R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")  

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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
                R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")  
  
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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
                R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")  
  
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
  
  boxplot(cbind(R[ix_org,ix[t]],R[ix_0,ix[t]],
                R[ix_1,ix[t]],R[ix_5,ix[t]],R[ix_10,ix[t]],
                R[ix_20,ix[t]],R[ix_30,ix[t]],R[ix_40,ix[t]],
                R[ix_60,ix[t]],R[ix_70,ix[t]]),col=colz,ylim=c(-10,10),main=group_nm)
  abline(h=0,col="gray")  
  
  

  
  
  
  
  
  
  
  
  
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
  