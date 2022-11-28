

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

{
#Load RMSE per cluster GUIDO 0% gaps
  sil_now <- rep(NA,12)
  repnums=1:10
  GapPercent1=0
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","RMSE","guido",
                              TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2= file.path(origin,"data_output","RMSE","guido",
                             TD_choice,GapPercent1,RepNum,paste0("RMSE_all_PlusCluster.csv"))
      if(file.exists(path_now1)&file.exists(path_now2)){
      rmse_0_general <- read.csv(file=path_now1)
      rmse_0_rawclust <- read.csv(file=path_now2)
    }
    print(dim(rmse_0_rawclust))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  
  # calculate per cluster an average
  rmseG_0_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
  rmseG_0_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
  rmseG_0_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
  rmseG_0_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
  rmseG_0_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
  rmseG_0_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)

  #Load RMSE per cluster RAINFOR 0% gaps
  sil_now <- rep(NA,12)
  repnums=1:10
  GapPercent1=0
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","RMSE","rainfor",
                           TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2= file.path(origin,"data_output","RMSE","rainfor",
                         TD_choice,GapPercent1,RepNum,paste0("RMSE_all_PlusCluster.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      rmse_0_general <- read.csv(file=path_now1)
      rmse_0_rawclust <- read.csv(file=path_now2)
    }
    print(dim(rmse_0_rawclust))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  
  # calculate per cluster an average
  rmseR_0_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
  rmseR_0_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
  rmseR_0_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
  rmseR_0_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
  rmseR_0_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
  rmseR_0_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)
  
  #70% gaps
  sil_now <- rep(NA,12)
  repnums=1:10
  GapPercent1=70
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","RMSE","guido",
                           TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2= file.path(origin,"data_output","RMSE","guido",
                         TD_choice,GapPercent1,RepNum,paste0("RMSE_all_PlusCluster.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      rmse_0_general <- read.csv(file=path_now1)
      rmse_0_rawclust <- read.csv(file=path_now2)
    }
    print(dim(rmse_0_rawclust))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  
  # calculate per cluster an average
  rmseG_70_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
  rmseG_70_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
  rmseG_70_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
  rmseG_70_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
  rmseG_70_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
  rmseG_70_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)

  #70% gaps
  sil_now <- rep(NA,12)
  repnums=1:10
  GapPercent1=70
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","RMSE","rainfor",
                           TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2= file.path(origin,"data_output","RMSE","rainfor",
                         TD_choice,GapPercent1,RepNum,paste0("RMSE_all_PlusCluster.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      rmse_0_general <- read.csv(file=path_now1)
      rmse_0_rawclust <- read.csv(file=path_now2)
    }
    print(dim(rmse_0_rawclust))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  
  # calculate per cluster an average
  rmseR_70_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
  rmseR_70_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
  rmseR_70_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
  rmseR_70_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
  rmseR_70_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
  rmseR_70_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)
  
  
# load the Silhouette data per cluster
# load the mean(sd) data per cluster within silhouettes I think
# load the number of observations somehowf
sil_now <- rep(NA,6)
repnums=1:3
for(RepNum in repnums){
  path_now=file.path(origin,"data_output","Sil","guido",
                     TD_choice,GapPercent1,
                     RepNum,paste0("all_clusters.csv"))
  if(file.exists(path_now)){sil_nowtmp <- read.csv(file=path_now,row.names = NULL)}
  head(sil_nowtmp)
  print(dim(sil_nowtmp))
 sil_now <- rbind(sil_now,sil_nowtmp)
}

  silG_0_spec <- sil_now[sil_now[,2]=="Species",3:ncol(sil_now)]
  silG_0_gen <- sil_now[sil_now[,2]=="Genus",3:ncol(sil_now)]
  silG_0_fam <- sil_now[sil_now[,2]=="Family",3:ncol(sil_now)]
  silG_0_pg <- sil_now[sil_now[,2]=="PG",3:ncol(sil_now)]
  silG_0_gf <- sil_now[sil_now[,2]=="GF",3:ncol(sil_now)]
  silG_0_pft <- sil_now[sil_now[,2]=="PFT",3:ncol(sil_now)]

  sil_now <- rep(NA,6)
  repnums=1:3
  for(RepNum in repnums){
    path_now = file.path(origin,"data_output","Sil","rainfor",
               TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))    
    if(file.exists(path_now)){
      sil_nowtmp <- read.csv(file=path_now,row.names = NULL)}
    head(sil_nowtmp)
    print(dim(sil_nowtmp))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  silR_0_spec <- sil_now[sil_now[,2]=="Species",3:ncol(sil_now)]
  silR_0_gen <- sil_now[sil_now[,2]=="Genus",3:ncol(sil_now)]
  silR_0_fam <- sil_now[sil_now[,2]=="Family",3:ncol(sil_now)]
  silR_0_pg <- sil_now[sil_now[,2]=="PG",3:ncol(sil_now)]
  silR_0_gf <- sil_now[sil_now[,2]=="GF",3:ncol(sil_now)]
  silR_0_pft <- sil_now[sil_now[,2]=="PFT",3:ncol(sil_now)]

  sil_now <- rep(NA,6)
  repnums=1:3
  for(RepNum in repnums){
    path_now = file.path(origin,"data_output","Sil","rainfor",
                         TD_choice,70,RepNum,paste0("all_clusters.csv"))   
    if(file.exists(path_now)){
      sil_nowtmp <- read.csv(file=path_now,row.names = NULL)}
    head(sil_nowtmp)
    print(dim(sil_nowtmp))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  
  silR_70_spec <- sil_now[sil_now[,2]=="Species",3:ncol(sil_now)]
  silR_70_gen <- sil_now[sil_now[,2]=="Genus",3:ncol(sil_now)]
  silR_70_fam <- sil_now[sil_now[,2]=="Family",3:ncol(sil_now)]
  silR_70_pg <- sil_now[sil_now[,2]=="PG",3:ncol(sil_now)]
  silR_70_gf <- sil_now[sil_now[,2]=="GF",3:ncol(sil_now)]
  silR_70_pft <- sil_now[sil_now[,2]=="PFT",3:ncol(sil_now)]

  sil_now <- rep(NA,6)
  repnums=1:3
  for(RepNum in repnums){
    path_now = file.path(origin,"data_output","Sil","guido",
                         TD_choice,70,RepNum,paste0("all_clusters.csv"))   
    if(file.exists(path_now)){
      sil_nowtmp <- read.csv(file=path_now,row.names = NULL)}
    head(sil_nowtmp)
    print(dim(sil_nowtmp))
    sil_now <- rbind(sil_now,sil_nowtmp)
  }
  silG_70_spec <- sil_now[sil_now[,2]=="Species",3:ncol(sil_now)]
  silG_70_gen <- sil_now[sil_now[,2]=="Genus",3:ncol(sil_now)]
  silG_70_fam <- sil_now[sil_now[,2]=="Family",3:ncol(sil_now)]
  silG_70_pg <- sil_now[sil_now[,2]=="PG",3:ncol(sil_now)]
  silG_70_gf <- sil_now[sil_now[,2]=="GF",3:ncol(sil_now)]
  silG_70_pft <- sil_now[sil_now[,2]=="PFT",3:ncol(sil_now)]

}



pdf(file=file.path(origin,"plots","figure_4","Silhouette_depends_on_orgSil.pdf"),
    width=8,height = 6)
{
  par(mfrow=c(3,4))
  plot(c(silR_0_spec[,3],silG_0_spec[,3]),
       c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
       pch=16,main="Species",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_spec[,2],silG_0_spec[,2]),
       c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
       pch=16,main="Species",ylim=c(-1.5,1.5),cex=.5,log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_70_spec[,2],silG_70_spec[,2]),
       c(silR_70_spec[,4],silG_70_spec[,4])+2,
       pch=16,main="Species",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_spec[,4],silG_70_spec[,4]),
       c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
       pch=16,main="Species",cex=.8,ylim=c(0,1.5),xlim = c(.3,.8),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  #------------------------------------
  # single data sets
  
  
  
  #---SPECIES
  x <- silR_70_spec$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('thistle1','red'))
  colzR <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsR <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsR <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  x <- silG_70_spec$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('lightblue2','blue'))
  colzG <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsG <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsG <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]

  rmseR_0_specM<- rmseR_0_spec[match(silR_0_spec$Cluster,rmseR_0_spec$Group.1),]
  rmseG_0_specM<- rmseG_0_spec[match(silG_0_spec$Cluster,rmseG_0_spec$Group.1),]
  rmseR_70_specM<- rmseR_70_spec[match(silR_70_spec$Cluster,rmseR_70_spec$Group.1),]
  rmseG_70_specM<- rmseG_70_spec[match(silG_70_spec$Cluster,rmseG_70_spec$Group.1),]
  dat_plot <- cbind(
    c(silR_0_spec[,2],silG_0_spec[,2]),
    log(c(silR_70_spec[,4],silG_70_spec[,4])+2),
    c(silR_0_spec[,3],silG_0_spec[,3]),
    c(silR_70_spec[,3],silG_70_spec[,3]),
    c((silR_70_spec[,3]-silR_0_spec[,3]),(silG_70_spec[,3]-silG_0_spec[,3])),
    c(rmseR_0_specM[,2],rmseG_0_specM[,2]),
    c(rmseR_0_specM[,3],rmseG_0_specM[,3]),
    c(rowMeans(rmseR_0_specM[,2:ncol(rmseR_0_specM)]),rowMeans(rmseG_0_specM[,2:ncol(rmseG_0_specM)])),
    c(rowMeans(rmseR_70_specM[,2:ncol(rmseR_70_specM)]),rowMeans(rmseG_70_specM[,2:ncol(rmseG_70_specM)]))
  )
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed",
                          "Silhouette_Change","RMSE_SLA","RMSE_Height","RMSE_0","RMSE70")
  dat_plot_spec <- dat_plot
  
  #---GENUS
  x <- silR_70_gen$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('thistle1','red'))
  colzR_gen <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsR_gen <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsR_gen<- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  x <- silG_70_gen$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('lightblue2','blue'))
  colzG_gen <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsG_gen <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsG_gen <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]

  rmseR_0_genM<- rmseR_0_gen[match(silR_0_gen$Cluster,rmseR_0_gen$Group.1),]
  rmseG_0_genM<- rmseG_0_gen[match(silG_0_gen$Cluster,rmseG_0_gen$Group.1),]
  rmseR_70_genM<- rmseR_70_gen[match(silR_70_gen$Cluster,rmseR_70_gen$Group.1),]
  rmseG_70_genM<- rmseG_70_gen[match(silG_70_gen$Cluster,rmseG_70_gen$Group.1),]
  
  dat_plot <- cbind(
    c(silR_0_gen[,2],silG_0_gen[,2]),
    log(c(silR_70_gen[,4],silG_70_gen[,4])+2),
    c(silR_0_gen[,3],silG_0_gen[,3]),
    c(silR_70_gen[,3],silG_70_gen[,3]),
    c((silR_70_gen[,3]-silR_0_gen[,3]),(silG_70_gen[,3]-silG_0_gen[,3])),
    c(rmseR_0_genM[,2],rmseG_0_genM[,2]),
    c(rmseR_0_genM[,3],rmseG_0_genM[,3]),
    c(rowMeans(rmseR_0_genM[,2:ncol(rmseR_0_genM)]),rowMeans(rmseG_0_genM[,2:ncol(rmseG_0_genM)])),
    c(rowMeans(rmseR_70_genM[,2:ncol(rmseR_70_genM)]),rowMeans(rmseG_70_genM[,2:ncol(rmseG_70_genM)]))
  )
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed",
                          "Silhouette_Change","RMSE_SLA","RMSE_Height","RMSE_0","RMSE70")
  dat_plot_gen <- dat_plot
  
  #---Family
  x <- silR_70_fam$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('thistle1','red'))
  colzR_fam <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsR_fam <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsR_fam<- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  x <- silG_70_fam$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('lightblue2','blue'))
  colzG_fam <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsG_fam <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsG_fam <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  
  rmseR_0_famM<- rmseR_0_fam[match(silR_0_fam$Cluster,rmseR_0_fam$Group.1),]
  rmseG_0_famM<- rmseG_0_fam[match(silG_0_fam$Cluster,rmseG_0_fam$Group.1),]
  rmseR_70_famM<- rmseR_70_fam[match(silR_70_fam$Cluster,rmseR_70_fam$Group.1),]
  rmseG_70_famM<- rmseG_70_fam[match(silG_70_fam$Cluster,rmseG_70_fam$Group.1),]
  
  dat_plot <- cbind(
    c(silR_0_fam[,2],silG_0_fam[,2]),
    log(c(silR_70_fam[,4],silG_70_fam[,4])+2),
    c(silR_0_fam[,3],silG_0_fam[,3]),
    c(silR_70_fam[,3],silG_70_fam[,3]),
    c((silR_70_fam[,3]-silR_0_fam[,3]),(silG_70_fam[,3]-silG_0_fam[,3])),
    c(rmseR_0_famM[,2],rmseG_0_famM[,2]),
    c(rmseR_0_famM[,3],rmseG_0_famM[,3]),
    c(rowMeans(rmseR_0_famM[,2:ncol(rmseR_0_famM)]),rowMeans(rmseG_0_famM[,2:ncol(rmseG_0_famM)])),
    c(rowMeans(rmseR_70_famM[,2:ncol(rmseR_70_famM)]),rowMeans(rmseG_70_famM[,2:ncol(rmseG_70_famM)]))
  )
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed",
                          "Silhouette_Change","RMSE_SLA","RMSE_Height","RMSE_0","RMSE70")
  dat_plot_fam <- dat_plot
  #---PG
  x <- silR_70_pg$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('thistle1','red'))
  colzR_pg <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsR_pg <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsR_pg<- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  x <- silG_70_pg$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('lightblue2','blue'))
  colzG_pg <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsG_pg <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsG_pg <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  
  rmseR_0_pgM<- rmseR_0_pg[match(silR_0_pg$Cluster,rmseR_0_pg$Group.1),]
  rmseG_0_pgM<- rmseG_0_pg[match(silG_0_pg$Cluster,rmseG_0_pg$Group.1),]
  rmseR_70_pgM<- rmseR_70_pg[match(silR_70_pg$Cluster,rmseR_70_pg$Group.1),]
  rmseG_70_pgM<- rmseG_70_pg[match(silG_70_pg$Cluster,rmseG_70_pg$Group.1),]
  
  dat_plot <- cbind(
    c(silR_0_pg[,2],silG_0_pg[,2]),
    log(c(silR_70_pg[,4],silG_70_pg[,4])+2),
    c(silR_0_pg[,3],silG_0_pg[,3]),
    c(silR_70_pg[,3],silG_70_pg[,3]),
    c((silR_70_pg[,3]-silR_0_pg[,3]),(silG_70_pg[,3]-silG_0_pg[,3])),
    c(rmseR_0_pgM[,2],rmseG_0_pgM[,2]),
    c(rmseR_0_pgM[,3],rmseG_0_pgM[,3]),
    c(rowMeans(rmseR_0_pgM[,2:ncol(rmseR_0_pgM)]),rowMeans(rmseG_0_pgM[,2:ncol(rmseG_0_pgM)])),
    c(rowMeans(rmseR_70_pgM[,2:ncol(rmseR_70_pgM)]),rowMeans(rmseG_70_pgM[,2:ncol(rmseG_70_pgM)]))
  )
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed",
                          "Silhouette_Change","RMSE_SLA","RMSE_Height","RMSE_0","RMSE70")
  dat_plot_pg <- dat_plot
  
  #---------------  
  #---Growth Form
  x <- silR_70_gf$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('thistle1','red'))
  colzR_gf <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsR_gf <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsR_gf<- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]
  x <- silG_70_gf$Cluster_size
  dat <- data.frame(x = x,y = x)
  rbPal <- colorRampPalette(c('lightblue2','blue'))
  colzG_gf <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
  pchsG_gf <- c(1,6,20,16,15)[as.numeric(cut(dat$y,breaks = 4))]
  cexsG_gf <- c(.4,.7,.8,.9,1)[as.numeric(cut(dat$y,breaks = 3))]

  rmseR_0_gfM<- rmseR_0_gf[match(silR_0_gf$Cluster,rmseR_0_gf$Group.1),]
  rmseG_0_gfM<- rmseG_0_gf[match(silG_0_gf$Cluster,rmseG_0_gf$Group.1),]
  rmseR_70_gfM<- rmseR_70_gf[match(silR_70_gf$Cluster,rmseR_70_gf$Group.1),]
  rmseG_70_gfM<- rmseG_70_gf[match(silG_70_gf$Cluster,rmseG_70_gf$Group.1),]
  
  dat_plot <- cbind(
    c(silR_0_gf[,2],silG_0_gf[,2]),
    log(c(silR_70_gf[,4],silG_70_gf[,4])+2),
    c(silR_0_gf[,3],silG_0_gf[,3]),
    c(silR_70_gf[,3],silG_70_gf[,3]),
    c((silR_70_gf[,3]-silR_0_gf[,3]),(silG_70_gf[,3]-silG_0_gf[,3])),
    c(rmseR_0_gfM[,2],rmseG_0_gfM[,2]),
    c(rmseR_0_gfM[,3],rmseG_0_gfM[,3]),
    c(rowMeans(rmseR_0_gfM[,2:ncol(rmseR_0_gfM)]),rowMeans(rmseG_0_gfM[,2:ncol(rmseG_0_gfM)])),
    c(rowMeans(rmseR_70_gfM[,2:ncol(rmseR_70_gfM)]),rowMeans(rmseG_70_gfM[,2:ncol(rmseG_70_gfM)]))
  )
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed",
                          "Silhouette_Change","RMSE_SLA","RMSE_Height","RMSE_0","RMSE70")
  dat_plot_gf <- dat_plot
  
  #---------------
  # DATA frame transform
  dat_plot_spec <- data.frame(dat_plot_spec)
  dat_plot_gen <- data.frame(dat_plot_gen)
  dat_plot_fam <- data.frame(dat_plot_fam)
  dat_plot_pg <- data.frame(dat_plot_pg)
  dat_plot_gf <- data.frame(dat_plot_gf)
  #---------------
  
  pairs(dat_plot_spec,col=c(colzR,colzG))
  pairs(dat_plot_gen)
  pairs(dat_plot_fam)
  pairs(dat_plot_pg)
  pairs(dat_plot_gf,col=c(colzR_gf,colzG_gf))
  
  require(FactoMineR)
  PCA(dat_plot_spec[complete.cases(dat_plot_spec),])
  PCA(dat_plot_gen[complete.cases(dat_plot_gen),])
  PCA(dat_plot_fam[complete.cases(dat_plot_fam),])
  PCA(dat_plot_pg[complete.cases(dat_plot_fam),])
  PCA(dat_plot[complete.cases(dat_plot),])
  
  
  plot(dat_plot$Silhouette_observed,dat_plot$Silhouette_Change,col="white")
  text(dat_plot$Silhouette_observed,dat_plot$Silhouette_Change,col=c(colzR,colzG),labels = dat_plot$Nb)
  
  plot(dat_plot$Nb,dat_plot$Silhouette_Change,col=c(colzR,colzG))

  plot(dat_plot$Silhouette_observed,dat_plot$Silhouette_Change,
       col=c(colzR,colzG),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Data set 1")
  abline(h=0,col="gray")
  
  
pdf(file=file.path(origin,"plots","Sil_Change_ClusterSize.pdf"))
{
  par(mfrow=c(2,2))
  plot(dat_plot_spec$Nb,dat_plot_spec$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Species,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_spec$Nb,dat_plot_spec$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Species,data set 2")
  abline(h=0,col="gray")
  
  plot(dat_plot_spec$Silhouette_observed,dat_plot_spec$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Species,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_spec$Silhouette_observed,dat_plot_spec$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Species,data set 2")
  abline(h=0,col="gray")
  
  # ---- GENUS
  plot(dat_plot_gen$Nb,dat_plot_gen$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_gen)),colzG_gen),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Genus,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_gen$Nb,dat_plot_gen$Silhouette_Change,
       col=c(colzR_gen,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_gen))),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Genus,data set 2")
  abline(h=0,col="gray")

  plot(dat_plot_gen$Silhouette_observed,dat_plot_gen$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_gen)),colzG_gen),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Genus,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_gen$Silhouette_observed,dat_plot_gen$Silhouette_Change,
       col=c(colzR_gen,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_gen))),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Genus,data set 2")
  abline(h=0,col="gray")
  
  # ---- Family
  plot(dat_plot_fam$Nb,dat_plot_fam$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_fam)),colzG_fam),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Family,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_fam$Nb,dat_plot_fam$Silhouette_Change,
       col=c(colzR_fam,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_fam))),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Family,data set 2")
  abline(h=0,col="gray")
  
  plot(dat_plot_fam$Silhouette_observed,dat_plot_fam$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_fam)),colzG_fam),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Family,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_fam$Silhouette_observed,dat_plot_fam$Silhouette_Change,
       col=c(colzR_fam,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_fam))),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Family,data set 2")
  abline(h=0,col="gray")
  
  # ---- PG
  plot(dat_plot_pg$Nb,dat_plot_pg$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_pg)),colzG_pg),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="PhyloG group,data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_pg$Nb,dat_plot_pg$Silhouette_Change,
       col=c(colzR_pg,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_pg))),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="PhyloG group,data set 2")
  abline(h=0,col="gray")
  
  plot(dat_plot_pg$Silhouette_observed,dat_plot_pg$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_pg)),colzG_pg),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="PhyloG group, data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_pg$Silhouette_observed,dat_plot_pg$Silhouette_Change,
       col=c(colzR_pg,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_pg))),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="PhyloG group, data set 2")
  abline(h=0,col="gray")
  
  # ---- GF
  plot(dat_plot_gf$Nb,dat_plot_gf$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_gf)),colzG_gf),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Growth form, data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_gf$Nb,dat_plot_gf$Silhouette_Change,
       col=c(colzR_gf,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_gf))),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)",main="Growth form, data set 2")
  abline(h=0,col="gray")
  
  plot(dat_plot_gf$Silhouette_observed,dat_plot_gf$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR_gf)),colzG_gf),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Growth form, data set 1")
  abline(h=0,col="gray")
  plot(dat_plot_gf$Silhouette_observed,dat_plot_gf$Silhouette_Change,
       col=c(colzR_gf,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG_gf))),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed",main="Growth form, data set 2")
  abline(h=0,col="gray")
  dev.off()
}

  pdf(file=file.path(origin,"plots","SD_info.pdf"))
  par(mfrow=c(3,2))
  plot(dat_plot$Silhouette_Change,dat_plot$SD,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylim=c(0.83,1),
       xlab="Silhouette 70%-observed",ylab="Confidence (sd)")
  plot(dat_plot$Silhouette_Change,dat_plot$SD,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylim=c(0.83,1),
       xlab="Silhouette 70%-observed",ylab="Confidence (sd)")
  
  plot(dat_plot$RMSE_0,dat_plot$SD,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylim=c(.8,1.2),
       xlab="Error (0%)",ylab="Confidence (sd)")
  plot(dat_plot$RMSE_0,dat_plot$SD,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylim=c(.8,1.2),
       xlab="Error (0%)",ylab="Confidence (sd)")
  
  plot(dat_plot$Nb,dat_plot$SD,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylim=c(.8,1.2),
       xlab="Cluster size (nb)",ylab="Confidence (sd)")
  plot(dat_plot$Nb,dat_plot$SD,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylim=c(.8,1.2),
       xlab="Cluster size (nb)",ylab="Confidence (sd)")
  dev.off()
  
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylab="Silhouette 70%-observed",xlab="Error (0% gaps)")
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylab="Silhouette 70%-observed",xlab="Error (0% gaps)")
  
  plot(dat_plot$RMSE70,dat_plot$Silhouette_70.gaps,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       xlim=c(0,5),
       ylab="Error (70%)",xlab="Silhouette (70% gaps)")
  plot(dat_plot$RMSE70,dat_plot$Silhouette_70.gaps,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       xlim=c(0,5),
       ylab="Error (70%)",xlab="Silhouette (70% gaps)")
  dev.off()

  
  pdf(file=file.path(origin,"plots","Error_AsIndicator.pdf"),height = 4,width=8)
  par(mfrow=c(1,2))
  plot(dat_plot$RMSE70,dat_plot$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       xlim=c(0,1.2),
       ylab="Silhouette 70%-observed",xlab="Error (70% gaps)")
  abline(v = mean(dat_plot$RMSE70[(length(colzR)+1):nrow(dat_plot)],na.rm = T),col="blue")
  plot(dat_plot$RMSE70,dat_plot$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       xlim=c(0.,1.2),
       ylab="Silhouette 70%-observed",xlab="Error (70% gaps)")
  abline(v = mean(dat_plot$RMSE70[1:length(colzG)],na.rm = T),col="red")
  dev.off()
  
  
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       xlim=c(0,1.2),
       ylab="Silhouette 70%-observed",xlab="RMSE (0% gaps)")
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       xlim=c(0.,1.2),
       ylab="Silhouette 70%-observed",xlab="RMSE (0% gaps)")
  
  
  plot(dat_plot$RMSE70,dat_plot$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       xlim=c(0,1.2),
       ylab="Silhouette 70%-observed",xlab="RMSE (70% gaps)")
  plot(dat_plot$RMSE70,dat_plot$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       xlim=c(0.,1.2),
       ylab="Silhouette 70%-observed",xlab="RMSE (70% gaps)")
  
  plot(dat_plot$RMSE70,dat_plot$Silhouette_70.gaps,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       xlim=c(0,1.2),
       ylab="Silhouette 70%",xlab="RMSE (70% gaps)")
  plot(dat_plot$RMSE70,dat_plot$Silhouette_70.gaps,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       xlim=c(0.,1.2),
       ylab="Silhouette 70%",xlab="RMSE (70% gaps)")
  
  plot(dat_plot$RMSE_0,dat_plot$SD,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylim=c(.8,1.1),
       ylab="Error (0%)",xlab="Confidence (sd)")
  plot(dat_plot$RMSE_0,dat_plot$SD,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylim=c(.8,1.1),
       ylab="Error (0%)",xlab="Confidence (sd)")
  
  plot(dat_plot$Nb,dat_plot$SD,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylim=c(.8,1.22),
       ylab="Error (0%)",xlab="Confidence (sd)")
  plot(dat_plot$Nb,dat_plot$SD,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylim=c(.8,1.22),
       ylab="Error (0%)",xlab="Confidence (sd)")
  dev.off()
  
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzR)),colzG),
       ylab="Silhouette 70%-observed",xlab="Error (0% gaps)")
  plot(dat_plot$RMSE_0,dat_plot$Silhouette_Change,
       col=c(colzR,rep(rgb(0, 0, 255, max = 255, alpha = 0, names = "white"),length(colzG))),
       ylab="Silhouette 70%-observed",xlab="Error (0% gaps)")

  dev.off()
  
  pdf(file=file.path(origin,"plots","Sil_Change_ClusterSize_2.pdf"))
  par(mfrow=c(2,2))
  plot(dat_plot$Nb,dat_plot$Silhouette_Change,
       col=c(colzR,colzG),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)")
  plot(dat_plot$Silhouette_observed,dat_plot$Silhouette_Change,col=c(colzR,colzG),
       ylab="Silhouette 70%-observed",xlab="Silhouette observed")
  plot(dat_plot$SD,dat_plot$Silhouette_Change,
       col=c(colzR,colzG),xlim=c(0.83,1),
       ylab="Silhouette 70%-observed",xlab="Cluster size (nb)")
  dev.off()
  
  dat_plot <- cbind(
    c(silR_0_spec[,2],silG_0_spec[,2]),
    log(c(silR_0_spec[,4],silG_0_spec[,4])+2),
    c(silR_0_spec[,3],silG_0_spec[,3]),
    c(silR_70_spec[,3],silG_70_spec[,3]),
    c((silR_70_spec[,3]-silR_0_spec[,3]),(silG_70_spec[,3]-silG_0_spec[,3])),
    c(rowMeans(rmseR_0_spec[,2:ncol(rmseR_0_spec)]),rowMeans(rmseG_0_spec[,2:ncol(rmseG_0_spec)])))
  colnames(dat_plot) <- c("Nb","SD","Silhouette_70%gaps","Silhouette_observed","Silhouette_Change")
  pairs(dat_plot,col=c(colzR,colzG))
  
  
  
  plot(c(silR_0_spec[,2],silG_0_spec[,2]),
       c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
       pch=16,main="Species",ylim=c(-1.5,1.5),cex=.5,log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_70_spec[,2],silG_70_spec[,2]),
       c(silR_70_spec[,4],silG_70_spec[,4])+2,
       pch=16,main="Species",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_spec[,4],silG_70_spec[,4]),
       c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
       pch=16,main="Species",cex=.8,ylim=c(0,1.5),xlim = c(.3,.8),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  #---------------------------------------------
  
  
  
  
  
  
  
  
  
  
  
  plot(c(silR_0_gen[,3],silG_0_gen[,3]),
       c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
       pch=16,main="Genus",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_fam[,3],silG_0_fam[,3]),
       c(silR_70_fam[,3]-silR_0_fam[,3],
         silG_70_fam[,3]-silG_0_fam[,3]),
       pch=16,main="Family",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_0_pg[,3],silG_0_pg[,3]),
       c(silR_70_pg[,3]-silR_0_pg[,3],
         silG_70_pg[,3]-silG_0_pg[,3]),
       pch=16,main="Phylogenetic Groups",xlim=c(-1,1),ylim=c(-1.5,1.5),
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_gf[,3],silG_0_gf[,3]),
       c(silR_70_gf[,3]-silR_0_gf[,3],
         silG_70_gf[,3]-silG_0_gf[,3]),
       main="Growth forms",xlim=c(-1,1),ylim=c(-1.5,1.5),
       pch=16,xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_pft[,3],silG_0_pft[,3]),
       c(silR_70_pft[,3]-silR_0_pft[,3],
         silG_70_pft[,3]-silG_0_pft[,3]),
       pch=16,main="PFTs",xlim=c(-1,1),ylim=c(-1.5,1.5),
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
}
dev.off()

pdf(file=file.path(origin,"plots","figure_4","Silhouette_ClusterSize.pdf"),
    width=8,height = 6)
{
  
  par(mfrow=c(2,3))
  
  plot(c(silR_0_gen[,2],silG_0_gen[,2]),
       c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
       pch=16,main="Genus",ylim=c(-1.5,1.5),cex=.5,log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_fam[,2],silG_0_fam[,2]),
       c(silR_70_fam[,3]-silR_0_fam[,3],
         silG_70_fam[,3]-silG_0_fam[,3]),
       pch=16,main="Family",ylim=c(-1.5,1.5),cex=.5,log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_0_pg[,2],silG_0_pg[,2]),
       c(silR_70_pg[,3]-silR_0_pg[,3],
         silG_70_pg[,3]-silG_0_pg[,3]),
       pch=16,main="Phylogenetic Groups",ylim=c(-1.5,1.5),log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_gf[,2],silG_0_gf[,2]),
       c(silR_70_gf[,3]-silR_0_gf[,3],
         silG_70_gf[,3]-silG_0_gf[,3]),
       main="Growth forms",ylim=c(-1.5,1.5),log="x",
       pch=16,xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_pft[,2],silG_0_pft[,2]),
       c(silR_70_pft[,3]-silR_0_pft[,3],
         silG_70_pft[,3]-silG_0_pft[,3]),
       pch=16,main="PFTs",ylim=c(-1.5,1.5),log="x",
       xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
}
dev.off()


pdf(file=file.path(origin,"plots","figure_4","ClusterSize_SD.pdf"),
    width=8,height = 6)
{
  par(mfrow=c(2,3))
 
  plot(c(silR_70_gen[,2],silG_70_gen[,2]),
       c(silR_70_gen[,4],silG_70_gen[,4])+2,
       pch=16,main="Genus",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_fam[,2],silG_70_fam[,2]),
       c(silR_70_fam[,4],silG_70_fam[,4])+2,
       pch=16,main="Family",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_70_pg[,2],silG_70_pg[,2]),
       c(silR_70_pg[,4],silG_70_pg[,4])+2,
       pch=16,main="Phylogenetic group",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_gf[,2],silG_70_gf[,2]),
       c(silR_70_gf[,4],silG_70_gf[,4]),
       pch=16,main="Growth form",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_pft[,2],silG_70_pft[,2]),
       c(silR_70_pft[,4],silG_70_pft[,4]),
       pch=16,main="PFTs",cex=.8,log="xy",
       xlab="Cluster size (nb)",ylab="SD confidence",
       col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
}
dev.off()


pdf(file=file.path(origin,"plots","figure_4","Silhouette_SD.pdf"),
    width=8,height = 6)
{
  
  plot(c(silR_70_gen[,4],silG_70_gen[,4]),
       c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
       pch=16,main="Genus",ylim=c(-1.5,1.5),cex=.8,xlim=c(.3,1),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_fam[,4],silG_70_fam[,4]),
       c(silR_70_fam[,3]-silR_0_fam[,3],
         silG_70_fam[,3]-silG_0_fam[,3]),
       pch=16,main="Family",ylim=c(-1.5,1.5),cex=.8,xlim=c(.3,1),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_70_pg[,4],silG_70_pg[,4]),
       c(silR_70_pg[,3]-silR_0_pg[,3],
         silG_70_pg[,3]-silG_0_pg[,3]),
       pch=16,main="Phylogenetic Groups",ylim=c(-1.5,1.5),xlim=c(.3,1),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_gf[,4],silG_70_gf[,4]),
       c(silR_70_gf[,3]-silR_0_gf[,3],
         silG_70_gf[,3]-silG_0_gf[,3]),
       main="Growth forms",ylim=c(-1.5,1.5),xlim=c(.3,1),
       pch=16,xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_70_pft[,4],silG_70_pft[,4]),
       c(silR_70_pft[,3]-silR_0_pft[,3],
         silG_70_pft[,3]-silG_0_pft[,3]),
       pch=16,main="PFTs",ylim=c(-1.5,1.5),xlim=c(.3,1),
       xlab="SD confidence",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
}
dev.off()











#-------------------------------------------------



  pdf(file=file.path(origin,"plots","figure_4","Silhouette_depends_on_orgSil.pdf"),
      width=8,height = 6)
  {
   par(mfrow=c(2,3))
    plot(c(silR_0_spec[,3],silG_0_spec[,3]),
    c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
    pch=16,main="Species",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
    xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
    col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)

  plot(c(silR_0_gen[,3],silG_0_gen[,3]),
       c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
       pch=16,main="Genus",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_fam[,3],silG_0_fam[,3]),
       c(silR_70_fam[,3]-silR_0_fam[,3],
         silG_70_fam[,3]-silG_0_fam[,3]),
       pch=16,main="Family",xlim=c(-1,1),ylim=c(-1.5,1.5),cex=.5,
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  
  plot(c(silR_0_pg[,3],silG_0_pg[,3]),
       c(silR_70_pg[,3]-silR_0_pg[,3],
         silG_70_pg[,3]-silG_0_pg[,3]),
       pch=16,main="Phylogenetic Groups",xlim=c(-1,1),ylim=c(-1.5,1.5),
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_gf[,3],silG_0_gf[,3]),
       c(silR_70_gf[,3]-silR_0_gf[,3],
         silG_70_gf[,3]-silG_0_gf[,3]),
       main="Growth forms",xlim=c(-1,1),ylim=c(-1.5,1.5),
       pch=16,xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  plot(c(silR_0_pft[,3],silG_0_pft[,3]),
       c(silR_70_pft[,3]-silR_0_pft[,3],
         silG_70_pft[,3]-silG_0_pft[,3]),
       pch=16,main="PFTs",xlim=c(-1,1),ylim=c(-1.5,1.5),
       xlab="Observed silhouette",ylab="70% gaps - observed silhouette",
       col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
  abline(h = 0,col="gray",lty=2)
  abline(v = 0,col="gray",lty=2)
  
  }
  dev.off()
  
  pdf(file=file.path(origin,"plots","figure_4","Silhouette_ClusterSize.pdf"),
      width=8,height = 6)
  {
    
    par(mfrow=c(2,3))
    plot(c(silR_0_spec[,2],silG_0_spec[,2]),
         c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
         pch=16,main="Species",ylim=c(-1.5,1.5),cex=.5,log="x",
         xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_0_gen[,2],silG_0_gen[,2]),
         c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
         pch=16,main="Genus",ylim=c(-1.5,1.5),cex=.5,log="x",
         xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_0_fam[,2],silG_0_fam[,2]),
         c(silR_70_fam[,3]-silR_0_fam[,3],
           silG_70_fam[,3]-silG_0_fam[,3]),
         pch=16,main="Family",ylim=c(-1.5,1.5),cex=.5,log="x",
         xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    
    plot(c(silR_0_pg[,2],silG_0_pg[,2]),
         c(silR_70_pg[,3]-silR_0_pg[,3],
           silG_70_pg[,3]-silG_0_pg[,3]),
         pch=16,main="Phylogenetic Groups",ylim=c(-1.5,1.5),log="x",
         xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_0_gf[,2],silG_0_gf[,2]),
         c(silR_70_gf[,3]-silR_0_gf[,3],
           silG_70_gf[,3]-silG_0_gf[,3]),
         main="Growth forms",ylim=c(-1.5,1.5),log="x",
         pch=16,xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_0_pft[,2],silG_0_pft[,2]),
         c(silR_70_pft[,3]-silR_0_pft[,3],
           silG_70_pft[,3]-silG_0_pft[,3]),
         pch=16,main="PFTs",ylim=c(-1.5,1.5),log="x",
         xlab="Cluster size (nb)",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
  }
  dev.off()
  
  
  pdf(file=file.path(origin,"plots","figure_4","ClusterSize_SD.pdf"),
      width=8,height = 6)
  {
    par(mfrow=c(2,3))
    plot(c(silR_70_spec[,2],silG_70_spec[,2]),
         c(silR_70_spec[,4],silG_70_spec[,4])+2,
         pch=16,main="Species",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_gen[,2],silG_70_gen[,2]),
         c(silR_70_gen[,4],silG_70_gen[,4])+2,
         pch=16,main="Genus",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_fam[,2],silG_70_fam[,2]),
         c(silR_70_fam[,4],silG_70_fam[,4])+2,
         pch=16,main="Family",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    
    plot(c(silR_70_pg[,2],silG_70_pg[,2]),
         c(silR_70_pg[,4],silG_70_pg[,4])+2,
         pch=16,main="Phylogenetic group",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_gf[,2],silG_70_gf[,2]),
         c(silR_70_gf[,4],silG_70_gf[,4]),
         pch=16,main="Growth form",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_pft[,2],silG_70_pft[,2]),
         c(silR_70_pft[,4],silG_70_pft[,4]),
         pch=16,main="PFTs",cex=.8,log="xy",
         xlab="Cluster size (nb)",ylab="SD confidence",
         col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
  }
  dev.off()

  
  pdf(file=file.path(origin,"plots","figure_4","Silhouette_SD.pdf"),
      width=8,height = 6)
  {
    par(mfrow=c(2,3))
    plot(c(silR_70_spec[,4],silG_70_spec[,4]),
         c(silR_70_spec[,3]-silR_0_spec[,3],silG_70_spec[,3]-silG_0_spec[,3]),
         pch=16,main="Species",cex=.8,ylim=c(0,1.5),xlim = c(.3,.8),
         xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_spec)),rep("blue",nrow(silG_70_spec))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_gen[,4],silG_70_gen[,4]),
         c(silR_70_gen[,3]-silR_0_gen[,3],silG_70_gen[,3]-silG_0_gen[,3]),
         pch=16,main="Genus",ylim=c(-1.5,1.5),cex=.8,xlim=c(.3,1),
         xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_gen)),rep("blue",nrow(silG_70_gen))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_fam[,4],silG_70_fam[,4]),
         c(silR_70_fam[,3]-silR_0_fam[,3],
           silG_70_fam[,3]-silG_0_fam[,3]),
         pch=16,main="Family",ylim=c(-1.5,1.5),cex=.8,xlim=c(.3,1),
         xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_fam)),rep("blue",nrow(silG_70_fam))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    
    plot(c(silR_70_pg[,4],silG_70_pg[,4]),
         c(silR_70_pg[,3]-silR_0_pg[,3],
           silG_70_pg[,3]-silG_0_pg[,3]),
         pch=16,main="Phylogenetic Groups",ylim=c(-1.5,1.5),xlim=c(.3,1),
         xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_pg)),rep("blue",nrow(silG_70_pg))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_gf[,4],silG_70_gf[,4]),
         c(silR_70_gf[,3]-silR_0_gf[,3],
           silG_70_gf[,3]-silG_0_gf[,3]),
         main="Growth forms",ylim=c(-1.5,1.5),xlim=c(.3,1),
         pch=16,xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_gf)),rep("blue",nrow(silG_70_gf))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
    plot(c(silR_70_pft[,4],silG_70_pft[,4]),
         c(silR_70_pft[,3]-silR_0_pft[,3],
           silG_70_pft[,3]-silG_0_pft[,3]),
         pch=16,main="PFTs",ylim=c(-1.5,1.5),xlim=c(.3,1),
         xlab="SD confidence",ylab="70% gaps - observed silhouette",
         col=c(rep("red",nrow(silR_70_pft)),rep("blue",nrow(silG_70_pft))))
    abline(h = 0,col="gray",lty=2)
    abline(v = 0,col="gray",lty=2)
    
  }
  dev.off()
  
  