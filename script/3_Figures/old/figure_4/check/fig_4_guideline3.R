

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
  #COVA
  {
    # load the coefficience of variance data per cluster GUIDO org
  cova_now <- rep(NA,7)
  repnums=1:10
  GapPercent1="org"
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "guido",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }

  covaG_0_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaG_0_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaG_0_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaG_0_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaG_0_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaG_0_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,8)
  repnums=1:10
  GapPercent1="org"
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }
  
  covaR_0_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaR_0_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaR_0_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaR_0_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaR_0_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaR_0_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  
  # load the coefficience of variance data per cluster GUIDO 20
  cova_now <- rep(NA,7)
  repnums=1:10
  GapPercent1=20
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "guido",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }
  
  covaG_20_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaG_20_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaG_20_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaG_20_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaG_20_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaG_20_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  
  # load the coefficience of variance data per cluster RAINFOR 20
  cova_now <- rep(NA,8)
  repnums=1:10
  GapPercent1=20
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }
  
  covaR_20_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaR_20_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaR_20_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaR_20_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaR_20_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaR_20_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  
  
  # load the coefficience of variance data per cluster GUIDO 70
  cova_now <- rep(NA,7)
  repnums=1:10
  GapPercent1=70
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "guido",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }
  
  covaG_70_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaG_70_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaG_70_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaG_70_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaG_70_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaG_70_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  
  # load the coefficience of variance data per cluster RAINFOR
  cova_now <- rep(NA,8)
  repnums=1:10
  GapPercent1=70
  for(RepNum in repnums){
    path_now1 <- file.path(origin,"data_output","CoVa",
                           "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
    path_now2=file.path(origin,"data_output","CoVa",
                        "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now1)&file.exists(path_now2)){
      cova_0_general <- read.csv(file=path_now1)
      cova_0_rawclust <- read.csv(file=path_now2)
      print(dim(cova_0_rawclust))
      cova_now <- rbind(cova_now,cova_0_rawclust)
    }
  }
  
  covaR_70_spec <- cova_now[cova_now[,1]=="Species",2:ncol(cova_now)]
  covaR_70_gen <- cova_now[cova_now[,1]=="Genus",2:ncol(cova_now)]
  covaR_70_fam <- cova_now[cova_now[,1]=="Family",2:ncol(cova_now)]
  covaR_70_pg <- cova_now[cova_now[,1]=="PG",2:ncol(cova_now)]
  covaR_70_gf <- cova_now[cova_now[,1]=="GF",2:ncol(cova_now)]
  covaR_70_pft <- cova_now[cova_now[,1]=="PFT",2:ncol(cova_now)]
  }

  #RMSE
  {
  #Load RMSE per cluster GUIDO 0% gaps
  rmse_trait <- rep(NA,12)
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
    rmse_trait <- rbind(rmse_trait,rmse_0_rawclust)
  }
  
  # calculate per cluster an average
  rmseG_0_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
  rmseG_0_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
  rmseG_0_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
  rmseG_0_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
  rmseG_0_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
  rmseG_0_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)

  #Load RMSE per cluster RAINFOR 0% gaps
  sil_now <- rep(NA,13)
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
    sil_now <- rbind(sil_now,rmse_0_rawclust)
  }
  
  # calculate per cluster an average
  rmseR_0_spec <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,2]),FUN=mean,na.rm=TRUE)
  rmseR_0_gen <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,3]),FUN=mean,na.rm=TRUE)
  rmseR_0_fam <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,4]),FUN=mean,na.rm=TRUE)
  rmseR_0_pg <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,5]),FUN=mean,na.rm=TRUE)
  rmseR_0_gf <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,6]),FUN=mean,na.rm=TRUE)
  rmseR_0_pft <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,7]),FUN=mean,na.rm=TRUE)
  
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
    sil_now <- rbind(sil_now,rmse_0_rawclust)
  }
  
  # calculate per cluster an average
  rmseG_70_spec <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,2]),FUN=mean,na.rm=TRUE)
  rmseG_70_gen <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,3]),FUN=mean,na.rm=TRUE)
  rmseG_70_fam <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,4]),FUN=mean,na.rm=TRUE)
  rmseG_70_pg <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,5]),FUN=mean,na.rm=TRUE)
  rmseG_70_gf <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,6]),FUN=mean,na.rm=TRUE)
  rmseG_70_pft <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,7]),FUN=mean,na.rm=TRUE)

  #70% gaps
  sil_now <- rep(NA,13)
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
    sil_now <- rbind(sil_now,rmse_0_rawclust)
  }
  
  # calculate per cluster an average
  rmseR_70_spec <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,2]),FUN=mean,na.rm=TRUE)
  rmseR_70_gen <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,3]),FUN=mean,na.rm=TRUE)
  rmseR_70_fam <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,4]),FUN=mean,na.rm=TRUE)
  rmseR_70_pg <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,5]),FUN=mean,na.rm=TRUE)
  rmseR_70_gf <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,6]),FUN=mean,na.rm=TRUE)
  rmseR_70_pft <- aggregate(sil_now[,8:ncol(sil_now)],by=list(sil_now[,7]),FUN=mean,na.rm=TRUE)
  }
  
# load the Silhouette data per cluster
  {
# load the mean(sd) data per cluster within silhouettes I think
# load the number of observations somehowf
  sil_now <- rep(NA,6)
  repnums=1:3
  GapPercent1="org"
  for(RepNum in repnums){
    path_now=file.path(origin,"data_output","Sil","guido",
                       TD_choice,GapPercent1,
                       RepNum,paste0("all_clusters.csv"))
    if(file.exists(path_now)){
      sil_nowtmp <- read.csv(file=path_now,row.names = NULL)
      head(sil_nowtmp)
      print(dim(sil_nowtmp))
      sil_now <- rbind(sil_now,sil_nowtmp)
    }
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
}




# Merge them together into one happy data frame per group
names(covaG_70_spec) <- paste0(names(covaG_70_spec),"_cova")
names(covaG_20_spec) <- paste0(names(covaG_20_spec),"_cova")
names(covaG_0_spec) <- paste0(names(covaG_0_spec),"_cova")
names(covaR_70_spec)<- paste0(names(covaR_70_spec),"_cova")
names(covaR_20_spec)<- paste0(names(covaR_20_spec),"_cova")
names(covaR_0_spec)<- paste0(names(covaR_0_spec),"_cova")
colnames(rmseG_0_spec)<- paste0(names(rmseG_0_spec),"_rmse")
colnames(rmseG_70_spec)<- paste0(names(rmseG_70_spec),"_rmse")
colnames(rmseR_0_spec)<- paste0(names(rmseR_0_spec),"_rmse")
colnames(rmseR_70_spec)<- paste0(names(rmseR_70_spec),"_rmse")

names(covaG_70_spec)[1] <- "Cluster"
names(covaG_20_spec)[1] <- "Cluster"
names(covaG_0_spec)[1] <- "Cluster"
names(covaR_70_spec)[1] <- "Cluster"
names(covaR_20_spec)[1] <- "Cluster"
names(covaR_0_spec)[1] <- "Cluster"
colnames(rmseG_0_spec)[1] <- "Cluster"
colnames(rmseG_70_spec)[1] <- "Cluster"
colnames(rmseR_0_spec)[1] <- "Cluster"
colnames(rmseR_70_spec)[1] <- "Cluster"

G0_spec <-  cbind(silG_0_spec,
                  covaG_0_spec[match(x = silG_0_spec$Cluster,table = covaG_0_spec$Cluster),2:ncol(covaG_0_spec)],
                  rmseG_0_spec[match(x = silG_0_spec$Cluster,table = rmseG_0_spec$Cluster),2:ncol(rmseG_0_spec)])
R0_spec <-  cbind(silR_0_spec,
                  covaR_0_spec[match(x = silR_0_spec$Cluster,table = covaR_0_spec$Cluster),2:ncol(covaR_0_spec)],
                  rmseR_0_spec[match(x = silR_0_spec$Cluster,table = rmseR_0_spec$Cluster),2:ncol(covaR_0_spec)])
G20_spec <-  cbind(covaG_20_spec[match(x = silG_0_spec$Cluster,table = covaG_20_spec$Cluster),
                                2:ncol(covaG_20_spec)])
R20_spec <-  cbind(covaR_20_spec[match(x = silR_0_spec$Cluster,table = covaR_20_spec$Cluster),
                                 2:ncol(covaR_20_spec)])
G70_spec <-  cbind(silG_70_spec,
                  covaG_70_spec[match(x = silG_70_spec$Cluster,table = covaG_70_spec$Cluster),2:ncol(covaG_0_spec)],
                  rmseG_70_spec[match(x = silG_70_spec$Cluster,table = rmseG_70_spec$Cluster),2:ncol(covaG_0_spec)])
R70_spec <-  cbind(silR_70_spec,
                  covaR_70_spec[match(x = silR_70_spec$Cluster,table = covaR_70_spec$Cluster),2:ncol(covaR_0_spec)],
                  rmseR_70_spec[match(x = silR_70_spec$Cluster,table = rmseR_70_spec$Cluster),2:ncol(covaR_0_spec)])
G70_spec <-  cbind(silG_70_spec,
                   covaG_70_spec[match(x = silG_70_spec$Cluster,table = covaG_70_spec$Cluster),2:ncol(covaG_0_spec)],
                   rmseG_70_spec[match(x = silG_70_spec$Cluster,table = rmseG_70_spec$Cluster),2:ncol(covaG_0_spec)])
R70_spec <-  cbind(silR_70_spec,
                   covaR_70_spec[match(x = silR_70_spec$Cluster,table = covaR_70_spec$Cluster),2:ncol(covaR_0_spec)],
                   rmseR_70_spec[match(x = silR_70_spec$Cluster,table = rmseR_70_spec$Cluster),2:ncol(covaR_0_spec)])


pdf(file = file.path(origin,"plots","Cova_Solution.pdf"))
{
covas <- grep(colnames(G0_spec),pattern = "cova")
t=1
par(mfrow=c(3,3))
for(t in 1:length(covas)){
  plot(G0_spec[,covas[t]],(G70_spec[,covas[t]]-G0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
       main=colnames(G0_spec)[covas[t]],col="blue",
       xlab="Cova observed",ylab="Cova 70% - observed")
    tR=which(colnames(R0_spec)==colnames(G0_spec)[covas[t]])
  if(sum(!is.na(R0_spec[,tR]))>3){
    points(R0_spec[,tR],(R70_spec[,tR]-R0_spec[,tR]),main=colnames(R0_spec)[tR],col="red")
  }
    abline(a = c(0,-1),col="gray")
}
covas <- grep(colnames(R0_spec),pattern = "cova")
t=1
for(t in 1:length(covas)){
  tR=which(colnames(G0_spec)==colnames(R0_spec)[covas[t]])
  if(sum(!is.na(G0_spec[,tR]))<3){
    plot(R0_spec[,covas[t]],(R70_spec[,covas[t]]-R0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
         main=colnames(R0_spec)[covas[t]],col="red",
         xlab="Cova observed",ylab="Cova 70% - observed")
    abline(a = c(0,-1),col="gray")
  }
}


covas <- grep(colnames(G0_spec),pattern = "cova")
t=2
par(mfrow=c(3,3))
for(t in 1:length(covas)){
  plot(G0_spec[,covas[t]],(G20_spec[,(covas-4)[t]]-G0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
       main=colnames(G0_spec)[covas[t]],col="blue",
       xlab="Cova observed",ylab="Cova 20% - observed")
  tR=which(colnames(R0_spec)==colnames(G0_spec)[covas[t]])
  if(sum(!is.na(R0_spec[,tR]))>3){
    points(R0_spec[,tR],(R20_spec[,tR]-R0_spec[,tR]),main=colnames(R0_spec)[tR],col="red")
  }
  abline(a = c(0,-1),col="gray")
}
covas <- grep(colnames(R0_spec),pattern = "cova")
t=1
for(t in 1:length(covas)){
  tR=which(colnames(G0_spec)==colnames(R0_spec)[covas[t]])
  if(sum(!is.na(G0_spec[,tR]))<3){
    plot(R0_spec[,covas[t]],(R20_spec[,(covas-4)[t]]-R0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
         main=colnames(R0_spec)[covas[t]],col="red",
         xlab="Cova observed",ylab="Cova 20% - observed")
    abline(a = c(0,-1),col="gray")
  }
}
dev.off()
}



pdf(file = file.path(origin,"plots","Cova_Trait_wise.pdf"))
{
  par(mfrow=c(2,2))
    boxplot(G0_spec[,grep(names(G0_spec),pattern = "cova")],ylim=c(-5,5),las=2,col="blue",
            main="Data set 1",ylab="Cova species",xaxt="n")
    axis(1,at = seq(1,to=length(grep(names(G0_spec),pattern = "cova")),by = 1),las=2,
         labels =gsub(names(G0_spec)[grep(names(G0_spec),pattern = "cova")],pattern = "_cova",replacement = ""))
    abline(h=0,col="gray",lty=2)
    boxplot(R0_spec[,grep(names(R0_spec),pattern = "cova")],ylim=c(-5,5),las=2,col="red",
            main="Data set 2",ylab="Cova species",xaxt="n")
    axis(1,at = seq(1,to=length(grep(names(R0_spec),pattern = "cova")),by = 1),las=2,
         labels =gsub(names(R0_spec)[grep(names(R0_spec),pattern = "cova")],pattern = "_cova",replacement = ""))
    abline(h=0,col="gray",lty=2)
    
    boxplot(G70_spec[,grep(names(G0_spec),pattern = "cova")]-G0_spec[,grep(names(G0_spec),pattern = "cova")],
            ylim=c(-5,5),las=2,col="blue",
            main="Data set 1",ylab="Deviation species",xaxt="n")
    axis(1,at = seq(1,to=length(grep(names(G0_spec),pattern = "cova")),by = 1),las=2,
    labels =gsub(names(G0_spec)[grep(names(G0_spec),pattern = "cova")],pattern = "_cova",replacement = ""))
    abline(h=0,col="gray",lty=2)
    boxplot(R70_spec[,grep(names(R0_spec),pattern = "cova")]-R0_spec[,grep(names(R0_spec),pattern = "cova")],
            ylim=c(-5,5),las=2,col="red",
            main="Data set 2",ylab="Deviation species",xaxt="n")
    axis(1,at = seq(1,to=length(grep(names(R0_spec),pattern = "cova")),by = 1),las=2,
         labels =gsub(names(R0_spec)[grep(names(R0_spec),pattern = "cova")],pattern = "_cova",replacement = ""))
    abline(h=0,col="gray",lty=2)
    dev.off()
}


pdf(file = file.path(origin,"plots","Cova_ClusterSize.pdf"))

  covas <- grep(colnames(G0_spec),pattern = "cova")
  t=1
  par(mfrow=c(3,3))
  for(t in 1:length(covas)){
    plot(G0_spec[,covas[t]],(G70_spec[,covas[t]]-G0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
         main=colnames(G0_spec)[covas[t]],col="blue",
         xlab="Cova observed",ylab="Cova 70% - observed")
    tR=which(colnames(R0_spec)==colnames(G0_spec)[covas[t]])
    if(sum(!is.na(R0_spec[,tR]))>3){
      points(R0_spec[,tR],(R70_spec[,tR]-R0_spec[,tR]),main=colnames(R0_spec)[tR],col="red")
    }
    abline(a = c(0,-1),col="gray")
  }
  covas <- grep(colnames(R0_spec),pattern = "cova")
  t=1
  for(t in 1:length(covas)){
    tR=which(colnames(G0_spec)==colnames(R0_spec)[covas[t]])
    if(sum(!is.na(G0_spec[,tR]))<3){
      plot(R0_spec[,covas[t]],(R70_spec[,covas[t]]-R0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
           main=colnames(R0_spec)[covas[t]],col="red",
           xlab="Cova observed",ylab="Cova 70% - observed")
      abline(a = c(0,-1),col="gray")
    }
  }
  
  
  covas <- grep(colnames(G0_spec),pattern = "cova")
  t=2
  par(mfrow=c(3,3))
  for(t in 1:length(covas)){
    plot(G0_spec[,covas[t]],(G20_spec[,(covas-4)[t]]-G0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
         main=colnames(G0_spec)[covas[t]],col="blue",
         xlab="Cova observed",ylab="Cova 20% - observed")
    tR=which(colnames(R0_spec)==colnames(G0_spec)[covas[t]])
    if(sum(!is.na(R0_spec[,tR]))>3){
      points(R0_spec[,tR],(R20_spec[,tR]-R0_spec[,tR]),main=colnames(R0_spec)[tR],col="red")
    }
    abline(a = c(0,-1),col="gray")
  }
  covas <- grep(colnames(R0_spec),pattern = "cova")
  t=1
  for(t in 1:length(covas)){
    tR=which(colnames(G0_spec)==colnames(R0_spec)[covas[t]])
    if(sum(!is.na(G0_spec[,tR]))<3){
      plot(R0_spec[,covas[t]],(R20_spec[,(covas-4)[t]]-R0_spec[,covas[t]]),xlim=c(-50,50),ylim=c(-50,50),
           main=colnames(R0_spec)[covas[t]],col="red",
           xlab="Cova observed",ylab="Cova 20% - observed")
      abline(a = c(0,-1),col="gray")
    }
  }
  dev.off()
  

t=2
par(mfrow=c(3,3))
covas <- grep(colnames(G0_spec),pattern = "cova")
for(t in 1:length(covas)){
  plot(G0_spec$Cluster_size,(G70_spec[,covas[t]]-G0_spec[,covas[t]]),xlim=c(0,10),ylim=c(-50,50),
       main=colnames(G0_spec)[covas[t]],col="blue",
       xlab="Cova observed",ylab="Cova 70% - observed")
  tR=which(colnames(R0_spec)==colnames(G0_spec)[covas[t]])
  if(sum(!is.na(R0_spec[,tR]))>3){
    points(R0_spec$Cluster_size,(R70_spec[,tR]-R0_spec[,tR]),main=colnames(R0_spec)[tR],col="red")
  }
}
covas <- grep(colnames(R0_spec),pattern = "cova")
t=1
for(t in 1:length(covas)){
  tR=which(colnames(G0_spec)==colnames(R0_spec)[covas[t]])
  if(sum(!is.na(G0_spec[,tR]))<3){
    plot(R0_spec$Cluster_size,(R70_spec[,covas[t]]-R0_spec[,covas[t]]),xlim=c(0,10),ylim=c(-50,50),
         main=colnames(R0_spec)[covas[t]],col="red",
         xlab="Cova observed",ylab="Cova 70% - observed")

  }
}

dev.off()








