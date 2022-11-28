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
  t_choices <- out$tsubs
  TDnos = out$TD_choices
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
  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  
  mean.fun <- function(input){if(sum(!is.na(input))==0){return(NA)}else{
      return(mean(input[!is.na(input)]))}}
  count.fun <- function(input){return(sum(!is.na(input)))}

  
  res <- matrix(NA,ncol=4,nrow=40)
  colnames(res) <- c("","TD","Envelope data", "TD2")
  
#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
repnums=3
t_choice="data"
if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}



#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data"){trait_names=trait_rainfor}
if(t_choice=="data_2"){trait_names=trait_guido}

colz=colz1[c(1,5,2,3,4)]


{
  # load functional group information
  GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
  PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")
  colnames(GF) <- c("rownames","ID","AccSpeciesName","GF")
  colnames(PFT) <- c("rownames","ID","AccSpeciesName","PFT")
  # load Envelope data
  # load TDenvelope
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  taxTD <- as.data.frame(read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                     "Obs_obs_TD","data","taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  taxTDGF <- merge(taxTD,GF,by = "ID",all.x = TRUE)
  taxTDGF <- taxTDGF[,c(1,2,3,4,5,8)]
  taxTDGFPFT <- merge(taxTDGF,PFT,by = "ID",all.x = TRUE)
  head(taxTDGFPFT)
  taxTD <- taxTDGFPFT[,c(1,2,3,4,5,6,9)]
  
  ID_TD <- taxTD[,1]
  head(taxTD)
  dim(taxTD)
  dim(TD)
  dim(TD_sparse)
  # load Envelope data
  # total trait data 
  EnvTot <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                             "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  taxEnv <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs","data","taxInfo.csv"),header=TRUE))
  taxEnv <- as.data.frame(read.table(file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                      "Obs_obs","data","taxInfo.csv"),
                                     sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  
  taxEnvGF <- merge(taxEnv,GF,by = "ID",all.x = TRUE)
  taxEnvGF <- taxEnvGF[,c(1,2,3,4,5,8)]
  taxEnvGFPFT <- merge(taxEnvGF,PFT,by = "ID",all.x = TRUE)
  head(taxEnvGFPFT)
  taxEnv <- taxEnvGFPFT[,c(1,2,3,4,5,6,9)]

  ID_env <- taxEnv[,1]
  dim(taxEnv)
  summary(EnvTot)
  Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
  head(Env)
  # chose only those Env with same species information
  Env_sp <- Env[taxEnv[,2]%in%taxTD[,2],]
  dim(Env_sp)
  Env_gen <- Env[taxEnv[,3]%in%taxTD[,3],]
  dim(Env_gen)
  Env_fam <- Env[taxEnv[,4]%in%taxTD[,4],]
  dim(Env_fam)
  Env_cl <- Env[taxEnv[,5]%in%taxTD[,5],]
  dim(Env_cl)
  
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],TDenv[,1])
  abline(0,1)
  plot(TD_sparse[,1],TD[,1])
  abline(0,1)
}


n=1
res[n,1] <- "Trait Nb"
if(t_choice=="data"){  res[n,2] <- ncol(TD)}
if(t_choice=="data_2"){  res[n,4] <- ncol(TD)}
res[n,3] <- ncol(EnvTot[,-which(colnames(EnvTot)=="ObservationID")])

n=n+1
res[n,1] <- "Trait names"
if(t_choice=="data"){  res[n,2] <- ncol(TD)}
if(t_choice=="data_2"){  res[n,4] <- ncol(TD)}
#res[n,2] <- colnames(TD)
#res[n,4] <- colnames(TD)
#res[3,3] <- colnames(EnvTot[,-which(colnames(EnvTot)=="ObservationID")])

#-----------------------------------------------------------------------------
# individualy
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "Individuals"
if(t_choice=="data"){res[n,2] <- nrow(taxTD)}
if(t_choice=="data_2"){res[n,4] <- nrow(taxTD)}
res[n,3] <- nrow(taxEnv)

#-----------------------------------------------------------------------------
# species
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "Species"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$Species))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$Species))}
res[n,3] <- length(unique(taxEnv$Species))
n=n+1
res[n,1] <- "Samples per Specie"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Species),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Species),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$Species),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for species"
if(t_choice=="data"){res[n,2] <- sum(taxTD$Species%in%taxEnv$Species[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$Species%in%taxEnv$Species[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for species"
a <- aggregate(x = EnvTot[!taxEnv$Species%in%taxTD$Species,],by = list(taxEnv$Species[!taxEnv$Species%in%taxTD$Species]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for species and traits"
a <- aggregate(x = EnvTot[!taxEnv$Species%in%taxTD$Species,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$Species[!taxEnv$Species%in%taxTD$Species]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}

#-----------------------------------------------------------------------------
# Genera
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "Genera"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$Genus))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$Genus))}
res[n,3] <- length(unique(taxEnv$Genus))
n=n+1
res[n,1] <- "Samples per Genus"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Genus),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Genus),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$Genus),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for Genera"
if(t_choice=="data"){res[n,2] <- sum(taxTD$Genus%in%taxEnv$Genus[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$Genus%in%taxEnv$Genus[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for Genus"
a <- aggregate(x = EnvTot[!taxEnv$Genus%in%taxTD$Genus,],by = list(taxEnv$Genus[!taxEnv$Genus%in%taxTD$Genus]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for Genus and traits"
a <- aggregate(x = EnvTot[!taxEnv$Genus%in%taxTD$Genus,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$Genus[!taxEnv$Genus%in%taxTD$Genus]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}

#-----------------------------------------------------------------------------
# Family
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "Family"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$Family))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$Family))}
res[n,3] <- length(unique(taxEnv$Family))
n=n+1
res[n,1] <- "Samples per Family"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Family),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Family),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$Family),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for Family"
if(t_choice=="data"){res[n,2] <- sum(taxTD$Family%in%taxEnv$Family[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$Family%in%taxEnv$Family[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for Family"
a <- aggregate(x = EnvTot[!taxEnv$Family%in%taxTD$Family,],by = list(taxEnv$Family[!taxEnv$Family%in%taxTD$Family]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for Family and traits"
a <- aggregate(x = EnvTot[!taxEnv$Family%in%taxTD$Family,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$Family[!taxEnv$Family%in%taxTD$Family]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}

#-----------------------------------------------------------------------------
# Clade
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "Clade"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$Clade))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$Clade))}
res[n,3] <- length(unique(taxEnv$Clade))
n=n+1
res[n,1] <- "Samples per Clade"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Clade),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$Clade),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$Clade),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for Clade"
if(t_choice=="data"){res[n,2] <- sum(taxTD$Clade%in%taxEnv$Clade[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$Clade%in%taxEnv$Clade[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for Clade"
a <- aggregate(x = EnvTot[!taxEnv$Clade%in%taxTD$Clade,],by = list(taxEnv$Clade[!taxEnv$Clade%in%taxTD$Clade]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for Clade and traits"
a <- aggregate(x = EnvTot[!taxEnv$Clade%in%taxTD$Clade,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$Clade[!taxEnv$Clade%in%taxTD$Clade]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}


#-----------------------------------------------------------------------------
# GF
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "GF"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$GF))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$GF))}
res[n,3] <- length(unique(taxEnv$GF))
n=n+1
res[n,1] <- "Samples per GF"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$GF),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$GF),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$GF),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for GFs"
if(t_choice=="data"){res[n,2] <- sum(taxTD$GF%in%taxEnv$GF[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$GF%in%taxEnv$GF[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for GF"
a <- aggregate(x = EnvTot[!taxEnv$GF%in%taxTD$GF,],by = list(taxEnv$GF[!taxEnv$GF%in%taxTD$GF]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for GF and traits"
a <- aggregate(x = EnvTot[!taxEnv$GF%in%taxTD$GF,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$GF[!taxEnv$GF%in%taxTD$GF]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}

#-----------------------------------------------------------------------------
# PFT
#-----------------------------------------------------------------------------
n=n+1
res[n,1] <- "PFT"
if(t_choice=="data"){res[n,2] <- length(unique(taxTD$PFT))}
if(t_choice=="data_2"){res[n,4] <- length(unique(taxTD$PFT))}
res[n,3] <- length(unique(taxEnv$PFT))
n=n+1
res[n,1] <- "Samples per PFT"
if(t_choice=="data"){res[n,2] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$PFT),FUN = count.fun)[,-1]))}
if(t_choice=="data_2"){res[n,4] <- mean(colMeans(aggregate(x = TD,by = list(taxTD$PFT),FUN = count.fun)[,-1]))}
a <- aggregate(x = EnvTot,by = list(taxEnv$PFT),FUN = count.fun)
res[n,3] <- mean(rowSums(a[,-1]))
n=n+1
res[n,1] <- "Additional info from Env for PFT"
if(t_choice=="data"){res[n,2] <- sum(taxTD$PFT%in%taxEnv$PFT[-which(taxEnv$ID%in%taxTD$ID)])}
if(t_choice=="data_2"){res[n,4] <- sum(taxTD$PFT%in%taxEnv$PFT[-which(taxEnv$ID%in%taxTD$ID)])}
n=n+1
res[n,1] <- "Additional info from Env for PFT"
a <- aggregate(x = EnvTot[!taxEnv$PFT%in%taxTD$PFT,],by = list(taxEnv$PFT[!taxEnv$PFT%in%taxTD$PFT]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}
n=n+1
res[n,1] <- "Additional info from Env for PFT and traits"
a <- aggregate(x = EnvTot[!taxEnv$PFT%in%taxTD$PFT,colnames(EnvTot)%in%colnames(TD)],by = list(taxEnv$PFT[!taxEnv$PFT%in%taxTD$PFT]),FUN = count.fun)
if(t_choice=="data"){res[n,2] <- mean(rowSums(a[,-1]))}
if(t_choice=="data_2"){res[n,4] <-  mean(rowSums(a[,-1]))}

n

res <- data.frame(res)
res$TD=round(as.numeric(res$TD),digits = 1)
res$Envelope.data=round(as.numeric(res$Envelope.data),digits = 1)
res$TD2=round(as.numeric(res$TD2),digits = 1)

require(xtable)
xtable(res,include.rownames=FALSE,digits = 1)
print(xtable(res,digits = 1), include.rownames=FALSE)

