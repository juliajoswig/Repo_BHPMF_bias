
#------------------------------------------------------------
# define path
#------------------------------------------------------------
setwd("/..")
origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git"
originData = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data"
list.files(file.path(origin,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","helper_scripts","fn_load_functions.R"))
load_functions(origin)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)
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

  #-------------------------------------------------------------------
  # chose trait data   
  #-------------------------------------------------------------------
  repnums=3
  t_choice="data_2"
  GapPercent=80
  RepNum=1
  
  
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(originData,"analyes","Point_wise","res.csv")
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data"){trait_names=trait_rainfor}
if(t_choice=="data_2"){trait_names=trait_guido}

colz=colz1[c(1,5,2,3,4)]


{
  # load Envelope data
  # load TDenvelope
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  taxTD <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_TD <- taxTD[,1]
  head(taxTD)
  dim(taxTD)
  dim(TD)
  dim(TD_sparse)
  # load Envelope data
  # total trait data 
  EnvTot <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                             "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  Envsparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  taxEnv <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs","data","taxInfo.csv"),header=TRUE))
  taxEnv <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                     "Obs_obs","data","taxInfo.csv"),
                                    sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
  ID_env <- taxEnv[,1]
  dim(taxEnv)
  summary(EnvTot)
  Env <- EnvTot[,colnames(EnvTot)%in%colnames(TD)]
  Env_sp <- Envsparse[,colnames(Envsparse)%in%colnames(TD)]
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  
  head(TDtd)
  head(TDenv)
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],Env[match(taxTD$ID,ID_env),1])
  abline(0,1)
  plot(TD[,1],Env_sp[match(taxTD$ID,ID_env),1])
  plot(TD[,1],TDenv[,1])
  abline(0,1)
  plot(TD_sparse[,1],TD[,1])
  abline(0,1)
}

#-------------------------------------------------------------------
# chose 
#-------------------------------------------------------------------
dev.off()
plot(TDtd[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDenv[,1])
abline(0,1)
plot(TD[,1],TDtd[,1])
abline(0,1)
length(ID_env)
length(ID_TD)

plot(Env[match(ID_TD,ID_env),5],TDtd[,5])
plot(Env[match(ID_TD,ID_env),1],TDenv[,1])
abline(0,1)
print("Yaaay loaded :)")

# create Ext
# cut OBS from Env to receive Ext(202212)
#Ext <- Env[!ID_env%in%match(ID_TD,ID_env),]

# table
m1 <- matrix(NA,ncol=6,nrow=length(trait_names))
m2 <- matrix(NA,ncol=6,nrow=length(trait_names))
m3 <- matrix(NA,ncol=6,nrow=length(trait_names))
m4 <- matrix(NA,ncol=6,nrow=length(trait_names))
colnames(m1) <- c("OBS","OBS_sparse","IMP_obs","TRY17","IMP_obsExt","TRY17_sparse")
colnames(m2) <-  c("OBS","OBS_sparse","IMP_obs","TRY17","IMP_obsExt","TRY17_sparse")
colnames(m3) <-  c("OBS","OBS_sparse","IMP_obs","TRY17","IMP_obsExt","TRY17_sparse")
colnames(m4) <-  c("OBS","OBS_sparse","IMP_obs","TRY17","IMP_obsExt","TRY17_sparse")
rownames(m1) <- trait_names
rownames(m2) <- trait_names
rownames(m3) <- trait_names
rownames(m4) <- trait_names

par(mfrow=c(1,1),mar=c(7,2,2,2))

t=1
for(t in 1:length(trait_names)){
  m1[t,1]<- quantile(TD[,colnames(TD)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,2]<- quantile(TD_sparse[,colnames(TD)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,3]<- quantile(TDtd[,colnames(TDtd)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,4] <- quantile(Env[,colnames(Env)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,5] <- quantile(TDenv[,colnames(TDenv)==trait_names[t]],probs = .5,na.rm = TRUE)
  m1[t,6] <- quantile(Env_sp[,colnames(Env_sp)==trait_names[t]],probs = .5,na.rm = TRUE)
  
  m2[t,1]<- max(x = TD[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m2[t,2]<- max(x = TD_sparse[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m2[t,3]<- max(x = TDtd[,colnames(TDtd)==trait_names[t]],na.rm = TRUE)
  m2[t,4] <-max(x = Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
  m2[t,5] <- max(x = TDenv[,colnames(TDenv)==trait_names[t]],na.rm = TRUE)
  m2[t,6] <- max(x = Env_sp[,colnames(Env_sp)==trait_names[t]],na.rm = TRUE)
  
  m3[t,1]<- min(x = TD[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m3[t,2]<- min(x = TD_sparse[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m3[t,3]<- min(x = TDtd[,colnames(TDtd)==trait_names[t]],na.rm = TRUE)
  m3[t,4] <-min(x = Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
  m3[t,5] <- min(x = TDenv[,colnames(TDenv)==trait_names[t]],na.rm = TRUE)
  m3[t,6] <- min(x = Env_sp[,colnames(Env_sp)==trait_names[t]],na.rm = TRUE)
  
  m4[t,1]<- mean(x = TD[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m4[t,2]<- mean(x = TD_sparse[,colnames(TD)==trait_names[t]],na.rm = TRUE)
  m4[t,3]<- mean(x = TDtd[,colnames(TDtd)==trait_names[t]],na.rm = TRUE)
  m4[t,4] <- mean(Env[,colnames(Env)==trait_names[t]],na.rm = TRUE)
  m4[t,5] <- mean(TDenv[,colnames(TDenv)==trait_names[t]],na.rm = TRUE)
  m4[t,6] <- mean(Env_sp[,colnames(Env_sp)==trait_names[t]],na.rm = TRUE)

}
  
getwd()
setwd(origin)
list.files(file.path("figures","Figure_2"))

head(m4)


m_now=m4[,c(1,2,3,4,6,5)]
#m_now=m4[,c(1,3,4,5)]
if(t_choice=="data"){
#pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_2_Mean_r.pdf"),width=18,height=5)
  pdf(file=file.path("figures","Figure_2","Fig_2_Mean_r.pdf"),width=18,height=5)
}

if(t_choice=="data_2"){
  #pdf(file=file.path(origin,"_2021","figures","figure_2","Fig_S2_Mean_r.pdf"),width=18,height=4)
  pdf(file=file.path("figures","Figure_2","Fig_S2_Mean_r.pdf"),width=18,height=4)
}
if(t_choice=="data_2"){
  colnames(m_now) <- c("OBS2","OBS2_sparse","IMP2_obs","TRY17","TRY17_2_sparse", "IMP2_obsExt")
}
{
  par(mfrow=c(1,length(trait_names)),mar=c(12.5,6,3,1))
  
  colz=c("#d1e5f0",#TD
         "#d1e5f0",#TDsparse
         "#ef8a62",#TDtd
         "#2166ac",#Env
         "#2166ac",#Env
         #"#67a9cf",#Env_sp
         "#b2182b",#TDenv
         "#b2182b"#TDenv
  )
  t=1
  for(t in 1:length(trait_names)){
    barplot(m_now[t,],col=colz,border = "white",ylab="",main=trait_names[t],
            cex.names = 1.8,las=2,yaxt="n",cex.main=3,cex.lab=2.5,
            cex.axis=2,width = c(1,0.5,1,1,0.5,1),ylim=c(0,max(m_now[t,])*1.2))
    axis(2,line = 2)
    axis(2,line = 3.5,at=max(m_now[t,])/2, labels = units[t],cex.axis=1.8,tick=FALSE)
    #abline(h=seq(0,max(m_now[t,]),length.out=6),col="gray",lty=2)
    #barplot(m_now[t,],col=colz,border = "white",add=TRUE)
  }
  
}

dev.off()
