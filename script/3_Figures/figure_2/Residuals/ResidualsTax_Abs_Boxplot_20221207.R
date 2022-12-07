# atm S7: Relative error to trait range for test data 2. Error (absolute) calculated
# as imputed - observed for test data imputed with and without the envelope (back-transformed data 


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
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  GapPercent=50
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3


#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
t_choice="data"
  
  if(t_choice=="data"){RepNum=1}
  if(t_choice=="data_2"){RepNum=2}
  
  ObsOrTD="Obs_obs_TD"
  Percent=80
  res <- read.csv(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))
  
  
  ObsOrTD="Obs_obs"
  Percent=80
  resENV <- read.csv(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))
  
  
  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  
  colnames(res)
# calc residuals
  # Species
  
  res$mean_spec_obs_zlog[which(res$Species==res$Species[121])]
  res$value_obs_zlog[which(res$Species==res$Species[121])]
    # get some numbers for the paper:
  table_now <- matrix(NA,ncol=8,nrow=length(trait_names))
  table_2 <- matrix(NA,ncol=8,nrow=length(trait_names))
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  t=1
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    ix_trait[is.na(ix_trait)] <- FALSE
    
    bxpl=res$value_pred[ix_trait]- res$value_obs[ix_trait]
    bxplENV=resENV$value_pred[ix_trait]- resENV$value_obs[ix_trait]
    
    table_now[t,1] <- quantile(abs(bxpl),probs = .25)
    table_now[t,2] <- quantile(abs(bxpl),probs = .5)
    table_now[t,3] <- quantile(abs(bxpl),probs = .75)
    table_now[t,4] <- max(abs(bxpl))
    table_now[t,5] <- quantile(abs(bxplENV),probs = .25)
    table_now[t,6] <- quantile(abs(bxplENV),probs = .5)
    table_now[t,7] <- quantile(abs(bxplENV),probs = .75)
    table_now[t,8] <- max(abs(bxplENV))
    
    print(quantile(abs(bxpl),probs = 1)<c(dist(range(abs(res$value_obs)[ix_trait]))))
    table_2[t,1] <- quantile(abs(bxpl),probs = .25)/c(dist(range(res$value_obs[ix_trait])))*100
    table_2[t,2] <- quantile(abs(bxpl),probs = .5)/c(dist(range(res$value_obs[ix_trait])))*100
    table_2[t,3] <- quantile(abs(bxpl),probs = .75)/c(dist(range(res$value_obs[ix_trait])))*100
    table_2[t,4] <- quantile(abs(bxpl),probs = 1)/c(dist(range(res$value_obs[ix_trait])))*100
    table_2[t,5] <- quantile(abs(bxplENV),probs = .25)/c(dist(range(resENV$value_obs[ix_trait])))*100
    table_2[t,6] <- quantile(abs(bxplENV),probs = .5)/c(dist(range(resENV$value_obs[ix_trait])))*100
    table_2[t,7] <- quantile(abs(bxplENV),probs = .75)/c(dist(range(resENV$value_obs[ix_trait])))*100
    table_2[t,8] <- quantile(abs(bxplENV),probs = 1)/c(dist(range(resENV$value_obs[ix_trait])))*100
  }   
  

   # par(mfrow=c(1,1),mar=c(12,4,1,1))
    pdf(file=file.path(origin,"figures","Figure_2",paste0(txs[tx],"_boxplots_",t_choice,".pdf")),width=8,height=2.3,pointsize = 9)
    par(mfrow=c(1,6),mar=c(12,4,1,1))
    t=1

    for(t in 1:ncol(TD)){
      head(bxplTDenv_TD)
      
      # (bxplTDenv_Envsparse),  (bxplTDenv_TD))[seq(from = t,to = ncol(TD)*5,by = ncol(TD))]
      dat_plot=cbind((bxplTDtd_TD),(bxplTDtd_TDsparse),(bxplTDtd_Envsparse),(bxplTDenv_Envsparse),(bxplTDenv_TD))[seq(from = t,to = ncol(TD)*5,by = ncol(TD))]
      colnames(dat_plot) =  c("IMP_obs - OBS","IMP_obs - OBS_sparse", 
                             "IMP_obs - TRY17_sparse", "IMP_obsExt - TRY17_sparse",  "IMP_obsExt - OBS")
      dat_plot <- abs(dat_plot)
      
      boxplot(dat_plot,main=trait_names[t],ylim=c(0,
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.95,na.rm=TRUE))),col=colz,las=2,ylab=paste("abs(zlog) distance"))
      rect(xleft = 0,ybottom = min(apply(dat_plot,MARGIN = 2,quantile,prob=0,na.rm=TRUE))-1,xright = 6,ytop = median(dat_plot[,1],na.rm = TRUE),col=colz2[1])
      rect(xleft = 0,ybottom = median(dat_plot[,1],na.rm = TRUE),xright = 6,ytop = max(apply(dat_plot,MARGIN = 2,quantile,prob=1,na.rm=TRUE))+1,col=colz2[3])
      abline(h = median(dat_plot[,1],na.rm = TRUE),col=colz2[8],lwd=1)
      abline(v = 3.5,col="white",lwd=1.5)
      boxplot(dat_plot,main=trait_names[t],ylim=c(min(apply(dat_plot,MARGIN = 2,quantile,prob=.01,na.rm=TRUE)),
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.99,na.rm=TRUE))),
              col=colz,las=2,ylab=paste("zlog(",units[t],")"),
              add=TRUE)
    }   
    dev.off()
    
  }
  
}


dev.off()



{
  # load Envelope data
  # load TDenvelope
  ObsOrTD="Obs_obs"
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
  TDtd <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                           "Obs_obs_TD","data","traitInfoTD_pred.csv"),header=TRUE))[,-c(1,2)]
  TD_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                 "Obs_obs_TD","data","taxInfo.csv"), sep=",")
  colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  head(TD_tax)
  colnames(TD_tax)
  dim(TD_tax)
  # load Envelope data
  # total trait data 
  Env <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                          "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  Env <- Env[,colnames(Env)%in%colnames(TD)]
  sum(colnames(Env)==colnames(TD))==ncol(Env)
  
  Env_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                  "Obs_obs","data","taxInfo.csv"), sep=",")
  colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  dim(Env_tax)==dim(Env)
  Env <- Env[Env_tax[,tx]%in%TD_tax[,tx],]
  Env_tax_c <- Env_tax[Env_tax[,tx]%in%TD_tax[,tx],]
  
  mtch <- match(TD_tax[,tx],Env_tax_c[,tx])
  print(sum(TD_tax[,tx]==Env_tax_c[mtch,tx])==length(TD_tax[,tx]))
  Env_tax_c<- Env_tax_c[mtch,]
  Env<- Env[mtch,]
  
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  # predicted 
  TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv
  head(TDtd)
  head(TDenv)
  head(TD_sparse)
  head(Env)
  par(mfrow=c(1,1))
  plot(TD[,1],TDtd[,1])
  abline(0,1)
  plot(TD[,1],Env[,1])      
  abline(0,1)
}