# atm S7: Relative error to trait range for test data 2. Error (absolute) calculated
# as imputed - observed for test data imputed with and without the envelope (back-transformed data 

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
zlog=FALSE
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}

if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}
tx=2 #(0=indiv,1=spec,2=genus,3=family,4=clades)
#for (tx in 0:4)
par(mfrow=c(4,4))
{
  tx=tx+1
  if(t_choice=="data"){trait_names=trait_rainfor}
  if(t_choice=="data_2"){trait_names=trait_guido}
  
  #-------------------------------------------------------------------
  # load trait data   
  #-------------------------------------------------------------------
  file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
  units <- c("mm2 mg-1","m","","","mm2")
  colz=colz1
  
  #load data
  {
    # load Envelope data
    # load TDenvelope
    ObsOrTD="Obs_obs"
    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
    
    # load TD data
    # total trait data 
    TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                           "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
    TD_sparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                  "Obs_obs_TD","data","traitInfoTD_obs.csv"),header=TRUE))[,-c(1,2)]
    TDtd <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                             "Obs_obs_TD","data","traitInfoTD_pred.csv"),header=TRUE))[,-c(1,2)]
    TD_tax <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                   "Obs_obs_TD","data","taxInfo.csv"), sep=",")
    colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
    head(TD_tax)
    colnames(TD_tax)
    dim(TD_tax)
    # load Envelope data
    # total trait data 
    Env <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
    summary(Env)
    Env <- Env[,colnames(Env)%in%colnames(TD)]
    Env_tax <- read.table(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                    "Obs_obs","data","taxInfo.csv"), sep=",")
    colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
    length(Env_tax[,tx]%in%TD_tax[,tx])
    sum(Env_tax[,tx]%in%TD_tax[,tx])
    dim(Env)
    Env <- Env[Env_tax[,tx]%in%TD_tax[,tx],]
    Env_tax_c <- Env_tax[Env_tax[,tx]%in%TD_tax[,tx],]
    length(unique(Env_tax_c[,tx]))
    length(unique(TD_tax[,tx]))
    
    list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
    # predicted 
    TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                          "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
    TDenv
    head(TDtd)
    head(TDenv)
    plot(TD[,1],TDtd[,1])
    abline(0,1)
    plot(TD[,1],TDenv[,1])
    abline(0,1)
  }
  
  trait_names=colnames(TD)
  
  input=c(1,2,3)
  mean.fun(input)
  rm("input")
  mean.fun <- function(input){
    if(sum(!is.na(input))==0){return(NA)}else{
      return(mean(input[!is.na(input)]))}
  }
  count.fun <- function(input){
    return(sum(!is.na(input)))
  }
  input=c(1:10,NA)
  rm("input")
  range.fun <- function(input){
    if(sum(!is.na(input))==0){return(NA)}else{
      return(max(input[!is.na(input)])-min(input[!is.na(input)]))
    }
  }
  
  
  # aggregate TD to taxon mean
  TD_tx <- aggregate(TD,by=list(TD_tax[,tx]),FUN = mean.fun)
  TD_sparse_tx <- aggregate(TD_sparse,by=list(TD_tax[,tx]),FUN = mean.fun)
  TDtd_tx <- aggregate(TDtd,by=list(TD_tax[,tx]),FUN = mean.fun)
  TDenv_tx <- aggregate(TDenv,by=list(TD_tax[,tx]),FUN = mean.fun)
  Env_sparse_tx <- aggregate(Env,by=list(Env_tax_c[,tx]),FUN = mean.fun)
  
  TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
  TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
  TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
  TDenv_tx <- TDenv_tx[,-which(colnames(TDenv_tx)=="Group.1")]
  Env_sparse_tx <- Env_sparse_tx[,-which(colnames(Env_sparse_tx)=="Group.1")]
  

# get some numbers for the paper:
  table_now <- matrix(NA,ncol=8,nrow=length(trait_names))
  table_2 <- matrix(NA,ncol=8,nrow=length(trait_names))
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  t=5
  
  for(t in 1:length(trait_names)){
    whichcol=which(colnames(TD_tx)==trait_names[t])
    
    print(dim(TD_tx))
    if(zlog){
      TDmn=mean(log(TD_tx[,whichcol]),na.rm = TRUE)
      TDsd=mean(log(TD_tx[,whichcol]),na.rm = TRUE)
      
      TDtd_now <- (log(TDtd_tx[,whichcol])-TDmn)/TDsd
      TD_now <- (log(TD_tx[,whichcol])-TDmn)/TDsd

      plot(TD_now,TDtd_now);abline(0,1)
      bxpl=TDtd_now-TD_now
      plot(abs(bxpl),ylim=c(0,1))
      
      TDenv_now <- (log(TDenv_tx[,whichcol])-TDmn)/TDsd
      bxplENV = TDenv_now - TD_now
      plot(TD_now,TDenv_now);abline(0,1)
      plot(abs(bxplENV),ylim=c(0,1))
      
    }else{
    bxpl=TDtd_tx[,whichcol]- TD_tx[,whichcol]
    bxplENV=TDenv_tx[,whichcol]- TD_tx[,whichcol]
    }
    table_now[t,1] <- quantile(abs(bxpl),probs = .25)
    table_now[t,2] <- quantile(abs(bxpl),probs = .5)
    table_now[t,3] <- quantile(abs(bxpl),probs = .75)
    table_now[t,4] <- max(abs(bxpl))
    table_now[t,5] <- quantile(abs(bxplENV),probs = .25)
    table_now[t,6] <- quantile(abs(bxplENV),probs = .5)
    table_now[t,7] <- quantile(abs(bxplENV),probs = .75)
    table_now[t,8] <- max(abs(bxplENV))
    
    print(quantile(abs(bxpl),probs = 1)< range.fun(abs(TD_tx[, whichcol])))
    table_2[t,1] <- quantile(abs(bxpl),probs = .25)/range.fun(TD_tx[, whichcol])*100
    table_2[t,2] <- quantile(abs(bxpl),probs = .5)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,3] <- quantile(abs(bxpl),probs = .75)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,4] <- quantile(abs(bxpl),probs = 1)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,5] <- quantile(abs(bxplENV),probs = .25)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,6] <- quantile(abs(bxplENV),probs = .5)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,7] <- quantile(abs(bxplENV),probs = .75)/ range.fun(TD_tx[, whichcol])*100
    table_2[t,8] <- quantile(abs(bxplENV),probs = 1)/ range.fun(TD_tx[,whichcol])*100
}   


  print(table_now[,c(4,8)])
  print(table_2[,c(4,8)])

#table_now <- table_now[c(-4,8)]
#table_2 <- table_2[c(-4,8)]
if(t_choice=="data"){
colnames(table_now) <- c("25th quantile TD","Median TD","75th quantile TD","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
rownames(table_now) <- trait_names
colnames(table_2) <- c("25th quantile TD","Median TD","75th quantile TD","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
}
if(t_choice=="data_2"){
  colnames(table_now) <- c("25th quantile TD2","Median TD2","75th quantile TD2","Max TD2",
                           "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
  rownames(table_now) <- trait_names
  colnames(table_2) <- c("25th quantile TD2","Median TD2","75th quantile TD2","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
}
rownames(table_2) <- trait_names
require(xtable)
}
#Absolute residuals of imputed - observed for test data imputed with and without the envelope (back-transformed data).
xtable(table_now,digits = 2)

#Residuals relative to trait range. Residuals (absolute) calculated as  imputed - observed for test data imputed with and without the envelope (back-transformed data).
xtable(table_2,digits = 2)

