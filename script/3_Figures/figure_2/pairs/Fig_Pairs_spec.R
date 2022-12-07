
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
repnums=3

colz2=c("#d53e4f","#f46d43","#fdae61",
        "#fee08b","#ffffbf","#e6f598","#abdda4",
        "#66c2a5","#3288bd")
#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
t_choice="data_2"
if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}
tx=0
for(tx in 0:4){
tx=tx+1

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(originData,"analyes","Point_wise","res.csv")
units <- c("mm2 mg-1","m","","","mm2")
colz=colz1

#load data
{
  # load Extended data
  # load TDenvelope
  ObsOrTD="Obs_obs"
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TDtd <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                         "Obs_obs_TD","data","traitInfo_pred.csv"),header=TRUE))[,-c(1,2)]
  TD_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                 "Obs_obs_TD","data","taxInfo.csv"), sep=",")
  colnames(TD_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  head(TD_tax)
  colnames(TD_tax)
  dim(TD_tax)
  # load Extended data
  # total trait data 
  Env <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                             "Obs_obs","data","traitInfo.csv"),header=TRUE))[,-1]
  summary(Env)
  Env <- Env[,colnames(Env)%in%colnames(TD)]
  summary(Env)
  Env_tax <- read.table(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","0"),
                                 "Obs_obs","data","taxInfo.csv"), sep=",")
  colnames(Env_tax) <- c("ObservationID","Species","Genus","Family","Clade")
  length(Env_tax[,tx]%in%TD_tax[,tx])
  sum(Env_tax[,tx]%in%TD_tax[,tx])
  dim(Env)
  Env <- Env[Env_tax[,tx]%in%TD_tax[,tx],]
  summary(Env)
  Env_tax_c <- Env_tax[Env_tax[,tx]%in%TD_tax[,tx],]
  length(unique(Env_tax_c[,tx]))
  length(unique(TD_tax[,tx]))
  
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),ObsOrTD,"data"))
  # predicted 
  TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
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

#input=c(1,2,3)
#mean.fun(input)
#rm("input")
mean.fun <- function(input){
  if(sum(!is.na(input))==0){return(NA)}else{
  return(mean(input[!is.na(input)]))}
}
count.fun <- function(input){
  return(sum(!is.na(input)))
}


Env[Env_tax_c$ObservationID%in%TD_tax$ObservationID,]

# aggregate TD to taxon mean
dim(TD)
dim(TD_tax[,tx])
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
head(Env_sparse_tx)

#zlog
{
  res <- list()
  res_zlog <- list()
  res$TD_tx <- TD_tx
  res$TD_sparse_tx <- TD_sparse_tx
  res$TDtd_tx <- TDtd_tx
  res$TDenv_tx <- TDenv_tx
  res$Env_sparse_tx <- Env_sparse_tx
  
  dat_now=res$TD_tx
  mn=apply(log(dat_now),MARGIN = 2,new.mean.fun)
  sd=apply(log(dat_now),MARGIN = 2,new.sd.fun)

  dat_sel=5
  for(dat_sel in 1:length(res)){
    dat_now=res[[dat_sel]]
    for(t in 1:ncol(dat_now)){
      dat_now[,t] <- (log(res[[dat_sel]][,t])-mn[t])/sd[t]
    }
    res_zlog[[dat_sel]] <- dat_now
  }
  
  
  TD_tx <- res_zlog[[1]]
  TD_sparse_tx <- res_zlog[[2]]
  TDtd_tx <- res_zlog[[3]]
  TDenv_tx <- res_zlog[[4]]
  Env_sparse_tx <- res_zlog[[5]]
  
}

head(Env_sparse_tx)
dim(TD)


plot(Env_sparse_tx$SLA,TD_sparse_tx$SLA)

plot_dat <- cbind(unlist(TD_tx),
                  unlist(TD_sparse_tx),
                  unlist(Env_sparse_tx),
                  unlist(TDtd_tx),
                  unlist(TDenv_tx))
colnames(plot_dat) <- c("TD","TD_sparse","ExTD_sparse","TDtd","TDext")
require(lattice)
require(ggplot2)
panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y,use="pairwise.complete.obs")) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y,use="pairwise.complete.obs") 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}
if(tx==1&t_choice=="data"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_indMean.pdf"))}
if(tx==2&t_choice=="data"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_specMean.pdf"))}
if(tx==3&t_choice=="data"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_genMean.pdf"))}
if(tx==4&t_choice=="data"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_famMean.pdf"))}
if(tx==5&t_choice=="data"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_cladMean.pdf"))}
if(tx==1&t_choice=="data_2"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_indMean_TD2.pdf"))}
if(tx==2&t_choice=="data_2"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_specMean_TD2.pdf"))}
if(tx==3&t_choice=="data_2"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_genMean_TD2.pdf"))}
if(tx==4&t_choice=="data_2"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_famMean_TD2.pdf"))}
if(tx==5&t_choice=="data_2"){pdf(file = file.path(origin,"figures","Figure_2","Pairs_cladMean_TD2.pdf"))}
if(t_choice=="data"){colnames(plot_dat) <- c("OBS", "OBSsparse","TRY17_sparse","IMPobs","IMPobsext")}
if(t_choice=="data_2"){colnames(plot_dat) <- c("OBS2", "OBS2sparse","TRY17_2_sparse","IMP2obs","IMP2obsext")}
pairs(plot_dat, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()
}

#require(FactoMineR)
#PCA(plot_dat)

