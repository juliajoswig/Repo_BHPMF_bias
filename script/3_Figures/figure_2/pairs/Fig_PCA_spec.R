
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
repnums=3

colz2=c("#d53e4f","#f46d43","#fdae61",
        "#fee08b","#ffffbf","#e6f598","#abdda4",
        "#66c2a5","#3288bd")
#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
RepNum=1
t_choice="data"


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
  length(Env_tax$Species%in%TD_tax$Species)
  sum(Env_tax$Species%in%TD_tax$Species)
  dim(Env)
  Env <- Env[Env_tax$Species%in%TD_tax$Species,]
  Env_tax_c <- Env_tax[Env_tax$Species%in%TD_tax$Species,]
  length(unique(Env_tax_c$Species))
  length(unique(TD_tax$Species))
  
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


TD_sparse
Env[Env_tax_c$ObservationID%in%TD_tax$ObservationID,]

# aggregate TD to taxon mean
tx=3
tx=tx+1
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

plot(Env_sparse_tx$SLA,TD_sparse_tx$SLA)

plot_dat <- cbind(unlist(TD_tx),
                  unlist(TD_sparse_tx),
                  unlist(Env_sparse_tx),
                  unlist(TDtd_tx),
                  unlist(TDenv_tx))
colnames(plot_dat) <- c("TD","TD_sparse","Env_sparse","TDtd","TDenv")
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
pdf(file = file.path(origin,"_2021","figures","Figure_2","Pairs.pdf"))
  pairs(plot_dat, lower.panel=panel.smooth, upper.panel=panel.cor)
dev.off()

#require(FactoMineR)
#PCA(plot_dat)
