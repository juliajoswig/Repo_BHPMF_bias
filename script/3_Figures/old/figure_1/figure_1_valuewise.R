
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
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs"
Percent=80
resENV <- read.csv(file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

colz=c("#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
colzENV="#b2182b"

colz=c("#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#b2182b")
res <- res[,colSums(!is.na(res))!=0]
trait_names=as.vector(unique(res$trait))
trait_names <- trait_names[!is.na(trait_names)]
missingness = unique(as.vector(res$missingness))
missingness <- missingness[!is.na(missingness)]
summary(res_now$value_obs)
m=1
t=2
w=1


units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")

pdf(file=file.path(origin,"_2021","figures","Figure_1","Fig_1_ValueWise.pdf"),width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4,2,2))
par(mfrow=c(1,6),mar=c(7,6,2,1))

limmin <- c(-30,-42,-.5,-25,-5,-2)
limmax <- c(30,42,.5,25,5,2)
ix_trait=res$trait==trait_names[t]
ix_trait[is.na(ix_trait)] <- FALSE
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  bxpl=res$value_pred[ix_trait]- res$value_obs[ix_trait]
  bxplENV=resENV$value_pred[ix_trait]- resENV$value_obs[ix_trait]
  dat_plot <- cbind(bxpl,bxplENV)
  colnames(dat_plot) <- c("Test data","Envelope")
  boxplot(dat_plot,
          ylim=c(quantile(bxpl,probs = .01),quantile(bxpl,probs = .99)),
          ylab=paste0("Error ",trait_names[t]," [",units[t],"]"),
          las=2,col=c(colz[t],colzENV),cex.axis=1.5,cex.lab=2.5,tick = FALSE,frame=FALSE)#width=1,
  abline(h=0)
}   
dev.off()

trait_names

pdf(file=file.path(origin,"_2021","figures","figure_1","Fig_1_ValueWiseGAPSonly.pdf"),width=10,height =8)
par(mfrow=c(1,1),mar=c(4,4,2,2))
par(mfrow=c(1,6),mar=c(0,6,0,0))

limmin <- c(-30,-42,-.5,-25,-5,-2)
limmax <- c(30,42,.5,25,5,2)
ix_trait=res$trait==trait_names[t]
ix_trait[is.na(ix_trait)] <- FALSE
t=2
for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$Gap
    ix_trait[is.na(ix_trait)] <- FALSE
     bxpl=res$value_pred[ix_trait]- res$value_obs[ix_trait]
    boxplot(bxpl,ylim=c(quantile(bxpl,probs = .05),quantile(bxpl,probs = .95)),ylab=paste0("Error ",trait_names[t]," [",units[t],"]"),width=1,las=2,col=colz[t],cex.axis=2,cex.lab=3,tick = FALSE,frame=FALSE)
    abline(h=0)
    
    }   
dev.off()


require(xtable)
# get some numbers for the paper:
table_now <- matrix()
ix_trait=res$trait==trait_names[t]
ix_trait[is.na(ix_trait)] <- FALSE
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  bxpl=res$value_pred[ix_trait]- res$value_obs[ix_trait]
  bxplENV=resENV$value_pred[ix_trait]- resENV$value_obs[ix_trait]
  quantile(abs(bxpl),probs = .25)
  quantile(abs(bxpl),probs = .5)
  quantile(abs(bxpl),probs = .75)
  quantile(abs(bxplENV),probs = .25)
  quantile(abs(bxplENV),probs = .5)
  quantile(abs(bxplENV),probs = .75)
}   

