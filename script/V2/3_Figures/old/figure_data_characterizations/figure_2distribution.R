
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
units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")



RepNum=1
t_choice="data"
Percent=0
# define dat path
list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))
list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),ObsOrTD,"data"))

ObsOrTD="Obs_obs_TD"
path_obs0TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfo_zlog.csv")
path_obs0TD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfo.csv")
path_taxTD <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","taxInfo.csv")
ObsOrTD="Obs_obs"
path_obs0ENV <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfo_zlog.csv")
path_obs0ENV <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","traitInfo.csv")
path_taxENV <- file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,"p_0",ObsOrTD,"data","taxInfo.csv")

# load data
traitInfo_obsTD <- as.data.frame(read.table(file = path_obs0TD, sep=",", dec=".",header=TRUE))[,-c(1,2)]
traitInfo_obsENV <- as.data.frame(read.table(file = path_obs0ENV, sep=",", dec=".",header=TRUE))[,-c(1,2)]

traitInfo_obsENV <- (log(traitInfo_obsENV)-apply(log(traitInfo_obsTD),2,mean,na.rm=TRUE))/apply(log(traitInfo_obsTD),2,sd,na.rm=TRUE)
traitInfo_obsTD <- (log(traitInfo_obsTD)-apply(log(traitInfo_obsTD),2,mean,na.rm=TRUE))/apply(log(traitInfo_obsTD),2,sd,na.rm=TRUE)

taxInfo_obsTD <- as.data.frame(read.table(file = path_taxTD, sep=",", dec=".",header=FALSE))
taxInfo_obsENV <- as.data.frame(read.table(file = path_taxENV, sep=",", dec=".",header=FALSE))
colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
colz=c("#d1e5f0","#b2182b","#b2182b","#b2182b","#b2182b","#b2182b","#b2182b")
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))
colnames(taxInfo_obsTD) <- c("ObsID","Species","Genus","Family","Clade")
colnames(taxInfo_obsENV) <- c("ObsID","Species","Genus","Family","Clade")


units=c("mm2 mg-1","m","mg mm-3","mg g-1","mg g-1","g m-2") # data 2
#units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")# data 1


pdf(file=file.path(origin,"_2021","figures","Figure_Supplementary","fig_distribution_o.pdf"),width=12,height=8)
par(mfrow=c(2,3),mar=c(5,6,3,1))
t=2
for(t in 1:length(traitInfo_obsTD)){
  dens_datENV <- traitInfo_obsENV[,colnames(traitInfo_obsENV)%in%colnames(traitInfo_obsTD)[t]]
  ix <- !taxInfo_obsENV$ObsID%in%taxInfo_obsTD$ObsID
  nme <- colnames(traitInfo_obsTD)[t]
  dens_datTD <- traitInfo_obsTD[,t]
  
  boxplot(dens_datTD,dens_datENV[ix],
          dens_datENV[ix&taxInfo_obsENV$Species%in%taxInfo_obsTD$Species],
          dens_datENV[ix&taxInfo_obsENV$Genus%in%taxInfo_obsTD$Genus],
          dens_datENV[ix&taxInfo_obsENV$Family%in%taxInfo_obsTD$Family],
          dens_datENV[ix&taxInfo_obsENV$Clade%in%taxInfo_obsTD$Clade],col=colz,xaxt="n",main=nme,cex.axis=2,cex.lab=2,
          cex.main=2,ylab="zlog")
  abline(h=median(dens_datTD,na.rm=TRUE))
  axis(1,at = 1:6,labels = c("Test data","Total envelope","Same species","Same genus","Same family","Same clade"),las=2)
  
}
dev.off()

