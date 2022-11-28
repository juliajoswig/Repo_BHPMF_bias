

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
Version_now="V1"
list.files(file.path(origin,"_2021","script","analysis",Version_now))
#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin = origin,Version_now = Version_now)

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
TD_choice="Obs_obs"
trait_sub="rainfor"

colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858","black")
colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000","black")

GapPercent=0
list.files(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                     paste0("p",GapPercent,"_",trait_sub),TD_choice,"data"))

trait_obs <- read.csv(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                                paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                                "traitInfo_TD.csv"))
dim(trait_obs)

TD_choice="Obs_obs_TD"
trait_obs_TD <- read.csv(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                                paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                                "traitInfo.csv"))

trait_obs <- read.csv(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                                paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                                "traitInfo.csv"))
taxInfo <- read.csv(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                               paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                               "taxInfo.csv"))
trait_pred <- read.csv(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum),
                             paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                             "mean.csv"),sep = "\t")
GF <- read.table(file.path(origin,"_2021","data","runs","META","GF_Obs.csv"),sep=",",dec=".")
PFT <- read.table(file.path(origin,"_2021","data","runs","META","PFT_Obs.csv"),sep=",",dec=".")

taxInfo$GF <- GF$GF
taxInfo$PFT <- PFT$PFT
head(trait_pred)

dim(taxInfo)
dim(trait_obs)
head(trait_obs)
summary(trait_obs)

count_fun <- function(x){
  return(sum(!is.na(x)))
}
#calc Cova Species
mean_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = sd,na.rm=TRUE)
spec_count_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = count_fun)
spec_count_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = count_fun)

spec_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
spec_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova Genus
mean_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = sd,na.rm=TRUE)
gen_count_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = count_fun)
gen_count_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = count_fun)

gen_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
gen_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova Family
mean_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = sd,na.rm=TRUE)
fam_count_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = count_fun)
fam_count_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = count_fun)

fam_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
fam_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova phylogroup
mean_obs <- aggregate(trait_obs,by = list(taxInfo$PhylogeneticGroup),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$PhylogeneticGroup),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$PhylogeneticGroup),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$PhylogeneticGroup),FUN = sd,na.rm=TRUE)
pg_count_obs <- aggregate(trait_obs,by = list(taxInfo$PhylogeneticGroup),FUN = count_fun)
pg_count_pred <- aggregate(trait_pred,by = list(taxInfo$PhylogeneticGroup),FUN = count_fun)

pg_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
pg_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova GF
mean_obs <- aggregate(trait_obs,by = list(taxInfo$GF),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$GF),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$GF),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$GF),FUN = sd,na.rm=TRUE)
gf_count_obs <- aggregate(trait_obs,by = list(taxInfo$GF),FUN = count_fun)
gf_count_pred <- aggregate(trait_pred,by = list(taxInfo$GF),FUN = count_fun)

gf_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
gf_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova PFT
mean_obs <- aggregate(trait_obs,by = list(taxInfo$PFT),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$PFT),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$PFT),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$PFT),FUN = sd,na.rm=TRUE)
pft_count_obs <- aggregate(trait_obs,by = list(taxInfo$PFT),FUN = count_fun)
pft_count_pred <- aggregate(trait_pred,by = list(taxInfo$PFT),FUN = count_fun)

pft_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
pft_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#----------------

tot_obs <- rbind(spec_obs,gen_obs,fam_obs,pg_obs,gf_obs,pft_obs)
tot_pred <- rbind(spec_pred,gen_pred,fam_pred,pg_pred,gf_pred,pft_pred)

spec_count <- spec_count_obs[,2:19]==1&spec_count_pred[,2:19]!=1
gen_count <- gen_count_obs[,2:19]==1&gen_count_pred[,2:19]!=1
fam_count <- fam_count_obs[,2:19]==1&fam_count_pred[,2:19]!=1
pg_count <- pg_count_obs[,2:19]==1&pg_count_pred[,2:19]!=1
gf_count <- gf_count_obs[,2:19]==1&gf_count_pred[,2:19]!=1
pft_count <- pft_count_obs[,2:19]==1&pft_count_pred[,2:19]!=1

#write.csv(spec_obs,file=file.path(origin, "data_output","CoVa","spec_obs.csv"))
#write.csv(spec_obs,file=file.path(origin, "data_output","CoVa","spec_obs.csv"))
tot_obs <- rbind(data.frame(group=rep("Species",nrow(spec_obs)),cluster=spec_count_obs$Group.1,spec_obs),
                 data.frame(group=rep("Genus",nrow(gen_obs)),cluster=gen_count_obs$Group.1,gen_obs),
                 data.frame(group=rep("Family",nrow(fam_obs)),cluster=fam_count_obs$Group.1,fam_obs),
                 data.frame(group=rep("PG",nrow(pg_obs)),cluster=pg_count_obs$Group.1,pg_obs),
                 data.frame(group=rep("GF",nrow(gf_obs)),cluster=gf_count_obs$Group.1,gf_obs),
                 data.frame(group=rep("PFT",nrow(pft_obs)),cluster=pft_count_obs$Group.1,pft_obs))
tot_pred <- rbind(data.frame(group=rep("Species",nrow(spec_pred)),cluster=spec_count_pred$Group.1,spec_pred),
                 data.frame(group=rep("Genus",nrow(gen_pred)),cluster=gen_count_pred$Group.1,gen_pred),
                 data.frame(group=rep("Family",nrow(fam_pred)),cluster=fam_count_pred$Group.1,fam_pred),
                 data.frame(group=rep("PG",nrow(pg_pred)),cluster=pg_count_pred$Group.1,pg_pred),
                 data.frame(group=rep("GF",nrow(gf_pred)),cluster=gf_count_pred$Group.1,gf_pred),
                 data.frame(group=rep("PFT",nrow(pft_pred)),cluster=pft_count_pred$Group.1,pft_pred))
write.csv(tot_obs,file=file.path(origin,"_2021","data", "data_output","CV","tot_obs.csv"))
write.csv(tot_pred,file=file.path(origin,"_2021","data", "data_output","CV","tot_pred.csv"))
head(tot_obs)
                    
