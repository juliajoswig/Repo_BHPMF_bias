

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
TD_choice="Obs_obs"
trait_sub="guido"

colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858","black")
colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000","black")

GapPercent=0
list.files(file.path(origin,"runs",paste0("Rep",RepNum),
                     paste0("p",GapPercent,"_",trait_sub),TD_choice,"data"))

trait_obs <- read.csv(file.path(origin,"runs",paste0("Rep",RepNum),
                                paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                                "traitInfo.csv"))
taxInfo <- read.csv(file.path(origin,"runs",paste0("Rep",RepNum),
                              paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                              "taxInfo.csv"))
trait_pred <- read.csv(file.path(origin,"runs",paste0("Rep",RepNum),
                                 paste0("p",GapPercent,"_",trait_sub),TD_choice,"data",
                                 "mean.csv"),sep = "\t")

head(trait_pred)

dim(taxInfo)
dim(trait_obs)
head(trait_obs)

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

tot_obs <- rbind(spec_obs,gen_obs,fam_obs,pg_obs)
tot_pred <- rbind(spec_pred,gen_pred,fam_pred,pg_pred)

spec_count <- spec_count_obs[,2:19]==1&spec_count_pred[,2:19]!=1
gen_count <- gen_count_obs[,2:19]==1&gen_count_pred[,2:19]!=1
fam_count <- fam_count_obs[,2:19]==1&fam_count_pred[,2:19]!=1
pg_count <- pg_count_obs[,2:19]==1&pg_count_pred[,2:19]!=1

# match the species_cova evnelope and species_cova subset
# plot Figure: 
# dist(cova_envelope_pred - cova_subset_pred) := difference in cova deviation for envelope and subset(predicted)
# vs. 
# dist(cova_envelope(observed) - cova_subset(observed)):=difference in cova deviation for envelope and subset(predicted)
# vs. 

#Figure: 
# dist(cova_subset_pred70% - cova_subset_obs) := deviation in cova for subset(predicted)
# vs. 
# cova_envelope(observed)

GapPercent1="org"
RepNum=1
TD_choice="Obs_obs_TD"
trait_sub="guido"

colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")
#COVA
# load the coefficience of variance data per cluster GUIDO org
cova_now <- rep(NA,19)
repnums=1:2
GapPercent=80
RepNum=1
TD_choice="Obs_obs_TD"
for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
  for(GapPercent in gappercents){
    if(GapPercent!=-1){GapPercent1=GapPercent}
    if(GapPercent==-1){GapPercent1="org"}
    for(RepNum in repnums){
      path_now2=file.path(origin,"data_output","CoVa",
                          "guido",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
      if(file.exists(path_now2)){
        print(GapPercent1)
        cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
        colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
        print(head(cova_0_rawclust))
        print(dim(cova_0_rawclust))
        cova_now <- rbind(cova_now,
                          cbind(cova_0_rawclust,rep(GapPercent1,nrow(cova_0_rawclust)),
                                rep(TD_choice,nrow(cova_0_rawclust))))
      }
    }
  }
}
print(dim(cova_now))
head(covaG)
colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
cova_now <- as.data.frame(cova_now)
covaG <- cova_now

# load the coefficience of variance data per cluster RAINFOR org
cova_now <- rep(NA,21)
repnums=1:2
GapPercent=80
TD_choice="Obs_obs_TD"
for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
  for(GapPercent in gappercents){
    if(GapPercent!=-1){GapPercent1=GapPercent}
    if(GapPercent==-1){GapPercent1="org"}
    for(RepNum in repnums){
      path_now1 <- file.path(origin,"data_output","CoVa",
                             "rainfor",TD_choice,GapPercent1,RepNum,paste0("all.csv"))
      path_now2=file.path(origin,"data_output","CoVa",
                          "rainfor",TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv"))
      if(file.exists(path_now2)){
        print(GapPercent1)
        cova_0_rawclust <- as.matrix(read.csv(file=path_now2))
        print(head(cova_0_rawclust))
        print(dim(cova_0_rawclust))
        cova_now <- rbind(cova_now,
                          cbind(cova_0_rawclust,
                                rep(GapPercent1,nrow(cova_0_rawclust)),
                                rep(TD_choice,nrow(cova_0_rawclust))))
        print(dim(cova_now))
      }
    }
  }
}
print(dim(cova_now))

colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
cova_now <- as.data.frame(cova_now)

covaR <- cova_now

#-------------------------------------
new.mean_fun <- function(input){
  out=mean(as.numeric(input),na.rm = TRUE)
  return(out)
}

R_now=covaR
ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent","TD_choice"))]);  mode(ag_now) <- "numeric"
R <- aggregate(x=ag_now,
               by=list(Group=R_now$Group,Cluster=R_now$Cluster,TD_choice=R_now$TD_choice,Gaps=R_now$GapPercent),FUN=new_mean_fun)

R_now=covaG
ag_now=as.matrix(R_now[,-which(colnames(R_now)%in%c("Group","Cluster","GapPercent"))]);  mode(ag_now) <- "numeric"
G <- aggregate(x=ag_now,
               by=list(Group=R_now$Group,Cluster=R_now$Cluster,Gaps=R_now$GapPercent),FUN=new_mean_fun)

