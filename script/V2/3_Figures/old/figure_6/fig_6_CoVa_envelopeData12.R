

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

TD_choice_TD=paste0(TD_choice,"_TD")
trait_TDpred <- read.csv(file.path(origin,"runs",paste0("Rep",RepNum),
                                 paste0("p",GapPercent,"_",trait_sub),TD_choice_TD,"data",
                                 "mean.csv"),sep = "\t")

#RESCALE!!!
dim(taxInfo)
dim(trait_obs)
# chose indiv with tax
head(trait_obs)

count_fun <- function(x){
  return(sum(!is.na(x)))
}

#RESCALE THE STUFF:




#calc Cova Species
taxInfo_spec <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = mean,na.rm=TRUE)
taxInfo_spec <- as.data.frame(taxInfo_spec);names(taxInfo_spec)[1] <- "Cluster"
mean_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = sd,na.rm=TRUE)
spec_count_obs <- aggregate(trait_obs,by = list(taxInfo$AccSpeciesName),FUN = count_fun)
spec_count_pred <- aggregate(trait_pred,by = list(taxInfo$AccSpeciesName),FUN = count_fun)
taxInfo_spec <- aggregate(taxInfo,by = list(taxInfo$AccSpeciesName),FUN = unique)

spec_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
spec_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova Genus
taxInfo_gen <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = mean,na.rm=TRUE)
taxInfo_gen <- as.data.frame(taxInfo_gen);names(taxInfo_gen)[1] <- "Cluster"
mean_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = sd,na.rm=TRUE)
gen_count_obs <- aggregate(trait_obs,by = list(taxInfo$Genus),FUN = count_fun)
gen_count_pred <- aggregate(trait_pred,by = list(taxInfo$Genus),FUN = count_fun)

gen_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
gen_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova Family
taxInfo_fam <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = mean,na.rm=TRUE)
taxInfo_fam <- as.data.frame(taxInfo_fam);names(taxInfo_fam)[1] <- "Cluster"
mean_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = mean,na.rm=TRUE)
sd_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = sd,na.rm=TRUE)
mean_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = mean,na.rm=TRUE)
sd_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = sd,na.rm=TRUE)
fam_count_obs <- aggregate(trait_obs,by = list(taxInfo$Family),FUN = count_fun)
fam_count_pred <- aggregate(trait_pred,by = list(taxInfo$Family),FUN = count_fun)

fam_obs = sd_obs[,2:ncol(sd_obs)]/mean_obs[,2:ncol(sd_obs)]
fam_pred = sd_pred[,2:ncol(sd_pred)]/mean_pred[,2:ncol(sd_pred)]

#calc Cova phylogroup
taxInfo_pg <- aggregate(trait_obs,by = list(taxInfo$PhylogeneticGroup),FUN = mean,na.rm=TRUE)
taxInfo_pg <- as.data.frame(taxInfo_pg);names(taxInfo_pg)[1] <- "Cluster"
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


#COVA of subsets
# load the coefficience of variance data per cluster GUIDO org
{
  cova_now <- rep(NA,19)
  repnums=1:2
  GapPercent=80
  RepNum=1
  TD_choice="Obs_obs_TD"
#  for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
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
 # }
  print(dim(cova_now))
  head(cova_now)
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  covaG <- cova_now
  colnames(covaG)[3:7] <- paste0(colnames(covaR)[3:7],"_covaG")
  
  # load the coefficience of variance data per cluster RAINFOR org
  cova_now <- rep(NA,21)
  repnums=1:2
  GapPercent=80
  TD_choice="Obs_obs_TD"
  #for(TD_choice in c("Obs_obs_TD","Obs_obs","Spec_spec_TD","Spec_spec")){
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
          colnames(cova_0_rawclust) <- gsub(colnames(cova_0_rawclust),pattern = ".1",replacement = "_gaps")
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
 
  # }
  print(dim(cova_now))
  
  colnames(cova_now)[ncol(cova_now)] <-"TD_choice" 
  colnames(cova_now)[(ncol(cova_now)-1)] <-"GapPercent" 
  cova_now <- as.data.frame(cova_now)
  
  covaR <- cova_now
  colnames(covaR)[3:8] <- paste0(colnames(covaR)[3:8],"_covaR")
}


  
pdf(file=file.path(origin,"figures","figure_Sx","NumberEnvelCova.pdf"))
par(mfrow=c(1,1))
tax_envelope <- taxInfo_spec
spec_count2 <- data.frame(Cluster=taxInfo_spec$Cluster, spec_count_obs[2:ncol(spec_count_obs)],spec_obs)

spec <- merge(covaR,spec_count2,by="Cluster")
spec <- spec[spec$Group=="Species",]
barplot(cbind(
  rbind(length(unique(spec$Cluster[!is.na(spec$SLA_obs)])),length(unique(as.vector(spec$Cluster)))-length(unique(spec$Cluster[!is.na(spec$SLA_obs)]))),
  rbind(length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0])),
        (length(unique(spec$Cluster[spec$SLA>1]))-length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))),
  names.arg = c("Data set 2","Envelope"),ylab="Nb of Clusters",ylim=c(0,500),
  main = "Species",col=colz[c(4,5,6)]
)
text(.7,330,labels = paste0("reps. within cluster == 1:"))
text(.7,280,labels = paste0("n=",length(unique(spec$Cluster))))
text(.7,110,labels = paste0("reps. within cluster > 1:"))
text(.7,80,labels = paste0("n=",length(unique(spec$Cluster[!is.na(spec$SLA_obs)]))))
text(1.9,110,labels = paste0("reps. within cluster == 1:"))
text(1.9,80,labels = paste0("n=",length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))
text(1.9,300,labels = paste0("reps. within cluster >= 1:"))
text(1.9,280,labels = paste0("(additionally)"))
text(1.9,250,labels = paste0("n=",length(unique(spec$Cluster[spec$SLA>1]))-
                               length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))

tax_envelope <- taxInfo_spec
spec_count2 <- data.frame(Cluster=taxInfo_spec$Cluster, spec_count_obs[2:ncol(spec_count_obs)],spec_obs)

spec <- merge(covaG,spec_count2,by="Cluster")
spec <- spec[spec$Group=="Species",]
head(spec)
barplot(cbind(
  rbind(length(unique(spec$Cluster[!is.na(spec$SLA_obs)])),#>1
        length(unique(as.vector(spec$Cluster)))-length(unique(spec$Cluster[!is.na(spec$SLA_obs)]))),
  rbind(length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0])),
        (length(unique(spec$Cluster[spec$SLA>1]))-
           length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))),
  names.arg = c("Data set 1","Envelope"),ylab="Nb of Clusters",ylim=c(0,300),
  main = "Species",col=colz1[c(4,5,6)]
)
text(.7,170,labels = paste0("reps. within cluster == 1:"))
text(.7,150,labels = paste0("n=",length(unique(spec$Cluster))))
text(.7,80,labels = paste0("reps. within cluster > 1:"))
text(.7,50,labels = paste0("n=",length(unique(spec$Cluster[!is.na(spec$SLA_obs)]))))
text(1.9,80,labels = paste0("reps. within cluster == 1:"))
text(1.9,50,labels = paste0("n=",length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))
text(1.9,155,labels = paste0("reps. within cluster >= 1:"))
text(1.9,140,labels = paste0("(additionally)"))
text(1.9,120,labels = paste0("n=",length(unique(spec$Cluster[spec$SLA>1]))-
                                length(unique(spec$Cluster[!is.na(spec$SLA_obs)&spec$SLA>0]))))
dev.off()

