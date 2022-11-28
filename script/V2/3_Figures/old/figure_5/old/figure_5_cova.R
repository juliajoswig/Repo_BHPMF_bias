

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
GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")

taxInfo$GF <- GF$GF
taxInfo$PFT <- PFT$PFT
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


dev.off()
pdf(file=file.path(origin,"figures","figure_5","figure_5_d.pdf"))
par(mfrow=c(2,2),mar=c(4,5,2,2))
if(trait_sub=="guido"){col_now=colz1}
if(trait_sub=="rainfor"){col_now=colz}
tr=1
for(tr in 2:ncol(spec_obs)){
  plot(cbind(spec_obs[,tr],spec_pred[,tr]-spec_obs[,tr]),
       main=paste0(colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",pch=16,
       xlim=c(-100,100),ylim=c(-100,100),xlab="Envelope cova obs",col=col_now[5])
  abline(0,-1)

  ix=spec_count_obs[,tr+1]==2
  points(cbind(spec_obs[ix,tr],spec_pred[ix,tr]-spec_obs[ix,tr]),col="red",pch=1)
}
  
  tr=1
  plot(cbind(gen_obs[,tr],gen_pred[,tr]-gen_obs[,tr]),
       main=paste0(colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
       xlim=c(-50,50),ylim=c(-50,50),xlab="Envelope cova obs",col=col_now[tr+1])
  abline(0,-1)
  for(tr in 2:ncol(gen_obs)){
    points(cbind(gen_obs[,tr],gen_pred[,tr]-gen_obs[,tr]),col=col_now[tr+1])
  }  
  
  plot(cbind(fam_obs[,tr],fam_pred[,tr]-fam_obs[,tr]),
       main=paste0(colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
       xlim=c(-50,50),ylim=c(-50,50),xlab="Envelope cova obs")
  abline(0,-1)
  
  plot(cbind(pg_obs[,tr],pg_pred[,tr]-pg_obs[,tr]),
       main=paste0(colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
       xlim=c(-50,50),ylim=c(-50,50),xlab="Envelope cova obs")
  abline(0,-1)

dev.off()



pdf(file=file.path(origin,"figures","figure_5","figure_5_e.pdf"),width = 12,height=4)
par(mfrow=c(1,6),mar=c(4,5,2,2))
col_env=c("#edf8e9","#bae4b3","#74c476","#238b45")
  
tr=1
for(tr in 1:ncol(spec_obs)){
  boxplot(cbind(spec_obs[,tr],spec_pred[,tr]-spec_obs[,tr]),
          main=paste0("Species envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-3,3),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")
  
  boxplot(cbind(gen_obs[,tr],gen_pred[,tr]-gen_obs[,tr]),
          main=paste0("Genera envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-3,3),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")
  
  boxplot(cbind(fam_obs[,tr],fam_pred[,tr]-fam_obs[,tr]),
          main=paste0("Families envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-10,10),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")
  
  boxplot(cbind(pg_obs[,tr],pg_pred[,tr]-pg_obs[,tr]),
          main=paste0("Phylo groups envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-10,10),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")

  
  boxplot(cbind(gf_obs[,tr],gf_pred[,tr]-gf_obs[,tr]),
          main=paste0("Growth forms envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-10,10),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")
  
  boxplot(cbind(pft_obs[,tr],pft_pred[,tr]-pft_obs[,tr]),
          main=paste0("PFTs envelope ",colnames(spec_obs)[tr]),ylab="Deviation Cova pred-obs",
          ylim=c(-10,10),xaxt="n",col=col_env)
  axis(1,at = 1:2,labels = c("obsered","predicted"),las=2)
  abline(h=0,col="red")
  
}
dev.off()



pdf(file=file.path(origin,"plots","figure_6","EnvelopeGuideCova.pdf"))
  par(mfrow=c(1,1),mar=c(4,5,2,2))
  for(tr in 1:ncol(cova_obs)){
    plot(cova_obs[,tr],cova_pred[,tr]-cova_obs[,tr],pch=16,cex.lab=2,cex.main=2,
         main=paste0(colnames(cova_obs)[tr]),xlab="Cova obs",ylab="Deviation Cova pred-obs",
         xlim=c(-40,40),ylim=c(-40,40))  
    abline(0,-1)
  }
dev.off()

pdf(file=file.path(origin,"plots","figure_6","EnvelopeGuideCova.pdf"))
dev.off()

cor_obs <- cor(trait_obs,use = "pairwise.complete.obs")
cor_pred <- cor(trait_pred,use = "pairwise.complete.obs")


i=1
j=2
pdf(file=file.path(origin,"plots","figure_6","EnvelopeGuideCorrel.pdf"),width=5,height=4)
par(mar=c(5,5,1,1))
  plot(1:10,col="white",ylim=c(-1,1),xlim=c(-1,1),cex.lab=1.2,
       xlab="Pearson observed",ylab="Pearson distance pred-obs")
  for(i in 1:ncol(cor_obs)){
    for(j in 1:ncol(cor_obs)){
      if(i!=j){
    points(x=cor_obs[i,j],y=c(cor_pred[i,j]-cor_obs[i,j]),pch=16)

  }}}
dev.off()


library(psych)
pdf(file=file.path(origin,"plots","figure_6","EnvelopeGuidePairs_obs.pdf"),width=10,height=10)
pairs.panels(trait_obs, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()
pdf(file=file.path(origin,"plots","figure_6","EnvelopeGuidePairs_pred.pdf"),width=10,height=10)
pairs.panels(trait_pred, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
dev.off()

pairs(cova_obs)
pairs(cova_pred)

tr=1
par(mfrow=c(1,2))
pie(sort(table(round(cova_obs[,tr],digits = 0)),decreasing = TRUE)/nrow(cova_pred)*100,main = "predicted")
pie(sort(table(round(cova_pred[,tr],digits = 0)),decreasing = TRUE)/nrow(cova_pred)*100,main = "observed")
boxplot(sort(table(round(cova_pred[,tr],digits = 0)),decreasing = TRUE)/nrow(cova_pred)*100,ylim=c(-2,2))
boxplot(sort(table(round(cova_obs[,tr],digits = 0)),decreasing = TRUE)/nrow(cova_pred)*100)

barplot(sort(table(round(cova_obs[,tr],digits = 0)),decreasing = TRUE))

plot(cova_obs$SLA,cova_pred$SLA)  

dist(cova_pred[,tr])
