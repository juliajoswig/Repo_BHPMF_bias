
# get per observed point value: 
# - error (pred-obs)
# - average distance to all other points observed 
# - Silhouette index for this group (or calculate somewhere else...?)  
# - ORIGINAL Silhouette index for this group (or calculate somewhere else...?)  
# - DEVIATION Silhouette index for this group (or calculate somewhere else...?)  
# - number of values available
# - original number of values available

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
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
TD_choices = out$TD_choices
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

GapPercent=50
RepNum=4
trait_sub="rainfor"
TD_choice="Obs_obs_TD"

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
# completely observed
if(TD_choice=="Obs_obs_TD"&trait_sub=="rainfor"){traitInfo_obs=out$rainforTD_observed;taxInfo=out$rainforTD_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs_TD}
if(TD_choice=="Obs_obs_TD"&trait_sub=="guido"){traitInfo_obs=out$guidoTD_observed;taxInfo=out$guidoTD_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs_TD}
if(TD_choice=="Obs_obs"&trait_sub=="rainfor"){traitInfo_obs=out$rainfor_observed;taxInfo=out$rainfor_tax;scaling_factor <- out$scaling_factors_rainfor_Obs_obs}
if(TD_choice=="Obs_obs"&trait_sub=="guido"){traitInfo_obs=out$guido_observed;taxInfo=out$guido_tax;scaling_factor <- out$scaling_factors_guido_Obs_obs}
#sparse observed
#  traitInfo_sparse <- load(file = file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","TestData_org.RData"))
traitInfo_sparse <- read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","traitInfo.csv"), sep=",", dec=".")
# load the output for trait predictions: mean.csv
traitInfo_pred_zlog <- as.matrix(read.table(file.path(origin,"_2021","data","runs",paste0("Rep",RepNum), 
                                                      paste0("p",GapPercent,"_",trait_sub),TD_choice,"data/mean.csv"),
                                            sep="\t", dec=".",header=TRUE))


#-------------------------------------------------------------------
# cut to test data only if necessary
#-------------------------------------------------------------------
if(trait_sub=="guido"){ traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
if(trait_sub=="rainfor"){traitInfoTD_pred_zlog <- traitInfo_pred_zlog[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  
if(trait_sub=="guido"){ traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$guidoTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_guido] }  
if(trait_sub=="rainfor"){traitInfo_obs_c <- traitInfo_obs[as.numeric(taxInfo[,1])%in%as.numeric(out$rainforTD_tax[,1]),colnames(traitInfo_obs)%in%out$trait_rainfor]}  

#-------------------------------------------------------------------
# back transformation necessary.
#-------------------------------------------------------------------

#-------------------------------------------------------------------
# back transformation necessary.
#-------------------------------------------------------------------
summary(traitInfo_sparse)
summary(traitInfo_obs)
summary(traitInfo_obs_c)
summary(traitInfoTD_pred_zlog)

plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_c[,1])
abline(0,1)

traitInfo_obs_zlog <- traitInfo_obs
for(i in 1:ncol(traitInfo_obs)){
  traitInfo_obs_zlog[,i] <- (log(traitInfo_obs[,i]) - scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),1])/
    scaling_factor[which(rownames(scaling_factor)==colnames(traitInfo_obs)[i]),2]
}

plot(traitInfoTD_pred_zlog[,1],traitInfo_obs_zlog[,1])
abline(0,1)

traitInfo_pred <- traitInfoTD_pred_zlog
for(i in 1:ncol(traitInfoTD_pred_zlog)){
  traitInfo_pred[,i] <- exp((traitInfoTD_pred_zlog[,i]*scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),2])+ 
                              scaling_factor[which(rownames(scaling_factor)==colnames(traitInfoTD_pred_zlog)[i]),1])
}

plot(traitInfo_pred[,1],traitInfo_obs_c[,1])
abline(0,1)


#---------------------------------------------------------------------------------------------------------------

summary(traitInfoTD_pred_zlog)
summary(traitInfo_obs_zlog)
summary(traitInfo_sparse)
traitInfoTD_pred_zlog_sparse <- traitInfoTD_pred_zlog
traitInfoTD_pred_zlog_sparse[is.na(traitInfo_sparse)] <- NA
summary(traitInfoTD_pred_zlog_sparse)

res <- as.data.frame(matrix(NA,ncol=36,nrow=ncol(traitInfo_obs_zlog)*nrow(traitInfo_obs_zlog)))
colnames(res) <- c("error","Gap","value_obs","value_pred",
                   "dist_spec","dist_gen","dist_fam","dist_clad",
                   "nb_spec","nb_gen","nb_fam","nb_clad",
                   "mean_spec","mean_gen","mean_fam","mean_clad",#input mean
                   "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap",
                   "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred",
                   "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs",
                   "trait","ID","Species","Genus","Family","Clade","GF","PFT")
t=1
i=20
n=1
for(t in 1:ncol(traitInfo_obs_zlog)){
  for(i in 1: nrow(traitInfo_obs_zlog)){

    res$trait[n] <- colnames(traitInfo_obs_zlog)[t]
    res$ID[n] <- as.numeric(as.vector(taxInfo$IDs[i]))
    res$Species[n] <- as.vector(taxInfo$AccSpeciesName[i])
    res$Genus[n] <- as.vector(taxInfo$Genus[i])
    res$Family[n] <- as.vector(taxInfo$Family[i])
    res$Clade[n] <- as.vector(taxInfo$PhylogeneticGroup[i])
    res$GF[n] <- as.vector(taxInfo$GF[i])
    res$PFT[n] <- as.vector(taxInfo$PFT[i])
    res$error[n] <- traitInfoTD_pred_zlog[i,t]-traitInfo_obs_zlog[i,t]
    res$Gap[n] <- is.na(traitInfo_sparse[i,t])
    res$value_obs[n] <- traitInfo_obs_zlog[i,t]
    res$value_pred[n] <- traitInfoTD_pred_zlog[i,t]
    
    #----------------------------------------------
    #Individuals within species
    #----------------------------------------------
    ix1 = taxInfo$AccSpeciesName==as.vector(taxInfo$AccSpeciesName[i])
    res$mean_spec[n] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input species mean
      res$nb_spec_gap[n] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res$nb_spec[n] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_ind[n] <- 0
      }else{
        res$dist_ind[n] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res$dist_AVspec_pred[n] <- abs(traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
        res$dist_AVspec_obs[n] <- abs(traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
      }
      #----------------------------------------------
      #Species within genera
      #----------------------------------------------
      ix1 = taxInfo$Genus==as.vector(taxInfo$Genus[i])
      res$mean_gen[n] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_gen_gap[n] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res$nb_gen[n] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_spec[n] <- 0
      }else{
        res$dist_spec[n] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res$dist_AVgen_pred[n] <- abs(traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
        res$dist_AVgen_obs[n] <- abs(traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
      }
      #----------------------------------------------
      #Genera within families
      #----------------------------------------------
      ix1 = taxInfo$Family==as.vector(taxInfo$Family[i])
      res$mean_fam[n] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_fam_gap[n] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res$nb_fam[n] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_fam[n] <- 0
      }else{
        res$dist_fam[n] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res$dist_AVfam_pred[n] <- abs(traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
        res$dist_AVfam_obs[n] <- abs(traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
      }
      #----------------------------------------------
      #Families within Clades
      #----------------------------------------------
      ix1 = taxInfo$PhylogeneticGroup==as.vector(taxInfo$PhylogeneticGroup[i])
      res$mean_clad[n] <- mean(traitInfoTD_pred_zlog_sparse[ix1],na.rm = TRUE) #input taxon mean
      res$nb_clad_gap[n] = sum(!is.na(traitInfo_sparse[ix1,t]))
      res$nb_clad[n] = sum(ix1)
      ix <- ix1
      ix[taxInfo$ID==as.vector(taxInfo$ID[i])] <- FALSE
      if(sum(ix)==0){
        res$dist_clad[n] <- 0
      }else{
        res$dist_clad[n] <- mean(as.matrix(
          dist(c(traitInfo_obs_zlog[i,t],
                 traitInfo_obs_zlog[ix,t])))[2:(sum(ix)+1),1])
        res$dist_AVclad_pred[n] <- abs(traitInfo_pred_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
        res$dist_AVclad_obs[n] <- abs(traitInfo_obs_zlog[i,t]-mean(traitInfo_obs_zlog[ix1,t]))
      }
      n=n+1
    }
  }

#-----------------------------
ix=res$Gap
#-----------------------------

# I am a point value and my deviation is attracted by:

plot(res$value_obs,res$mean_spec,pch=16,xlim=c(-4,4),cex=.5)
plot(res$value_pred,res$mean_spec,pch=16,xlim=c(-4,4),cex=.5)
plot(res$value_obs,res$mean_gen,pch=16,xlim=c(-4,4),cex=.5)
plot(res$value_pred,res$mean_gen,pch=16,xlim=c(-4,4),cex=.5)
plot(res$dist_spec,res$mean_spec,pch=16,xlim=c(-4,4),cex=.5)

# how to best display "choice" of value where to head 
library("GGally")
ggparcoord(
  iris,
  columns = 1:4, groupColumn = 5, order = "anyClass",# the class input
  showPoints = TRUE, 
  title = "Where the points gets data from",
  alphaLines = 0.3
) + 
  theme_bw() +
  theme(legend.position = "top")

res_plot <- res[,colnames(res)%in%c("value_obs","value_pred","mean_spec","mean_gen","mean_fam","mean_clad",
                                    "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs")]
res_plot$value_dev <- res$value_pred-res$value_obs

res_plot <- res[,colnames(res)%in%c("value_obs","value_pred","mean_spec","mean_gen","mean_fam","mean_clad")]
#res_plot$Gap <- rep(NA,length(res_plot$mean_spec))
#res_plot$Gap[res$Gap] <- "Gap"
#res_plot$Gap[!res$Gap] <- "No Gap"
res_plot <- res_plot[complete.cases(res_plot),]
dim(res_plot)
head(res_plot)
# ------


# Library
#install.packages("fmsb")
library(fmsb)

# Create data: note in High school for Jonathan:
data <- as.data.frame(matrix( sample( 2:20 , 10 , replace=T) , ncol=10))
colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding", "data-viz" , "french" , "physic", "statistic", "sport" )

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data <- rbind(rep(20,10) , rep(0,10) , data)

# Check your data, it has to look like this!
# head(data)

# Custom the radarChart !
radarchart( data  , axistype=1 , 
            
            #custom polygon
            pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 , 
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            
            #custom labels
            vlcex=0.8 
)

# Custom the radarChart !
radarchart( res_plot[1:3,1:6]  , axistype=1 , 
            
            #custom polygon
            pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 , 
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            
            #custom labels
            vlcex=0.8 
)



require(FactoMineR)
#res_pca=res[,colnames(res)%in%c("error","dist_ind","dist_spec","dist_fam","dist_clad",
#                                                   "nb_spec","nb_gen","nb_fam","nb_clad",
#                                                   "nb_spec_gap","nb_gen_gap","nb_fam_gap","nb_clad_gap",
#                                                   "dist_AVspec_pred","dist_AVgen_pred","dist_AVfam_pred","dist_AVclad_pred",
#                                                   "dist_AVspec_obs","dist_AVgen_obs","dist_AVfam_obs","dist_AVclad_obs")]
res_pca <- data.frame(error=res$error,
                      dist_ind=res$dist_ind,dist_spec=res$dist_spec,dist_fam=res$dist_fam,dist_clad=res$dist_clad,
                      nb_spec=res$nb_spec,nb_gen=res$nb_gen,nb_fam=res$nb_fam,nb_clad=res$nb_clad,
                      nb_spec_gap=res$nb_spec_gap,nb_gen_gap=res$nb_gen_gap,nb_fam_gap=res$nb_fam_gap,nb_clad_gap=res$nb_clad_gap,
                      dist_AVind_dev=res$dist_AVspec_pred-res$dist_AVspec_obs,dist_AVspec_dev=res$dist_AVgen_pred,dist_AVfam_dev=res$dist_AVfam_pred,dist_AVclad_dev=res$dist_AVclad_pred,
                      dist_AVspec_pred=res$dist_AVspec_pred,dist_AVgen_pred=res$dist_AVgen_pred-res$dist_AVgen_obs,
                      dist_AVfam_pred=res$dist_AVfam_pred-res$dist_AVfam_obs,dist_AVclad_pred=res$dist_AVclad_pred-res$dist_AVclad_obs,
                      dist_AVspec_obs=res$dist_AVspec_obs,dist_AVgen_obs=res$dist_AVgen_obs,dist_AVfam_obs=res$dist_AVfam_obs,dist_AVclad_obs=res$dist_AVclad_obs)
#res_pca[is.na(res_pca)] <- 0
res_pca <- res_pca[complete.cases(res_pca[is.na(res_pca)])] <- 0
pdf(file=file.path(origin,"_2021","figures","pca.pdf"),height=15, width=15)
  PCA(res_pca)
dev.off()

#-----
# div av spec pred
# div av spec_dev
# dist avfam pred
plot(res_pca$dist_AVgen_pred,res_pca$dist_AVspec_dev)
plot(res_pca$dist_AVgen_pred,res_pca$dist_AVfam_pred)
dev.off()
res_pca <- data.frame(error=res$error,
                      dist_ind=res$dist_ind,dist_spec=res$dist_spec,dist_fam=res$dist_fam,dist_clad=res$dist_clad,
                      # nb_spec=res$nb_spec,nb_gen=res$nb_gen,nb_fam=res$nb_fam,nb_clad=res$nb_clad,
                      nb_spec_gap=res$nb_spec_gap,nb_gen_gap=res$nb_gen_gap,nb_fam_gap=res$nb_fam_gap,nb_clad_gap=res$nb_clad_gap,
                      dist_AVind_dev=res$dist_AVspec_pred-res$dist_AVspec_obs,
                      dist_AVspec_dev=res$dist_AVgen_pred-res$dist_AVgen_obs,
                      dist_AVfam_dev=res$dist_AVfam_pred-res$dist_AVfam_obs,
                      dist_AVclad_dev=res$dist_AVclad_pred-res$dist_AVclad_obs,
                      )

pdf(file=file.path(origin,"_2021","figures","pairs.pdf"),height=15, width=15)
  pairs(as.matrix(res_pca[,c(grep(colnames(res_pca),pattern = ),
                             )]))
dev.off()
# dist AVind obs
# dist_ind

# dist avspec obs
# dist spec

# dist AV fam obs
# dist fam

# dist_fam
#dist AVfam obs
# dist AVspec_dev

# distAVfam
# distAV fam dev

# we cannot calculate any deviaion for n=1
ix_ind=res$nb_spec_gap==2
dat_plot <- data.frame(deviation=res$dist_AVspec_pred[ix_ind] - res$dist_ind[ix_ind],error=abs(res$error)[ix_ind])
plot(dat_plot,col=colz[ix_ind],pch=16)
ix_ind=res$nb_spec_gap==3
dat_plot <- data.frame(deviation=res$dist_AVspec_pred[ix_ind] - res$dist_ind[ix_ind],error=abs(res$error)[ix_ind])
plot(dat_plot,col=colz[ix_ind],pch=16)


dat_plot <- data.frame(deviation=res$dist_AVspec_pred[ix] - res$dist_ind[ix],error=abs(res$error)[ix])
plot(dat_plot,col=colz[ix]);abline(0,-1)
lm_now <- lm(error~deviation,data = dat_plot)
abline(lm_now$coefficients,col="red")
cor(dat_plot[complete.cases(dat_plot),])
smr_spec <- summary(lm_now);print(smr_spec)
residuals_spec <- abs(abs(dat_plot[,1])-dat_plot[,2])

# chose only high
ix_resS <- residuals_spec>quantile(residuals_spec,probs = .5,na.rm = TRUE)
ix=res$Gap&ix_resS[res$Gap]
plot(dat_plot,col=colz[ix]);abline(0,-1)
dat_plot <- data.frame(deviation=res$dist_AVgen_pred[ix] - res$dist_spec[ix],error=abs(res$error)[ix])
plot(dat_plot,col=colz[ix]);abline(0,-1)
lm_now <- lm(error~deviation,data = dat_plot)
cor(dat_plot[complete.cases(dat_plot),])
smr_fam <- summary(lm_now);print(smr_fam)
residuals_fam <- smr_fam$residuals

dat_plot <- data.frame(deviation=res$dist_AVfam_pred[ix] - res$dist_fam[ix],error=abs(res$error)[ix])
plot(dat_plot,col=colz[ix]);abline(0,-1)
lm_now <- lm(error~deviation,data = dat_plot)
cor(dat_plot[complete.cases(dat_plot),])
smr_fam <- summary(lm_now);print(smr_fam)

dat_plot <- data.frame(deviation=res$dist_AVfam_pred[ix] - res$dist_fam[ix],error=abs(res$error)[ix])
plot(dat_plot,col=colz[ix]);abline(0,-1)
lm_now <- lm(error~deviation,data = dat_plot)
cor(dat_plot[complete.cases(dat_plot),])
smr_clad <- summary(lm_now);print(smr_clad)
smr_spec$adj.r.squared
smr_fam$adj.r.squared
smr_fam$adj.r.squared
smr_clad$adj.r.squared


plot(abs(res$dist_AVspec_pred[ix] - res$dist_ind[ix])-abs(res$error)[ix])
cor(abs(res$dist_AVspec_pred[ix] - res$dist_ind[ix])-abs(res$error)[ix])

boxplot(cbind(
  abs(res$dist_AVspec_pred[ix] - res$dist_ind[ix])-abs(res$error)[ix],
  abs(res$dist_AVgen_pred[ix] - res$dist_spec[ix])-abs(res$error)[ix],
  abs(res$dist_AVfam_pred[ix] - res$dist_fam[ix])-abs(res$error)[ix],
  abs(res$dist_AVfam_pred[ix] - res$dist_fam[ix])-abs(res$error)[ix]),ylim=c(-1,1),
  ylab="Distance of 1:1 line",xlab="",xaxt="n")
axis(side = 1,at = 1:4,labels = c("Species","Genera","Families","Clades"),las=2)
abline(h=0,col="gray")


Imp_groups <-  cbind(
  (res$dist_AVspec_pred[ix] - res$dist_ind[ix])-res$error[ix],
  (res$dist_AVgen_pred[ix] - res$dist_spec[ix])-res$error[ix],
  (res$dist_AVfam_pred[ix] - res$dist_fam[ix])-res$error[ix],
  (res$dist_AVclad_pred[ix] - res$dist_clad[ix])-res$error[ix])
boxplot(Imp_groups,ylim=c(-2,2),xaxt="n",ylab="Deviation - error")
axis(side = 1,at = 1:4,labels = c("Species","Genera","Families","Clades"),las=2)
abline(h=0,col="gray")





par(mfrow=c(1,1))
x <- res$nb_spec_gap[ix]
dat <- data.frame(x = x,y = x^(.5))
rbPal <- colorRampPalette(c('#ca0020',"#f4a582","#92c5de","#0571b0"))
colz_ieNB <- rbPal(20)[as.numeric(cut(dat$y,breaks = 20))]
plot(x,x,col=colz_ieNB,pch=16)

par(mfrow=c(1,1))
x <- res$dist_ind[ix]
dat <- data.frame(x = x,y = x^(.5))
rbPal <- colorRampPalette(c('red',"orange","blue","darkblue"))
rbPal <- colorRampPalette(c('#ca0020',"#f4a582","#92c5de","#0571b0")[4:1])
colz_ie <- rbPal(20)[as.numeric(cut(dat$y,breaks = 20))]
plot(x,x,col=colz_ie)
pchs=rep(16,nrow(res))
pchs[x>quantile(x,probs = .95)] <- 15 
pchs[x<quantile(x,probs = .50)] <- 25 
plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(species)",ylab = "pred-obs (absolute)")

  
pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_deviation.pdf"),width = 10,height = 14)
par(mfrow=c(2,2))
{
  colz=colz_ie
  plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(species)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("+ +","+","-","- -"),title = "Distance to individuals of same species",pch=16,col=c('#ca0020',"#f4a582","#92c5de","#0571b0"))
  plot(res$dist_AVgen_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("+ +","+","-","- -"),title = "Distance to individuals of same species",pch=16,col=c('#ca0020',"#f4a582","#92c5de","#0571b0"))
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("+ +","+","-","- -"),title = "Distance to individuals of same species",pch=16,col=c('#ca0020',"#f4a582","#92c5de","#0571b0"))
  plot(res$dist_AVclad_pred[ix] - res$dist_clad[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("+ +","+","-","- -"),title = "Distance to individuals of same species",pch=16,col=c('#ca0020',"#f4a582","#92c5de","#0571b0"))
}
dev.off()


par(mfrow=c(1,1))
xi <- round(res$nb_spec_gap[ix]/max(res$nb_spec_gap[ix]),digits = 5)
xs <- round(res$nb_gen_gap[ix]/max(res$nb_gen_gap[ix]),digits = 5)
xg <- round(res$nb_fam_gap[ix]/max(res$nb_fam_gap[ix]),digits = 5)

xi <- round(res$dist_ind/max(res$dist_ind),digits = 3)
xs <- round(res$dist_spec/max(res$dist_spec),digits = 3)
xg <- round(res$dist_fam/max(res$dist_fam),digits = 3)
xc <- round(res$dist_clad/max(res$dist_clad),digits = 3)

# kick out those values which are too dissimilar
xi[xi>quantile(xi,probs = .50)] <- .8
xs[xs>quantile(xs,probs = .50)] <- .8
xg[xg>quantile(xg,probs = .50)] <- .8
xc[xc>quantile(xc,probs = .50)] <- .8
xi[xi>quantile(xi,probs = .90)] <- .95
xs[xs>quantile(xs,probs = .90)] <- .95
xg[xg>quantile(xg,probs = .90)] <- .95
xc[xc>quantile(xc,probs = .90)] <- .95
summary(xi)

pchs=xc
pchs[xc<quantile(xc,probs = .50)&xc>quantile(xc,probs = .30)] <- 15
pchs[xc<quantile(xc,probs = .30)] <- 11
pchs[xc>=quantile(xc,probs = .50)] <- 16

plot(1:10,1:10,col=rgb(.8,.8,.8),pch=16)
rgb_dissimilarity <- cbind(xi,xs,xg)
colz_rgb_reverse=rgb(1-rgb_dissimilarity,alpha = 1)
colz=colz_rgb_reverse
plot(res$dist_AVspec_pred - res$dist_ind,res$error,pch=pchs,
     cex=.5,col=colz,xlab = "Deviation of average distance to ind(species)",ylab = "pred-obs (absolute)")
plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Deviation of average distance to ind(species)",ylab = "pred-obs (absolute)")
abline(h=0,v=0)
plot(res$nb_spec[ix]+rnorm(sum(ix),sd = .2),res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Number of co-values (same species)",ylab = "pred-obs (absolute)")

plot(res$dist_AVgen_pred[ix] - res$dist_spec[ix],res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Deviation of average distance to ind(genera)",ylab = "pred-obs (absolute)")
abline(h=0,v=0)
plot(res$nb_gen[ix]+rnorm(sum(ix),sd = .2),res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Number of co-values (same genus)",ylab = "pred-obs (absolute)")

plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Deviation of average distance to ind(families)",ylab = "pred-obs (absolute)")
abline(h=0,v=0)
plot(res$nb_fam[ix]+rnorm(sum(ix),sd = .2),res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Number of co-values (same family)",ylab = "pred-obs (absolute)")

plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],res$error[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Deviation of average distance to ind(families)",ylab = "pred-obs (absolute)")
plot(res$nb_clad[ix]+rnorm(sum(ix),sd = .2),abs(res$error)[ix],pch=pchs[ix],
     cex=.5,col=colz[ix],xlab = "Number of co-values (same clade)",ylab = "pred-obs (absolute)")

rgb_input <- cbind(xi,xs,xg)
colz_rgb=rgb(rgb_input,alpha = .5)
colz=colz_rgb
plot(1:length(xi),rnorm(length(xi)),col=colz_rgb,pch=16)
pchs=rep(16,nrow(res))


dev.off()
pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_dissimilarity.pdf"),width = 14,height = 14)
par(mfrow=c(2,2))
colz=colz_rgb_reverse[ix]
pchs=pchs[ix]
{
  plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  abline(0,1)
  abline(0,-1)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  
  plot(res$dist_AVgen_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-genera",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-families",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVclad_pred[ix] - res$dist_clad[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-clades",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
}
dev.off()



pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_dissimilarity.pdf"),width = 14,height = 14)
par(mfrow=c(2,2))
colz=colz_rgb_reverse[ix]
pchs=pchs[ix]
{
  i=150
  
  
  plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  abline(0,1)
  abline(0,-1)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  
  plot(res$dist_AVgen_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-genera",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-families",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVclad_pred[ix] - res$dist_clad[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-clades",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
}
dev.off()





pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_similarity.pdf"),width = 10,height = 14)
par(mfrow=c(2,2))
colz=colz_rgb
{
  plot(res$dist_AVspec_pred[ix] - res$dist_ind[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of average distance to ind(species)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "High similarity to indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  plot(res$dist_AVgen_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "High similarity to indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "High similarity to indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,cex=.7,col=colz,xlab = "Deviation of average distance to ind(gen)",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "High similarity to indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
}
dev.off()

pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_RGB.pdf"),width = 10,height = 10)
par(mfrow=c(4,4),mar=c(4,4,0,0))
{
  #------------------------
  #Ind within Species
  #------------------------
  plot(res$dist_ind[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_species (mean)",ylab="pred-obs (absolute)")
  legend(x = 0,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "High similarity to...",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  plot(res$nb_spec_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVspec_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVspec_obs[ix],res$dist_AVspec_pred[ix],pch=16,cex=1,col=colz,ylab="Predicted av distance",xlab="Observed av distance")
  #------------------------
  #Species within Genera
  #------------------------
  plot(res$dist_spec[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_genus (mean)",ylab="pred-obs (absolute)")
  plot(res$nb_gen_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVgen_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVgen_obs[ix],res$dist_AVgen_pred[ix],pch=16,cex=1,col=colz,ylab="Predicted av distance",xlab="Observed av distance")
  #------------------------
  #Genera within Families
  #------------------------
  plot(res$dist_fam[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_family (mean)",ylab="pred-obs (absolute)")
  plot(res$nb_fam_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVfam_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVfam_obs[ix],res$dist_AVfam_pred[ix],pch=16,cex=1,col=colz,ylab="Predicted av distance",xlab="Observed av distance")
  #------------------------
  #Familes within Clades
  #------------------------
  plot(res$dist_clad[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_clade (mean)",ylab="pred-obs (absolute)")
  plot(res$nb_clad_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  plot(res$dist_AVclad_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Pred distance to other individuals_clade (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVclad_obs[ix],res$dist_AVclad_pred[ix],pch=16,cex=1,col=colz,ylab="Predicted av distance",xlab="Observed av distance")
}
dev.off()

x <- res$nb_spec
dat <- data.frame(x = x,y = log(x))
rbPal <- colorRampPalette(c('red','blue'))
colz_i <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]

pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor.pdf"),width = 7,height = 10)
colz=colz_ie
#colz=colz_rgb
par(mfrow=c(2,2))
{  
  plot(res$dist_ind[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_species (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_spec[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_genus (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_fam[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_clad[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  
  par(mfrow=c(2,2))
  plot(res$nb_spec[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Nb of other individuals_species",ylab="pred-obs (absolute)")
  plot(res$nb_gen[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Nb of other individuals_genus",ylab="pred-obs (absolute)")
  plot(res$nb_fam[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Nb of other individuals_fam",ylab="pred-obs (absolute)")
  plot(res$nb_clad[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Nb of other individuals_fam",ylab="pred-obs (absolute)")
  
  par(mfrow=c(2,2))
  plot(res$dist_AVspec_pred[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_species (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVgen_pred[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_genus (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVfam_pred[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVclad_pred[ix],abs(res$error)[ix],pch=16,cex=.3,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  
  
  par(mfrow=c(2,2))
  plot(res$dist_AVspec_pred[ix]-res$dist_AVspec_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_species (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVgen_pred[ix]-res$dist_AVgen_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_genus (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVfam_pred[ix]-res$dist_AVfam_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVclad_pred[ix]-res$dist_AVclad_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_fam (mean)",ylab="pred-obs (absolute)")
  
  par(mfrow=c(2,2))
  plot(res$dist_AVspec_obs[ix] ,res$dist_AVspec_pred[ix],pch=16,cex=1,col=colz,xlab="Observed distance to other individuals_species (mean)",ylab="Pred distance to other individuals_spec (mean)")
  plot(res$dist_AVgen_obs[ix] ,res$dist_AVgen_pred[ix],pch=16,cex=1,col=colz,xlab="Observed distance to other individuals_genus (mean)",ylab="Pred distance to other individuals_gen (mean)")
  plot(res$dist_AVfam_obs[ix] ,res$dist_AVfam_pred[ix],pch=16,cex=1,col=colz,xlab="Observed distance to other individuals_fam (mean)",ylab="Pred distance to other individuals_fam (mean)")
  plot(res$dist_AVclad_obs[ix] ,res$dist_AVclad_pred[ix],pch=16,cex=1,col=colz,xlab="Observed distance to other individuals_clade (mean)",ylab="Pred distance to other individuals_clade (mean)")
}
  dev.off()
  
boxplot(cbind(res$dist_AVspec_obs[ix],res$dist_AVspec_pred[ix]))

colz=colz_rgb
plot(res$nb_spec_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Nb of individuals_species (mean)",ylab="pred-obs (absolute)")
par(mfrow=c(1,1))
plot(res$dist_AVspec_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")

pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor2.pdf"),width = 7,height = 10)
par(mfrow=c(4,3),mar=c(4,4,0,0))
{
  #------------------------
  #Ind within Species
  #------------------------
  plot(res$dist_ind[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_species (mean)",ylab="pred-obs (absolute)")
 # plot(res$nb_spec[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  plot(res$nb_spec_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Nb of individuals_species (mean)",ylab="pred-obs (absolute)")
  plot(res$dist_AVspec_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")

  x <- res$nb_gen
  dat <- data.frame(x = x,y = log(x))
  rbPal <- colorRampPalette(c('yellow','green'))
  colz_s <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
#  colz=colz_ie
  
  #------------------------
  #Species within Genera
  #------------------------
  plot(res$dist_spec[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_genus (mean)",ylab="pred-obs (absolute)")
  #plot(res$nb_gen[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  plot(res$nb_gen_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Nb of individuals_genera (mean)",ylab="pred-obs (absolute)")
  #plot(res$nb_gen,res$error,pch=16,cex=1,col=colz)
  #plot(res$nb_gen,res$dist_spec,pch=16,cex=1,col=colz)
  plot(res$dist_AVgen_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  #plot(res$dist_AVgen_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  
  #------------------------
  #Genera within Families
  #------------------------
  x <- res$nb_fam
  dat <- data.frame(x = x,y = log(x))
  rbPal <- colorRampPalette(c('yellow','green'))
  colz_g <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
#  colz=colz_ie
  plot(res$dist_fam[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_family (mean)",ylab="pred-obs (absolute)")
  #plot(res$nb_fam[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  plot(res$nb_fam_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Nb of individuals_family (mean)",ylab="pred-obs (absolute)")
  #plot(res$nb_fam[ix],res$dist_fam[ix],pch=16,cex=.5,col=colz)
  plot(res$dist_AVfam_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,ylab="pred-obs (absolute)")
  #plot(res$dist_AVfam_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  #Species within famera
  x <- res$nb_clad
  
  dat <- data.frame(x = x,y = log(x))
  rbPal <- colorRampPalette(c('yellow','green'))
  colz_f <- rbPal(10)[as.numeric(cut(dat$y,breaks = 10))]
#  colz=colz_ie
  plot(res$dist_clad[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Distance to other individuals_clade (mean)",ylab="pred-obs (absolute)")
  plot(res$nb_clad_gap[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Nb of individuals_clade (mean)",ylab="pred-obs (absolute)")
  #plot(res$nb_fam,res$dist_fam,pch=16,cex=.5,col=colz)
  plot(res$dist_AVclad_pred[ix],abs(res$error)[ix],pch=16,cex=1,col=colz,xlab="Pred distance to other individuals_clade (mean)",ylab="pred-obs (absolute)")
  #plot(res$dist_AVfam_obs[ix],abs(res$error)[ix],pch=16,cex=1,col=colz)
  }
dev.off()
# get per observed point value: 
# - error (pred-obs)
# - average distance to all other points observed 
# - Silhouette index for this group (or calculate somewhere else...?)  
# - ORIGINAL Silhouette index for this group (or calculate somewhere else...?)  
# - DEVIATION Silhouette index for this group (or calculate somewhere else...?)  
# - number of values available
# - original number of values available





















