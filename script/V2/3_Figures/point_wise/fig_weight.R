
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

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")

res <- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  m=2
pdf(file=file.path(origin,"_2021","figures","barplot.pdf"))
par(mfrow=c(1,1),mar=c(12,4,2,2))
for(m in 1:length(missingness)){
for(t in 1:length(trait_names)){
  
  res_now=res[res$trait==trait_names[t],]
  res_now <- res_now[res_now$missingness==missingness[m],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
  dim(res_now)

# calc dev = pred-obs OR use "error"
# plot 

#plot(res_now$mean_spec,abs(res_now$value_pred-res_now$value_obs))

#dat_plot=cbind(abs(res_now$value_pred-res_now$mean_spec)-abs(res_now$value_obs-res_now$mean_spec),
#               abs(res_now$value_pred-res_now$mean_gen)-abs(res_now$value_obs-res_now$mean_gen),
#               abs(res_now$value_pred-res_now$mean_fam)-abs(res_now$value_obs-res_now$mean_fam),
#               abs(res_now$value_pred-res_now$mean_clad)-abs(res_now$value_obs-res_now$mean_clad))

#boxplot(dat_plot)
#abline(h=0)

dat_plot=cbind(res_now$error,res_now$value_obs,res_now$value_pred,res_now$mean_spec,res_now$mean_gen,
               res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT,
               res_now[,grep(colnames(res_now),pattern = "_lm")])
#require(FactoMineR)
#PCA(dat_plot[complete.cases(dat_plot),])
dat_plot=cbind(res_now$value_obs,res_now$value_pred,res_now$mean_spec,res_now$mean_gen,
               res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT,
               res_now[,grep(colnames(res_now),pattern = "_lm")])


cor_out <- cor(dat_plot,use = "pairwise.complete.obs")[,1:2]
rownames(cor_out) <- c("value_obs","value_pred","mean_spec","mean_gen",
                       "mean_fam","mean_clad","mean_GF","mean_PFT",colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
barplot(c(median(res_now$value_obs-res_now$value_pred,na.rm = TRUE),cor_out[3:nrow(cor_out),2]-cor_out[3:nrow(cor_out),1]),las=2,ylim=c(-.3,.3),
        main=paste0(trait_names[t],"_",missingness[m]),col=c("gray",colz1),ylab="cor_pred - cor_obs")
abline(h=seq(-.5,to=.5,by=.1),col="gray",lty=2)

}
}
dev.off()


pairs(dat_plot)


barplot(c(median(res_now$value_pred,na.rm = TRUE),median(res_now$value_obs,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_spec,na.rm = TRUE),median(res_now$value_obs-res_now$mean_spec,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_gen,na.rm = TRUE),median(res_now$value_obs-res_now$mean_gen,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_fam,na.rm = TRUE),median(res_now$value_obs-res_now$mean_fam,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_clad,na.rm = TRUE),median(res_now$value_obs-res_now$mean_clad,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_GF,na.rm = TRUE),median(res_now$value_obs-res_now$mean_GF,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_PFT,na.rm = TRUE),median(res_now$value_obs-res_now$mean_PFT,na.rm = TRUE)))

boxplot(cbind((res_now$value_pred-res_now$mean_spec),(res_now$value_obs-res_now$mean_spec)))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_gen),(res_now$value_obs-res_now$mean_gen))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_fam),(res_now$value_obs-res_now$mean_fam))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_clad),(res_now$value_obs-res_now$mean_clad))
abline(h=0)

plot(res_now$mean_gen,abs(res_now$error))
plot(abs(res_now$error),res_now$mean_gen)
plot(res_now$mean_gen,abs(res_now$error))
plot(res_now$mean_fam,res_now$error)
plot(res_now$mean_clad,res_now$error)
plot(res_now$mean_GF,res_now$error)
plot(res_now$mean_PFT,res_now$error)
# all lms
# plot(res_now$error,res_now$lm)

dat_cor <- cbind(res_now$error,res_now$mean_spec,res_now$mean_gen,res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT)
dat_cor <- cbind(res_now$error,res_now$mean_spec,res_now$mean_gen,res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT)
pairs(dat_cor)
cor(dat_cor)


