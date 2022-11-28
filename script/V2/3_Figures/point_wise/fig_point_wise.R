
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



#-----------------------------
ix=res$Gap
#-----------------------------

xs <- round(res$nb_spec_gap[ix]/max(res$nb_spec_gap[ix]),digits = 5)
xg <- round(res$nb_gen_gap[ix]/max(res$nb_gen_gap[ix]),digits = 5)
xf <- round(res$nb_fam_gap[ix]/max(res$nb_fam_gap[ix]),digits = 5)

xs <- round(res$dist_spec/max(res$dist_spec),digits = 3)
xg <- round(res$dist_gen/max(res$dist_gen),digits = 3)
xf <- round(res$dist_fam/max(res$dist_fam),digits = 3)
xc <- round(res$dist_fam/max(res$dist_fam),digits = 3)

# kick out those values which are too dissimilar
xs[xs>quantile(xs,probs = .50)] <- .8
xg[xg>quantile(xg,probs = .50)] <- .8
xf[xf>quantile(xf,probs = .50)] <- .8
xc[xc>quantile(xc,probs = .50)] <- .8

xs[xs>quantile(xs,probs = .90)] <- .95
xg[xg>quantile(xg,probs = .90)] <- .95
xc[xf>quantile(xf,probs = .90)] <- .95
xc[xc>quantile(xc,probs = .90)] <- .95
summary(xs)

pchs=xc
pchs[xc<quantile(xc,probs = .50)&xc>quantile(xc,probs = .30)] <- 15
pchs[xc<quantile(xc,probs = .30)] <- 11
pchs[xc>=quantile(xc,probs = .50)] <- 16

plot(1:10,1:10,col=rgb(.8,.8,.8),pch=16)
rgb_dissimilarity <- cbind(xs,xg,xf)
colz_rgb_reverse=rgb(1-rgb_dissimilarity,alpha = 1)
colz=colz_rgb_reverse


par(mfrow=c(1,1))

rgb_input <- cbind(xs,xg,xf)
colz_rgb=rgb(rgb_input,alpha = .5)
plot(1:length(xs),rnorm(length(xs)),col=colz_rgb,pch=16)
pchs=rep(16,nrow(res))



pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_dissimilarity.pdf"),width = 14,height = 14)
par(mfrow=c(2,2))
colz=colz_rgb_reverse[ix]
pchs=pchs[ix]
{
  plot(res$dist_AVspec_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,
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
  
  plot(res$dist_AVgen_pred[ix] - res$dist_gen[ix],abs(res$error)[ix],pch=pchs,
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

#--------------------------------------------------------------------------------
# PCA
#--------------------------------------------------------------------------------

res_pca <- data.frame(error=res$error,
                      dist_spec=res$dist_spec,dist_gen=res$dist_gen,dist_fam=res$dist_fam,dist_clad=res$dist_clad,
                      # nb_ind=res$nb_ind,nb_spec=res$nb_spec,nb_gen=res$nb_gen,nb_fam=res$nb_fam,
                      nb_spec_gap=res$nb_spec_gap,nb_gen_gap=res$nb_gen_gap,nb_fam_gap=res$nb_fam_gap,nb_clad_gap=res$nb_clad_gap,
                      dist_AVspec_dev=res$dist_AVgen_pred-res$dist_AVspec_obs,
                      dist_AVgen_dev=res$dist_AVfam_pred-res$dist_AVgen_obs,
                      dist_AVfam_dev=res$dist_AVclad_pred-res$dist_AVfam_obs,
                      dist_AVclad_dev=res$dist_AVclad_pred-res$dist_AVclad_obs
)

pdf(file=file.path(origin,"_2021","figures","pairs.pdf"),height=15, width=15)
pairs(as.matrix(res_pca[,c(grep(colnames(res_pca),pattern = ),
)]))
dev.off()

