

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

#------------------------------------------------------------
res_matrix_name="res_20201126"#"res_20201112"
res <- read.table(file.path(origin,"runs","META",paste0(res_matrix_name,".csv")),sep=",",dec=".")

for(tc in 1:4){
  
  TD_choice <- TD_choices[tc]
  pdf(file=file.path(origin,"plots","figure_6",paste0("Deviation_correlation",
                                                      TD_choice <- TD_choices[tc],".pdf")),width=10,height=4)
  
  par(mfrow=c(1,2),mar=c(5,5,1,1))
  
for(td in 1:2){
  trait_choice <- tsubs[td]
  if(trait_choice=="rainfor"){col_now=colz}
  if(trait_choice=="guido"){col_now=colz1}
  
  res_sub <- res[res$TraitChoice==trait_choice&res$Obs_or_Spec==TD_choice,c(4,grep(colnames(res),pattern = "cor_"))]
  res_sub <- res_sub[,-grep(colnames(res_sub),pattern = "gappy")]
  res_sub <- res_sub[,colSums(!is.na(res_sub))!=0]
  
  res_sub <- as.matrix(res_sub)
  mode(res_sub) <- "numeric"
  res_obs <- colMeans(res_sub[res_sub[,1]==-1,],na.rm=TRUE)
  res_obs2 <- res_sub
  
  for(i in 2:ncol(res_obs2)){
    res_obs2[,i] <- rep(res_obs[i],nrow(res_sub))
  }
  res_dist <- res_sub-res_obs2

  plot(1:10,col="white",ylab="Deviation pred-observed",xlab="observed correl",xlim=c(-1,1),ylim=c(-1,1))

  i=2
  for(i in 2:ncol(res_dist)){
    print(i)
    dat_now <- cbind(res_sub[,1],res_obs2[,i],res_dist[,i])
    for(j in 1:10){
      now=c(0,1,5,10,20,30,40,50,60,70)[j]
      cexs=c(3,2.5,2,1.5,1,.8,.6,.5,.4,.3)
      ix=dat_now[,1]==now
      ix[is.na(ix)] <- FALSE
      points(dat_now[ix,2:3],col=col_now[j],pch=16,cex=cexs[j])
    }
  }
  

}
  dev.off()
}


#aussortiert:   
plot(1:10,col="white",ylab="predicted correl",xlab="observed correl",xlim=c(-1,1),ylim=c(-1,1))
abline(0,1)
i=2
for(i in 2:ncol(res_dist)){
  print(i)
  dat_now <- cbind(res_sub[,1],res_obs2[,i],res_sub[,i])
  for(j in 1:10){
    now=c(0,1,5,10,20,30,40,50,60,70)[j]
    cexs=c(3,2.5,2,1.5,1,.8,.6,.5,.4,.3)
    ix=dat_now[,1]==now
    ix[is.na(ix)] <- FALSE
    points(dat_now[ix,2:3],col=col_now[j],pch=16,cex=cexs[j])
  }
}