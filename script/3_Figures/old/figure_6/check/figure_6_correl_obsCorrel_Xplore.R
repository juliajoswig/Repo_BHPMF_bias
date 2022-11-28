## figure correlation deviation against observed correlation


  colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
  colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")
  
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
  trait_guido =out$trait_guido
  trait_rainfor =out$trait_rainfor
  colz1 =out$colz1
  colz2 =out$colz2
  Indeces <- c("rmse_man")

  
  res_matrix_name="res_20201126"
  res <- read.table(file.path(origin,"runs","META",paste0(res_matrix_name,".csv")),sep=",",dec=".")
  res <- as.data.frame(res)
  
  ts=1
  tsub=tsubs[ts]
  whichcols_TD <- res$Obs_or_Spec=="Obs_obs_TD"&res$TraitChoice==tsub
  whichcols <- res$Obs_or_Spec=="Obs_obs"&res$TraitChoice==tsub
  res_sub <- res[whichcols,c(1:4,grep(colnames(res),pattern = "cor_"))]
  res_TD <- res[whichcols_TD,c(1:4,grep(colnames(res),pattern = "cor_"))]
  colnames(res_sub)  

  res_sub <- as.matrix(res_sub)
  mode(res_sub) <- "numeric"
  res_sub_obs <- colMeans(res_sub[res_sub[,4]==-1,grep(colnames(res_sub),pattern = "cor_")],na.rm = TRUE)
  res_dist <- res_sub
  i=5
  for(i in 5:ncol(res_dist)){
    res_dist[,i] <- res_sub[,i] - res_sub_obs[(i-4)]
  }  
  
  res_dist <- as.data.frame(res_dist)
  res_dist <- res_dist[,colSums(!is.na(res_dist))!=0]
  res_sub <- as.matrix(res_sub)
  mode(res_sub) <- "numeric"
  i=5
  par(mfrow=c(2,4))
  plot(res_dist[,c(i+5)],res_dist[,c(i)],
       ylim = c(-.2,.2),xlab="Observed change",ylab="pred change")
  plot(res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i+5])],
       res_dist[,c(i)],
       ylim = c(-.2,.2),xlab="Observed cor",ylab="pred change")
  plot(res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i])],
       res_dist[,c(i+5)],
       ylim = c(-.2,.2),ylab="Observed cor",xlab="pred change")
  plot(res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i+5])],
       res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i])],
       xlab="Observed cor",ylab="Pred cor")
  boxplot(res_dist[,c(i)],ylim=c(-.6,.5),main="Pred distance to -1")  
  boxplot(res_dist[,c(i+5)],ylim=c(-.6,.5),main="Obs distance to -1")  
  boxplot(res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i])],ylim=c(-.6,.5),main="Pred")  
  boxplot(res_sub[,which(colnames(res_sub)%in%colnames(res_dist)[i+5])],ylim=c(-.6,.5),main="Obs")  
  
  plot(res_dist$GapPercent,res_dist$cor_SLA_PlantHeight)
  plot(res_dist$GapPercent,res_dist$cor_SLA_PlantHeight_gappy)
  
  plot(res_dist$cor_SLA_PlantHeight,res_dist$cor_SLA_PlantHeight_gappy,
       xlim = c(-1,1),ylim = c(-1,1))
  boxplot(res_dist$cor_SLA_PlantHeight,ylim=c(-.2,.2))  
  boxplot(res_dist$cor_SLA_PlantHeight_gappy,ylim=c(-.2,.2))  
  
  