# Julia Joswig 
# 20221003, last change 202212


#------------------------------------------------------------
# define path
#------------------------------------------------------------
setwd("/..")
origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git"
originData = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data"
list.files(file.path(origin,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","helper_scripts","fn_load_functions.R"))
load_functions(origin)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)
  t_choices <- out$t_choices
  TDnos = out$TDnos
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
  colz=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
         "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
         "#fbb4ae", "#b3cde3", "#ccebc5","#decbe4", "#fed9a6", "#ffffcc","#e5d8bd","#fddaec")
  
  
  
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
g=1
# start 20221003 ##############################################
RepNum = 1
t_choice="data"
t_choice="data_2"
res_matrix_name="res_20201020"#"res_20201112"
res_matrix_name= "res_20210303"
res <- read.table(file.path(originData,"analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
res <- as.data.frame(res)

correl.cols = grep(colnames(res),pattern = "cor_")#correl.cols#c(23,25:33)
Correls <- res[,correl.cols]
#colz = rainbow(ncol(Correls))
colnames(Correls)
Correls <- Correls[,colSums(!is.na(Correls))!=0]
# end 20221003 ##############################################

res_TDsparse <- list()
res_TD <- list()
res_TDtd <- list()
res_TDext <- list()

# load 
Percent <- gappercents[g]
{
  # load Envelope data
  # load TDenvelope
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),"Obs_obs_TD","data"))
  list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  dim(TD)
  dim(TD_sparse)
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
}

  
  td=1
  ts=1
#  for(td in 1:4){
#   for(ts in 1:2){
    
  #t_choice = "data_2"# see above

    
pdf(file=file.path(origin,"figures","Figure_3","Figure_Correl_review.pdf"),width=7,height=7)
{
  layout(mat = matrix(c(2, 1, 0, 3), 
                      nrow = 2, 
                      ncol = 2),
         heights = c(1, 2),    # Heights of the two rows
         widths = c(2, 1))     # Widths of the two columns
  layout(mat = matrix(c(2, 1,1,1, 
                        2, 1,1,1,
                        2, 1,1,1,
                        0, 3,3,3), 
                      nrow = 4, 
                      ncol = 4),
         heights = c(1, 2),    # Heights of the two rows
         widths = c(2, 1))     # Widths of the two columns
  
  
  TDno = "Obs_obs_TD"
  t_choice = "data"
  
  
  ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
  Correls2 <- as.matrix(Correls[ix_now,])
  Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
  
  Percent <- as.numeric(res$GapPercent)[ix_now]
  ltys=c(rep(1,9),rep(2,9))
  
  Correl_now <- cbind(Percent,abs(Correls2))
  colnames(Correl_now)
  whichcols=which(!duplicated(Correl_now[1,]))
  Correl_now <- Correl_now[,whichcols]
  #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
  colnames(Correl_now)
  
  Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
  Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
  
  Correl_80 <- Correl_now[Percent==80,2:ncol(Correl_now)]
  Correl_0 <- Correl_now[Percent==0,2:ncol(Correl_now)]
  CorrelTD_80 <- Correl_80
  CorrelTD_0 <- Correl_0
  #------------------------------------------------------------------
  TDno = "Obs_obs_TD"
  t_choice = "data_2"
  
  
  ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
  Correls2 <- as.matrix(Correls[ix_now,])
  Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
  
  Percent <- as.numeric(res$GapPercent)[ix_now]
  ltys=c(rep(1,9),rep(2,9))
  
  Correl_now <- cbind(Percent,abs(Correls2))
  colnames(Correl_now)
  whichcols=which(!duplicated(Correl_now[1,]))
  Correl_now <- Correl_now[,whichcols]
  #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
  colnames(Correl_now)
  
  Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
  Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
  
  CorrelTD2_80 <- Correl_now[Percent==80,2:ncol(Correl_now)]
  CorrelTD2_0 <- Correl_now[Percent==0,2:ncol(Correl_now)]
  
  par(mar = c(7, 7, 0, 0))
#  plot(1:10,type="n",frame=FALSE,xaxt="n", yaxt="n" ,ylab="" ,xlab="")
  
  plot(CorrelTD_0,CorrelTD_80,xlim=c(0,1),ylim=c(0,1),
       col=colz[8],pch=16,cex=2,cex.lab=3,
       ylab="Imputed Pearson corr.",
       xlab="Observed Pearson corr.")

  
#  par(mfrow=c(1,2))
#  boxplot(c(CorrelTD_0),ylim=c(0,1))
#  boxplot(c(CorrelTD_80),ylim=c(0,1))
  
  points(CorrelTD2_0,CorrelTD2_80,
       col=colz[10],pch=15,cex=2)

  abline(v=median(c(CorrelTD2_0)),col=colz[10])
  abline(v=median(c(CorrelTD_0)),col=colz[8])
  
  abline(h=median(c(CorrelTD2_80)),col=colz[10])
  abline(h=median(c(CorrelTD_80)),col=colz[8])
  
  
  legend(.1,.95,legend = c("OBS - IMPobs","OBS2 - IMP2obs"),
         pch = c(16,15),col = c(colz[8],colz[10]),
         cex=2, bty = "n")
  
  abline(0,1)
  
  ttest_out <- t.test(c(CorrelTD_0),c(CorrelTD_80))
  ttest2_out <- t.test(c(CorrelTD2_0),c(CorrelTD2_80))
  text(.69,.5, "p-val. OBS",cex=2)
  text(.935,.5, paste(round(ttest_out$p.value,digits = 3),"***"),cex=2)
  text(.67,.4, "p-val. OBS2",cex=2)
  text(.935,.4, paste(round(ttest2_out$p.value,digits = 3),"n.s."),cex=2)
  
  par(mar = c(0, 7, 0, 0))
  boxplot(list(OBS=c(CorrelTD_0),OBS2=c(CorrelTD2_0)),las=2, 
          horizontal=TRUE,
          ylim=c(0,1),col=colz[c(8,10)],xaxt="n",frame=FALSE,
          cex.axis=1.2)
  par(mar = c(7, 0, 0, 0))
  boxplot(list(IMPobsExt=c(CorrelTD_80),IMP2obsExt=c(CorrelTD2_80)),las=2, 
          ylim=c(0,1),col=colz[c(8,10)],yaxt="n",frame=FALSE,
          cex.axis=1.2)
  
#----------------------------------
#---------------------------------
  
  
  TDno = "Obs_obs"
  t_choice = "data"
  
  TDno = "Obs_obs"
  t_choice = "data"
  
  
  ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
  Correls2 <- as.matrix(Correls[ix_now,])
  Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
  
  Percent <- as.numeric(res$GapPercent)[ix_now]
  ltys=c(rep(1,9),rep(2,9))
  
  Correl_now <- cbind(Percent,abs(Correls2))
  colnames(Correl_now)
  whichcols=which(!duplicated(Correl_now[1,]))
  Correl_now <- Correl_now[,whichcols]
  #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
  colnames(Correl_now)
  
  Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
  Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
  
  Correl_80 <- Correl_now[Percent==80,2:ncol(Correl_now)]
  Correl_0 <- Correl_now[Percent==0,2:ncol(Correl_now)]
  CorrelTD_80 <- Correl_80
  CorrelTD_0 <- Correl_0
  
  #------------------------------------------------------------------
  TDno = "Obs_obs"
  t_choice = "data_2"
  
  
  ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
  Correls2 <- as.matrix(Correls[ix_now,])
  Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
  
  Percent <- as.numeric(res$GapPercent)[ix_now]
  ltys=c(rep(1,9),rep(2,9))
  
  Correl_now <- cbind(Percent,abs(Correls2))
  colnames(Correl_now)
  whichcols=which(!duplicated(Correl_now[1,]))
  Correl_now <- Correl_now[,whichcols]
  #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
  colnames(Correl_now)
  
  Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
  Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
  
  CorrelTD2_80 <- Correl_now[Percent==80,2:ncol(Correl_now)]
  CorrelTD2_0 <- Correl_now[Percent==0,2:ncol(Correl_now)]
  
  par(mar = c(7, 7, 0, 0))
  #  plot(1:10,type="n",frame=FALSE,xaxt="n", yaxt="n" ,ylab="" ,xlab="")
  
  plot(CorrelTD_0,CorrelTD_80,xlim=c(0,1),ylim=c(0,1),
       col=colz[8],pch=16,cex=2,cex.lab=3,
       ylab="Imputed Pearson corr.",
       xlab="Observed Pearson corr.")
  
  
  points(CorrelTD2_0,CorrelTD2_80,
         col=colz[10],pch=15,cex=2)
  abline(v=median(c(CorrelTD2_0)),col=colz[10])
  abline(v=median(c(CorrelTD_0)),col=colz[8])
  
  abline(h=median(c(CorrelTD2_80)),col=colz[10])
  abline(h=median(c(CorrelTD_80)),col=colz[8])
  
  
  legend(.1,.95,
         legend = c("OBS - IMPobsExt","OBS2 - IMP2obsExt"),pch = c(16,15),col = c(colz[8],colz[10]),
         cex=2,bty = "n")
  
  abline(0,1)
  ttest_out <- t.test(c(CorrelTD_0),c(CorrelTD_80))
  ttest2_out <- t.test(c(CorrelTD2_0),c(CorrelTD2_80))
  text(.69,.5, "p-val. OBS",cex=2)
  text(.935,.5, paste(round(ttest_out$p.value,digits = 3),"n.s."),cex=2)
  text(.67,.4, "p-val. OBS2",cex=2)
  text(.935,.4, paste(round(ttest2_out$p.value,digits = 3),"n.s."),cex=2)
  
  par(mar = c(0, 7, 0, 0))
  boxplot(list(OBS=c(CorrelTD_0),OBS2=c(CorrelTD2_0)),las=2, 
          horizontal=TRUE,
          ylim=c(0,1),col=colz[c(8,10)],xaxt="n",frame=FALSE,
          cex.axis=1.2)
  par(mar = c(7, 0, 0, 0))
  boxplot(list(IMPobsExt=c(CorrelTD_80),IMP2obsExt=c(CorrelTD2_80)),las=2, 
          ylim=c(0,1),col=colz[c(8,10)],yaxt="n",frame=FALSE,
          cex.axis=1.2)
  

  #---------------------------------
}
    
  dev.off()
  
  
  if(t_choice=="data"&TDno=="Obs_obs_TD"){pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_review_TDtd1.pdf"),width=6,height=8)}
  if(t_choice=="data_2"&TDno=="Obs_obs_TD"){pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_S_review_TDtd2.pdf"),width=6,height=8)}
  if(t_choice=="data"&TDno=="Obs_obs"){pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_review_TDex1.pdf"),width=6,height=8)}
  if(t_choice=="data_2"&TDno=="Obs_obs"){pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_S_review_TDex2.pdf"),width=6,height=8)}
  
    
    TDno = "Obs_obs_TD"
    TDno="Obs_obs"
    t_choice = "data_2"
    
    
    ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
    Correls2 <- as.matrix(Correls[ix_now,])
    Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
    
    Percent <- as.numeric(res$GapPercent)[ix_now]
    ltys=c(rep(1,9),rep(2,9))
    
    Correl_now <- cbind(Percent,abs(Correls2))
    colnames(Correl_now)
    whichcols=which(!duplicated(Correl_now[1,]))
    Correl_now <- Correl_now[,whichcols]
    #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
    colnames(Correl_now)
    
    Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
    Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
    
    Correl_80 <- Correl_now[Percent==80,2:ncol(Correl_now)]
    Correl_0 <- Correl_now[Percent==0,2:ncol(Correl_now)]
    plot(Correl_0,Correl_80,xlim=c(0,1),ylim=c(0,1),
         col="#a6cee3",pch=16,cex=2,
         ylab="Imputed Pearson corr.",
         xlab="Observed Pearson corr.")
    abline(0,1)
    
    dev.off()
    