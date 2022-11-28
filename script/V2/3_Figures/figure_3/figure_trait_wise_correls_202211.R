# Julia Joswig 
# 20221003

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
  # start 20221003 ##############################################
  origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/"
  # end 20221003 ##############################################
}


Version_now="V2"
list.files(file.path(origin,"_2021","script",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
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
res <- read.table(file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
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
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),"Obs_obs_TD","data"))
  list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
  # load TD data
  # total trait data 
  TD <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                         "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  TD_sparse <- as.data.frame(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
  dim(TD)
  dim(TD_sparse)
  # predicted 
  TDtd <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                       "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
  TDenv <- as.matrix(read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                        "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
}

  
  td=1
  ts=1
#  for(td in 1:4){
#   for(ts in 1:2){
    
  t_choice = "data_2"
  if(t_choice=="data"){
    trait_names=trait_rainfor
    clz=8}
    if(t_choice=="data_2"){ 
    trait_names=trait_guido
    clz=10}
  print(trait_names[tr])
  tr=4
  for(tr in 1:length(trait_names)){
pdf(file=file.path(origin,"_2021","figures","Figure_4","trait_wise",paste0("Figure_review_",t_choice,"_",trait_names[tr],".pdf")),
    width=7,height=7)
  {
  layout(mat = matrix(c(2, 1,1,1, 
                        2, 1,1,1,
                        2, 1,1,1,
                        0, 3,3,3), 
                      nrow = 4, 
                      ncol = 4),
         heights = c(1, 2),    # Heights of the two rows
         widths = c(2, 1))     # Widths of the two columns
  
  
  TDno = "Obs_obs_TD"

  
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

  par(mar = c(7, 7, 0, 0))
  ix=grep(pattern = trait_names[tr],x = colnames(CorrelTD_0))
  if(trait_names[tr]=="LeafN"){
    ix2 <- grep(x = colnames(CorrelTD_0),pattern = "LeafNArea")
    ix3 <- ix[!(ix%in%ix2)]
    ix4 <- ix2[grep(colnames(CorrelTD_0)[ix2],pattern = "LeafN_")]
    colnames(CorrelTD_0)[c(ix3,ix4)]
    ix=c(ix3,ix4)
    }
  plot(c(CorrelTD_0[,ix]),CorrelTD_80[,ix],xlim=c(0,1),ylim=c(0,1),
       col=colz[c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))],pch=16,cex=2,cex.lab=3,
       ylab="Imputed Pearson corr.",
       xlab="Observed Pearson corr.")

  legend(x=.4,y = 1,colnames(CorrelTD_0)[ix],
         col=colz[1:length(ix)],cex=1.5,lwd=4,lty=ltys,
         bty = "n")   
  
  abline(v=median(c(CorrelTD_0[,ix])),col=colz[clz])
  abline(h=median(c(CorrelTD_80[,ix])),col=colz[clz])
  
#if(t_choice=="data"){
#  legend(.1,.95,legend = c(paste0("TD-TDtd ",trait_names[tr])),
#         pch = c(16,15),col = c(colz[clz]),
#         cex=2, bty = "n")
#}else{
#  legend(.1,.95,legend = c(paste0("TD2-TD2td ",trait_names[tr])),
#         pch = c(16,15),col = c(colz[clz]),
#         cex=2, bty = "n")
#} 
  
  abline(0,1)
  
  ttest_out <- t.test(c(CorrelTD_0[,ix]),c(CorrelTD_80[,ix]))
  text(.71,.5, "p-val.",cex=2)
  if(ttest_out$p.value<=0.05){star="*"}
  if(ttest_out$p.value<0.01){star="**"}
  if(ttest_out$p.value<0.005){star="***"}
  if(ttest_out$p.value>0.1){star="n.s."}
  text(.935,.5, paste(round(ttest_out$p.value,digits = 3),star),cex=2)

  par(mar = c(0, 7, 0, 0))
  boxplot(data.frame(TD=c(CorrelTD_0[,ix])),las=2, 
          horizontal=TRUE,
          ylim=c(0,1),col=colz[clz],xaxt="n",frame=FALSE,
          cex.axis=2)
  par(mar = c(7, 0, 0, 0))
  boxplot(list(TDtd=c(CorrelTD_80[,ix])),las=2, 
          ylim=c(0,1),col=colz[clz],yaxt="n",frame=FALSE,
          cex.axis=2)
  
#----------------------------------
#---------------------------------
  TDno = "Obs_obs"
  
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
  par(mar = c(7, 7, 0, 0))
  plot(c(CorrelTD_0[,ix]),CorrelTD_80[,ix],xlim=c(0,1),ylim=c(0,1),
       col=colz[c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))],pch=16,cex=2,cex.lab=3,
       ylab="Imputed Pearson corr.",
       xlab="Observed Pearson corr.")
  
  legend(x=.4,y = 1,colnames(CorrelTD_0)[ix],
         col=colz[1:length(ix)],cex=1.5,lwd=4,lty=ltys,
         bty = "n")   
  
  
  abline(v=median(c(CorrelTD_0[,ix])),col=colz[clz])
  abline(h=median(c(CorrelTD_80[,ix])),col=colz[clz])
  
  
  abline(0,1)
  ttest_out <- t.test(c(CorrelTD_0[,ix]),c(CorrelTD_80[,ix]))
  text(.71,.5, "p-val. TD",cex=2)
  if(ttest_out$p.value<0.05){star="*"}
  if(ttest_out$p.value<0.005){star="**"}
  if(ttest_out$p.value<0.0005){star="***"}
  if(ttest_out$p.value>=0.05){star="n.s."}
  
  text(.935,.5, paste(round(ttest_out$p.value,digits = 3),star),cex=2)

  par(mar = c(0, 7, 0, 0))
  boxplot(list(TD=c(CorrelTD_0[,ix])),las=2, 
          horizontal=TRUE,
          ylim=c(0,1),col=colz[clz],xaxt="n",frame=FALSE,
          cex.axis=2)
  par(mar = c(7, 0, 0, 0))
  boxplot(list(TDext=c(CorrelTD_80[,ix])),las=2, 
          ylim=c(0,1),col=colz[clz],yaxt="n",frame=FALSE,
          cex.axis=2)
  
  #----------------------------------
  #---------------------------------
}
    
  dev.off()
  }
  