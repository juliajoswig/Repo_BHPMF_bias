
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
gappercents=c(1,5,10,20,30,40,50,60,70,80)
colz = rainbow(ncol(Correls))
colz=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
       "#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6",
       "#fbb4ae", "#b3cde3", "#ccebc5","#decbe4", "#fed9a6", "#ffffcc","#e5d8bd","#fddaec")


res_matrix_name="res_20201020"#"res_20201112"
res_matrix_name= "res_20210303"
res <- read.table(file.path(originData,"analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
res <- as.data.frame(res)

  correl.cols = grep(colnames(res),pattern = "cor_")#correl.cols#c(23,25:33)
  Correls <- res[,correl.cols]
  colnames(Correls)
  Correls <- Correls[,colSums(!is.na(Correls))!=0]

  
  list.files(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                       "Obs_obs",t_choice))
  
  EnvTot <- as.data.frame(read.csv(file.path(origin,"_2021",t_choice,"_runs",paste0("Rep_",1),t_choice,paste0("p_","0"),
                                             "Obs_obs",t_choice,"traitInfo_obs_REzlog.csv"),header=TRUE))[,-1]
  
  trait_names=trait_rainfor
  Correl_Env <- cor(EnvTot[,which(colnames(EnvTot)%in%trait_names)],use = "pairwise.complete.obs") 
  
  td=1
  ts=1
#  for(td in 1:4){
#   for(ts in 1:2){
    

    

    dev.off()
    pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_xyPearson.pdf"),width=15,height=5)
    
    layout(matrix(c(1,1,2,2,3),nrow = 1))
    par(mar=c(7,7,5,1))
    #    par(mfrow=c(1,2),mar=c(7,7,1,1))
    
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
    cexs=c(0,seq(.6,2.5,length.out = 14))
#    Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
#   Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
    i=2
    plot(Correl_now[1,-1],Correl_now[i,-1],pch=16,col=colz[1:ncol(Correl_agg)],cex=cexs[i])
    abline(0,1,lwd=2,col="gray")
    i=3
    for(i in 3:(ncol(Correl_now)-1)){
      points(Correl_now[1,-1],Correl_now[i,-1],col=colz[1:ncol(Correl_agg)],pch=16,cex=cexs[i])
    }
    
    
    
    #---------------------------------------------------------------
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
    cexs=c(0,seq(.6,2.5,length.out = 14))
    #    Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
    #   Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
    i=2
    plot(Correl_now[1,-1],Correl_now[i,-1],pch=16,col=colz[1:ncol(Correl_agg)],cex=cexs[i])
    abline(0,1,lwd=2,col="gray")
    i=3
    for(i in 3:(ncol(Correl_now)-1)){
      points(Correl_now[1,-1],Correl_now[i,-1],col=colz[1:ncol(Correl_agg)],pch=16,cex=cexs[i])
    }
    
    
    
    par(mar=c(0,0,0,0))
    plot(x=1:10,y=1:10,col="white",xaxt="n",yaxt="n",xlab="",ylab="",frame=FALSE)
    legend(x=1,y = 9,colnames(Correl_agg)[2:ncol(Correl_agg)],col=colz[1:ncol(Correl_agg)],cex=1.5,lwd=4,lty=ltys,
           bty = "n")    
    
    
    dev.off()
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_3","Figure_3_Corsd.pdf"),width=6,height=8)
    par(mfrow=c(1,1),mar=c(7,7,1,1))
    
    Correl_now <- cbind(Percent,Correls2)
    colnames(Correl_now)
    whichcols=which(!duplicated(Correl_now[1,]))
    Correl_now <- Correl_now[,whichcols]
    #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
    colnames(Correl_now)
    
    Correl_now <- aggregate(x=Correl_now,by=list(Percent),FUN=sd,na.rm=TRUE)
    Correl_now <- Correl_now[,c(1,3:ncol(Correl_now))]
    i=1
    plot(Correl_now[,c(1,(i+1))],col="white",cex.axis=2,
         xlab="Missingness [% gaps]",ylab="Pearson correlation coeff (sd)",xlim=c(0,120),ylim=c(0,1),cex.lab=3)
    lines(Correl_now[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])
    i=2
    for(i in 2:(ncol(Correl_now)-1)){
      points(Correl_now[,c(1,(i+1))],col="white",pch=16,ylim=c(-0,1.1))
      lines(Correl_now[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])      
    }
    
    abline(v=seq(1,80,by=20),col="gray")
    
    colnames(Correl_now) <- gsub(colnames(Correl_now),pattern = "cor_",replacement = "")
    legend(x=50,y = 1,colnames(Correl_now)[2:ncol(Correl_now)],col=colz[1:ncol(Correl_now)],cex=.8,lwd=3,lty=ltys,
           bty = "n",
           border = "white")    
    
    dev.off()
    
    
    
    # -----------------------------------------------------------------------------------------------------
    
    res_matrix_name="res_20201020"#"res_20201112"
    res_matrix_name= "res_20210303"
    res <- read.table(file.path(originData,"analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
    res <- as.data.frame(res)
    
    correl.cols = grep(colnames(res),pattern = "corzlog_")#correl.cols#c(23,25:33)
    Correls <- res[,correl.cols]
    colnames(Correls)
    Correls <- Correls[,colSums(!is.na(Correls))!=0]
    
    
    td=1
    ts=1
    #  for(td in 1:4){
    #   for(ts in 1:2){
    
    TDno = "Obs_obs_TD"
    t_choice = "data"
    
    ix_now=res$Obs_or_Spec==TDno & res$TraitChoice == t_choice
    Correls2 <- as.matrix(Correls[ix_now,])
    Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
    
    Percent <- as.numeric(res$GapPercent)[ix_now]
    ltys=c(rep(1,9),rep(2,9))
    
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_3","Figure_3_PearsonGaps_zlog.pdf"),width=6,height=8)
    par(mfrow=c(1,1),mar=c(7,7,1,1))
    
    Correl_now <- cbind(Percent,abs(Correls2))
    colnames(Correl_now)
    whichcols=which(!duplicated(Correl_now[1,]))
    Correl_now <- Correl_now[,whichcols]
    #Correl_now <- Correl_now[,-c(2,12,13,17,18,19,22,23,24,25,27,28,29,30,31)]
    colnames(Correl_now)
    
    Correl_agg <- aggregate(x=Correl_now,by=list(Percent),FUN=median,na.rm=TRUE)
    Correl_agg <- Correl_agg[,2:ncol(Correl_agg)]
    i=1
    plot(Correl_agg[,c(1,(i+1))],col="white",cex.axis=2,
         xlab="Missingness [% gaps]",ylab="Pearson correlation coeff",xlim=c(0,80),ylim=c(0,1.2),cex.lab=2.8)
    points(Correl_agg[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])
    i=2
    for(i in 2:(ncol(Correl_agg)-1)){
      points(Correl_agg[,c(1,(i+1))],col=colz[i],pch=16,ylim=c(-0,1.1))
      lines(Correl_agg[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])      
    }
    
    abline(v=seq(1,81,by=20),col="gray",lty=2)
    abline(h=seq(0,1,by=.20),col="gray",lty=2)
    rect(xleft = 80,ybottom = 0,xright = 125,ytop = 1,col="white",border = NA)
    colnames(Correl_agg) <- gsub(colnames(Correl_agg),pattern = "corzlog_",replacement = "")
    legend(x=35,y = 1.25,colnames(Correl_agg)[2:ncol(Correl_agg)],col=colz[1:ncol(Correl_agg)],cex=.8,lwd=2,lty=ltys,
           bty = "n",
           border = "white")    
    
    dev.off()
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_3","Figure_3_Corsd_zlog.pdf"),width=6,height=8)
    par(mfrow=c(1,1),mar=c(7,7,1,1))
    
    Correl_now <- cbind(Percent,Correls2)
    colnames(Correl_now)
    whichcols=which(!duplicated(Correl_now[1,]))
    Correl_now <- Correl_now[,whichcols]
    colnames(Correl_now)
    
    Correl_now <- aggregate(x=Correl_now,by=list(Percent),FUN=sd,na.rm=TRUE)
    Correl_now <- Correl_now[,c(1,3:ncol(Correl_now))]
    i=1
    plot(Correl_now[,c(1,(i+1))],col="white",cex.axis=2,
         xlab="Missingness [% gaps]",ylab="Pearson correlation coeff (sd)",xlim=c(0,120),ylim=c(0,1),cex.lab=3)
    lines(Correl_now[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])
    i=2
    for(i in 2:(ncol(Correl_now)-1)){
      points(Correl_now[,c(1,(i+1))],col="white",pch=16,ylim=c(-0,1.1))
      lines(Correl_now[,c(1,(i+1))],col=colz[i],lwd=6,lty=ltys[i])      
    }
    
    abline(v=seq(1,80,by=20),col="gray")
    
    colnames(Correl_now) <- gsub(colnames(Correl_now),pattern = "corzlog_",replacement = "")
    legend(x=50,y = 1,colnames(Correl_now)[2:ncol(Correl_now)],col=colz[1:ncol(Correl_now)],cex=1,lwd=3,lty=ltys,
           bty = "n",
           border = "white")    
    
    dev.off()
    