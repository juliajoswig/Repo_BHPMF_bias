
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
t_choices <- out$tsubs
TDnos = out$TD_choices
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
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)

res_matrix_name="res_20201020"
#res_o=res
res_matrix_name="res_20210303"
res <- read.csv(file=file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")))
summary(res)
list.files(file.path(origin,"runs","META"))
#colz=rainbow(7)#c("red","blue","green","magenta","turquoise","orange","yellow")#
#colz=c("#d73027","#f46d43","#fdae61","#fee08b","#d9ef8b","#a6d96a","#66bd63","#1a9850")
#colz=c("#d7191c","#2c7bb6","#fdae61","#abd9e9")#"#ffffbf",
#colz1=c("#f0f9e8","#bae4bc","#7bccc4","#43a2ca","#0868ac")#"#ffffbf",
#colz2=c("#fef0d9","#fdcc8a","#fc8d59","#e34a33","#b30000")#"#ffffbf",
output_term="2020"

# RMSE increases with gap-size
# RMSE different for traits
# RMSE different for data set

#--------------------------------------
# RMSE increases with gap-size
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="RMSE"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="data_2"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse)
head(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
dat_plot <- dat_plot[,colSums(!is.na(dat_plot))!=0]
dat_plot <- dat_plot[complete.cases(dat_plot),]
#chose envelope guido:

#chose envelope rainfor
ix2=res$TraitChoice=="data_2"&res$Obs_or_Spec=="Obs_obs"
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent2) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

dat_plotENV <- aggregate(dat_rmse2,by=list(Percent2),FUN=mean,na.rm = TRUE)
dat_plotENV <- dat_plotENV[,colSums(!is.na(dat_plotENV))!=0]
dat_plotENV <- dat_plotENV[complete.cases(dat_plotENV),]
dat_plotENV <- dat_plotENV[dat_plotENV[,1]!=0,]

plot(dat_plotENV[,2],dat_plot2[,2])
cbind(dat_plotENV[,2],dat_plot[,2],dat_plot2[,2])

#chose data (rainfor)
ix3=res[,colnames(res)=="TraitChoice"]=="data"&res[,colnames(res)=="Obs_or_Spec"]=="Obs_obs_TD"
Percent2 <- res[ix3,colnames(res)=="GapPercent"]
mode(Percent3) <- "numeric"

dat_rmse3 <- res[ix3,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse3)

dat_rmse3 <- as.matrix(dat_rmse3)
mode(dat_rmse3)="numeric"
colnames(dat_rmse3)

dat_plot2 <- aggregate(dat_rmse3,by=list(Percent3),FUN=mean,na.rm = TRUE)
dat_plot2 <- dat_plot2[,colSums(!is.na(dat_plot2))!=0]
dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]
dat_plot2 <- dat_plot2[dat_plot2[,1]!=0,]

cbind(dat_plot[,2],dat_plot2[,2])

#chose envelope rainfor
ix4=res$TraitChoice=="data"&res$Obs_or_Spec=="Obs_obs"
Percent4 <- res[ix4,colnames(res)=="GapPercent"]
mode(Percent4) <- "numeric"

dat_rmse4 <- res[ix4,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse4)

dat_rmse4 <- as.matrix(dat_rmse3)
mode(dat_rmse4)="numeric"
colnames(dat_rmse4)

dat_plot2ENV <- aggregate(dat_rmse4,by=list(Percent4),FUN=mean,na.rm = TRUE)
dat_plot2ENV <- dat_plot2ENV[,colSums(!is.na(dat_plot2ENV))!=0]
dat_plot2ENV <- dat_plot2ENV[complete.cases(dat_plot2ENV),]
dat_plot2ENV <- dat_plot2ENV[dat_plot2ENV[,1]!=0,]

plot(dat_plot2ENV[,2],dat_plot2[,2])
cbind(dat_plotENV[,2],dat_plot2ENV[,2],dat_plot[,2],dat_plot2[,2])



trait_names=trait_rainfor
#--------------------------------------
# data 1
#--------------------------------------
colz=c("#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#b2182b")

dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"]
dat_plot2ENV[,colnames(dat_plot2)=="RMSE_zlog_Total"]
dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"]
dat_plot2ENV[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"]

  pdf(file=file.path(origin,"_2021","Figures","Figure_1","Fig_1_RMSE.pdf"),width=15,height=10)
  {
  layout(matrix(c(1,1,2,2,3,
                        4,4,5,5,6),nrow=2,ncol=5,byrow=TRUE))
  #layout.show(nf)
  #par(mar=c(7,7,2,2),mfrow=c(1,1))
  par(mar=c(8,10,2,2))
  i=2
      
      #------------------------------TDtd Total= all values
      plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
           xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
          xaxt="n",yaxt="n", main="Total values",cex.main=2.5)
      axis(1,line = 2.5,cex.axis=2)
      axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
      axis(2,line = 2.5,cex.axis=2)
      axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
      abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
      colnames(dat_plot2)
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
      
      i=1
      for(i in 1:length(trait_names)){
        print(trait_names[i])
        lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
      }  
      
      i=2
      # TDtd
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="white",lwd=5.5)
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
      
      
      #------------------------------TDenv Total = all values
      plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
           xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
           xaxt="n",yaxt="n", main="Total values (imputed with envelope)",cex.main=2.5)
      axis(1,line = 2.5,cex.axis=2)
      axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
      axis(2,line = 2.5,cex.axis=2)
      axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
      abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_zlog_Total"],col="black",lwd=5)
      
      i=1
      for(i in 1:length(trait_names)){
        print(trait_names[i])
        lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
      }  
      
      i=2
      # TDenvelope
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_zlog_Total"],col="white",lwd=5.5)
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_zlog_Total"],col="black",lwd=5)
      par(mar=c(0,0,0,0))     
      plot(x = 1:10,y = 1:10,xaxt="n",yaxt="n",xlab="",ylab="",col="white",frame = FALSE)
      level_names = c(trait_names,"Mean")
      legend(1, 10, level_names, col = c(colz[1:length(trait_names)],"black"),bty = "n",
             text.col = "black", lty = c(rep(1,length(trait_names)),1), lwd=4,cex=3.15,
             merge = TRUE, bg = "white")  
      
      #-----------------------------------------------------------------------------------
      #------------------------------TDtd gap= only missing values
      par(mar=c(8,10,2,2))
      plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
           xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
           xaxt="n",yaxt="n", main="Gaps",cex.main=2.5)
      axis(1,line = 2.5,cex.axis=2)
      axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
      axis(2,line = 2.5,cex.axis=2)
      axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
      abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
      colnames(dat_plot2)
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
      
      i=1
      for(i in 1:length(trait_names)){
        print(trait_names[i])
        lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],lwd=5)
      }  
      
      i=2
      # TDtd
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"],col="white",lwd=5.5)
      lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
      
      
      #------------------------------TDenv  gap= only missing values
      plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
           xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
           xaxt="n",yaxt="n", main="Gaps (imputed with envelope)",cex.main=2.5)
      axis(1,line = 2.5,cex.axis=2)
      axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
      axis(2,line = 2.5,cex.axis=2)
      axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
      abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
      colnames(dat_plot2ENV)
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
      
      i=1
      for(i in 1:length(trait_names)){
        print(trait_names[i])
        lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],lwd=5)
      }  
      
      i=2
      # TDenvelope
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_gap_zlog_Total"],col="white",lwd=5.5)
      lines(dat_plot2ENV[,1],dat_plot2ENV[,colnames(dat_plot2ENV)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
  }
  dev.off()
  

  
  
  
  
  trait_names = trait_guido
 
  pdf(file=file.path(origin,"_2021","Figures","Figure_S1","Fig_1data2_RMSE_zlog.pdf"),width=15,height=10)
  {
  layout(matrix(c(1,1,2,2,3,
                  4,4,5,5,6),nrow=2,ncol=5,byrow=TRUE))
  #layout.show(nf)
  #par(mar=c(7,7,2,2),mfrow=c(1,1))
  par(mar=c(8,10,2,2))
  i=2
  
  #------------------------------TDtd Total= all values
  plot(dat_plot[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
       xaxt="n",yaxt="n", main="Total values",cex.main=2.5)
  axis(1,line = 2.5,cex.axis=2)
  axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
  axis(2,line = 2.5,cex.axis=2)
  axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
  abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
  
  i=1
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plot[,1],dat_plot[,colnames(dat_plot)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
  }  
  
  i=2
  # TDtd
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot)=="RMSE_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot)=="RMSE_zlog_Total"],col="black",lwd=5)
  
  
  #------------------------------TDenv Total = all values
  plot(dat_plot[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
       xaxt="n",yaxt="n", main="Total values (imputed with envelope)",cex.main=2.5)
  axis(1,line = 2.5,cex.axis=2)
  axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
  axis(2,line = 2.5,cex.axis=2)
  axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
  abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col="black",lwd=5)
  
  i=1
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
  }  
  
  i=2
  # TDenvelope
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col="black",lwd=5)
  
  par(mar=c(0,0,0,0))     
  plot(x = 1:10,y = 1:10,xaxt="n",yaxt="n",xlab="",ylab="",col="white",frame = FALSE)
  level_names = c(trait_names,"Mean")
  legend(1, 10, level_names, col = c(colz[1:length(trait_names)],"black"),bty = "n",
         text.col = "black", lty = c(rep(1,length(trait_names)),1), lwd=4,cex=3.15,
         merge = TRUE, bg = "white")  
  
  #-----------------------------------------------------------------------------------
  #------------------------------TDtd gap= only missing values
  par(mar=c(8,10,2,2))
  plot(dat_plot[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
       xaxt="n",yaxt="n", main="Gaps",cex.main=2.5)
  axis(1,line = 2.5,cex.axis=2)
  axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
  axis(2,line = 2.5,cex.axis=2)
  axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
  abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
  colnames(dat_plot)
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
  
  i=3
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plot[,1],dat_plot[,colnames(dat_plot)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],lwd=5)
  }  
  
  i=2
  # TDtd
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot)=="RMSE_gap_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plot[,1],dat_plot[,colnames(dat_plot)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
  
  
  #------------------------------TDenv  gap= only missing values
  plot(dat_plot[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE,
       xaxt="n",yaxt="n", main="Gaps (imputed with envelope)",cex.main=2.5)
  axis(1,line = 2.5,cex.axis=2)
  axis(1,line=5,at=c(40),labels = "Missingness [% gaps]",tick=FALSE,cex.axis=2)
  axis(2,line = 2.5,cex.axis=2)
  axis(2,line=5,at=.3,labels = "RMSE",tick=FALSE,cex.axis=3)
  abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
  colnames(dat_plotENV)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_gap_zlog_Total"],col="black",lwd=5)
  
  i=1
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],lwd=5)
  }  
  
  i=2
  # TDenvelope
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_gap_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_gap_zlog_Total"],col="black",lwd=5)

  }
  dev.off()
  

  
  
  
  
  
  
  
  #old: 
  
  pdf(file=file.path(origin,"_2021","Figures","Figure_1","Fig_1_RMSE.pdf"),width=10,height=8)
  par(mfrow=c(1,2),mar=c(7,7,2,2))
  i=2
  plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="Missingness [% gaps]",ylab="RMSE",xlim=c(0,110),ylim=c(0,.6),frame=FALSE)
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  colnames(dat_plot2)
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
  
  i=1
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
  }  
  
  i=2
  # TDtd
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
  # TDenvelope
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col="white",lwd=5.5)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col=colz[7],lwd=5)
  
  rect(xleft = 82,ybottom = 0,xright = 109,ytop = .95,col = "white",border = FALSE)
  level_names = c(trait_names,"","Average","Total","Missing only","With envelope")
  legend(83, .6, level_names, col = c(colz[1:length(trait_names)],"white","black","gray","gray",colz[7]),bty = "n",
         text.col = "black", lty = c(rep(1,length(trait_names)),1,1,1,3,1), lwd=4,cex=1.3,
         merge = TRUE, bg = "white")  
  
  #-----------------------------------------------------------------------------------
  plot(dat_plot2[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,
       xlab="Missingness [% gaps]",ylab="RMSE",xlim=c(0,110),ylim=c(0,.6),frame=FALSE)
  abline(h=seq(-1.2,1.2,by=.2),col="gray",lty=2)
  colnames(dat_plot2)
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"],col="black",lwd=3,lty=5)
  
  i=1
  for(i in 1:length(trait_names)){
    print(trait_names[i])
    lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=5)
    lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],lwd=3,lty=5)
  }  
  
  i=2
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_zlog_Total"],col="black",lwd=5)
  lines(dat_plot2[,1],dat_plot2[,colnames(dat_plot2)=="RMSE_gap_zlog_Total"],col="black",lwd=3,lty=5)
  
  
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_zlog_Total"],col=colz[7],lwd=5)
  lines(dat_plotENV[,1],dat_plotENV[,colnames(dat_plotENV)=="RMSE_gap_zlog_Total"],col=colz[7],lwd=3,lty=5)
  
  rect(xleft = 82,ybottom = 0,xright = 109,ytop = .95,col = "white",border = FALSE)
  level_names = c(trait_names,"","Average","Total","Missing only","With envelope")
  legend(83, .6, level_names, col = c(colz[1:length(trait_names)],"white","black","gray","gray",colz[7]),bty = "n",
         text.col = "black", lty = c(rep(1,length(trait_names)),1,1,1,3,1), lwd=4,cex=1.3,
         merge = TRUE, bg = "white")  
  dev.off()
  
  