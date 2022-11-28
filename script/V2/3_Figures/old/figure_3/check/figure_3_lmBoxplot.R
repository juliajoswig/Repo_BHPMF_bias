
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
  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  Percent=80
  res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  head(res)
  
res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res$value_obs)
m=1
t=2
w=1
colz <- list()
colzTRAITS <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")

  {
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  colz_zlog = colz_now[grep(x = colnames(res)[colz_now],pattern = "zlog")]
  colz_now = colz_now[-grep(x = colnames(res)[colz_now],pattern = "zlog")]
  colnames(res[,colz_zlog])
  colnames(res[,colz_now])
  bxpl_SLA=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_SLA_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_zlog],res$value_pred_zlog-res[ix_trait,colz_zlog])
  
  t=2
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_H=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_H_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  
  t=3
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_SSD=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_SSD_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  
  t=4
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_N=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_N_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  
  t=5
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_P=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_P_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  
  t=6
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_Nar=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  bxpl_Nar_zlog=cbind(res$value_obs_zlog-res[ix_trait,colz_now],res$value_pred_zlog-res[ix_trait,colz_now])
  }
  
  colnames(bxpl_SLA) <- gsub(colnames(bxpl_SLA),pattern = paste0(trait_names[1],"_from_"),replacement = "")
  colnames(bxpl_H) <- gsub(colnames(bxpl_H),pattern = paste0(trait_names[2],"_from_"),replacement = "")
  colnames(bxpl_SSD) <- gsub(colnames(bxpl_SSD),pattern = paste0(trait_names[3],"_from_"),replacement = "")
  colnames(bxpl_N) <- gsub(colnames(bxpl_N),pattern = paste0(trait_names[4],"_from_"),replacement = "")
  colnames(bxpl_P) <- gsub(colnames(bxpl_P),pattern = paste0(trait_names[5],"_from_"),replacement = "")
  colnames(bxpl_Nar) <- gsub(colnames(bxpl_Nar),pattern = paste0(trait_names[6],"_from_"),replacement = "")
  
  colz <- list()
  colnames(bxpl_SLA) <- gsub(colnames(bxpl_SLA),pattern = "_lm",replacement = "")
  colz[[1]] <- rep(colzTRAITS[match(unique(colnames(bxpl_SLA)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(bxpl_H) <- gsub(colnames(bxpl_H),pattern = "_lm",replacement = "")
  colz[[2]] <- rep(colzTRAITS[match(unique(colnames(bxpl_H)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(bxpl_SSD) <- gsub(colnames(bxpl_SSD),pattern = "_lm",replacement = "")
  colz[[3]] <- rep(colzTRAITS[match(unique(colnames(bxpl_SSD)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(bxpl_N) <- gsub(colnames(bxpl_N),pattern = "_lm",replacement = "")
  colz[[4]] <- rep(colzTRAITS[match(unique(colnames(bxpl_N)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(bxpl_P) <- gsub(colnames(bxpl_P),pattern = "_lm",replacement = "")
  colz[[5]] <- rep(colzTRAITS[match(unique(colnames(bxpl_P)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colnames(bxpl_Nar) <- gsub(colnames(bxpl_Nar),pattern = "_lm",replacement = "")
  colz[[6]] <- rep(colzTRAITS[match(unique(colnames(bxpl_Nar)),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  
  par(mfrow=c(1,1),mar=c(10,0,3,0))
  
   pdf(file=file.path(origin,"_2021","figures","Figure_3",paste0("Figure_3_lm",missingness[m],"2.pdf")),width=20,height=8)
  {
    par(mar=c(2,0,4,0))
    layout(matrix(c(1:7,1:7,1:7,1:7,8,rep(9,6)), 
                  nrow = 5, ncol = 7, byrow = TRUE))
    ytop=4
    
    plot(1:10,col="white",ylim=c(0,ytop),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    axis(side = 2,line = -17,cex.axis=2)
    axis(side = 2,at = seq(0,4,by=2),line = -10,labels = c("","Distance",""),tick = FALSE,cex.axis=4)
    axis(side = 2,at = seq(0,4,by=2),line = -14,labels = c("","to lm",""),tick = FALSE,cex.axis=4)
    
    t=1
    for(t in 1:6){
      if(t==1){bxpl=bxpl_SLA}
      if(t==2){bxpl=bxpl_H}
      if(t==3){bxpl=bxpl_SSD}
      if(t==4){bxpl=bxpl_N}
      if(t==5){bxpl=bxpl_P}
      if(t==6){bxpl=bxpl_Nar}
    
      boxplot(abs(bxpl[c(seq(1,ncol(bxpl),by=5),
                       seq(2,ncol(bxpl),by=5),
                       seq(3,ncol(bxpl),by=5),
                       seq(4,ncol(bxpl),by=5),
                       seq(5,ncol(bxpl),by=5))]),main=trait_names[t],cex.main=3,
            ylim=c(0,ytop),ylab="CV pred - obs",col=colz[[t]],xaxt="n",
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    yBot_now=apply(abs(bxpl)[,1:5],2,FUN = median,na.rm = TRUE)
      i=1
      rect(xleft = 0,ybottom = 0,xright = 2.5+(i-1)*2,ytop = yBot_now[i],col=colz2[6],border = NA)
      rect(xleft = 0,ybottom = yBot_now[i],xright = 2.5+(i-1)*2,ytop = ytop,col=colz2[10],border = NA)
      i=2
      rect(xleft = 2.5,ybottom = 0,xright = 2.5+(i-1)*2,ytop = yBot_now[i],col=colz2[6],border = NA)
      rect(xleft = 2.5,ybottom = yBot_now[i],xright = 2.5+(i-1)*2,ytop = ytop,col=colz2[10],border = NA)
      i=3
      rect(xleft = 4.5,ybottom = 0,xright = 2.5+(i-1)*2,ytop = yBot_now[i],col=colz2[6],border = NA)
      rect(xleft = 4.5,ybottom = yBot_now[i],xright = 2.5+(i-1)*2,ytop = ytop,col=colz2[10],border = NA)
      i=4
      rect(xleft = 6.5,ybottom = 0,xright = 2.5+(i-1)*2,ytop = yBot_now[i],col=colz2[6],border = NA)
      rect(xleft = 6.5,ybottom = yBot_now[i],xright = 2.5+(i-1)*2,ytop = ytop,col=colz2[10],border = NA)
      i=5
      rect(xleft = 8.5,ybottom = 0,xright = 11,ytop = yBot_now[i],col=colz2[6],border = NA)
      rect(xleft = 8.5,ybottom = yBot_now[i],xright = 11,ytop = ytop,col=colz2[10],border = NA)
    boxplot(abs(bxpl[c(seq(1,ncol(bxpl),by=5),
                       seq(2,ncol(bxpl),by=5),
                       seq(3,ncol(bxpl),by=5),
                       seq(4,ncol(bxpl),by=5),
                       seq(5,ncol(bxpl),by=5))]),
            ylim=c(0,ytop),ylab="CV pred - obs",col=colz[[t]],xaxt="n",
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE,add=TRUE)
    axis(side = 1,at = seq(from=1,to = 10,by = 1),labels = c("obs","pred","obs","pred","obs","pred","obs","pred","obs","pred"),las=2,cex.axis=1.6)
    #axis(side = 1,line = 5, at = seq(from=1.5,to = 10,by = 2),las=2,labels=colnames(bxpl)[seq(from=1,to = ncol(bxpl),by = 2)],tick = FALSE,cex.axis=2)
    abline(v=seq(.5,11,2),col="white",lty=1,lwd=3)
    
    }
    
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    rect(xleft = 2,ybottom = 0,xright = 10,ytop = 10,col = "lightgray",border = FALSE)
    legend(3, 11, "SLA", col = colzTRAITS[1], pch = 15, cex = 3,bty = "n")
    legend(3.9,11, "Height", col = colzTRAITS[2], pch = 15, cex = 3,bty = "n")
    legend(5,11, "SSD", col = colzTRAITS[3], pch = 15, cex = 3,bty = "n")
    legend(6,11, "LeafN", col = colzTRAITS[4], pch = 15, cex = 3,bty = "n")
    legend(7,11, "LeafP", col = colzTRAITS[5], pch = 15, cex = 3,bty = "n")
    legend(8,11, "LeafNArea", col = colzTRAITS[6], pch = 15, cex = 3,bty = "n")
  }
  dev.off()
  
  par(mar=c(0,0,0,0))
  layout(matrix(c(1:7,rep(8,7)), 
                nrow = 2, ncol = 7, byrow = TRUE))
  
  layout(matrix(c(1,1,2,3,4,4), 
                nrow = 3, ncol = 2, byrow = TRUE))
  plot(1,main=1)
  plot(2,main=2)
  plot(3,main=3)
  plot(4,main=4)
  
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_SLA = abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])

  t=2
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_H=abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])

  t=3
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_SSD=abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])

  t=4
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_N=abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])

  t=5
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_P=abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])

  t=6
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_Nar=abs(res$value_pred_zlog-res[ix_trait,colz_now])-
    abs(res$value_obs_zlog-res[ix_trait,colz_now])
  
  colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
  bxpl_SLA <- bxpl_SLA[,colSums(!is.na(bxpl_SLA))!=0]
  bxpl_H <- bxpl_H[,colSums(!is.na(bxpl_H))!=0]
  bxpl_N <- bxpl_N[,colSums(!is.na(bxpl_N))!=0]
  bxpl_Nar <- bxpl_Nar[,colSums(!is.na(bxpl_Nar))!=0]
  bxpl_SSD <- bxpl_SSD[,colSums(!is.na(bxpl_SSD))!=0]
  bxpl_P <- bxpl_P[,colSums(!is.na(bxpl_P))!=0]
  
  colnames(bxpl_SLA) <- gsub(colnames(bxpl_SLA),pattern = paste0(trait_names[1],"_from_"),replacement = "")
  colnames(bxpl_H) <- gsub(colnames(bxpl_H),pattern = paste0(trait_names[2],"_from_"),replacement = "")
  colnames(bxpl_SSD) <- gsub(colnames(bxpl_SSD),pattern = paste0(trait_names[3],"_from_"),replacement = "")
  colnames(bxpl_N) <- gsub(colnames(bxpl_N),pattern = paste0(trait_names[4],"_from_"),replacement = "")
  colnames(bxpl_P) <- gsub(colnames(bxpl_P),pattern = paste0(trait_names[5],"_from_"),replacement = "")
  colnames(bxpl_Nar) <- gsub(colnames(bxpl_Nar),pattern = paste0(trait_names[6],"_from_"),replacement = "")
  colnames(bxpl_SLA) <- gsub(colnames(bxpl_SLA),pattern = "_lm",replacement = "")
  colnames(bxpl_H) <- gsub(colnames(bxpl_H),pattern = "_lm",replacement = "")
  colnames(bxpl_SSD) <- gsub(colnames(bxpl_SSD),pattern = "_lm",replacement = "")
  colnames(bxpl_N) <- gsub(colnames(bxpl_N),pattern = "_lm",replacement = "")
  colnames(bxpl_P) <- gsub(colnames(bxpl_P),pattern = "_lm",replacement = "")
  colnames(bxpl_Nar) <- gsub(colnames(bxpl_Nar),pattern = "_lm",replacement = "")
  
  pdf(file=file.path(origin,"_2021","figures","Figure_3",paste0("Figure_3_lm",missingness[m],".pdf")),width=10,height=8)
  {
    par(mfrow=c(1,1))
    par(mfrow=c(1,8),mar=c(10,0,0,1))
    
    plot(1:10,col="white",ylim=c(-.5,.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    axis(side = 2,line = -7,cex.axis=2)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -5.5,labels = c("","","pred - obs",""),tick = FALSE,cex.axis=3)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -3.25,labels = c("","","Distance",""),tick = FALSE,cex.axis=3)
    
    bxpl=bxpl_SLA
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "SLA",cex=3)
    
    bxpl=bxpl_H
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "Height",cex=3)
    
    bxpl=bxpl_SSD
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "SSD",cex=3)
    
    
    bxpl=bxpl_N
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "LeafN",cex=3)
    
    bxpl=bxpl_P
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "LeafP",cex=3)
    
    bxpl=bxpl_Nar
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0,lwd=2)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    rect(xleft = 0,ybottom = .4,xright = ncol(bxpl)+2,ytop =.6,col = "white",border = FALSE)
    text(x=3,y=.45,labels = "LeafNArea",cex=2)
    
    
    plot(1:10,col="white",ylim=c(-.5,.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = 7,ytop =.4,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = 7,ytop =0,col = colz[4],border = FALSE)
    text(2.5,.25,labels="less",col="black",srt=90,cex=3)
    text(5,.25,labels="weight",col="black",srt=90,cex=3)
    text(2.5,-.25,labels="more",col="black",srt=90,cex=3)
    text(5,-.25,labels="weight",col="black",srt=90,cex=3)
  }
  dev.off()

  

  
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl_SLA <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_SLA=cbind(bxpl_SLA,
                 abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                   abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  t=2
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  print(sum(ix_trait))
  bxpl_H <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_H=cbind(bxpl_H,
               abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                 abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  t=3
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl_SSD <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_SSD=cbind(bxpl_SSD,
                 abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                   abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  t=4
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl_N <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_N=cbind(bxpl_N,
               abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                 abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  t=5
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl_P <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_P=cbind(bxpl_P,
               abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                 abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  t=6
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]&res$Gap
  ix_trait[is.na(ix_trait)] <- FALSE
  bxpl_Nar <- matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(sum(ix_trait))
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from"))
  bxpl_Nar=cbind(bxpl_Nar,
                 abs(res[ix_trait,colz_now]-res$value_pred_zlog)-
                   abs(res[ix_trait,colz_now]-res$value_obs_zlog))
  
  
  colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
  bxpl_SLA <- bxpl_SLA[,colSums(!is.na(bxpl_SLA))!=0]
  bxpl_H <- bxpl_H[,colSums(!is.na(bxpl_H))!=0]
  bxpl_N <- bxpl_N[,colSums(!is.na(bxpl_N))!=0]
  bxpl_Nar <- bxpl_Nar[,colSums(!is.na(bxpl_Nar))!=0]
  bxpl_SSD <- bxpl_SSD[,colSums(!is.na(bxpl_SSD))!=0]
  bxpl_P <- bxpl_P[,colSums(!is.na(bxpl_P))!=0]
  
  colnames(bxpl_SLA) <- gsub(colnames(bxpl_SLA),pattern = paste0(trait_names[1],"_from_"),replacement = "")
  colnames(bxpl_H) <- gsub(colnames(bxpl_H),pattern = paste0(trait_names[2],"_from_"),replacement = "")
  colnames(bxpl_SSD) <- gsub(colnames(bxpl_SSD),pattern = paste0(trait_names[3],"_from_"),replacement = "")
  colnames(bxpl_N) <- gsub(colnames(bxpl_N),pattern = paste0(trait_names[4],"_from_"),replacement = "")
  colnames(bxpl_P) <- gsub(colnames(bxpl_P),pattern = paste0(trait_names[5],"_from_"),replacement = "")
  colnames(bxpl_Nar) <- gsub(colnames(bxpl_Nar),pattern = paste0(trait_names[6],"_from_"),replacement = "")
  
  pdf(file=file.path(origin,"_2021","figures","Figure_3",paste0("Figure_3_lm",missingness[m],"ONLYGAPS.pdf")),width=20,height=8)
  {
    par(mfrow=c(1,1))
    par(mfrow=c(1,8),mar=c(10,0,0,1))
    
    plot(1:10,col="white",ylim=c(-.5,.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    axis(side = 2,line = -7,cex.axis=2)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -5.5,labels = c("","","pred - obs",""),tick = FALSE,cex.axis=3)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -3.25,labels = c("","","Distance",""),tick = FALSE,cex.axis=3)
    
    bxpl=bxpl_SLA
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "SLA",cex=3)
    
    bxpl=bxpl_H
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "Height",cex=3)
    
    bxpl=bxpl_SSD
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "SSD",cex=3)
    
    
    bxpl=bxpl_N
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "LeafN",cex=3)
    
    bxpl=bxpl_P
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "LeafP",cex=3)
    
    bxpl=bxpl_Nar
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl)+2,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl)+2,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    boxplot((bxpl),ylim=c(-.5,.5),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    text(x=3,y=2.1,labels = "LeafNArea",cex=3)
    
    
    plot(1:10,col="white",ylim=c(-.5,.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = 7,ytop =2,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = 7,ytop =0,col = colz[4],border = FALSE)
    abline(h=0)
    text(2.5,.5,labels="less",col="black",srt=90,cex=3)
    text(5,.5,labels="weight",col="black",srt=90,cex=3)
    text(2.5,-.5,labels="more",col="black",srt=90,cex=3)
    text(5,-.5,labels="weight",col="black",srt=90,cex=3)
  }
  dev.off()
  
}
boxplot(bxpl_sp,bxpl_gen,bxpl_fam,bxpl_clad,bxpl_GF,bxpl_PFT,rep(NA,length(bxpl_PFT)),
        ylim=c(-5,5),ylab="CV pred - obs",
        col=colz,
        las=2,cex.axis=1.5,cex.lab=2,lwd=2)
axis(side = 1,at = 1:7,labels =  c("Species","Genus","Family","Clade","Growth Form","PFT",""),las=2,tick = FALSE)
abline(h=0,col="gray",lty=2,lwd=2)
rect(xleft = 6.7,ybottom = 0,xright = 10,ytop = 10,col = colz[1],border = FALSE)
rect(xleft = 6.7,ybottom = -10,xright = 10,ytop = 0,col = colz[6],border = FALSE)
text(7.2,2.5,labels="more clustering",col="white",srt=90,cex=2)
text(7.2,-2.5,labels="less clustering",col="white",srt=90,cex=2)

#bxpl <- cbind(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op)
colnames(bxpl) <- c("Species_pred","Species_pred","Genus_obs","Genus_pred",
                    "Family_obs","Family_pred","Clade_obs","Clade_pred","GF_obs","GF_pred",
                    "PFT_obs","PFT_pred")
boxplot(bxpl_sp_op,ylim=c(-5,5),col=c(colz[1],colz[1]))
boxplot(bxpl_gen_op,ylim=c(-5,5),col=c(colz[2],colz[2]))
boxplot(bxpl_fam_op,ylim=c(-5,5),col=c(colz[3],colz[3]))
boxplot(bxpl_clad_op,ylim=c(-5,5),col=c(colz[4],colz[4]))
boxplot(bxpl_GF_op,ylim=c(-5,5),col=c(colz[4],colz[4]))
boxplot(bxpl_PFT_op,ylim=c(-5,5),col=c(colz[4],colz[4]))

boxplot(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op,
        ylim=c(0,15),ylab="CV",
        #        col=c(colz[1],colz[1],colz[2],colz[2],colz[3],colz[3],colz[4],colz[4],colz[5],colz[5],colz[6],colz[6]),
        las=2,cex.axis=1,cex.lab=2)
abline(h=0)
abline(v=seq(.5,to=50,by=2),col="gray",lwd=3)
dev.off()
}




for(m in 1:length(missingness)){
  pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_CV",missingness[m],".pdf")),width=6,height=8)
  par(mfrow=c(1,1),mar=c(10,6,2,2))
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  sum(ix_trait)
  
  bxpl_sp=NA
  bxpl_gen=NA#rep(NA,sum(ix_trait))
  bxpl_fam=NA#rep(NA,sum(ix_trait))
  bxpl_clad=NA#rep(NA,sum(ix_trait))
  bxpl_GF=NA#rep(NA,sum(ix_trait))
  bxpl_PFT=NA#rep(NA,sum(ix_trait))
  bxpl_sp_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Species)[ix_trait]),ncol=2)
  bxpl_gen_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Genus)[ix_trait]),ncol=2)
  bxpl_fam_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Family)[ix_trait]),ncol=2)
  bxpl_clad_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Clade)[ix_trait]),ncol=2)
  bxpl_GF_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$GF)[ix_trait]),ncol=2)
  bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$PFT)[ix_trait]),ncol=2)
  
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    sum(ix_trait)
    
    bxpl_sp=c(bxpl_sp,
              abs(res$CV_spec_pred[ix_trait&!duplicated(res$Species)])-
                abs(res$CV_spec_obs[ix_trait&!duplicated(res$Species)]))
    bxpl_gen=c(bxpl_gen,
               abs(res$CV_gen_pred[ix_trait&!duplicated(res$Genus)])-
                 abs(res$CV_gen_obs[ix_trait&!duplicated(res$Genus)]))
    bxpl_fam=c(bxpl_fam,
               abs(res$CV_fam_pred[ix_trait&!duplicated(res$Family)])-
                 abs(res$CV_fam_obs[ix_trait&!duplicated(res$Family)]))
    bxpl_clad=c(bxpl_clad,
                abs(res$CV_clad_pred[ix_trait&!duplicated(res$Clade)]) -
                  abs(res$CV_clad_obs[ix_trait&!duplicated(res$Clade)]))
    bxpl_GF=c(bxpl_GF,
              abs(res$CV_GF_pred[ix_trait&!duplicated(res$GF)])-
                abs(res$CV_GF_obs[ix_trait&!duplicated(res$GF)]))
    bxpl_PFT=c(bxpl_PFT,
               abs(res$CV_PFT_pred[ix_trait&!duplicated(res$PFT)])-
                 abs(res$CV_PFT_obs[ix_trait&!duplicated(res$PFT)]))
    
    bxpl_sp_op=rbind(bxpl_sp_op,
                     cbind(res$CV_spec_obs[ix_trait&!duplicated(res$Species)],
                           res$CV_spec_pred[ix_trait&!duplicated(res$Species)]))
    bxpl_gen_op=rbind(bxpl_gen_op,
                      cbind(res$CV_gen_obs[ix_trait&!duplicated(res$Genus)],
                            res$CV_gen_pred[ix_trait&!duplicated(res$Genus)]))
    bxpl_fam_op=rbind(bxpl_fam_op,
                      cbind(res$CV_fam_obs[ix_trait&!duplicated(res$Family)],
                            res$CV_fam_pred[ix_trait&!duplicated(res$Family)]))
    bxpl_clad_op=rbind(bxpl_clad_op,
                       cbind(res$CV_clad_obs[ix_trait&!duplicated(res$Clade)],
                             res$CV_clad_pred[ix_trait&!duplicated(res$Clade)]))
    bxpl_GF_op=rbind(bxpl_GF_op,
                     cbind(res$CV_GF_obs[ix_trait&!duplicated(res$GF)],
                           res$CV_GF_pred[ix_trait&!duplicated(res$GF)]))
    bxpl_PFT_op=rbind(bxpl_PFT_op,
                      cbind(res$CV_PFT_obs[ix_trait&!duplicated(res$PFT)],
                            res$CV_PFT_pred[ix_trait&!duplicated(res$PFT)]))
    
  }   
  
  #bxpl <- cbind(abs(bxpl_sp),abs(bxpl_gen),abs(bxpl_fam),abs(bxpl_clad),bxpl_GF,bxpl_PFT,rep(NA,length(bxpl_PFT)))
  #dim(bxpl_sp)
  #colnames(bxpl) <- c("Species","Genus","Family","Clade","Growth Form","PFT","")
  boxplot(bxpl_sp,bxpl_gen,bxpl_fam,bxpl_clad,bxpl_GF,bxpl_PFT,rep(NA,length(bxpl_PFT)),
          ylim=c(-5,5),ylab="CV pred - obs",
          col=colz,
          las=2,cex.axis=1.5,cex.lab=2,lwd=2)
  axis(side = 1,at = 1:7,labels =  c("Species","Genus","Family","Clade","Growth Form","PFT",""),las=2,tick = FALSE)
  abline(h=0,col="gray",lty=2,lwd=2)
  rect(xleft = 6.7,ybottom = 0,xright = 10,ytop = 10,col = colz[1],border = FALSE)
  rect(xleft = 6.7,ybottom = -10,xright = 10,ytop = 0,col = colz[6],border = FALSE)
  text(7.2,2.5,labels="more clustering",col="white",srt=90,cex=2)
  text(7.2,-2.5,labels="less clustering",col="white",srt=90,cex=2)
  
  #bxpl <- cbind(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op)
  colnames(bxpl) <- c("Species_pred","Species_pred","Genus_obs","Genus_pred",
                      "Family_obs","Family_pred","Clade_obs","Clade_pred","GF_obs","GF_pred",
                      "PFT_obs","PFT_pred")
  boxplot(bxpl_sp_op,ylim=c(-5,5),col=c(colz[1],colz[1]))
  boxplot(bxpl_gen_op,ylim=c(-5,5),col=c(colz[2],colz[2]))
  boxplot(bxpl_fam_op,ylim=c(-5,5),col=c(colz[3],colz[3]))
  boxplot(bxpl_clad_op,ylim=c(-5,5),col=c(colz[4],colz[4]))
  boxplot(bxpl_GF_op,ylim=c(-5,5),col=c(colz[4],colz[4]))
  boxplot(bxpl_PFT_op,ylim=c(-5,5),col=c(colz[4],colz[4]))
  
  boxplot(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op,
          ylim=c(0,15),ylab="CV",
          #        col=c(colz[1],colz[1],colz[2],colz[2],colz[3],colz[3],colz[4],colz[4],colz[5],colz[5],colz[6],colz[6]),
          las=2,cex.axis=1,cex.lab=2)
  abline(h=0)
  abline(v=seq(.5,to=50,by=2),col="gray",lwd=3)
  dev.off()
}

