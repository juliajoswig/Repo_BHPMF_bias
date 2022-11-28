

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
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)


GapPercent=50
RepNum=1

t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
res <- res[,colSums(!is.na(res))!=0]
trait_names=as.vector(unique(res$trait))
trait_names <- trait_names[!is.na(trait_names)]
missingness = unique(as.vector(res$missingness))
missingness <- missingness[!is.na(missingness)]
m=1
t=2
w=1

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
colz2=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
        "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
colzBoTop=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
        "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")

res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res$value_obs)
m=1
t=2
w=1
#colz <- list()
colzTRAITS <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")

bxpl_tot <- list()
bxpl_tot_op <- list()
colz <- list()
t=2
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from_"))
  colz_now = colz_now[-grep(x = colnames(res)[colz_now],pattern = "zlog")]
  colz_obs = colz_now[grep(x = colnames(res)[colz_now],pattern = "obs")]
  colz_pred = colz_now[grep(x = colnames(res)[colz_now],pattern = "pred")]
  colnames(res)[colz_obs]
  colnames(res)[colz_pred]
  colnames(res)[colz_now]
  obs  =    res$value_obs - res[ix_trait,colz_obs]
  pred =    res$value_pred - res[ix_trait,colz_pred]
  colnames(obs) <- gsub(colnames(obs),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(obs) <- gsub(colnames(obs),pattern = "_lm_obs",replacement = "")
  colnames(obs) <- paste0(colnames(obs),"_obs")
  colnames(pred) <- gsub(colnames(pred),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(pred) <- gsub(colnames(pred),pattern = "_lm_pred",replacement = "")
  colnames(pred) <- paste0(colnames(pred),"_pred")
  bxpl_tot[[t]]  <- pred-obs
  bxpl_tot_op[[t]]<- cbind(obs,pred)

  names(bxpl_tot)[t] <- trait_names[t]
  names(bxpl_tot_op)[t] <- trait_names[t]
  
#  colz[[t]] <- rep(colzTRAITS[match(unique(gsub(colnames(bxpl_tot[[t]])[1:5],pattern = "_obs",replacement = "")),trait_names)],2)[c(1,6,2,7,3,8,4,9,5,10)]
  colz[[t]] <- rep(colzTRAITS[match(unique(gsub(colnames(bxpl_tot[[t]])[1:5],pattern = "_pred",replacement = "")),trait_names)],2)
}
  
    

  par(mfrow=c(1,1),mar=c(8,4,3,2))
  
  pdf(file=file.path(origin,"_2021","figures","Figure_3",paste0("Figure_3_lm",missingness[m],"2.pdf")),width=18,height=5)
  {
    par(mar=c(10,3,4,0),mfrow=c(1,6))
    
    ytop <- list()
    ybot <- list()
    
    ytop[[1]]=c(40,65,40,65,5)
    ybot[[1]]=c(0,0,15,-0.1,-0.1)
    ytop[[2]]=c(20,65,25,65,2)
    ybot[[2]]=c(-0.1,-0.1,20,-0.1,-0.1)
    ytop[[3]]=c(1000,2500,300,50,10)
    ybot[[3]]=c(-0.1,-0.1,20,-0.1,-0.1)
    ytop[[4]]=c(12,55,.6,50,2)
    ybot[[4]]=c(-0.1,-0.1,-.1,-0.1,-0.1)
    ytop[[5]]=c(120,400,.6,200,2)
    ybot[[5]]=c(-0.1,-0.1,-.1,40,-0.1)
    ytop[[6]]=c(200,400,.6,70,2)
    ybot[[6]]=c(-0.1,-0.1,-.1,20,-0.1)
    
    
    t=6
    t2=2
    for(t in 1:6){
      plot(1:10,col="white",ylim=c(-2,2),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      #axis(side = 2,line = -17,cex.axis=2)
      axis(side = 2,at=seq(-2,2,by=2),line = -6,labels = c("",trait_names[t],""),tick = FALSE,cex.axis=4)
      axis(side = 2,at=seq(-2,2,by=2),line = -10,labels = c("","Distance",""),tick = FALSE,cex.axis=4)
      axis(side = 2,at = seq(-2,2,by=2),line = -14,labels = c("","to lm",""),tick = FALSE,cex.axis=4)
      
      for(t2 in 1:5){
      bxpl=abs(bxpl_tot_op[[t]][,c(1,6)+(t2-1)])
      boxplot(bxpl,cex.main=3,
              ylab="",col=colz[[t]][t2],xaxt="n",ylim=c(ybot[[t]][t2],ytop[[t]][t2]),
              las=2,cex.axis=1.5,cex.lab=1,lwd=2,frame=FALSE)
      yBot_now=apply(bxpl,2,FUN = median,na.rm = TRUE)
      axis(side = 1,at = seq(from=1,to = 2,by = 1),labels = colnames(bxpl),las=2,cex.axis=1.6)
      rect(xleft = .5,xright = 2.5,ybottom = ybot[[t]][t2],ytop = median(bxpl[,1]),col=colzBoTop[5])
      rect(xleft = .5,xright = 2.5,ybottom = median(bxpl[,1]),ytop = ytop[[t]][t2],col=colzBoTop[9])
      boxplot(bxpl,cex.main=3,
              ylab="",col=colz[[t]],xaxt="n",#ylim=c(ybot,ytop)
              las=2,cex.axis=1.5,cex.lab=1,lwd=2,frame=FALSE,add=TRUE)
      }
    }
    
  }
  dev.off()
  

  
  for(t in 1:6){
    bxpl=bxpl_tot[[t]]
    bxpl=bxpl_tot_op[[t]]
    boxplot(bxpl[,c(1,2)],main=trait_names[t],cex.main=3,
            ylim=c(ybot,ytop),ylab="CV pred - obs",col=colz[[t]],xaxt="n",
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    yBot_now=apply(abs(bxpl)[,1:5],2,FUN = median,na.rm = TRUE)
    axis(side = 1,at = seq(from=1,to = 5,by = 1),labels = gsub(colnames(bxpl_tot[[t]]),pattern = "_pred",replacement = ""),las=2,cex.axis=1.6)
    abline(h=0,col="gray",lty=1,lwd=3)
    boxplot(bxpl,main=trait_names[t],cex.main=3,
            ylim=c(ybot,ytop),ylab="CV pred - obs",col=colz[[t]],xaxt="n",
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE,add=TRUE)
    
  }
  
  
  
  
  for(t in 1:6){
    bxpl=bxpl_tot[[t]]
    
    boxplot(abs(bxpl[c(seq(1,ncol(bxpl),by=5),
                       seq(2,ncol(bxpl),by=5),
                       seq(3,ncol(bxpl),by=5),
                       seq(4,ncol(bxpl),by=5),
                       seq(5,ncol(bxpl),by=5))]),main=trait_names[t],cex.main=3,
            ylim=c(0,ytop),ylab="CV pred - obs",col=colz[[t]],xaxt="n",
            las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
    yBot_now=apply(abs(bxpl)[,1:5],2,FUN = median,na.rm = TRUE)
    i=1
    rect(xleft = 0.5,ybottom = 0,xright = 2.5+(i-1)*2,ytop = yBot_now[i],col=colz2[6],border = NA)
    rect(xleft = 0.5,ybottom = yBot_now[i],xright = 2.5+(i-1)*2,ytop = ytop,col=colz2[10],border = NA)
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
    rect(xleft = 8.5,ybottom = 0,xright = 10.5,ytop = yBot_now[i],col=colz2[6],border = NA)
    rect(xleft = 8.5,ybottom = yBot_now[i],xright = 10.5,ytop = ytop,col=colz2[10],border = NA)
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

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
