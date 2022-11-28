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
list.files(file.path(origin,"script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin = origin,Version_now = "V1")

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
ObsOrTDs =  out$ObsOrTDs
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido =out$trait_guido
trait_rainfor =out$trait_rainfor
colz1 = out$colz1
colz2 = out$colz2
gappercents= c("1","5","10","20","30","40","50","60")
repnums=1:3

GapPercent=50
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
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

colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")


  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res$value_obs)
  m=1
  t=2
  w=1

    par(mfrow=c(1,1),mar=c(10,6,2,2))
    t=1
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    sum(ix_trait)
    
    bxpl_sp=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Species)[ix_trait]),ncol=1)
    bxpl_gen=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Genus)[ix_trait]),ncol=1)#rep(NA,sum(ix_trait))
    bxpl_fam=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Family)[ix_trait]),ncol=1)
    bxpl_clad=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Clade)[ix_trait]),ncol=1)
    bxpl_GF=matrix(NA,nrow=sum(ix_trait&!duplicated(res$GF)[ix_trait]),ncol=1)
    bxpl_PFT=matrix(NA,nrow=sum(ix_trait&!duplicated(res$PFT)[ix_trait]),ncol=1)
    bxpl_sp_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Species)[ix_trait]),ncol=1)
    bxpl_gen_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Genus)[ix_trait]),ncol=1)
    bxpl_fam_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Family)[ix_trait]),ncol=1)
    bxpl_clad_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$Clade)[ix_trait]),ncol=1)
    bxpl_GF_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$GF)[ix_trait]),ncol=1)
    bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res$PFT)[ix_trait]),ncol=1)
    print(head(bxpl_sp))
    
    for(t in 1:length(trait_names)){
      ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
      ix_trait[is.na(ix_trait)] <- FALSE
      print(sum(ix_trait))
      
      bxpl_sp=cbind(bxpl_sp,
                    res$CV_spec_pred[ix_trait&!duplicated(res$Species[ix_trait])]-
                      res$CV_spec_obs[ix_trait&!duplicated(res$Species[ix_trait])])
      bxpl_gen=cbind(bxpl_gen,
                 res$CV_gen_pred[ix_trait&!duplicated(res$Genus[ix_trait])]-
                   res$CV_gen_obs[ix_trait&!duplicated(res$Genus[ix_trait])])
      bxpl_fam=cbind(bxpl_fam,
                 res$CV_fam_pred[ix_trait&!duplicated(res$Family[ix_trait])]-
                   res$CV_fam_obs[ix_trait&!duplicated(res$Family[ix_trait])])
      bxpl_clad=cbind(bxpl_clad,
                 res$CV_clad_pred[ix_trait&!duplicated(res$Clade[ix_trait])] -
                 res$CV_clad_obs[ix_trait&!duplicated(res$Clade[ix_trait])])
      bxpl_GF=cbind(bxpl_GF,
                res$CV_GF_pred[ix_trait&!duplicated(res$GF[ix_trait])]-
                res$CV_GF_obs[ix_trait&!duplicated(res$GF[ix_trait])])
      bxpl_PFT=cbind(bxpl_PFT,
                 res$CV_PFT_pred[ix_trait&!duplicated(res$PFT[ix_trait])]-
                 res$CV_PFT_obs[ix_trait&!duplicated(res$PFT[ix_trait])])
      
      bxpl_sp_op=cbind(bxpl_sp_op,
                       cbind(res$CV_spec_obs[ix_trait&!duplicated(res$Species[ix_trait])],
                             res$CV_spec_pred[ix_trait&!duplicated(res$Species[ix_trait])]))
      bxpl_gen_op=cbind(bxpl_gen_op,
                        cbind(res$CV_gen_obs[ix_trait&!duplicated(res$Genus[ix_trait])],
                              res$CV_gen_pred[ix_trait&!duplicated(res$Genus[ix_trait])]))
      bxpl_fam_op=cbind(bxpl_fam_op,
                        cbind(res$CV_fam_obs[ix_trait&!duplicated(res$Family[ix_trait])],
                              res$CV_fam_pred[ix_trait&!duplicated(res$Family[ix_trait])]))
      bxpl_clad_op=cbind(bxpl_clad_op,
                         cbind(res$CV_clad_obs[ix_trait&!duplicated(res$Clade[ix_trait])],
                               res$CV_clad_pred[ix_trait&!duplicated(res$Clade[ix_trait])]))
      bxpl_GF_op=cbind(bxpl_GF_op,
                       cbind(res$CV_GF_obs[ix_trait&!duplicated(res$GF[ix_trait])],
                             res$CV_GF_pred[ix_trait&!duplicated(res$GF[ix_trait])]))
      bxpl_PFT_op=cbind(bxpl_PFT_op,
                        cbind(res$CV_PFT_obs[ix_trait&!duplicated(res$PFT[ix_trait])],
                              res$CV_PFT_pred[ix_trait&!duplicated(res$PFT[ix_trait])]))
      
    }   
    
    bxpl_sp <- bxpl_sp[,colSums(!is.na(bxpl_sp))!=0]
    bxpl_gen <- bxpl_gen[,colSums(!is.na(bxpl_gen))!=0]
    bxpl_fam <- bxpl_fam[,colSums(!is.na(bxpl_fam))!=0]
    bxpl_clad <- bxpl_clad[,colSums(!is.na(bxpl_clad))!=0]
    bxpl_GF <- bxpl_GF[,colSums(!is.na(bxpl_GF))!=0]
    bxpl_PFT <- bxpl_PFT[,colSums(!is.na(bxpl_PFT))!=0]
    colnames(bxpl_sp) <- trait_names
    colnames(bxpl_gen) <- trait_names
    colnames(bxpl_fam) <- trait_names
    colnames(bxpl_clad) <- trait_names
    colnames(bxpl_GF) <- trait_names
    colnames(bxpl_PFT) <- trait_names
    
    colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
           "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
    bxpl_sp_op <- bxpl_sp_op[,colSums(!is.na(bxpl_sp_op))!=0]
    bxpl_gen_op <- bxpl_gen_op[,colSums(!is.na(bxpl_gen_op))!=0]
    bxpl_fam_op <- bxpl_fam_op[,colSums(!is.na(bxpl_fam_op))!=0]
    bxpl_clad_op <- bxpl_clad_op[,colSums(!is.na(bxpl_clad_op))!=0]
    bxpl_GF_op <- bxpl_GF_op[,colSums(!is.na(bxpl_GF_op))!=0]
    bxpl_PFT_op <- bxpl_PFT_op[,colSums(!is.na(bxpl_PFT_op))!=0]
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_CV",missingness[m],"2.pdf")),width=20,height=8)
    
    par(mar=c(4,0,4,0))
    layout(matrix(c(1:7,1:7,1:7,1:7,8,rep(9,6)), 
                  nrow = 5, ncol = 7, byrow = TRUE))
    {
      ymax=1.4
      plot(1:10,col="white",ylim=c(0,ymax),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -17,cex.axis=2)
      axis(side = 2,at = .75,line = -14,labels = "CV",tick = FALSE,cex.axis=4)

      names_now=c("Species","Genus","Family","Clades","GF","PFT")
      bxpl=bxpl_sp_op
      t=1
      for(t in 1:6){
        if(t==1){bxpl=bxpl_sp_op}
        if(t==2){bxpl=bxpl_gen_op}
        if(t==3){bxpl=bxpl_fam_op}
        if(t==4){bxpl=bxpl_clad_op}
        if(t==5){bxpl=bxpl_GF_op}
        if(t==6){bxpl=bxpl_PFT_op}
        dat_now=bxpl[,c(seq(1,to=ncol(bxpl),by = 6),
                        seq(2,to=ncol(bxpl),by = 6),
                        seq(3,to=ncol(bxpl),by = 6),
                        seq(4,to=ncol(bxpl),by = 6),
                        seq(5,to=ncol(bxpl),by = 6),
                        seq(6,to=ncol(bxpl),by = 6))]
        boxplot(dat_now,ylim=c(0,ymax),ylab="CV pred - obs",col=colz,main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
        axis(1,at = 1:ncol(bxpl),labels = rep(c("obs","pred"),6),las=2,cex.axis=2)
        yBot_now=apply(dat_now[,seq(1,ncol(dat_now),2)],2,FUN = median,na.rm = TRUE)
        i=1
        rect(xleft = .5,ybottom = 0,xright = 2.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = .5,ybottom = yBot_now[i],xright = 2.5,ytop = ymax,col=colz[10],border = NA)
        i=2
        rect(xleft = 2.5,ybottom = 0,xright = 4.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = 2.5,ybottom = yBot_now[i],xright = 4.5,ytop = ymax,col=colz[10],border = NA)
        i=3
        rect(xleft = 4.5,ybottom = 0,xright = 6.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = 4.5,ybottom = yBot_now[i],xright = 6.5,ytop = ymax,col=colz[10],border = NA)
        i=4
        rect(xleft = 6.5,ybottom = 0,xright = 8.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = 6.5,ybottom = yBot_now[i],xright = 8.5,ytop = ymax,col=colz[10],border = NA)
        i=5
        rect(xleft = 8.5,ybottom = 0,xright = 10.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = 8.5,ybottom = yBot_now[i],xright = 10.5,ytop = ymax,col=colz[10],border = NA)
        i=6
        rect(xleft = 10.5,ybottom = 0,xright = 12.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = 10.5,ybottom = yBot_now[i],xright = 12.5,ytop = ymax,col=colz[10],border = NA)
        
        boxplot(dat_now,ylim=c(0,1),ylab="CV pred - obs",col=colz,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n",add=TRUE)
        abline(v=seq(from = .5,to = 12.5,by = 2),lwd=4,col="white")
#        text(x=6.5,y=ymax-.5,labels = names_now[t],cex=3)
        
      }
 
    
      plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
      plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
      rect(xleft = 2,ybottom = 0,xright = 10,ytop = 12,col = "lightgray",border = FALSE)
      legend(3, 14, "SLA", col = colz[1], pch = 15, cex = 3,bty = "n")
      legend(3.9,14, "Height", col = colz[3], pch = 15, cex = 3,bty = "n")
      legend(5, 14, "SSD", col = colz[5], pch = 15, cex = 3,bty = "n")
      legend(6, 14, "LeafN", col = colz[7], pch = 15, cex = 3,bty = "n")
      legend(7, 14, "LeafP", col = colz[9], pch = 15, cex = 3,bty = "n")
      legend(8, 14, "LeafNArea", col = colz[11], pch = 15, cex = 3,bty = "n")
    }
    
    dev.off()
    
    #-----------------------------------------------------------------------------------------------------
    colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_CV",missingness[m],"3.pdf")),width=10,height=8)

    par(mar=c(2,0,2,0))
    layout(matrix(c(1:8,1:8,1:8,1:8,9,rep(10,6)), 
                nrow = 5, ncol = 8, byrow = TRUE))
  {
    plot(1:10,col="white",ylim=c(-.75,.75),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    axis(side = 2,line = -7,cex.axis=2)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -5.5,labels = c("","","pred - obs",""),tick = FALSE,cex.axis=3)
    axis(side = 2,at = seq(-1,.5,by=.5),line = -3.25,labels = c("","","CV",""),tick = FALSE,cex.axis=3)
    
    boxplot((bxpl_sp),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_sp),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
    abline(h=median((bxpl_sp),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_sp),na.rm = TRUE),lwd=3)
    text(x=3.5,y=0.55,labels = "Species",cex=3)
    
    boxplot((bxpl_gen),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_gen),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE,xaxt="n")
    abline(h=median((bxpl_gen),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_gen),na.rm = TRUE),lwd=3)
    text(x=3.5,y=0.55,labels = "Genus",cex=3)

    boxplot((bxpl_fam),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_fam),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE,xaxt="n")
    abline(h=median((bxpl_fam),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_fam),na.rm = TRUE),lwd=3)
    text(x=3.5,y=0.55,labels = "Family",cex=3)
    
    boxplot((bxpl_clad),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_clad),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE,xaxt="n")
    abline(h=median((bxpl_clad),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_clad),na.rm = TRUE),lwd=3)
    text(x=3.5,y=0.55,labels = "Clade",cex=3)

    boxplot((bxpl_GF),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_GF),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE,xaxt="n")
    abline(h=median((bxpl_GF),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_GF),na.rm = TRUE),lwd=3)
    text(x=3.5,y=0.55,labels = "GF",cex=3)
    
    boxplot((bxpl_PFT),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
    rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
    boxplot((bxpl_PFT),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
            las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE,xaxt="n")
    abline(h=median((bxpl_PFT),na.rm = TRUE),lwd=4,col="white")
    abline(h=median((bxpl_PFT),na.rm = TRUE),lwd=3)
    text(x=3.5,y=.55,labels = "PFT",cex=3)
    

    plot(1:10,col="white",ylim=c(-.75,.75),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    rect(xleft = 0,ybottom = 0,xright = 7,ytop =.5,col = colz[5],border = FALSE)
    rect(xleft = 0,ybottom = -5,xright = 7,ytop =0,col = colz[4],border = FALSE)
    text(2.5,.25,labels="less",col="black",srt=90,cex=3)
    text(5,.25,labels="clustering",col="black",srt=90,cex=3)
    text(2.5,-.25,labels="more",col="black",srt=90,cex=3)
    text(5,-.25,labels="clustering",col="black",srt=90,cex=3)
    }
    
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    rect(xleft = 2,ybottom = 0,xright = 10,ytop = 10,col = "lightgray",border = FALSE)
    legend(3, 10, "SLA", col = colz[1], pch = 15, cex = 3,bty = "n")
    legend(3.9, 10, "Height", col = colz[3], pch = 15, cex = 3,bty = "n")
    legend(5, 10, "SSD", col = colz[5], pch = 15, cex = 3,bty = "n")
    legend(6, 10, "LeafN", col = colz[7], pch = 15, cex = 3,bty = "n")
    legend(7, 10, "LeafP", col = colz[9], pch = 15, cex = 3,bty = "n")
    legend(8,10, "LeafNArea", col = colz[11], pch = 15, cex = 3,bty = "n")
    
    dev.off()
    


  
  
  
  
    
  
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
  boxplot(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op,rep(NA,length(bxpl_PFT)),
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
  
 