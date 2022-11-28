
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


#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
  file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
  colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
  colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")

res <- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res2.csv"))

  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res$value_obs)
  m=1
  t=2
  w=1

  
  
  for(m in 1:length(missingness)){
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
    
      ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
      ix_trait[is.na(ix_trait)] <- FALSE
      print(sum(ix_trait))
      
      bxpl_sp=cbind(bxpl_sp,
                    res$Sil_spec_pred[ix_trait&!duplicated(res$Species[ix_trait])]-
                      res$Sil_spec_obs[ix_trait&!duplicated(res$Species[ix_trait])])
      bxpl_gen=cbind(bxpl_gen,
                 res$Sil_gen_pred[ix_trait&!duplicated(res$Genus[ix_trait])]-
                   res$Sil_gen_obs[ix_trait&!duplicated(res$Genus[ix_trait])])
      bxpl_fam=cbind(bxpl_fam,
                 res$Sil_fam_pred[ix_trait&!duplicated(res$Family[ix_trait])]-
                   res$Sil_fam_obs[ix_trait&!duplicated(res$Family[ix_trait])])
      bxpl_clad=cbind(bxpl_clad,
                 res$Sil_clad_pred[ix_trait&!duplicated(res$Clade[ix_trait])] -
                 res$Sil_clad_obs[ix_trait&!duplicated(res$Clade[ix_trait])])
      bxpl_GF=cbind(bxpl_GF,
                res$Sil_GF_pred[ix_trait&!duplicated(res$GF[ix_trait])]-
                res$Sil_GF_obs[ix_trait&!duplicated(res$GF[ix_trait])])
      bxpl_PFT=cbind(bxpl_PFT,
                 res$Sil_PFT_pred[ix_trait&!duplicated(res$PFT[ix_trait])]-
                 res$Sil_PFT_obs[ix_trait&!duplicated(res$PFT[ix_trait])])
      
      bxpl_sp_op=cbind(bxpl_sp_op,
                       res$Sil_spec_obs[ix_trait&!duplicated(res$Species[ix_trait])],
                             res$Sil_spec_pred[ix_trait&!duplicated(res$Species[ix_trait])])
      bxpl_gen_op=cbind(bxpl_gen_op,
                        res$Sil_gen_obs[ix_trait&!duplicated(res$Genus[ix_trait])],
                              res$Sil_gen_pred[ix_trait&!duplicated(res$Genus[ix_trait])])
      bxpl_fam_op=cbind(bxpl_fam_op,
                        res$Sil_fam_obs[ix_trait&!duplicated(res$Family[ix_trait])],
                              res$Sil_fam_pred[ix_trait&!duplicated(res$Family[ix_trait])])
      bxpl_clad_op=cbind(bxpl_clad_op,
                        res$Sil_clad_obs[ix_trait&!duplicated(res$Clade[ix_trait])],
                        res$Sil_clad_pred[ix_trait&!duplicated(res$Clade[ix_trait])])
      bxpl_GF_op=cbind(bxpl_GF_op,
                      res$Sil_GF_obs[ix_trait&!duplicated(res$GF[ix_trait])],
                      res$Sil_GF_pred[ix_trait&!duplicated(res$GF[ix_trait])])
      bxpl_PFT_op=cbind(bxpl_PFT_op,
                        res$Sil_PFT_obs[ix_trait&!duplicated(res$PFT[ix_trait])],
                        res$Sil_PFT_pred[ix_trait&!duplicated(res$PFT[ix_trait])])
}
      
    colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
    bxpl_sp <- bxpl_sp[,colSums(!is.na(bxpl_sp))!=0]
    bxpl_gen <- bxpl_gen[,colSums(!is.na(bxpl_gen))!=0]
    bxpl_fam <- bxpl_fam[,colSums(!is.na(bxpl_fam))!=0]
    bxpl_clad <- bxpl_clad[,colSums(!is.na(bxpl_clad))!=0]
    bxpl_GF <- bxpl_GF[,colSums(!is.na(bxpl_GF))!=0]
    bxpl_PFT <- bxpl_PFT[,colSums(!is.na(bxpl_PFT))!=0]
    bxpl_sp_op <- bxpl_sp_op[,colSums(!is.na(bxpl_sp_op))!=0]
    bxpl_gen_op <- bxpl_gen_op[,colSums(!is.na(bxpl_gen_op))!=0]
    bxpl_fam_op <- bxpl_fam_op[,colSums(!is.na(bxpl_fam_op))!=0]
    bxpl_clad_op <- bxpl_clad_op[,colSums(!is.na(bxpl_clad_op))!=0]
    bxpl_GF_op <- bxpl_GF_op[,colSums(!is.na(bxpl_GF_op))!=0]
    bxpl_PFT_op <- bxpl_PFT_op[,colSums(!is.na(bxpl_PFT_op))!=0]
    
    colnames(bxpl_sp_op) <- rep(c("obs","pred"),1)
    colnames(bxpl_gen_op) <- rep(c("obs","pred"),1)
    colnames(bxpl_fam_op) <- rep(c("obs","pred"),1)
    colnames(bxpl_clad_op) <- rep(c("obs","pred"),1)
    colnames(bxpl_GF_op) <- rep(c("obs","pred"),1)
    colnames(bxpl_PFT_op) <- rep(c("obs","pred"),1)
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Sil",missingness[m],"2.pdf")),width=10,height=8)
    {
      par(mfrow=c(1,1))
      par(mfrow=c(1,7),mar=c(10,0,0,1))
      
      plot(1:10,col="white",ylim=c(-1,1),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -7,cex.axis=2)
      axis(side = 2,at = seq(-1,1,by=1),line = -5,labels = c("","Silhouette",""),tick = FALSE,cex.axis=4)
      
      boxplot(bxpl_sp_op,ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=4,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_sp_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot(bxpl_sp_op,ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,add=TRUE,xaxt="n",yaxt="n",
              las=2,lwd=2,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
#      abline(h=median((bxpl_sp),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "Species",cex=3)
      
      boxplot((bxpl_gen_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_gen_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot((bxpl_gen_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,xaxt="n",yaxt="n",add=TRUE,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
#            abline(h=median((bxpl_gen),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "Genus",cex=3)
      
      boxplot((bxpl_fam_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_fam_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot((bxpl_fam_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
      #      abline(h=median((bxpl_fam),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "Family",cex=3)
      
      boxplot((bxpl_clad_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_clad_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot((bxpl_clad_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,yaxt="n",xaxt="n",
              las=2,cex.axis=3,cex.lab=2,lwd=2,add=TRUE,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
      #      abline(h=median((bxpl_clad),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "Clade",cex=3)
      
      boxplot((bxpl_GF_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_GF_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot((bxpl_GF_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,xaxt="n",yaxt="n",
              las=2,cex.axis=3,cex.lab=2,lwd=2,add=TRUE,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
      #      abline(h=median((bxpl_GF),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "GF",cex=3)
      
      boxplot((bxpl_PFT_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      ynow=median(bxpl_PFT_op[,1],na.rm = TRUE)
      rect(xleft = 0,ybottom = ynow,xright = 3,ytop =1,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 3,ytop =ynow,col = colz[5],border = FALSE)
      boxplot((bxpl_PFT_op),ylim=c(-1,1.1),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=3,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
      abline(h=ynow,lwd=3,col="black")
      #      abline(h=median((bxpl_PFT),na.rm = TRUE),lwd=3)
      text(x=1.5,y=1.1,labels = "PFT",cex=3)
      
      }
    dev.off()
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Sil",missingness[m],".pdf")),width=10,height=8)
    {
      par(mfrow=c(1,1))
      par(mfrow=c(1,8),mar=c(10,0,0,1))
      
      plot(1:10,col="white",ylim=c(-1,1.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -7,cex.axis=2)
      axis(side = 2,at = seq(-1,1,by=1),line = -5.5,labels = c("","pred - obs",""),tick = FALSE,cex.axis=3)
      axis(side = 2,at = seq(-1,1,by=1),line = -3.25,labels = c("","Silhouette",""),tick = FALSE,cex.axis=3)
      
      boxplot((bxpl_sp),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_sp),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_sp),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_sp),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "Species",cex=3)
      
      boxplot((bxpl_gen),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_gen),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_gen),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_gen),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "Genus",cex=3)
      
      boxplot((bxpl_fam),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_fam),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_fam),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_fam),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "Family",cex=3)
      
      boxplot((bxpl_clad),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_clad),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_clad),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_clad),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "Clade",cex=3)
      
      boxplot((bxpl_GF),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_GF),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_GF),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_GF),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "GF",cex=3)
      
      boxplot((bxpl_PFT),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 2,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 2,ytop =0,col = colz[5],border = FALSE)
      boxplot((bxpl_PFT),ylim=c(-1,1.5),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_PFT),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_PFT),na.rm = TRUE),lwd=3)
      text(x=1,y=1.4,labels = "PFT",cex=3)
      
      
      plot(1:10,col="white",ylim=c(-1,1.5),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = 10,ytop =1.3,col = colz[4],border = FALSE)
      rect(xleft = 0,ybottom = -1,xright = 10,ytop =0,col = colz[5],border = FALSE)
      text(3,.75,labels="more",col="black",srt=90,cex=3)
      text(6,.75,labels="clustering",col="black",srt=90,cex=3)
      text(3,-.5,labels="less",col="black",srt=90,cex=3)
      text(6,-.5,labels="clustering",col="black",srt=90,cex=3)
    }
    dev.off()
    
    
    
  }