
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

  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  m=1
  t=2
  w=1

  
  
  for(m in 1:length(missingness)){
    par(mfrow=c(1,1),mar=c(10,6,2,2))
    t=1
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    sum(ix_trait)
    
    bxpl_sp=matrix(NA,nrow=sum(ix_trait[ix_trait]),ncol=1)
    bxpl_gen=matrix(NA,nrow=sum(ix_trait),ncol=1)#rep(NA,sum(ix_trait))
    bxpl_fam=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_clad=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_GF=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_PFT=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_sp_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_gen_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_fam_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_clad_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_GF_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
    print(head(bxpl_sp))
    
    t=1
    for(t in 1:length(trait_names)){
      ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
      ix_trait[is.na(ix_trait)] <- FALSE
      print(sum(ix_trait))
      
      plot(res$dist_spec_obs[ix_trait],res$dist_spec_pred[ix_trait])
      plot(res$dist_spec_obs_zlog[ix_trait],res$dist_spec_pred_zlog[ix_trait])
      plot(res$dist_gen_obs[ix_trait],res$dist_gen_pred[ix_trait])
      plot(res$dist_gen_obs[ix_trait],res$dist_gen_pred[ix_trait])
      mean_now=res$mean_gen_obs[ix_trait]

          bxpl_sp=cbind(bxpl_sp,
                         (res$dist_spec_pred_zlog[ix_trait])-
                           (res$dist_spec_obs_zlog[ix_trait]))
          bxpl_sp_op=cbind(bxpl_sp_op,
                            res$dist_spec_obs_zlog[ix_trait],
                           res$dist_spec_pred_zlog[ix_trait])
          

          bxpl_gen=cbind(bxpl_gen,
                     (res$dist_gen_pred_zlog[ix_trait])-
                       (res$dist_gen_obs_zlog[ix_trait]))
          bxpl_gen_op=cbind(bxpl_gen_op,
                        res$dist_gen_obs_zlog[ix_trait],
                        res$dist_gen_pred_zlog[ix_trait])
      
      bxpl_fam=cbind(bxpl_fam,
                     (res$dist_fam_pred_zlog[ix_trait])-
                       (res$dist_fam_obs_zlog[ix_trait]))
      bxpl_fam_op=cbind(bxpl_fam_op,
                        (res$dist_fam_pred_zlog[ix_trait]),
                        (res$dist_fam_obs_zlog[ix_trait]))
      
      bxpl_clad=cbind(bxpl_clad,
                     (res$dist_clad_pred_zlog[ix_trait])-
                       (res$dist_clad_obs_zlog[ix_trait]))
      bxpl_clad_op=cbind(bxpl_clad_op,
                        res$dist_clad_obs_zlog[ix_trait],
                        res$dist_clad_pred_zlog[ix_trait])
      
      bxpl_GF=cbind(bxpl_GF,
                     (res$dist_GF_pred_zlog[ix_trait])-
                       (res$dist_GF_obs_zlog[ix_trait]))
      bxpl_GF_op=cbind(bxpl_GF_op,
                        res$dist_GF_obs_zlog[ix_trait],
                       res$dist_GF_pred_zlog[ix_trait])
      
      bxpl_PFT=cbind(bxpl_PFT,
                     (res$dist_PFT_pred_zlog[ix_trait])-
                       (res$dist_PFT_obs_zlog[ix_trait]))
      bxpl_PFT_op=cbind(bxpl_PFT_op,
                        res$dist_PFT_obs_zlog[ix_trait],res$dist_PFT_pred_zlog[ix_trait])
  
    }
    
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
    colnames(bxpl_sp) <- trait_names
    colnames(bxpl_gen) <- trait_names
    colnames(bxpl_fam) <- trait_names
    colnames(bxpl_clad) <- trait_names
    colnames(bxpl_GF) <- trait_names
    colnames(bxpl_PFT) <- trait_names
    colnames(bxpl_sp_op) <- rep(c("obs","pred"),6)
    colnames(bxpl_gen_op) <- rep(c("obs","pred"),6)
    colnames(bxpl_fam_op) <- rep(c("obs","pred"),6)
    colnames(bxpl_clad_op) <- rep(c("obs","pred"),6)
    colnames(bxpl_GF_op) <- rep(c("obs","pred"),6)
    colnames(bxpl_PFT_op) <- rep(c("obs","pred"),6)
    
    colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
           "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Dist",missingness[m],".pdf")),width=20,height=8)
    {
      par(mar=c(4,0,4,0))
      layout(matrix(c(1:7,1:7,1:7,1:7,8,rep(9,6)), 
                    nrow = 5, ncol = 7, byrow = TRUE))
      
      ytop=1.25
      plot(1:10,col="white",ylim=c(0,ytop),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -17,cex.axis=2)
      axis(side = 2,at = .6,line = -10,labels = "Distance to",tick = FALSE,cex.axis=4)
      axis(side = 2,at = .6,line = -14,labels = "co-cluster individuals",tick = FALSE,cex.axis=4)
      
      names_now=c("Species","Genus","Family","Clades","GF","PFT")
      for(t in 1:6){
        c1=6
        c2=10
        if(t==1){bxpl <- bxpl_sp_op}
        if(t==2){bxpl <- bxpl_gen_op}
        if(t==3){bxpl <- bxpl_fam_op}
        if(t==4){bxpl <- bxpl_clad_op}
        if(t==5){bxpl <- bxpl_GF_op}
        if(t==6){bxpl <- bxpl_PFT_op}
        
        boxplot((bxpl),ylim=c(0,ytop),ylab="CV pred - obs",col=colz,cex.main=5,
                las=2,cex.axis=2,cex.lab=1,lwd=2,yaxt="n",frame=FALSE,main=names_now[t])
        yBot_now=apply(bxpl[,seq(from=1,to = 12,by = 2)],2,FUN = median,na.rm = TRUE)
        i=1
        rect(xleft = .5,ybottom = 0,xright = 2.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = .5,ybottom = yBot_now[i],xright = 2.5,ytop = ytop,col=colz[c2],border = NA)
        i=2
        rect(xleft = 2.5,ybottom = 0,xright = 4.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = 2.5,ybottom = yBot_now[i],xright = 4.5,ytop = ytop,col=colz[c2],border = NA)
        i=3
        rect(xleft = 4.5,ybottom = 0,xright = 6.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = 4.5,ybottom = yBot_now[i],xright = 6.5,ytop = ytop,col=colz[c2],border = NA)
        i=4
        rect(xleft = 6.5,ybottom = 0,xright = 8.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = 6.5,ybottom = yBot_now[i],xright = 8.5,ytop = ytop,col=colz[c2],border = NA)
        i=5
        rect(xleft = 8.5,ybottom = 0,xright = 10.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = 8.5,ybottom = yBot_now[i],xright = 10.5,ytop = ytop,col=colz[c2],border = NA)
        i=6
        rect(xleft = 10.5,ybottom = 0,xright = 12.5,ytop = yBot_now[i],col=colz[c1],border = NA)
        rect(xleft = 10.5,ybottom = yBot_now[i],xright = 12.5,ytop = ytop,col=colz[c2],border = NA)
        abline(v=seq(.5,12.5,2),col="white",lty=1,lwd=3)
       # axis(side = 1,line = 3, at = seq(from=1.5,to = 12,by = 2),las=2,labels=trait_names,tick = FALSE,cex.axis=2)
        boxplot(bxpl,ylim=c(0,ytop),ylab="CV pred - obs",col=colz,xaxt="n",xlab="",
                las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE,add=TRUE)
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
    
    colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Dist",missingness[m],"2.pdf")),width=10,height=8)
    {
      par(mfrow=c(1,1))
      par(mfrow=c(1,8),mar=c(10,0,0,1))
      
      plot(1:10,col="white",ylim=c(-.75,.75),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -7,cex.axis=2)
      axis(side = 2,at = seq(-1,.5,by=.5),line = -5.5,labels = c("","","pred - obs",""),tick = FALSE,cex.axis=3)
      axis(side = 2,at = seq(-1,.5,by=.5),line = -3.25,labels = c("","","Distance",""),tick = FALSE,cex.axis=3)
      
      boxplot((bxpl_sp),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=1,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_sp),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",xaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_sp),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_sp),na.rm = TRUE),lwd=3)
      text(x=3.5,y=0.55,labels = "Species",cex=3)
      
      boxplot((bxpl_gen),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_gen),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_gen),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_gen),na.rm = TRUE),lwd=3)
      text(x=3.5,y=0.55,labels = "Genus",cex=3)
      
      boxplot((bxpl_fam),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_fam),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_fam),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_fam),na.rm = TRUE),lwd=3)
      text(x=3.5,y=0.55,labels = "Family",cex=3)
      
      boxplot((bxpl_clad),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_clad),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_clad),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_clad),na.rm = TRUE),lwd=3)
      text(x=3.5,y=0.55,labels = "Clade",cex=3)
      
      boxplot((bxpl_GF),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_GF),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
      abline(h=median((bxpl_GF),na.rm = TRUE),lwd=4,col="white")
      abline(h=median((bxpl_GF),na.rm = TRUE),lwd=3)
      text(x=3.5,y=0.55,labels = "GF",cex=3)
      
      boxplot((bxpl_PFT),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE)
      rect(xleft = 0,ybottom = 0,xright = ncol(bxpl_clad)+2,ytop =.5,col = colz[5],border = FALSE)
      rect(xleft = 0,ybottom = -5,xright = ncol(bxpl_clad)+2,ytop =0,col = colz[4],border = FALSE)
      boxplot((bxpl_PFT),ylim=c(-.75,.75),ylab="CV pred - obs",col=colz,
              las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",add=TRUE,frame=FALSE)
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
    dev.off()
  }
  
  
  
  
  
  
  
  # NO TRANSFORMATION: 

for(m in 1:length(missingness)){
  par(mfrow=c(1,1),mar=c(10,6,2,2))
  t=1
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  sum(ix_trait)
  
  bxpl_sp=matrix(NA,nrow=sum(ix_trait[ix_trait]),ncol=1)
  bxpl_gen=matrix(NA,nrow=sum(ix_trait),ncol=1)#rep(NA,sum(ix_trait))
  bxpl_fam=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_clad=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_GF=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_PFT=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_sp_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_gen_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_fam_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_clad_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_GF_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait),ncol=1)
  print(head(bxpl_sp))
  
  t=1
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    print(sum(ix_trait))
    
    mean_now=res$mean_gen_obs[ix_trait]
    
    bxpl_sp=cbind(bxpl_sp,
                  (res$dist_spec_pred[ix_trait])-
                    (res$dist_spec_obs[ix_trait]))
    bxpl_sp_op=cbind(bxpl_sp_op,
                     res$dist_spec_pred[ix_trait],
                     res$dist_spec_obs[ix_trait])
    
    
    bxpl_gen=cbind(bxpl_gen,
                   (res$dist_gen_pred[ix_trait])-
                     (res$dist_gen_obs[ix_trait]))
    bxpl_gen_op=cbind(bxpl_gen_op,
                      res$dist_gen_pred[ix_trait],
                      res$dist_gen_obs[ix_trait])
    
    bxpl_fam=cbind(bxpl_fam,
                   (res$dist_fam_pred[ix_trait])-
                     (res$dist_fam_obs[ix_trait]))
    bxpl_fam_op=cbind(bxpl_fam_op,
                      (res$dist_fam_pred[ix_trait]),
                      (res$dist_fam_obs[ix_trait]))
    
    bxpl_clad=cbind(bxpl_clad,
                    (res$dist_clad_pred[ix_trait])-
                      (res$dist_clad_obs[ix_trait]))
    bxpl_clad_op=cbind(bxpl_clad_op,
                       res$dist_clad_pred[ix_trait],
                       res$dist_clad_obs[ix_trait])
    
    bxpl_GF=cbind(bxpl_GF,
                  (res$dist_GF_pred[ix_trait])-
                    (res$dist_GF_obs[ix_trait]))
    bxpl_GF_op=cbind(bxpl_GF_op,
                     res$dist_GF_pred[ix_trait],
                     res$dist_GF_obs[ix_trait])
    
    bxpl_PFT=cbind(bxpl_PFT,
                   (res$dist_PFT_pred[ix_trait])-
                     (res$dist_PFT_obs[ix_trait]))
    bxpl_PFT_op=cbind(bxpl_PFT_op,
                      res$dist_PFT_pred[ix_trait],res$dist_PFT_obs[ix_trait])
    
  }
  
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
  colnames(bxpl_sp) <- trait_names
  colnames(bxpl_gen) <- trait_names
  colnames(bxpl_fam) <- trait_names
  colnames(bxpl_clad) <- trait_names
  colnames(bxpl_GF) <- trait_names
  colnames(bxpl_PFT) <- trait_names
  colnames(bxpl_sp_op) <- rep(c("pred","obs"),6)
  colnames(bxpl_gen_op) <- rep(c("pred","obs"),6)
  colnames(bxpl_fam_op) <- rep(c("pred","obs"),6)
  colnames(bxpl_clad_op) <- rep(c("pred","obs"),6)
  colnames(bxpl_GF_op) <- rep(c("pred","obs"),6)
  colnames(bxpl_PFT_op) <- rep(c("pred","obs"),6)
  
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  par(mfrow=c(1,1))
  pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Dist",missingness[m],"notrans.pdf")),width=20,height=8)
  {
    par(mar=c(4,5,4,0))
    layout(matrix(c(1:6,1:6,1:6,1:6,7,rep(8,5)), 
                  nrow = 5, ncol = 6, byrow = TRUE))
    
    limmin <- c(0,0,-.01,0,0,0)
#    limmax <- c(50,40,.6,40,5,5)
    limmax <- c(20,30,.4,20,2.5,2)
    names_now=c("Species","Genus","Family","Clades","GF","PFT")
    units_names <- c("","Meters","m3","","","")
    t=1
    for(t in 1:6){
      c1=6
      c2=10
      t2=seq(1,12,2)[t]
      t2=c(t2+1,t2)
      bxpl <- cbind(bxpl_sp_op[,t2],bxpl_gen_op[,t2],bxpl_fam_op[,t2],bxpl_clad_op[,t2],bxpl_GF_op[,t2],bxpl_PFT_op[,t2])
      ytop=limmax[t]  
      boxplot((bxpl),ylim=c(limmin[t],limmax[t]),ylab=units_names[t],col=colz,cex.main=4,cex.lab=3,
              las=2,cex.axis=2,lwd=2,frame=FALSE,main=trait_names[t])
      yBot_now=apply(bxpl[,seq(from=1,to = 12,by = 2)],2,FUN = median,na.rm = TRUE)
      i=1
      rect(xleft = .5,ybottom = 0,xright = 2.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = .5,ybottom = yBot_now[i],xright = 2.5,ytop = ytop,col=colz[c2],border = NA)
      i=2
      rect(xleft = 2.5,ybottom = 0,xright = 4.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = 2.5,ybottom = yBot_now[i],xright = 4.5,ytop = ytop,col=colz[c2],border = NA)
      i=3
      rect(xleft = 4.5,ybottom = 0,xright = 6.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = 4.5,ybottom = yBot_now[i],xright = 6.5,ytop = ytop,col=colz[c2],border = NA)
      i=4
      rect(xleft = 6.5,ybottom = 0,xright = 8.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = 6.5,ybottom = yBot_now[i],xright = 8.5,ytop = ytop,col=colz[c2],border = NA)
      i=5
      rect(xleft = 8.5,ybottom = 0,xright = 10.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = 8.5,ybottom = yBot_now[i],xright = 10.5,ytop = ytop,col=colz[c2],border = NA)
      i=6
      rect(xleft = 10.5,ybottom = 0,xright = 12.5,ytop = yBot_now[i],col=colz[c1],border = NA)
      rect(xleft = 10.5,ybottom = yBot_now[i],xright = 12.5,ytop = ytop,col=colz[c2],border = NA)
      abline(v=seq(.5,12.5,2),col="white",lty=1,lwd=3)
      # axis(side = 1,line = 3, at = seq(from=1.5,to = 12,by = 2),las=2,labels=trait_names,tick = FALSE,cex.axis=2)
      boxplot((bxpl),ylim=c(limmin[t],limmax[t]),ylab="CV pred - obs",col=colz,cex.main=5,
              las=2,cex.axis=2,cex.lab=1,lwd=2,frame=FALSE,main=trait_names[t],add=TRUE)
    }  
    
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    rect(xleft = 2,ybottom = 0,xright = 10,ytop = 12,col = "lightgray",border = FALSE)
    legend(2.5, 14, names_now[1], col = colz[1], pch = 15, cex = 3,bty = "n")
    legend(3.8,14,names_now[2], col = colz[3], pch = 15, cex = 3,bty = "n")
    legend(5, 14, names_now[3], col = colz[5], pch = 15, cex = 3,bty = "n")
    legend(6.3, 14, names_now[4], col = colz[7], pch = 15, cex = 3,bty = "n")
    legend(7.5, 14, names_now[5], col = colz[9], pch = 15, cex = 3,bty = "n")
    legend(8.5, 14, names_now[6], col = colz[11], pch = 15, cex = 3,bty = "n")
    
  }
  dev.off()
  
  colz=c("#b2182b","#d6604d","#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac")
  
}

