
# get per observed point value: 
# - error (pred-obs)
# - average distance to all other points observed 
# - Silhouette index for this group (or calculate somewhere else...?)  
# - ORIGINAL Silhouette index for this group (or calculate somewhere else...?)  
# - DEVIATION Silhouette index for this group (or calculate somewhere else...?)  
# - number of values available
# - original number of values available

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
list.files(file.path(origin,"_2021","script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
TD_choices = out$TD_choices
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

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

res <- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))

  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res_now$value_obs)
  m=1
  t=2
  w=1
  
  
  for(m in 1:length(missingness)){
  pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Distance",missingness[m],".pdf")),width=6,height=8)
  par(mfrow=c(1,1),mar=c(10,6,2,2))
  sum(ix_trait)
  
  bxpl_sp=NA
  bxpl_gen=NA#rep(NA,sum(ix_trait))
  bxpl_fam=NA#rep(NA,sum(ix_trait))
  bxpl_clad=NA#rep(NA,sum(ix_trait))
  bxpl_GF=NA#rep(NA,sum(ix_trait))
  bxpl_PFT=NA#rep(NA,sum(ix_trait))
  bxpl_sp_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  bxpl_gen_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  bxpl_fam_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  bxpl_clad_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  bxpl_GF_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait),ncol=2)
  
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    
    bxpl_sp=c(bxpl_sp,-res$mean_spec[ix_trait]-res$value_obs[ix_trait]+
                                      res$mean_spec[ix_trait]-res$value_pred[ix_trait])
    bxpl_gen=c(bxpl_gen,-res$mean_gen[ix_trait]-res$value_obs[ix_trait]+
                                        res$mean_gen[ix_trait]-res$value_pred[ix_trait])
    bxpl_fam=c(bxpl_fam,res$mean_fam[ix_trait]-res$value_obs[ix_trait]+
                                        res$mean_fam[ix_trait]-res$value_pred[ix_trait])
    bxpl_clad=c(bxpl_clad,-res$mean_clad[ix_trait]-res$value_obs[ix_trait]+
                                          res$mean_clad[ix_trait]-res$value_pred[ix_trait])
    bxpl_GF=c(bxpl_GF,-res$mean_GF[ix_trait]-res$value_obs[ix_trait]+
                                      res$mean_GF[ix_trait]-res$value_pred[ix_trait])
    bxpl_PFT=c(bxpl_PFT,-res$mean_PFT[ix_trait]-res$value_obs[ix_trait]+
                                        res$mean_PFT[ix_trait]-res$value_pred[ix_trait])
    
    bxpl_sp_op=rbind(bxpl_sp_op,cbind(res$mean_spec[ix_trait]-res$value_obs[ix_trait],
                                      res$mean_spec[ix_trait]-res$value_pred[ix_trait]))
    bxpl_gen_op=rbind(bxpl_gen_op,cbind(res$mean_gen[ix_trait]-res$value_obs[ix_trait],
                                        res$mean_gen[ix_trait]-res$value_pred[ix_trait]))
    bxpl_fam_op=rbind(bxpl_fam_op,cbind(res$mean_fam[ix_trait]-res$value_obs[ix_trait],
                                        res$mean_fam[ix_trait]-res$value_pred[ix_trait]))
    bxpl_clad_op=rbind(bxpl_clad_op,cbind(res$mean_clad[ix_trait]-res$value_obs[ix_trait],
                                          res$mean_clad[ix_trait]-res$value_pred[ix_trait]))
    bxpl_GF_op=rbind(bxpl_GF_op,cbind(res$mean_GF[ix_trait]-res$value_obs[ix_trait],
                                      res$mean_GF[ix_trait]-res$value_pred[ix_trait]))
    bxpl_PFT_op=rbind(bxpl_PFT_op,cbind(res$mean_PFT[ix_trait]-res$value_obs[ix_trait],
                                        res$mean_PFT[ix_trait]-res$value_pred[ix_trait]))
    
#    boxplot(res$dist_spec_obs[ix_trait],res$dist_spec_pred[ix_trait],
#            res$dist_gen_obs[ix_trait],res$dist_gen_pred[ix_trait],
#            res$dist_fam_obs[ix_trait],res$dist_fam_pred[ix_trait],
#            res$dist_clad_obs[ix_trait],res$dist_clad_pred[ix_trait],
#            col=colz)
#    abline(v=seq(.5,to=50,by=6),col="gray",lwd=3)
  }   
length(bxpl_fam)
  bxpl <- cbind(bxpl_sp,bxpl_gen,bxpl_fam,bxpl_clad,bxpl_GF,bxpl_PFT,rep(NA,length(bxpl_PFT)))
  colnames(bxpl) <- c("Species","Genus","Family","Clade","Growth Form","PFT","")
  boxplot(bxpl,ylim=c(-5,5),ylab="Distance to pred_mean - obs_mean",
          col=colz,
          las=2,cex.axis=1.5,cex.lab=2,lwd=2)
  abline(h=0,col="gray",lty=2,lwd=2)
  rect(xleft = 6.7,ybottom = 0,xright = 10,ytop = 10,col = colz[6],border = FALSE)
  rect(xleft = 6.7,ybottom = -10,xright = 10,ytop = 0,col = colz[1],border = FALSE)
  text(7.2,2.5,labels="less clustering",col="white",srt=90,cex=2)
  text(7.2,-2.5,labels="more clustering",col="white",srt=90,cex=2)
  
  bxpl <- cbind(bxpl_sp_op,bxpl_gen_op,bxpl_fam_op,bxpl_clad_op,bxpl_GF_op,bxpl_PFT_op)
  colnames(bxpl) <- c("Species_pred","Species_pred","Genus_obs","Genus_pred",
                      "Family_obs","Family_pred","Clade_obs","Clade_pred","GF_obs","GF_pred",
                      "PFT_obs","PFT_pred")
  boxplot(abs(bxpl),ylim=c(0,3),ylab="Distance to cluster mean",
          col=c(colz[1],colz[1],colz[2],colz[2],colz[3],colz[3],colz[4],colz[4],colz[5],colz[5],colz[6],colz[6]),
          las=2,cex.axis=2.5,cex.lab=2)
  abline(h=0)
  abline(v=seq(.5,to=50,by=2),col="gray",lwd=3)
  dev.off()
  }
  
  pdf(file=file.path(origin,"_2021","figures","Figure_2","Figure_2_DistBPT.pdf"),width=16,height=4)
  par(mfrow=c(1,1),mar=c(7,4,2,2))
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  bxpl_sp=rep(NA,sum(ix_trait))
  bxpl_gen=rep(NA,sum(ix_trait))
  bxpl_fam=rep(NA,sum(ix_trait))
  bxpl_clad=rep(NA,sum(ix_trait))
  bxpl_GF=rep(NA,sum(ix_trait))
  bxpl_PFT=rep(NA,sum(ix_trait))
  
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    ix_trait[is.na(ix_trait)] <- FALSE
    
    bxpl_sp=cbind(bxpl_sp,res$dist_spec_pred[ix_trait]- res$dist_spec_obs[ix_trait])
    bxpl_gen=cbind(bxpl_gen,res$dist_gen_pred[ix_trait]- res$dist_gen_obs[ix_trait])
    bxpl_fam=cbind(bxpl_fam,res$dist_fam_pred[ix_trait]- res$dist_fam_obs[ix_trait])
    bxpl_clad=cbind(bxpl_clad,res$dist_clad_pred[ix_trait]- res$dist_clad_obs[ix_trait])
    bxpl_GF=cbind(bxpl_GF,res$dist_GF_pred[ix_trait]- res$dist_GF_obs[ix_trait])
    bxpl_PFT=cbind(bxpl_PFT,res$dist_PFT_pred[ix_trait]- res$dist_PFT_obs[ix_trait])
    
  }   
  colnames(bxpl_sp) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  colnames(bxpl_gen) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  colnames(bxpl_fam) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  colnames(bxpl_clad) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  colnames(bxpl_GF) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  colnames(bxpl_PFT) <- c("","Species","Genus","Family","Clade","Growth Form","PFT")
  
  boxplot(cbind(bxpl_sp,bxpl_gen,bxpl_fam,bxpl_clad,bxpl_GF,bxpl_PFT),ylim=c(-1,.5),
          col=c(rep(colz[1],7),rep(colz[2],7),rep(colz[3],7),rep(colz[4],7),rep(colz[5],7),rep(colz[6],7)),
          las=2,cex.axis=1.5)
  abline(h=0)
 # abline(v=seq(.5,to=50,by=6),col="gray",lwd=3)
  dev.off()
  
  
  bxpl_sp   <- bxpl_sp[,colSums(!is.na(bxpl_sp))!=0]
  bxpl_gen  <- bxpl_gen[,colSums(!is.na(bxpl_gen))!=0]
  bxpl_fam  <- bxpl_fam[,colSums(!is.na(bxpl_fam))!=0]
  bxpl_clad <- bxpl_clad[,colSums(!is.na(bxpl_clad))!=0]
  bxpl_GF   <- bxpl_GF[,colSums(!is.na(bxpl_GF))!=0]
  bxpl_PFT  <- bxpl_PFT[,colSums(!is.na(bxpl_PFT))!=0]
  
  
  pdf(file=file.path(origin,"_2021","figures","Figure_2","Figure_2_DistBPT.pdf"))
  par(mfrow=c(1,1),mar=c(4,5,2,2))
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  for(t in 1:length(trait_names)){
    ix_trait=res$trait==trait_names[t]
    ix_trait[is.na(ix_trait)] <- FALSE
    
    bxpl=cbind(res$dist_spec_obs[ix_trait]- res$dist_spec_pred[ix_trait],
               res$dist_gen_obs[ix_trait]- res$dist_gen_pred[ix_trait],
               res$dist_fam_obs[ix_trait]- res$dist_fam_pred[ix_trait],
               res$dist_clad_obs[ix_trait]- res$dist_clad_pred[ix_trait],
               res$dist_GF_obs[ix_trait]- res$dist_GF_pred[ix_trait],
               res$dist_PFT_obs[ix_trait]- res$dist_PFT_pred[ix_trait])
    colnames(bxpl) <- c("Species","Genus","Family","Clades","GF","PFT")
    boxplot(bxpl,ylim=c(-.5,1.5),col=colz,las=2,main=trait_names[t])
    abline(h=0)
  }   
    dev.off()
  
  pdf(file=file.path(origin,"_2021","figures","Pies.pdf"))
  par(mfrow=c(1,3))
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    res_lm <- res_now[,grep(colnames(res_now),pattern = "_lm")]
    res_now <- res_now[,-grep(colnames(res_now),pattern = "_lm")]
    res_lm <- res_lm[,grep(colnames(res_lm),pattern = paste0(trait_names[t],"_from"))]
    res_now <- cbind(res_now,res_lm)
    error_now = res_now$value_pred-res_now$value_obs
    target_value=.25
    part_high = sum((res_now$value_pred-res_now$value_obs>target_value),na.rm = TRUE)
    part_low=sum((res_now$value_pred-res_now$value_obs<(target_value*-1)),na.rm = TRUE)
    part_ok=sum((res_now$value_pred-res_now$value_obs>(target_value*-1)&res_now$value_pred-res_now$value_obs<target_value),na.rm = TRUE)
    
    #-----------------------------------------------
    pie_input  <- c(part_high,part_low,part_ok)
    pie(pie_input,labels = c("+","-","="),main=trait_names[t])
    #-----------------------------------------------   
    
    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    res_now <- res_now[rowSums(!is.na(res_now))!=0,]
    
    weights=c("mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",
              colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
    spider_1=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_1) <- c("",weights)
    spider_obs=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_obs) <- c("",weights)
    spider_pred=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_pred) <- c("",weights)
    w=1
    for(w in 1:length(weights)){
      # deviation towards weight y/n
      obs = res_now$value_obs-res_now[,colnames(res_now)%in%weights[w]]      
      pred = res_now$value_pred-res_now[,colnames(res_now)%in%weights[w]]
      
      spider_obs[,colnames(spider_obs)==weights[w]] <- obs
      spider_pred[,colnames(spider_pred)==weights[w]] <- pred
      
      res_yn <- cbind(abs(pred),abs(obs))
      spider_1[,colnames(spider_1)==weights[w]] <- res_yn[,1]-res_yn[,2]
    }
    spider_obs <- as.data.frame(spider_obs[,colSums(!is.na(spider_obs))!=0])
    spider_pred <- as.data.frame(spider_pred[,colSums(!is.na(spider_pred))!=0])
    spider_1 <- as.data.frame(spider_1[,colSums(!is.na(spider_1))!=0])
    spider_1$abs_error <- abs(res_now$value_pred-res_now$value_obs)*-1
    
    target_value=.25
    part_high = sum((spider_1$mean_spec>target_value),na.rm = TRUE)
    part_low=sum((spider_1$mean_spec<(target_value*-1)),na.rm = TRUE)
    part_ok=sum((spider_1$mean_spec>(target_value*-1)&spider_1$mean_spec<target_value),na.rm = TRUE)
    
    #-----------------------------------------------
    pie_input  <- c(part_high,part_low,part_ok)
    pie(pie_input,labels = c("+","-","="),main=paste0(trait_names[t]),"_taxonomy")
    #-----------------------------------------------   
    
    target_value=.25
    part_high = sum((spider_1$mean_GF>target_value),na.rm = TRUE)
    part_low=sum((spider_1$mean_GF<(target_value*-1)),na.rm = TRUE)
    part_ok=sum((spider_1$mean_GF>(target_value*-1)&spider_1$mean_GF<target_value),na.rm = TRUE)
    
    #-----------------------------------------------
    pie_input  <- c(part_high,part_low,part_ok)
    pie(pie_input,labels = c("+","-","="),main=paste0(trait_names[t],"_function"))
    #-----------------------------------------------   
    
       
  }
  dev.off()
  
  pdf(file=file.path(origin,"_2021","figures","Proper_spider.pdf"))
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    res_lm <- res_now[,grep(colnames(res_now),pattern = "_lm")]
    res_now <- res_now[,-grep(colnames(res_now),pattern = "_lm")]
    res_lm <- res_lm[,grep(colnames(res_lm),pattern = paste0(trait_names[t],"_from"))]
    res_now <- cbind(res_now,res_lm)
    colnames(res_now)
    
    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    res_now <- res_now[rowSums(!is.na(res_now))!=0,]
    
    weights=c("mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",
              colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
    spider_1=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_1) <- c("",weights)
    spider_obs=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_obs) <- c("",weights)
    spider_pred=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_pred) <- c("",weights)
    w=1
    for(w in 1:length(weights)){
      # deviation towards weight y/n
      obs = res_now$value_obs-res_now[,colnames(res_now)%in%weights[w]]      
      pred = res_now$value_pred-res_now[,colnames(res_now)%in%weights[w]]
      
      spider_obs[,colnames(spider_obs)==weights[w]] <- obs
      spider_pred[,colnames(spider_pred)==weights[w]] <- pred
      
      res_yn <- cbind(abs(pred),abs(obs))
      spider_1[,colnames(spider_1)==weights[w]] <- res_yn[,1]-res_yn[,2]
    }
    spider_obs <- as.data.frame(spider_obs[,colSums(!is.na(spider_obs))!=0])
    spider_pred <- as.data.frame(spider_pred[,colSums(!is.na(spider_pred))!=0])
    spider_1 <- as.data.frame(spider_1[,colSums(!is.na(spider_1))!=0])
    spider_1$abs_error <- abs(res_now$value_pred-res_now$value_obs)*-1
    
#    par(mfrow=c(1,1),mar=c(10,2,2,2))
#    boxplot(spider_1,ylim=c(-2,2),las=2)
#    abline(h=0)
#    boxplot(spider_1[res_now$dist_spec_obs>quantile(res_now$dist_spec_obs,probs = .9,na.rm = TRUE),],ylim=c(-2,2),las=2)
#    abline(h=0)
#    boxplot(spider_1[res_now$dist_spec_obs<quantile(res_now$dist_spec_obs,probs = .1,na.rm = TRUE),],ylim=c(-2,2),las=2)
#    abline(h=0)
    #install.packages("fmsb")  
    library(fmsb)
    par(mfrow=c(1,1),mar=c(2,2,2,2))
    
    summary(spider_obs)
    data1.0 <- spider_1
    data1 <- spider_1
    #    data1 <- data.frame(Taxonomy=data1.0[,which(apply(data1.0[,1:4],2,mean,na.rm=TRUE)==max(apply(data1.0[,1:4],2,mean,na.rm=TRUE)))],
#                        Function=data1.0$mean_GF,
#                        Trait=data1.0[,which(apply(data1.0[,7:ncol(data1.0)],2,mean,na.rm=TRUE)==max(apply(data1.0[,7:ncol(data1.0)],2,mean,na.rm=TRUE)))])
    info=res_now
#    hist(data1$Taxonomy)
#    hist(data2$Taxonomy)
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    
    # colours
    {    
      colz=rep(rgb(.9,.9,.9,alpha = .01),nrow(data1))
      
      ix=res_now$dist_spec_obs>quantile(res_now$dist_spec_obs,probs = .9,na.rm = TRUE)
      ix[is.na(ix)] <- FALSE
      data1 <- rbind(data1, apply(data1[ix,],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
      colz <-     c(colz,"orange")
      ix=res_now$dist_spec_obs<quantile(res_now$dist_spec_obs,probs = .1,na.rm = TRUE)
      ix[is.na(ix)] <- FALSE
      data1 <- rbind(data1, apply(data1[ix,],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
      colz <-     c(colz,"lightblue")

      ix=res_now$dist_spec_obs>quantile(res_now$dist_spec_obs,probs = .9,na.rm = TRUE)&res_now$Gap
      ix[is.na(ix)] <- FALSE
      data1 <- rbind(data1, apply(data1[ix,],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
      colz <-     c(colz,"red")
      ix=res_now$dist_spec_obs<quantile(res_now$dist_spec_obs,probs = .1,na.rm = TRUE)&res_now$Gap
      ix[is.na(ix)] <- FALSE
      data1 <- rbind(data1, apply(data1[ix,],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
      colz <-     c(colz,"blue")
      
#      ix=res_now$dist_spec_obs>quantile(res_now$dist_spec_obs,probs = .75,na.rm = TRUE)
#      ix[is.na(ix)] <- FALSE
#      data1 <- rbind(data1, data1[ix,][sample(1:nrow(data1[ix,]),10),])
#      colz <-     c(colz,rep(rgb(.9,.9,.9,alpha = .5),10))
      
      # data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.1,na.rm=TRUE))
      # colz <-     c(colz,colorset1[1])
      # data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.25,na.rm=TRUE))
      # colz <-     c(colz,colorset1[2])
      # data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
      # colz <-     c(colz,colorset1[3])
      # data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.9,na.rm=TRUE))
      # colz <-     c(colz,colorset1[4])
      # data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = mean,na.rm=TRUE))
      # colz <-     c(colz,"red")
      # data1 <- rbind(data1, rep(0,ncol(data1)))
      # colz <-     c(colz,"black")
    }    
    
    # for data 1 = predicted - observed
    data_plot <- data1
    data_plot[data_plot<-1.5] <- -1.5
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = paste0(trait_names[t],"_from_"),replacement = "")
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = "_lm",replacement = "")
    data_plot <- data_plot[rowSums(!is.na(data_plot))!=0,]
    data_plot <- as.data.frame(data_plot)
    data_now <- rbind(rep(-1.25,10) , rep(0,10) , data_plot[nrow(res_now):nrow(data_plot),])
    radarchart(data_now,pcol = colz[nrow(res_now):nrow(data_plot)],axistype = 2,plty = 1)
    
    
  }
  dev.off()
  
  
  pdf(file=file.path(origin,"_2021","figures","Triangle_spider.pdf"))
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    res_lm <- res_now[,grep(colnames(res_now),pattern = "_lm")]
    res_now <- res_now[,-grep(colnames(res_now),pattern = "_lm")]
    res_lm <- res_lm[,grep(colnames(res_lm),pattern = paste0(trait_names[t],"_from"))]
    res_now <- cbind(res_now,res_lm)
    colnames(res_now)
    
    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    res_now <- res_now[rowSums(!is.na(res_now))!=0,]
    
    weights=c("mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",
              colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
    spider_1=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_1) <- c("",weights)
    spider_obs=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_obs) <- c("",weights)
    spider_pred=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_pred) <- c("",weights)
    w=1
    for(w in 1:length(weights)){
      # deviation towards weight y/n
      obs = res_now$value_obs-res_now[,colnames(res_now)%in%weights[w]]      
      pred = res_now$value_pred-res_now[,colnames(res_now)%in%weights[w]]
      
      spider_obs[,colnames(spider_obs)==weights[w]] <- obs
      spider_pred[,colnames(spider_pred)==weights[w]] <- pred
      
      res_yn <- cbind(abs(pred),abs(obs))
      spider_1[,colnames(spider_1)==weights[w]] <- res_yn[,1]-res_yn[,2]
    }
    spider_obs <- as.data.frame(spider_obs[,colSums(!is.na(spider_obs))!=0])
    spider_pred <- as.data.frame(spider_pred[,colSums(!is.na(spider_pred))!=0])
    spider_1 <- as.data.frame(spider_1[,colSums(!is.na(spider_1))!=0])
    
    #install.packages("fmsb")  
    library(fmsb)
    par(mfrow=c(1,2),mar=c(2,2,2,2))
    
    summary(spider_obs)
    data1.0 <- abs(spider_obs)
    data2.0 <- abs(spider_pred)
    data1 <- data.frame(Taxonomy=data1.0[,which(apply(data1.0[,1:4],2,mean,na.rm=TRUE)==max(apply(data1.0[,1:4],2,mean,na.rm=TRUE)))],
                   Function=data1.0$mean_GF,
                   Trait=data1.0[,which(apply(data1.0[,7:ncol(data1.0)],2,mean,na.rm=TRUE)==max(apply(data1.0[,7:ncol(data1.0)],2,mean,na.rm=TRUE)))])
    data2 <- data.frame(Taxonomy=data2.0[,which(apply(data2.0[,1:4],2,mean,na.rm=TRUE)==max(apply(data2.0[,1:4],2,mean,na.rm=TRUE)))],
                        Function=data2.0$mean_GF,
                        Trait=data2.0[,which(apply(data2.0[,7:ncol(data2.0)],2,mean,na.rm=TRUE)==max(apply(data2.0[,7:ncol(data2.0)],2,mean,na.rm=TRUE)))])
    info=res_now
    hist(data1$Taxonomy)
    hist(data2$Taxonomy)
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    
    # colours
    {    
    colz=rep(rgb(.9,.9,.9,alpha = .01),nrow(data1))
    data1 <- rbind(data1, apply(data1[info$nb_spec>5&info$nb_gen>5,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,"orange")
    data1 <- rbind(data1, apply(data1[info$nb_spec==1&info$nb_gen==1,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,"blue")
    
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.1,na.rm=TRUE))
    colz <-     c(colz,colorset1[1])
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.25,na.rm=TRUE))
    colz <-     c(colz,colorset1[2])
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,colorset1[3])
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.9,na.rm=TRUE))
    colz <-     c(colz,colorset1[4])
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = mean,na.rm=TRUE))
    colz <-     c(colz,"red")
    data1 <- rbind(data1, rep(0,ncol(data1)))
    colz <-     c(colz,"black")
    }    
    {    
      colz2=rep(rgb(.9,.9,.9,alpha = .01),nrow(data2))
      data2 <- rbind(data2, apply(data2[info$nb_spec>5&info$nb_gen>5,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
      colz2 <-     c(colz2,"orange")
      data2 <- rbind(data2, apply(data2[info$nb_spec==1&info$nb_gen==1,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
      colz2 <-     c(colz2,"blue")
      
      data2 <- rbind(data2, apply(data2,MARGIN = 2,FUN = quantile,probs=.1,na.rm=TRUE))
      colz2 <-     c(colz2,colorset1[1])
      data2 <- rbind(data2, apply(data2,MARGIN = 2,FUN = quantile,probs=.25,na.rm=TRUE))
      colz2 <-     c(colz2,colorset1[2])
      data2 <- rbind(data2, apply(data2,MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
      colz2 <-     c(colz2,colorset1[3])
      data2 <- rbind(data2, apply(data2,MARGIN = 2,FUN = quantile,probs=.9,na.rm=TRUE))
      colz2 <-     c(colz2,colorset1[4])
      data2 <- rbind(data2, apply(data2,MARGIN = 2,FUN = mean,na.rm=TRUE))
      colz2 <-     c(colz2,"red")
      data2 <- rbind(data2, rep(0,ncol(data2)))
      colz2 <-     c(colz2,"black")
    }  
    
    # for data 1 = observed
    data_plot <- data1
    data_plot[data_plot>2] <- 2
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = paste0(trait_names[t],"_from_"),replacement = "")
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = "_lm",replacement = "")
    data_plot <- data_plot[rowSums(!is.na(data_plot))!=0,]
    data_plot <- as.data.frame(data_plot)
    data_now <- rbind(rep(0,10) , rep(2,10) , data_plot)
    radarchart(data_now,pcol = colz)
    
    # for data 2 = predicted
    data_plot <- data2
    data_plot[data_plot>2] <- 2
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = paste0(trait_names[t],"_from_"),replacement = "")
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = "_lm",replacement = "")
    data_plot <- data_plot[rowSums(!is.na(data_plot))!=0,]
    data_plot <- as.data.frame(data_plot)
    data_now <- rbind(rep(0,10) , rep(2,10) , data_plot)
    radarchart(data_now,pcol = colz2)
    
  }
  dev.off()
  
  
  pdf(file=file.path(origin,"_2021","figures","Spiderplot.pdf"))
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    res_lm <- res_now[,grep(colnames(res_now),pattern = "_lm")]
    summary(res_lm)
    res_now <- res_now[,-grep(colnames(res_now),pattern = "_lm")]
    res_lm <- res_lm[,grep(colnames(res_lm),pattern = paste0(trait_names[t],"_from"))]
    res_now <- cbind(res_now,res_lm)
    colnames(res_now)

    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    res_now <- res_now[rowSums(!is.na(res_now))!=0,]
    dim(res_now)
    weights=c("mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",
              colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
    spider_1=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_1) <- c("",weights)
    
    for(w in 1:length(weights)){
    # deviation towards weight y/n
        res_yn <- cbind(res_now$value_obs-res_now[,colnames(res_now)%in%weights[w]],res_now$value_pred-res_now[,colnames(res_now)%in%weights[w]])
        res_yn <- abs(res_yn)
        spider_1[,colnames(spider_1)==weights[w]] <- res_yn[,2]-res_yn[,1]
    }
    
    #install.packages("fmsb")  
    library(fmsb)
    dim(spider_1)
    data <- as.data.frame(spider_1)
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    par(mfrow=c(1,1),mar=c(2,2,2,2))
    data1 <- data[,!colnames(data)%in%"V1"]
#    data2 <- data1#rep(0,ncol(data1))
    
    colz=rep(rgb(.9,.9,.9,alpha = .01),nrow(data1))
    data1 <- rbind(data1, apply(data1[info$nb_spec>5&info$nb_gen>5,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,"orange")
    data1 <- rbind(data1, apply(data1[info$nb_spec==1&info$nb_gen==1,],MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,"blue")
    
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.25,na.rm=TRUE))
    colz <-     c(colz,"gray")
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.75,na.rm=TRUE))
    colz <-     c(colz,"gray")
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = mean,na.rm=TRUE))
    colz <-     c(colz,"black")
    data1 <- rbind(data1, rep(0,ncol(data1)))
    colz <-     c(colz,"red")

    # for nb_spec=1
    data_plot <- data1
    data_plot[data_plot>0.5] <- .5
    data_plot[data_plot<(-.6)] <- -.6
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = paste0(trait_names[t],"_from_"),replacement = "")
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = "_lm",replacement = "")
    data_plot <- data_plot[rowSums(!is.na(data_plot))!=0,]
    data_plot <- as.data.frame(data_plot)
    data_now <- rbind(rep(-.5,10) , rep(.5,10) , data_plot)
    radarchart(data_now,pcol = colz)
  
    }
  dev.off()
  
  
  res <- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))
  
  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res_now$value_obs)
  m=1
  t=1
  w=1
  pdf(file=file.path(origin,"_2021","figures","Spiderplot.pdf"))
  for(t in 1:length(trait_names)){
    res_now=res[res$trait==trait_names[t],]
    res_lm <- res_now[,grep(colnames(res_now),pattern = "_lm")]
    summary(res_lm)
    res_now <- res_now[,-grep(colnames(res_now),pattern = "_lm")]
    res_lm <- res_lm[,grep(colnames(res_lm),pattern = paste0(trait_names[t],"_from"))]
    res_now <- cbind(res_now,res_lm)
    colnames(res_now)
    
    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    res_now <- res_now[rowSums(!is.na(res_now))!=0,]
    dim(res_now)
    weights=c("mean_spec","mean_gen","mean_fam","mean_clad","mean_GF","mean_PFT",
              colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
    spider_1=matrix(NA,ncol=1+length(weights),nrow=nrow(res_now))
    colnames(spider_1) <- c("",weights)
    
    for(w in 1:length(weights)){
      # deviation towards weight y/n
      res_yn <- cbind(res_now$value_obs-res_now[,colnames(res_now)%in%weights[w]],res_now$value_pred-res_now[,colnames(res_now)%in%weights[w]])
      res_yn <- abs(res_yn)
      spider_1[,colnames(spider_1)==weights[w]] <- res_yn[,2]-res_yn[,1]
    }
    spider_1 <- cbind(spider_1,abs(res_now$value_pred-res_now$value_obs)*-1)
    colnames(spider_1)[ncol(spider_1)] <- "value"
    
    require(FactoMineR)
    PCA(data1[complete.cases(data1),])
    
    #install.packages("fmsb")  
    library(fmsb)
    dim(spider_1)
    data <- as.data.frame(spider_1)
    # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
    par(mfrow=c(1,1),mar=c(2,2,2,2))
    data1 <- data[,!colnames(data)%in%"V1"]
    
    colz=rep(rgb(.9,.9,.9,alpha = .01),nrow(data1))
    
#    data1 <- rbind(data1, data1[which(info$nb_spec>20&info$nb_gen>20),])
#    colz <-     c(colz,rep("orange",sum(info$nb_spec>20&info$nb_gen>20)))
    data1 <- rbind(data1, apply(data1[which(info$nb_spec>20&info$nb_gen>20),],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
    colz <-     c(colz,rep("orange",sum(info$nb_spec>20&info$nb_gen>20)))
    
#    data1 <- rbind(data1, data1[which(info$nb_spec==1&info$nb_gen==1),])
#    colz <-     c(colz,rep("blue",sum(info$nb_spec==1&info$nb_gen==1)))
    data1 <- rbind(data1, apply(data1[which(info$nb_spec==1&info$nb_gen==1),],MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
    colz <-     c(colz,rep("blue",sum(info$nb_spec==1&info$nb_gen==1)))
    
    data1 <- rbind(data1, apply(data1,MARGIN = 2,FUN = quantile,probs=.5,na.rm=TRUE))
    colz <-     c(colz,"gray")

    data1 <- rbind(data1, rep(0,ncol(data1)))
    colz <-     c(colz,"red")
    
    length(colz)
    dim(data1)
    
    # for nb_spec=1
    data_plot <- data1
    data_plot[data_plot>0.5] <- .5
    data_plot[data_plot<(-.6)] <- -.6
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = paste0(trait_names[t],"_from_"),replacement = "")
    colnames(data_plot) <- gsub(colnames(data_plot),pattern = "_lm",replacement = "")
    data_plot <- data_plot[rowSums(!is.na(data_plot))!=0,]
    data_plot <- as.data.frame(data_plot)
    data_now <- rbind(rep(-.5,10) , rep(.5,10) , data_plot)
    radarchart(data_now,pcol = colz)
    
  }
  dev.off()
  
  
  
  
  
  # Create data: note in High school for Jonathan:
  data <- as.data.frame(matrix( sample( 2:20 , 10 , replace=T) , ncol=10))
  colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding", "data-viz" , "french" , "physic", "statistic", "sport" )
  
  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
  data <- rbind(rep(20,10) , rep(0,10) , data)
  
  # Check your data, it has to look like this!
  # head(data)
  
  # Custom the radarChart !
  radarchart( data  , axistype=1 , 
              
              #custom polygon
              pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 , 
              
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
              
              #custom labels
              vlcex=0.8 
  )
  
  
  
  
  
  
  
  
  
pdf(file=file.path(origin,"_2021","figures","barplot.pdf"))
par(mfrow=c(1,1),mar=c(12,4,2,2))
for(m in 1:length(missingness)){
for(t in 1:length(trait_names)){
  
  res_now=res[res$trait==trait_names[t],]
  res_now <- res_now[res_now$missingness==missingness[m],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
  dim(res_now)

# calc dev = pred-obs OR use "error"
# plot 

#plot(res_now$mean_spec,abs(res_now$value_pred-res_now$value_obs))

#dat_plot=cbind(abs(res_now$value_pred-res_now$mean_spec)-abs(res_now$value_obs-res_now$mean_spec),
#               abs(res_now$value_pred-res_now$mean_gen)-abs(res_now$value_obs-res_now$mean_gen),
#               abs(res_now$value_pred-res_now$mean_fam)-abs(res_now$value_obs-res_now$mean_fam),
#               abs(res_now$value_pred-res_now$mean_clad)-abs(res_now$value_obs-res_now$mean_clad))

#boxplot(dat_plot)
#abline(h=0)

dat_plot=cbind(res_now$error,res_now$value_obs,res_now$value_pred,res_now$mean_spec,res_now$mean_gen,
               res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT,
               res_now[,grep(colnames(res_now),pattern = "_lm")])
#require(FactoMineR)
#PCA(dat_plot[complete.cases(dat_plot),])
dat_plot=cbind(res_now$value_obs,res_now$value_pred,res_now$mean_spec,res_now$mean_gen,
               res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT,
               res_now[,grep(colnames(res_now),pattern = "_lm")])


cor_out <- cor(dat_plot,use = "pairwise.complete.obs")[,1:2]
rownames(cor_out) <- c("value_obs","value_pred","mean_spec","mean_gen",
                       "mean_fam","mean_clad","mean_GF","mean_PFT",colnames(res_now)[grep(colnames(res_now),pattern = "_lm")])
barplot(c(median(res_now$value_obs-res_now$value_pred,na.rm = TRUE),cor_out[3:nrow(cor_out),2]-cor_out[3:nrow(cor_out),1]),las=2,ylim=c(-.3,.3),
        main=paste0(trait_names[t],"_",missingness[m]),col=c("gray",colz1),ylab="cor_pred - cor_obs")
abline(h=seq(-.5,to=.5,by=.1),col="gray",lty=2)

}
}
dev.off()


pairs(dat_plot)


barplot(c(median(res_now$value_pred,na.rm = TRUE),median(res_now$value_obs,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_spec,na.rm = TRUE),median(res_now$value_obs-res_now$mean_spec,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_gen,na.rm = TRUE),median(res_now$value_obs-res_now$mean_gen,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_fam,na.rm = TRUE),median(res_now$value_obs-res_now$mean_fam,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_clad,na.rm = TRUE),median(res_now$value_obs-res_now$mean_clad,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_GF,na.rm = TRUE),median(res_now$value_obs-res_now$mean_GF,na.rm = TRUE)))
barplot(c(median(res_now$value_pred-res_now$mean_PFT,na.rm = TRUE),median(res_now$value_obs-res_now$mean_PFT,na.rm = TRUE)))

boxplot(cbind((res_now$value_pred-res_now$mean_spec),(res_now$value_obs-res_now$mean_spec)))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_gen),(res_now$value_obs-res_now$mean_gen))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_fam),(res_now$value_obs-res_now$mean_fam))
abline(h=0)
boxplot((res_now$value_pred-res_now$mean_clad),(res_now$value_obs-res_now$mean_clad))
abline(h=0)

plot(res_now$mean_gen,abs(res_now$error))
plot(abs(res_now$error),res_now$mean_gen)
plot(res_now$mean_gen,abs(res_now$error))
plot(res_now$mean_fam,res_now$error)
plot(res_now$mean_clad,res_now$error)
plot(res_now$mean_GF,res_now$error)
plot(res_now$mean_PFT,res_now$error)
# all lms
# plot(res_now$error,res_now$lm)

dat_cor <- cbind(res_now$error,res_now$mean_spec,res_now$mean_gen,res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT)
dat_cor <- cbind(res_now$error,res_now$mean_spec,res_now$mean_gen,res_now$mean_fam,res_now$mean_clad,res_now$mean_GF,res_now$mean_PFT)
pairs(dat_cor)
cor(dat_cor)


