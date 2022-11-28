
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

res <- read.csv(file=file.path(origin,"_2021","data","analyes","Point_wise","res.csv"))

#cut down res for now:
trait_names=unique(res$trait)
trait_names <- as.character(trait_names[!is.na(trait_names)])
ix_res=!is.na(res$trait)
res <- res[ix_res,]

colz=c(
  "#f7f4f9",
  "#e7e1ef",
  "#d4b9da",
  "#c994c7",
  
  "#df65b0",
  "#e7298a",
  "#ce1256",
  "#980043",
  "#67001f")

colz2=c(
  "#f7f4f9",#tax
  "#f7f4f9",#tax
  "#e7e1ef",#tax
  "#e7e1ef",#tax
  "#d4b9da",#tax
  "#d4b9da",#tax
  "#c994c7",#tax
  "#c994c7",#tax
  "#addd8e",#fun
  "#addd8e",#fun
  "#31a354",#fun
  "#31a354",#fun
  "#f7fbff",
  "#f7fbff",
  "#deebf7",
  "#deebf7",
  "#c6dbef",
  "#c6dbef",
  "#9ecae1",
  "#9ecae1",
  "#6baed6",
  "#6baed6",
  "#4292c6",
  "#4292c6",
  "#2171b5",
  "#2171b5",
  "#08519c",
  "#08519c",
  "#08306b",
  "#08306b"
)
colz1=c(
  "#f7f4f9",#tax
  "#e7e1ef",#tax
  "#d4b9da",#tax
  "#c994c7",#tax
  "#addd8e",#fun
  "#31a354",#fun
  "#f7fbff",
  "#deebf7",
  "#c6dbef",
  "#9ecae1",
  "#6baed6",
  "#4292c6",
  "#2171b5",
  "#08519c",
  "#08306b"
)

t=1
pdf(file=file.path(origin,"_2021","figures","Weight_strength1.pdf"),width = 14,height = 10)
par(mfrow=c(1,1),mar=c(10,4,2,2))
for(t in 1:length(trait_names)){
  res_now <- res[res$trait==trait_names[t],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
  # gibt es deviations?
  data_t1_pred=res_now[,c(grep(colnames(res_now),pattern = "mean"),
                          grep(colnames(res_now),pattern = "lm"))]
  data_t1_act=res_now[,c(grep(colnames(res_now),pattern = "value_obs"),
                         grep(colnames(res_now),pattern = "value_pred"))]
  cn=1
  nms=NA
  bxpl=rep(NA,nrow(data_t1_pred))
  for(cn in 1:ncol(data_t1_pred)){
    bxpl <- cbind(cbind(bxpl,cbind(abs(data_t1_act[,1])-abs(data_t1_pred[,cn]),
                                   abs(data_t1_act[,2])-abs(data_t1_pred[,cn]))
    ))
    nms=c(nms,paste0("dist(obs,", colnames(data_t1_pred)[cn],")"),
          paste0("dist(pred,",colnames(data_t1_pred)[cn],")"))
  }
  colnames(bxpl) <- nms
  bxpl <- bxpl[,colSums(!is.na(bxpl))!=0]
  boxplot(bxpl,ylab="distance to weight",las=2,ylim=c(-2,2),main=trait_names[t],col=colz2)
  abline(h=0,col="gray")
  abline(h=seq(-5,to=5,by = .5),col="gray",lty=3)
  abline(v = seq(0.52,to = 60,by = 2),col="gray",lty=1,lwd=2)
  boxplot(bxpl,ylab="distance to weight",las=2,ylim=c(-2,2),col=colz2,add=TRUE)
}

dev.off()

pdf(file=file.path(origin,"_2021","figures","Weight_strength1.pdf"),width = 14,height = 10)
par(mfrow=c(1,1),mar=c(10,4,2,2))
for(t in 1:length(trait_names)){
  res_now <- res[res$trait==trait_names[t],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
  # gibt es deviations?
  data_t1_pred=res_now[,c(grep(colnames(res_now),pattern = "mean"),
                          grep(colnames(res_now),pattern = "lm"))]
  data_t1_act=res_now[,c(grep(colnames(res_now),pattern = "value_obs"),
                          grep(colnames(res_now),pattern = "value_pred"))]
  cn=1
  nms=NA
  bxpl=rep(NA,nrow(data_t1_pred))
  for(cn in 1:ncol(data_t1_pred)){
    bxpl <- cbind(cbind(bxpl,cbind(abs(data_t1_act[,1])-abs(data_t1_pred[,cn]),
                                   abs(data_t1_act[,2])-abs(data_t1_pred[,cn]))
                        ))
    nms=c(nms,paste0("dist(obs,", colnames(data_t1_pred)[cn],")"),
              paste0("dist(pred,",colnames(data_t1_pred)[cn],")"))
  }
  colnames(bxpl) <- nms
  bxpl <- bxpl[,colSums(!is.na(bxpl))!=0]
  boxplot(bxpl,ylab="distance to weight",las=2,ylim=c(-2,2),main=trait_names[t],col=colz2)
  abline(h=0,col="gray")
  abline(h=seq(-5,to=5,by = .5),col="gray",lty=3)
  abline(v = seq(0.52,to = 60,by = 2),col="gray",lty=1,lwd=2)
  boxplot(bxpl,ylab="distance to weight",las=2,ylim=c(-2,2),col=colz2,add=TRUE)
}

dev.off()
    
  pdf(file=file.path(origin,"_2021","figures","Weight_strength.pdf"),width = 14,height = 10)
  par(mfrow=c(4,6),mar=c(4,4,2,2))
  par(mfrow=c(1,1),mar=c(4,4,2,2))
  for(t in 1:length(trait_names)){
    res_now <- res[res$trait==trait_names[t],]
    res_now <- res_now[,colSums(!is.na(res_now))!=0]
    # gibt es deviations?
    data_t1_pred=res_now[,c(grep(colnames(res_now),pattern = "mean"),
                            grep(colnames(res_now),pattern = "lm"))]
    data_t1_act=res_now[,c(grep(colnames(res_now),pattern = "value_obs"),
                           grep(colnames(res_now),pattern = "value_pred"))]
    cn=1
    nms=NA
    bxpl=rep(NA,nrow(data_t1_pred))
    for(cn in 1:ncol(data_t1_pred)){
      bxpl <- cbind(cbind(bxpl,cbind(abs(data_t1_act[,1])-abs(data_t1_pred[,cn]),
                                     abs(data_t1_act[,2])-abs(data_t1_pred[,cn]))
      ))
      nms=c(nms,paste0("dist(obs,", colnames(data_t1_pred)[cn],")"),
            paste0("dist(pred,",colnames(data_t1_pred)[cn],")"))
    }
    colnames(bxpl) <- nms
    bxpl <- bxpl[,colSums(!is.na(bxpl))!=0]
    
    boxplot(bxpl[,seq(from=2,to=ncol(bxpl),by=2)]-bxpl[,seq(1,ncol(bxpl),2)],las=2,main=trait_names[t],ylim=c(-.5,.2),col=colz1)
    abline(h=0,col="gray")
    
    bxpl_median <- apply(bxpl,MARGIN = 2,FUN = median,na.rm=TRUE)
    barplot(bxpl_median[seq(from=2,to=length(bxpl_median),by=2)]-bxpl_median[seq(1,length(bxpl_median),2)],las=2,main=trait_names[t],ylim=c(-.2,.1))
    bxplt_dist<-bxpl_median[seq(from=2,to=length(bxpl_median),by=2)]-bxpl_median[seq(3,length(bxpl_median),2)] 
#  barplot(c(c(dist(bxpl_median[seq(from=3,to=16,by=12)])),
#            c(dist(bxpl_median[seq(from=5,to=24,by=12)])),
#            c(dist(bxpl_median[seq(from=7,to=24,by=12)])),
#            c(dist(bxpl_median[seq(from=9,to=30,by=12)])),
#            c(dist(bxpl_median[seq(from=11,to=30,by=12)])),
#            c(dist(bxpl_median[seq(from=13,to=30,by=12)]))))
  axis(1,at = 1:6,labels =   colnames(bxpl)[c(3,5,7,9,11,13)],las=2)
  bxpl_median
#  barplot(bxpl_median[seq(from=2,to=16,by=11)],las=2)
  
  cn=17
#  for(cn in 1:ncol(data_t1_pred)){
#    plot(abs(data_t1_pred[,cn]),abs(data_t1_act[,1]),
#         ylab=colnames(data_t1_act)[1],
#         xlab=colnames(data_t1_pred)[cn],pch=16,cex=.5,
#         main=paste0(colnames(data_t1_pred)[cn],"_",trait_names[t]),xlim = c(0,3),ylim = c(0,3))
#    for(i in 1:nrow(data_t1_pred)){
#      lines(x = cbind(abs(data_t1_pred[i,cn]),abs(data_t1_pred[i,cn])),
#            y= abs(data_t1_act[i,1:2]),col="gray",pch=16,cex=1)
#    }
#    points(abs(data_t1_pred[,cn]),abs(data_t1_act[,2]),col="green",pch=16,cex=.5)
#    abline(0,1)
#  }
  
}
dev.off()

pdf(file=file.path(origin,"_2021","figures","Pred_Obs2.pdf"),width = 14,height = 10)
par(mfrow=c(4,6),mar=c(4,4,2,2))
for(t in 1:length(trait_names)){
  res_now <- res[res$trait==trait_names[t],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
  # gibt es deviations?
  data_t1_o=res_now[,c(grep(colnames(res_now),pattern = "obs"))]
  data_t1_pred=res_now[,c(grep(colnames(res_now),pattern = "pred"))]
  data_t1_dev=res_now[,c(grep(colnames(res_now),pattern = "RESdev"))]

  for(cn in 1:ncol(data_t1_o)){
  plot(data_t1_o[,cn],data_t1_pred[,cn],main=colnames(data_t1_o)[cn],xlim = c(-5,5),ylim = c(-5,5))
    abline(0,1)
  }
}

par(mfrow=c(4,6),mar=c(4,4,2,2))
  for(cn in 1:ncol(data_t1_o)){
#    plot(data_t1_o[,cn],data_t1_pred[,cn],main=colnames(data_t1_o)[cn],xlim = c(-5,5),ylim = c(-5,5))
#    abline(0,1)
    plot(data_t1_pred[,cn]-data_t1_o[,cn],
         res_now$error,
         main=colnames(data_t1_o)[cn],
         xlim = c(-5,5),ylim = c(-5,5))
    txt_now <- cor(data_t1_pred[,cn]-data_t1_o[,cn],res_now$error,use = "pairwise.complete.obs")
    text(-2,2,labels = round(txt_now,digits = 2),col="blue",cex=2)
  }
dev.off()
  
  
  bx=rep(NA,nrow(data_t1_o))
  for(cn in 1:ncol(data_t1_o)){
    bx=  cbind(bx,-(data_t1_o[,cn]-data_t1_pred[,cn])-res_now$error)
  }
  par(mfrow=c(1,1),mar=c(10,2,2,2))
  boxplot(bx,ylim=c(-2,2),xaxt="n")
  axis(side = 1,at = seq(1,to = ncol(bx)*1,by = 1),
       labels = c("",gsub(colnames(data_t1_o),pattern = "obs",replacement = "")),
       las=2)
  boxplot(bx,ylim=c(-4,4),xaxt="n")
  axis(side = 1,at = seq(1,to = ncol(bx)*1,by = 1),
       labels = c("",gsub(colnames(data_t1_o),pattern = "obs",replacement = "")),
       las=2)
  
  # Library
#  install.packages("fmsb")
  library(fmsb)
  data <- as.data.frame(matrix( sample( 2:20 , 10 , replace=T) , ncol=10))
  colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding", "data-viz" , "french" , "physic", "statistic", "sport" )
  
  colnames(bx) <- c("",colnames(data_t1_o))
  bx2 <- bx[,2:ncol(bx)]
  bx2 <- bx2[,-grep(colnames(bx2),pattern = "RM")]
  data <- as.data.frame(abs(bx2))
  data[data>2] <- 2
  names(data) <- gsub(colnames(bx2),pattern = "obs",replacement = "")
  # To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
  #data <- rbind(rep(2,10) , rep(0,10) , data)

  par(mfrow=c(1,1),mar=c(2,2,2,2))
  data1 <- data
  data1 <- rbind(data,apply(data,MARGIN =  2,FUN=mean,na.rm=TRUE))
  colz=rep(rgb(.9,.9,.9,alpha = .01),nrow(data1))
  colz[nrow(data1)] <- "red"
  data_now <- rbind(rep(2,10) , rep(0,10) , data1)
  radarchart(data_now ,pcol = colz)
  #radarchart(data1[nrow(data1),],pcol = "red")
  dim(data1)
  
  pca_input=bx[,2:ncol(bx)]
  colnames(pca_input) <- gsub(colnames(data_t1_o),pattern = "obs",replacement = "")
  require(FactoMineR)
  pca_out <- PCA(pca_input[complete.cases(pca_input),])
  
  pchs=rep(1,nrow(res_now))
  pchs[res_now$nb_spec==1|res_now$nb_gen==1|res_now$nb_fam==1|res_now$nb_clad==1] <- 16
  colz=rep("gray",nrow(res_now))
  colz[res_now$nb_spec==1] <- "red" 
  colz[res_now$nb_gen==1] <- "orange" 
  colz[res_now$nb_fam==1] <- "lightblue" 
  colz[res_now$nb_clad==1] <- "blue" 
  par(mfrow=c(1,1),mar=c(4,4,2,2))
  plot(pca_out$ind$coord[,1:2],col=colz,pch=pchs,xlim=c(-5,5),ylim=c(-5,5))
  # pairs(cbind(data_t1_o,data_t1_dev))
  # pairs(cbind(res_now$error,data_t1_dev))
  
  # jo
  #ja.
  i=2
  cors=NA
  name_cors=NA
  for(i in 1:ncol(data_t1_o)){
    cors <- c(cors,cor(cbind(data_t1_o[,i],data_t1_dev[,i]),use = "pairwise.complete.obs")[1,2])
    #  name_cors <- c(name_cors,gsub(colnames(data_t1_o)[i],pattern = "_RESobs",replacement = ""))
    name_cors <- c(name_cors,gsub(colnames(data_t1_o)[i],pattern = "_RESpred",replacement = ""))
  }
  
  # taxonomic data
  data_now=res_now[grep(colnames(res_now),pattern = "dist_AV")]
  data_now=data_now[-grep(colnames(res_now),pattern = "RM")]
  head(data_now)
  data_now_o=data_now[c(grep(colnames(data_now),pattern = "obs"))]
  data_now_dev=data_now[c(grep(colnames(data_now),pattern = "pred"))]-data_now_o
  colnames(data_now_dev) <- gsub(colnames(data_now_dev),pattern = "_pred",replacement = "_dev")
  
  for(i in 1:ncol(data_now_o)){
    cors <- c(cors,cor(cbind(data_now_o[,i],data_now_dev[,i]),use = "pairwise.complete.obs")[1,2])
    name_cors <- c(name_cors,gsub(colnames(data_now_o)[i],pattern = "_obs",replacement = ""))
  }
  par(mar=c(10,4,2,1))
  barplot(abs(cors),ylab="Pearson cor deviation-observaration",main=trait_names[t])
  axis(1,seq(1,length(cors)*1.15,by = 1.15),labels = name_cors,las=2,tick = FALSE)
}
dev.off()



t=1
pdf(file=file.path(origin,"_2021","figures","tmp.pdf"))
for(t in 1:length(trait_names)){
  res_now <- res[res$trait==trait_names[t],]
  res_now <- res_now[,colSums(!is.na(res_now))!=0]
# gibt es deviations?
  data_t1_o=res_now[,c(grep(colnames(res_now),pattern = "RESobs"))]
  data_t1_pred=res_now[,c(grep(colnames(res_now),pattern = "RESpred"))]
  data_t1_dev=res_now[,c(grep(colnames(res_now),pattern = "RESdev"))]
  plot(data_t1_o$LeafNArea_from_SLA_RESobs,data_t1_pred$LeafNArea_from_SLA_RESpred)
  
# pairs(cbind(data_t1_o,data_t1_dev))
# pairs(cbind(res_now$error,data_t1_dev))

# jo
#ja.
i=2
cors=NA
name_cors=NA
for(i in 1:ncol(data_t1_o)){
  cors <- c(cors,cor(cbind(data_t1_o[,i],data_t1_dev[,i]),use = "pairwise.complete.obs")[1,2])
#  name_cors <- c(name_cors,gsub(colnames(data_t1_o)[i],pattern = "_RESobs",replacement = ""))
  name_cors <- c(name_cors,gsub(colnames(data_t1_o)[i],pattern = "_RESpred",replacement = ""))
}

# taxonomic data
data_now=res_now[grep(colnames(res_now),pattern = "dist_AV")]
data_now=data_now[-grep(colnames(res_now),pattern = "RM")]
head(data_now)
data_now_o=data_now[c(grep(colnames(data_now),pattern = "obs"))]
data_now_dev=data_now[c(grep(colnames(data_now),pattern = "pred"))]-data_now_o
colnames(data_now_dev) <- gsub(colnames(data_now_dev),pattern = "_pred",replacement = "_dev")

for(i in 1:ncol(data_now_o)){
  cors <- c(cors,cor(cbind(data_now_o[,i],data_now_dev[,i]),use = "pairwise.complete.obs")[1,2])
  name_cors <- c(name_cors,gsub(colnames(data_now_o)[i],pattern = "_obs",replacement = ""))
}
par(mar=c(10,4,2,1))
  barplot(abs(cors),ylab="Pearson cor deviation-observaration",main=trait_names[t])
  axis(1,seq(1,length(cors)*1.15,by = 1.15),labels = name_cors,las=2,tick = FALSE)
}
dev.off()

cor_m <- as.matrix(cor(cbind(data_now_o,data_now_dev),use = "pairwise.complete.obs"))
nms <- grep(colnames(cor_m),pattern = "obs")
cor_mc <- cor_m[grep(colnames(cor_m),pattern = "dev"),grep(colnames(cor_m),pattern = "obs")]


#Zusammenhänge?

  xs <- round(res$dist_spec/max(res$dist_spec,na.rm = TRUE),digits = 3)
  xg <- round(res$dist_gen/max(res$dist_gen,na.rm = TRUE),digits = 3)
  xf <- round(res$dist_fam/max(res$dist_fam,na.rm = TRUE),digits = 3)
  xc <- round(res$dist_clad/max(res$dist_clad,na.rm = TRUE),digits = 3)
  xc2 <- round((res$dist_fam+res$dist_clad)/max((res$dist_fam+res$dist_clad)),digits = 3)
  
  # kick out those values which have any tax info
  c_DIST=rep(0,nrow(res))
  c_DIST[xs>quantile(xs,probs = .90)] <- 
    c_DIST[xs>quantile(xs,probs = .90)]+.25
  c_DIST[xg>quantile(xg,probs = .90)] <- 
    c_DIST[xg>quantile(xg,probs = .90)]+.25
  c_DIST[xf>quantile(xf,probs = .90,na.rm = TRUE)] <- 
    c_DIST[xf>quantile(xf,probs = .90,na.rm = TRUE)]+.25
  c_DIST[xc>quantile(xc,probs = .90,na.rm = TRUE)] <- 
    c_DIST[xc>quantile(xc,probs = .90,na.rm = TRUE)]+.25
  #-----
  perc=.60
  c_DIST_sp=rep(0,nrow(res))
  c_DIST_sp[xs>quantile(xs,probs = perc)] <- 
    c_DIST_sp[xs>quantile(xs,probs = perc)]+1
  c_DIST_gen=rep(0,nrow(res))
  c_DIST_gen[xg>quantile(xg,probs = perc)] <- 
    c_DIST_gen[xg>quantile(xg,probs = perc)]+1
  c_DIST_fc=rep(0,nrow(res))
  c_DIST_fc[xc2>quantile(xc2,probs = perc,na.rm = TRUE)] <- 
    c_DIST_fc[xc2>quantile(xc2,probs = perc,na.rm = TRUE)]+1
  
# xs[xs>quantile(xs,probs = .90,na.rm = TRUE)] <- 0
# xg[xg>quantile(xg,probs = .90,na.rm = TRUE)] <- 0
# xf[xf>quantile(xf,probs = .90,na.rm = TRUE)] <- 0
# xc[xc>quantile(xc,probs = .90,na.rm = TRUE)] <- 0
  
#  xs[xs<=quantile(xs,probs = .90,na.rm = TRUE)] <- 1
#  xg[xg<=quantile(xg,probs = .90,na.rm = TRUE)] <- 1
#  xf[xf<=quantile(xf,probs = .90,na.rm = TRUE)] <- 1
#  xc[xc<=quantile(xc,probs = .90,na.rm = TRUE)] <- 1
  
  rgb_input <- cbind(xs,xg,xf)
  rgb_input <- cbind(c_DIST_sp,c_DIST_gen,c_DIST_fc)
  #rgb_input[is.na(rgb_input)] <- 0
  colz_rgbinv=rgb((1-rgb_input),alpha = .9)
  colz_rgb=rgb(rgb_input,alpha = .9)
  
  data_t1 <- data_t1_a
  ix=colSums(!is.na(data_t1))!=0
  pairs(data_t1[,ix],col=colz_rgbinv,pch=16)
  pairs(data_t1[,ix],col=colz_rgb,pch=16)
  














Dev=res$dist_AVspec_pred[ix] - res$dist_spec[ix]
Dev_pred = (-1)*(res$dist_AVspec_pred[ix] - res$dist_spec[ix])
importance_Spec <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVgen_pred[ix] - res$dist_gen[ix])
importance_Gen <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVfam_pred[ix] - res$dist_fam[ix])
importance_Fam <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVclad_pred[ix] - res$dist_clad[ix])
importance_Clad <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])

Dev=res$dist_AVspecRM_pred[ix] - res$dist_specRM[ix]
Dev_pred = (-1)*(res$dist_AVspecRM_pred[ix] - res$dist_specRM[ix])
importance_SpecRM <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVgenRM_pred[ix] - res$dist_genRM[ix])
importance_GenRM <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVfamRM_pred[ix] - res$dist_famRM[ix])
importance_FamRM <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])
Dev_pred = (-1)*(res$dist_AVcladRM_pred[ix] - res$dist_cladRM[ix])
importance_CladRM <- max(abs(Dev_pred-abs(res$error)[ix]),na.rm = TRUE)-abs(Dev_pred-abs(res$error)[ix])

pdf(file=file.path(origin,"_2021","figures","Value_wise_importance_pairs.pdf"),width = 14,height = 14)
  pairs(cbind(importance_Spec,importance_Gen,importance_Fam,importance_Clad),col=colz[ix],pch=16)
dev.off()
require(FactoMineR)
PCA(cbind(importance_Spec,importance_Gen,importance_Fam,importance_Clad))

pdf(file=file.path(origin,"_2021","figures","Value_wise_importanceRM_pairs.pdf"),width = 14,height = 14)
  pairs(cbind(importance_SpecRM,importance_GenRM,importance_FamRM,importance_CladRM),col=colz[ix],pch=16)
dev.off()

pdf(file=file.path(origin,"_2021","figures","Value_wise_importance.pdf"))
  boxplot(cbind(importance_Spec,importance_Gen,importance_Fam,importance_Clad),col=c('#ca0020',"#f4a582","#92c5de","#0571b0"),
            xaxt="n",ylab="Attractor importance",xlab="Deviation origin")
  axis(1,at = 1:4,labels = c("Species","Genus","Family","Clade"),las=2)
dev.off()

pdf(file=file.path(origin,"_2021","figures","Value_wise_importanceRM.pdf"))
  par(mfrow=c(1,1))
  boxplot(cbind(importance_Spec,importance_SpecRM,importance_Gen,importance_GenRM,
                importance_Fam,importance_FamRM,importance_Clad,importance_CladRM),
          col=c('#ca0020','#ca0020',"#f4a582","#f4a582","#92c5de","#92c5de","#0571b0","#0571b0"),
          xaxt="n",ylab="Attractor importance",xlab="Deviation origin")
  axis(1,at = 1:8,labels = c("Species","randomized Sp","Genus","randomized gen","Family","randomized Fam","Clade","randomized Clade"),las=2)
dev.off()

x <- importance_Spec
x <- importance_Gen
dat <- data.frame(x = x,y = x^(.5))
rbPal <- colorRampPalette(c('#ca0020',"#f4a582","#92c5de","#0571b0"))
colz_sp <- rbPal(20)[as.numeric(cut(dat$y,breaks = 20))]
par(mfrow=c(1,1))
plot(Dev,abs(res$error)[ix],pch=pchs,xlim=c(3,-3),
     cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
plot(Dev,abs(res$error)[ix],pch=pchs,xlim=c(3,-3),
     cex=.7,col=colz_sp,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
plot(Dev_pred,abs(res$error)[ix],pch=pchs,
     cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
error=abs(res$error)[ix]

dev.off()
plot(error,Dev,col=colz_rgb,pch=16)
plot(error,Dev)

# No species information | Species information
# No genus information | genus information
# No fam information | fam information


pdf(file=file.path(origin,"_2021","figures","Value_wise_rainfor_dissimilarity.pdf"),width = 14,height = 14)
par(mfrow=c(2,2))
colz=colz_rgb_reverse[ix]
pchs=pchs[ix]
{
  plot(res$dist_AVspec_pred[ix] - res$dist_spec[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-species",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  abline(0,1)
  abline(0,-1)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  
  plot(res$dist_AVgen_pred[ix] - res$dist_gen[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-genera",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVfam_pred[ix] - res$dist_fam[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-families",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
  
  plot(res$dist_AVclad_pred[ix] - res$dist_clad[ix],abs(res$error)[ix],pch=pchs,
       cex=.7,col=colz,xlab = "Deviation of single values' distance to their co-clades",ylab = "pred-obs (absolute)")
  abline(v = 0,col="gray",lty=2)
  abline(h = 0,col="gray",lty=2)
  text(1.5,.5,labels = "no go")
  text(1.5,.3,labels = "area")
  text(-3.3,3.3,labels = "BHPMF")
  text(-3,3,labels = "target")
  text(-2.7,2.7,labels = "relationship")
  legend(x = -3,y = 6.5,legend = c("species","species&genus","genus","genus&family","family","none","all groups"),
         title = "Different from indiv of",pch=16,col=c(rgb(1,0,0),rgb(1,1,0),rgb(0,1,0),rgb(0,1,1),rgb(0,0,1),rgb(0,0,0),rgb(1,1,1)))
  abline(0,1)
  abline(0,-1)
}
dev.off()

#--------------------------------------------------------------------------------
# PCA
#--------------------------------------------------------------------------------

res_pca <- data.frame(error=res$error,
                      dist_spec=res$dist_spec,dist_gen=res$dist_gen,dist_fam=res$dist_fam,dist_clad=res$dist_clad,
                      # nb_ind=res$nb_ind,nb_spec=res$nb_spec,nb_gen=res$nb_gen,nb_fam=res$nb_fam,
                      nb_spec_gap=res$nb_spec_gap,nb_gen_gap=res$nb_gen_gap,nb_fam_gap=res$nb_fam_gap,nb_clad_gap=res$nb_clad_gap,
                      dist_AVspec_dev=res$dist_AVgen_pred-res$dist_AVspec_obs,
                      dist_AVgen_dev=res$dist_AVfam_pred-res$dist_AVgen_obs,
                      dist_AVfam_dev=res$dist_AVclad_pred-res$dist_AVfam_obs,
                      dist_AVclad_dev=res$dist_AVclad_pred-res$dist_AVclad_obs
)

pdf(file=file.path(origin,"_2021","figures","pairs.pdf"),height=15, width=15)
pairs(as.matrix(res_pca[,c(grep(colnames(res_pca),pattern = ),
)]))
dev.off()

