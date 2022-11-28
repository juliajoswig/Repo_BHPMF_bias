

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
t_choice="data_2"
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
plot_vals_obs <- list()
plot_vals_pred <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  plot_vals_pred[[t]] <- res$value_pred_zlog[ix_trait]
  plot_vals_obs[[t]] <- res$value_obs_zlog[ix_trait]
  colz_now = grep(x = colnames(res),pattern = paste0(trait_names[t],"_from_"))
  colz_now = colz_now[grep(x = colnames(res)[colz_now],pattern = "zlog")]
  colz_obs = colz_now[grep(x = colnames(res)[colz_now],pattern = "obs")]
  colz_pred = colz_now[grep(x = colnames(res)[colz_now],pattern = "pred")]
  colnames(res)[colz_obs]
  colnames(res)[colz_pred]
  colnames(res)[colz_now]
  obs  =    res$value_obs_zlog - res[ix_trait,colz_obs]
  pred =    res$value_pred_zlog - res[ix_trait,colz_pred]
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
  colz[[t]] <- rep(colzTRAITS[match(unique(gsub(colnames(bxpl_tot[[t]])[1:4],pattern = "_zlog_pred",replacement = "")),trait_names)],2)[c(1,1,2,2,3,3,4,4,5,5)]
}

colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))


  par(mfrow=c(1,1),mar=c(8,4,3,2))
  
  pdf(file=file.path(origin,"_2021","figures","Figure_3",paste0("Figure_S3_lm",missingness[m],"_2_zlog.pdf")),width=20,height=8)
  {
    par(mar=c(10,3,4,0),mfrow=c(1,5))
    par(mar=c(4,0,4,0))
    layout(mat = matrix(c(1:6,
                          1:6,
                          1:6,
                          1:6,
                          7,rep(8,5)),nrow = 5,byrow = TRUE))
    
    ytop=c(4,4,4,4,4,4)
    ybot=c(0,0,0,0,0,0)

    plot(1:10,col="white",ylim=c(-2,2),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
    #axis(side = 2,line = -17,cex.axis=2)
    axis(side = 2,at=seq(-2,2,by=2),line = -10,labels = c("","Distance",""),tick = FALSE,cex.axis=4)
    axis(side = 2,at = seq(-2,2,by=2),line = -14,labels = c("",paste0("to lm"),""),tick = FALSE,cex.axis=4)

    t=1
    t2=2
    for(t in 1:5){

      t2=1
      bxpl=rep(NA,nrow(bxpl_tot_op[[t]]))
      for(t2 in 1:4){
        bxpl=cbind(bxpl,abs(bxpl_tot_op[[t]][,c(1,5)+(t2-1)]))
      }
      bxpl <- bxpl[,-1]
      mnn=gsub(colnames(bxpl),pattern = "_zlog",replacement = "")
      mnn=gsub(mnn,pattern = "_obs",replacement = "")
      mnn=gsub(mnn,pattern = "_pred",replacement = "")
      boxplot(bxpl,cex.main=3,
                ylab="",col=c(colz_solid[1],colz_solid[2]),xaxt="n",ylim=c(ybot[t2],ytop[t2]),
                las=2,cex.axis=1.5,cex.lab=1,lwd=2,frame=FALSE,main=trait_names[t],xlab="")
      axis(side = 1,at = seq(from=1,to = 8,by = 1),labels = rep(c("obs","pred"),4),las=2,cex.axis=1.6)
      
      for(t2 in 1:4){
        yBot_now=median(bxpl[,(t2*2-1)],na.rm = TRUE)
        rect(xleft = .5+(2*t2-2),xright = 2.5+(2*t2-2),ybottom = ybot[t2],ytop = yBot_now,col=colzBoTop[5],border = NA)
        rect(xleft = .5+(2*t2-2),xright = 2.5+(2*t2-2),ybottom = yBot_now,ytop = ytop[t2],col=colzBoTop[9],border = NA)
        rect(xleft = .5+(2*t2-2),xright = 2.5+(2*t2-2),ybottom = yBot_now,ytop = ytop[t2],col=colzBoTop[9],border = NA)
        #text(x=1.5,y = ytop[[t]][t2],labels=trait_names[t2],cex=3)
      }
        boxplot(bxpl,cex.main=3,
                ylab="",col=colz[[t]],xaxt="n",#ylim=c(ybot,ytop)
                las=2,cex.axis=1.5,cex.lab=1,lwd=2,frame=FALSE,add=TRUE)
    }
    
    
    
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    plot(1:12,xaxt="n",yaxt="n",frame=FALSE,col="white",xlab="",ylab="")
    rect(xleft = 2.5,ybottom = 0,xright = 10.5,ytop = 12,col = "lightgray",border = FALSE)
    legend(3, 14, "SLA", col = colzTRAITS[1], pch = 15, cex = 3,bty = "n")
    legend(3.9,14, "Height", col = colzTRAITS[3], pch = 15, cex = 3,bty = "n")
    legend(5, 14, "SSD", col = colzTRAITS[5], pch = 15, cex = 3,bty = "n")
    legend(6, 14, "LeafN", col = colzTRAITS[7], pch = 15, cex = 3,bty = "n")
    legend(7, 14, "LeafP", col = colzTRAITS[9], pch = 15, cex = 3,bty = "n")
    legend(8, 14, "LeafNArea", col = colzTRAITS[11], pch = 15, cex = 3,bty = "n")
    
    }
  dev.off()
  
 