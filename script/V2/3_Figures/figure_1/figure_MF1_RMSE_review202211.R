
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
  # start 20221003 ##############################################
  origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/"
  # end 20221003 ##############################################
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
# Test if RMSE increases with gap-size
#--------------------------------------

sel_ObSpe1="Obs_obs_TD" # TD2td 
pattern_now="RMSE"

#--------------------------------------
# select TD2td RMSE results of all repetitions
#--------------------------------------
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

plot_TD2td <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
plot_TD2td <- plot_TD2td[,colSums(!is.na(plot_TD2td))!=0]
plot_TD2td <- plot_TD2td[complete.cases(plot_TD2td),]

#--------------------------------------
# select TDtd RMSE results of all repetitions
#--------------------------------------
#chose data (rainfor)
ix2=res[,colnames(res)=="TraitChoice"]=="data"&res[,colnames(res)=="Obs_or_Spec"]=="Obs_obs_TD"
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent2) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

plot_TDtd <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
plot_TDtd <- plot_TDtd[,colSums(!is.na(plot_TDtd))!=0]
plot_TDtd <- plot_TDtd[complete.cases(plot_TDtd),]
plot_TDtd <- plot_TDtd[plot_TDtd[,1]!=0,]

cbind(plot_TD2td[,2],plot_TDtd[,2])
#--------------------------------------
# select TDenv RMSE results of all repetitions
#--------------------------------------
#chose envelope rainfor
ix3=res$TraitChoice=="data"&res$Obs_or_Spec=="Obs_obs"
Percent3 <- res[ix3,colnames(res)=="GapPercent"]
mode(Percent3) <- "numeric"

dat_rmse3 <- res[ix3,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse3)

dat_rmse3 <- as.matrix(dat_rmse3)
mode(dat_rmse3)="numeric"
colnames(dat_rmse3)

plot_TDext <- aggregate(dat_rmse3,by=list(Percent3),FUN=mean,na.rm = TRUE)
plot_TDext <- plot_TDext[,colSums(!is.na(plot_TDext))!=0]
plot_TDext <- plot_TDext[complete.cases(plot_TDext),]
plot_TDext <- plot_TDext[plot_TDext[,1]!=0,]

#--------------------------------------
# select TD2env RMSE results of all repetitions
#--------------------------------------
#chose envelope rainfor
ix4=res$TraitChoice=="data_2"&res$Obs_or_Spec=="Obs_obs"
Percent4 <- res[ix4,colnames(res)=="GapPercent"]
mode(Percent4) <- "numeric"

dat_rmse4 <- res[ix4,grep(colnames(res),pattern = pattern_now)]
colnames(dat_rmse4)

dat_rmse4 <- as.matrix(dat_rmse4)
mode(dat_rmse4)="numeric"
colnames(dat_rmse4)

plot_TD2ext <- aggregate(dat_rmse4,by=list(Percent4),FUN=mean,na.rm = TRUE)
plot_TD2ext <- plot_TD2ext[,colSums(!is.na(plot_TD2ext))!=0]
plot_TD2ext <- plot_TD2ext[complete.cases(plot_TD2ext),]
plot_TD2ext <- plot_TD2ext[plot_TD2ext[,1]!=0,]



tst <- data.frame(TDext=plot_TDext[,2],
                  TDtd=plot_TDtd[,2],
                  TD2ext=plot_TD2ext[,2],
                  TD2td=plot_TD2td[,2])
# check general output
plot(tst$TDext,tst$TDtd,xlim=c(0,1),ylim=c(0,1));abline(0,1)
plot(tst$TD2ext,tst$TD2td,xlim=c(0,1),ylim=c(0,1));abline(0,1)
plot(tst$TDtd,tst$TD2td,xlim=c(0,1),ylim=c(0,1));abline(0,1)




#--------------------------------------
# data 1
#--------------------------------------
trait_names=trait_rainfor
colz=c("#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#b2182b")

plot_TDtd[,colnames(plot_TDtd)=="RMSE_zlog_Total"]
plot_TDext[,colnames(plot_TDtd)=="RMSE_zlog_Total"]
plot_TDtd[,colnames(plot_TDtd)=="RMSE_gap_zlog_Total"]
plot_TDext[,colnames(plot_TDtd)=="RMSE_gap_zlog_Total"]

list.files(file.path(origin,"_2021","figures","Figure_1"))

  pdf(file=file.path(origin,"_2021","figures","Figure_1","Fig_1_RMSE_zlog_.pdf"),width=8,height=7)
  layout(matrix(c(1,1,2,2,3,
                  4,4,5,5,6),nrow=2,ncol=5,byrow=TRUE))
  par(mar=c(7,7,1,1))
      i=2
      # Total data TDtd
      {
      plot(plot_TDtd[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,xaxt="n",yaxt="n",
           xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE)
        axis(2,line = 3,at = c(0,.3,.6),labels = c("","RMSE",""),tick = FALSE,cex.axis=2)
        axis(2,line = 1,at = seq(0,.6,.1),labels = seq(0,.6,.1),cex.axis=2)
        
        axis(1,line = 3.5,at = c(40),labels = "Missingness [% gaps]",cex.axis=2,tick=F)
        axis(1,line = 1,at = seq(0,80,20),labels = rep("",length(seq(0,80,20))),cex.axis=2)
        axis(1,line = 1.5,at = seq(0,80,20),labels = seq(0,80,20),cex.axis=2,tick=FALSE)
        
        #        abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
        rect(ybottom = 0,ytop = .05,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .1,ytop = .15,xleft =-1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .2,ytop = .25,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .3,ytop = .35,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .4,ytop = .45,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .5,ytop = .55,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        colnames(plot_TDtd)
        text(10,.58,"TDtd",cex=2)
        
        lines(plot_TDtd[,1],plot_TDtd[,colnames(plot_TDtd)=="RMSE_zlog_Total"],col="black",lwd=3)
        
      i=1
      for(i in 1:length(trait_names)){
        print(trait_names[i])
        lines(plot_TDtd[,1],plot_TDtd[,colnames(plot_TDtd)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=3)
      }  
      
      i=2
      lines(plot_TDtd[,1],plot_TDtd[,colnames(plot_TDtd)=="RMSE_zlog_Total"],col="black",lwd=3)
      
      }
      # Total TDext
      {
        i=2
        # Total data TDext
        plot(plot_TDext[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,xaxt="n",yaxt="n",
             xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE)
        axis(2,line = 3,at = c(0,.3,.6),labels = c("","RMSE",""),tick = FALSE,cex.axis=2)
        axis(2,line = 1,at = seq(0,.6,.1),labels = seq(0,.6,.1),cex.axis=2)
        
        axis(1,line = 3.5,at = c(40),labels = "Missingness [% gaps]",cex.axis=2,tick=F)
        axis(1,line = 1,at = seq(0,80,20),labels = rep("",length(seq(0,80,20))),cex.axis=2)
        axis(1,line = 1.5,at = seq(0,80,20),labels = seq(0,80,20),cex.axis=2,tick=FALSE)
        
        #        abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
        rect(ybottom = 0,ytop = .05,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .1,ytop = .15,xleft =-1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .2,ytop = .25,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .3,ytop = .35,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .4,ytop = .45,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .5,ytop = .55,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        
        text(10,.58,"TDext",cex=2)
        colnames(plot_TDext)
        lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)=="RMSE_zlog_Total"],col="black",lwd=3)
        
        i=1
        for(i in 1:length(trait_names)){
          print(trait_names[i])
          lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)==paste0("RMSE_zlog_",trait_names[i])],col=colz[i],lwd=3)
        }  
        
        i=2
        lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)=="RMSE_zlog_Total"],col="black",lwd=3)
        
      }
      # ----------------------------------------------------------------------------
      par(mar=c(0,0,0,0))
      plot(0:10,0:10,col="white",xaxt="n",yaxt="n",ylab="",xlab="",frame=FALSE)
      level_names = c(trait_names,"","Average")
      legend(1, 10, level_names, col = c(colz[1:length(trait_names)],"white","black"),bty = "n",
             text.col = "black", lty = c(rep(1,length(trait_names)),1,1,1), lwd=4,cex=1.4,
             merge = TRUE, bg = "white")  
      
      
      # ----------------------------------------------------------------------------
      par(mar=c(7,7,1,1))
      #Gaps TDtd
      {
        i=2
        plot(plot_TDtd[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,xaxt="n",yaxt="n",
             xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE)
        axis(2,line = 3,at = c(0,.3,.6),labels = c("","RMSE",""),tick = FALSE,cex.axis=2)
        axis(2,line = 1,at = seq(0,.6,.1),labels = seq(0,.6,.1),cex.axis=2)
        
        axis(1,line = 3.5,at = c(40),labels = "Missingness [% gaps]",cex.axis=2,tick=F)
        axis(1,line = 1,at = seq(0,80,20),labels = rep("",length(seq(0,80,20))),cex.axis=2)
        axis(1,line = 1.5,at = seq(0,80,20),labels = seq(0,80,20),cex.axis=2,tick=FALSE)
        
        
        #        abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
        rect(ybottom = 0,ytop = .05,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .1,ytop = .15,xleft =-1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .2,ytop = .25,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .3,ytop = .35,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .4,ytop = .45,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .5,ytop = .55,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        text(17,.58,"TDtd gaps",cex=2)
        
        colnames(plot_TDtd)
        
        i=1
        for(i in 1:length(trait_names)){
          print(trait_names[i])
          lines(plot_TDtd[,1],plot_TDtd[,colnames(plot_TDtd)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],
                lwd=3,lty=1)
        }  
        
        i=2
        lines(plot_TDtd[,1],plot_TDtd[,colnames(plot_TDtd)=="RMSE_gap_zlog_Total"],col="black",lwd=3,lty=1)
        
        #lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)=="RMSE_gap_zlog_Total"],col=colz[7],lwd=3,lty=5)
        
        
      }
      #Gaps TDext
      {
        i=2
        plot(plot_TDext[,c(1,i)],col=colz[i],pch=16,cex=.1,cex.lab=3,cex.axis=3,xaxt="n",yaxt="n",
             xlab="",ylab="",xlim=c(0,80),ylim=c(0,.6),frame=FALSE)
        axis(2,line = 3,at = c(0,.3,.6),labels = c("","RMSE",""),tick = FALSE,cex.axis=2)
        axis(2,line = 1,at = seq(0,.6,.1),labels = seq(0,.6,.1),cex.axis=2)
        
        axis(1,line = 3.5,at = c(40),labels = "Missingness [% gaps]",cex.axis=2,tick=F)
        axis(1,line = 1,at = seq(0,80,20),labels = rep("",length(seq(0,80,20))),cex.axis=2)
        axis(1,line = 1.5,at = seq(0,80,20),labels = seq(0,80,20),cex.axis=2,tick=FALSE)
        
       
#        abline(h=seq(-1.2,1.2,by=.1),col="gray",lty=2)
        rect(ybottom = 0,ytop = .05,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .1,ytop = .15,xleft =-1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .2,ytop = .25,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .3,ytop = .35,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .4,ytop = .45,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        rect(ybottom = .5,ytop = .55,xleft = -1,xright = 100,col="#EBEBEB",border = NA)
        colnames(plot_TDext)
        text(20,.58,"TDext gaps",cex=2) 
        i=1
        for(i in 1:length(trait_names)){
          print(trait_names[i])
          lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)==paste0("RMSE_gap_zlog_",trait_names[i])],col=colz[i],
                lwd=3,lty=1)
        }  
        
        i=2
        lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)=="RMSE_gap_zlog_Total"],col="black",lwd=3,lty=1)
        
        #lines(plot_TDext[,1],plot_TDext[,colnames(plot_TDext)=="RMSE_gap_zlog_Total"],col=colz[7],lwd=3,lty=5)
        
      }
  dev.off()
  