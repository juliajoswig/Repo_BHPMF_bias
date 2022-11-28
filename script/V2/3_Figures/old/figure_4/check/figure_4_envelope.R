figure_envelope <- function(){

  colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858","black","black","black","black")
  colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000","black","black","black","black")
  col_back=c("#7b3294","#c2a5cf","#a6dba0","#008837")
  col_back=c("#f4a582","#fddbc7","#f7f7f7","#d1e5f0")
  col_back=c("#ffffe5","#f7fcb9","#d9f0a3","#addd8e")
  col_back=c("#e5f5e0","#f7fcf5","#c7e9c0","#a1d99b")
  
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
  Version_now="V3"
  list.files(file.path(origin,"script","analysis",Version_now))
  
  #------------------------------------------------------------
  # load some functions
  #------------------------------------------------------------
  source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
  load_functions()
  
  #------------------------------------------------------------
  # define data set approaches/choices
  #------------------------------------------------------------
  out <- choices()
  tsubs <- out$tsubs
  TD_choices =  out$TD_choices
  repnums = out$repnums
  gappercents = out$gappercents
  whichDataSet = out$whichDataSet
  ObsSpec = out$ObsSpec
  obsspec = ObsSpec
  preparation = out$preparation
  trait_guido =out$trait_guido
  trait_rainfor =out$trait_rainfor
#  colz1 =out$colz1
#  colz2 =out$colz2
  Indeces <- c("rmse_man")

  
  res_matrix_name="res_20201020"#"res_20201112"
  res_matrix_name="res_20201126"
  res <- read.table(file.path(origin,"runs","META",paste0(res_matrix_name,".csv")),sep=",",dec=".")
  res <- as.data.frame(res)
  head(res)
  i=1
  j=1
  groups=c("Spec","Gen","Fam","PG","GF","PFT")
  g=1
  for(i in c(1,2)){
  ObsSpec_nowTD <- TD_choices[i]
  ObsSpec_now <- TD_choices[i+2]
  
  pdf(file=file.path(origin,"figures","figure_4",paste0("figure_4_d",ObsSpec_now,".pdf")),width=10,height=3.5)
  par(mfrow=c(2,6),mar=c(4,4,2,1))
  
  for(j in 1:2){
    test_data_now <- tsubs[j]
    if(test_data_now=="guido"){col_now=colz1[5]}
    if(test_data_now=="rainfor"){col_now=colz[5]}
    ## distance to observed test data vs distance to observed envelope data
  res_TD=res[res$Obs_or_Spec==ObsSpec_nowTD&res$TraitChoice==test_data_now,]
  res_evlp=res[res$Obs_or_Spec==ObsSpec_now&res$TraitChoice==test_data_now,]

## 
  for(g in 1:length(groups)){
  dat_TD <- res_TD[,c(4,grep(colnames(res_TD),pattern = paste0("Sil_",groups[g])))]   
  dat_obs_TD <- dat_TD
  dat_obs_1 <- colMeans(dat_TD[dat_TD$GapPercent==-1,],na.rm=TRUE)
  for(i in 1:length(dat_obs_1)){dat_obs_TD[,i] <- rep(dat_obs_1[i],nrow(dat_obs_TD))}
  dat_TD_d <- dat_TD - dat_obs_TD

  dat_evlp <- res_evlp[,c(4,grep(colnames(res_evlp),pattern = paste0("Sil_",groups[g])))]   
  dat_obs_evlp <- dat_evlp
  dat_obs_1 <- colMeans(dat_evlp[dat_evlp$GapPercent==-1,],na.rm=TRUE)
  for(i in 1:length(dat_obs_1)){dat_obs_evlp[,i] <- rep(dat_obs_1[i],nrow(dat_obs_evlp))}
  dat_evlp_d <- dat_evlp - dat_obs_evlp
  
  
  if(test_data_now=="guido"){test_data_name="data 1"}
  if(test_data_now=="rainfor"){test_data_name="data 2"}
  
  plot(dat_TD_d[,c(3)],dat_evlp_d[,c(3)],xlim=c(-.4,.4),ylim=c(-.4,.4),
       ylab="Envelope sil. deviation",xlab="Subset sil. deviation",
       pch=16,main=paste0(groups[g]," ",test_data_name))
  rect(xleft = -2,ybottom = -2,xright = 0,ytop = 0,col=col_back[2],border = col_back[2])
  rect(xleft = 0,ybottom = -2,xright = 2,ytop = 0,col=col_back[1],border = col_back[1])
  rect(xleft = -2,ybottom = 0,xright = 0,ytop = 2,col=col_back[3],border = col_back[3])
  rect(xleft = 0,ybottom = 0,xright = 2,ytop = 2,col=col_back[4],border = col_back[4])
  abline(0,1,col="white",lwd=3)
  points(dat_TD_d[,c(3)],dat_evlp_d[,c(3)],xlim=c(-1,1),ylim=c(-1,1),pch=16,col="black",cex=1.6)
  points(dat_TD_d[,c(3)],dat_evlp_d[,c(3)],xlim=c(-1,1),ylim=c(-1,1),pch=16,col=col_now,cex=1)
}  
  }
  dev.off()
 }
  
  
  pdf(file=file.path(origin,"figures","figure_4","figure_4_e.pdf"),width=8,height=4)
  par(mfrow=c(1,2))
  i=1
  for(i in 1){
    ObsSpec_nowTD <- TD_choices[i]
    ObsSpec_now <- TD_choices[i+2]
    for(j in 1:2){
      test_data_now <- tsubs[j]
      ## distance to observed test data vs distance to observed envelope data
      res_TD=res[res$Obs_or_Spec==ObsSpec_nowTD&res$TraitChoice==test_data_now,]
      res_evlp=res[res$Obs_or_Spec==ObsSpec_now&res$TraitChoice==test_data_now,]
      
      ## 
        dat_TD <- res_TD[,c(4,grep(colnames(res_TD),pattern = "cor_"))]   
        dat_obs_TD <- dat_TD
        dat_obs_1 <- colMeans(dat_TD[dat_TD$GapPercent==-1,],na.rm=TRUE)
        for(k in 1:length(dat_obs_1)){dat_obs_TD[,k] <- rep(dat_obs_1[k],nrow(dat_obs_TD))}
        dat_TD_d <- dat_TD - dat_obs_TD
        
        dat_evlp <- res_evlp[,c(4,grep(colnames(res_evlp),pattern = "cor_"))]   
        dat_obs_evlp <- dat_evlp
        dat_obs_1 <- colMeans(dat_evlp[dat_evlp$GapPercent==-1,],na.rm=TRUE)
        for(k in 1:length(dat_obs_1)){dat_obs_evlp[,k] <- rep(dat_obs_1[k],nrow(dat_obs_evlp))}
        dat_evlp_d <- dat_evlp - dat_obs_evlp
    
        if(test_data_now=="guido"){test_data_name="data 1";col_now=colz1}
        if(test_data_now=="rainfor"){test_data_name="data 2";col_now=colz}
        
        plot(1:10,xlim=c(-.5,.5),ylim=c(-.5,.5),
             main=paste0("Correlations ","and envelope"," ",test_data_name),col="white",
             ylab="Envelope deviation",xlab="Subset data deviation")
        rect(xleft = -2,ybottom = -2,xright = 0,ytop = 0,col=col_back[2],border = col_back[2])
        rect(xleft = 0,ybottom = -2,xright = 2,ytop = 0,col=col_back[1],border = col_back[1])
        rect(xleft = -2,ybottom = 0,xright = 0,ytop = 2,col=col_back[3],border = col_back[3])
        rect(xleft = 0,ybottom = 0,xright = 2,ytop = 2,col=col_back[4],border = col_back[4])
        abline(0,1,col="white",lwd=3)
        cn=2
        for(cn in 2:ncol(dat_TD_d)){
          tmp=cbind(dat_TD_d[,cn],dat_evlp_d[,cn])
          if(sum(complete.cases(tmp))>10){
            points(dat_TD_d[,cn],dat_evlp_d[,cn],pch=16,main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1)
            points(dat_TD_d[,cn],dat_evlp_d[,cn],pch=16,main=paste0(groups[g],ObsSpec_now,test_data_now),cex=.8,col=col_now)
          }
      }  
    }
  }
  dev.off()
  
  i=1
  for(i in 1:2){
    ObsSpec_nowTD <- TD_choices[i]
    ObsSpec_now <- TD_choices[i+2]
pdf(file=file.path(origin,"figures","figure_4",paste0("figure_4_f_",ObsSpec_now,".pdf")),width=8,height=4)
par(mfrow=c(1,2))
    j=2
    for(j in 1:2){
      test_data_now <- tsubs[j]
      ## distance to observed test data vs distance to observed envelope data
      res_TD=res[res$Obs_or_Spec==ObsSpec_nowTD&res$TraitChoice==test_data_now,]
      res_evlp=res[res$Obs_or_Spec==ObsSpec_now&res$TraitChoice==test_data_now,]
      
      ## 
      dat_TD <- res_TD[,c(4,grep(colnames(res_TD),pattern = "RMSE_"))]   
      dat_TD <- dat_TD[,-grep(colnames(dat_TD),pattern = "gap")]   
      dat_obs_TD <- dat_TD
      dat_TD_d <- dat_obs_TD
      
      dat_evlp <- res_evlp[,c(4,grep(colnames(res_evlp),pattern = "RMSE_"))]   
      dat_evlp <- dat_evlp[,-grep(colnames(dat_evlp),pattern = "gap")]   
      dat_evlp_d <- dat_evlp
      
      if(test_data_now=="guido"){test_data_name="data 1";col_now=colz1}
      if(test_data_now=="rainfor"){test_data_name="data 2";col_now=colz}
      
      plot(1:10,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
           main=paste0("RMSE ","and envelope"," ",test_data_name),col="white",
           ylab="Envelope deviation",xlab="Subset data deviation")
      rect(xleft = -2,ybottom = -2,xright = 0,ytop = 0,col=col_back[2],border = col_back[2])
      rect(xleft = 0,ybottom = -2,xright = 2,ytop = 0,col=col_back[1],border = col_back[1])
      rect(xleft = -2,ybottom = 0,xright = 0,ytop = 2,col=col_back[3],border = col_back[3])
      rect(xleft = 0,ybottom = 0,xright = 2,ytop = 2,col=col_back[4],border = col_back[4])
      abline(0,1,col="white",lwd=3)
      cn=3
      for(cn in 2:length(gappercents)){
        tmp_td=dat_TD_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_TD_d)]
        tmp_e=dat_evlp_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_evlp_d)]
        summary(tmp_e)
        ct=1
        for(ct in 1:ncol(tmp_td)){
          tmp_trait <- cbind(tmp_td[,ct],tmp_e[,ct])
        if(sum(complete.cases(tmp_trait))>0){
          points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                 main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1.5)
          points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                 main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1,col=col_now[cn])
        }
      }  
    }
  }
  dev.off()
  }

  
  
  
  for(i in 1:2){
    ObsSpec_nowTD <- TD_choices[i]
    ObsSpec_now <- TD_choices[i+2]
    pdf(file=file.path(origin,"figures","figure_4",paste0("figure_4_g_",ObsSpec_now,".pdf")),
        width=10,height=3)
    par(mfrow=c(1,4))
    j=2
    for(j in 1:2){
      test_data_now <- tsubs[j]
      ## distance to observed test data vs distance to observed envelope data
      res_TD=res[res$Obs_or_Spec==ObsSpec_nowTD&res$TraitChoice==test_data_now,]
      res_evlp=res[res$Obs_or_Spec==ObsSpec_now&res$TraitChoice==test_data_now,]
      
      ## 
      dat_TD <- res_TD[,c(4,grep(colnames(res_TD),pattern = "RMSE_"))]   
      dat_TD <- dat_TD[,-grep(colnames(dat_TD),pattern = "gap")]   
      dat_obs_TD <- dat_TD
      dat_TD_d <- dat_obs_TD
      
      dat_evlp <- res_evlp[,c(4,grep(colnames(res_evlp),pattern = "RMSE_"))]   
      dat_evlp <- dat_evlp[,-grep(colnames(dat_evlp),pattern = "gap")]   
      dat_evlp_d <- dat_evlp
      
      if(test_data_now=="guido"){test_data_name="data 1";col_now=colz1}
      if(test_data_now=="rainfor"){test_data_name="data 2";col_now=colz}
      
      plot(1:10,xlim=c(-1.3,1.3),ylim=c(-1.3,1.3),
           main=paste0("RMSE all envelope"," ",test_data_name),col="white",
           ylab="Envelope data",xlab="Subset data")
      rect(xleft = -2,ybottom = -2,xright = 0,ytop = 0,col=col_back[2],border = col_back[2])
      rect(xleft = 0,ybottom = -2,xright = 2,ytop = 0,col=col_back[1],border = col_back[1])
      rect(xleft = -2,ybottom = 0,xright = 0,ytop = 2,col=col_back[3],border = col_back[3])
      rect(xleft = 0,ybottom = 0,xright = 2,ytop = 2,col=col_back[4],border = col_back[4])
      abline(0,1,col="white",lwd=3)
      cn=3
      for(cn in 2:length(gappercents)){
        tmp_td=dat_TD_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_TD_d)]
        tmp_e=dat_evlp_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_evlp_d)]
        summary(tmp_e)
        ct=1
        for(ct in 1:ncol(tmp_td)){
          tmp_trait <- cbind(tmp_td[,ct],tmp_e[,ct])
          if(sum(complete.cases(tmp_trait))>0){
            points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                   main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1.5)
            points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                   main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1,col=col_now[cn])
          }
        }  
      }
      
      
      test_data_now <- tsubs[j]
      ## distance to observed test data vs distance to observed envelope data
      res_TD=res[res$Obs_or_Spec==ObsSpec_nowTD&res$TraitChoice==test_data_now,]
      res_evlp=res[res$Obs_or_Spec==ObsSpec_now&res$TraitChoice==test_data_now,]
      
      ## 
      dat_TD <- res_TD[,c(4,grep(colnames(res_TD),pattern = "RMSE_gap"))]   
      dat_obs_TD <- dat_TD
      dat_TD_d <- dat_obs_TD
      
      dat_evlp <- res_evlp[,c(4,grep(colnames(res_evlp),pattern = "RMSE_gap"))]   
      dat_evlp_d <- dat_evlp
      
      if(test_data_now=="guido"){test_data_name="data 1";col_now=colz1}
      if(test_data_now=="rainfor"){test_data_name="data 2";col_now=colz}
      
      plot(1:10,xlim=c(-1.5,1.5),ylim=c(-1.5,1.5),
           main=paste0("RMSE filled envelope"," ",test_data_name),col="white",
           ylab="Envelope data",xlab="Subset data")
      rect(xleft = -2,ybottom = -2,xright = 0,ytop = 0,col=col_back[2],border = col_back[2])
      rect(xleft = 0,ybottom = -2,xright = 2,ytop = 0,col=col_back[1],border = col_back[1])
      rect(xleft = -2,ybottom = 0,xright = 0,ytop = 2,col=col_back[3],border = col_back[3])
      rect(xleft = 0,ybottom = 0,xright = 2,ytop = 2,col=col_back[4],border = col_back[4])
      abline(0,1,col="white",lwd=3)
      cn=3
      for(cn in 2:length(gappercents)){
        tmp_td=dat_TD_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_TD_d)]
        tmp_e=dat_evlp_d[dat_TD_d[,1]==gappercents[cn],2:ncol(dat_evlp_d)]
        summary(tmp_e)
        ct=1
        for(ct in 1:ncol(tmp_td)){
          tmp_trait <- cbind(tmp_td[,ct],tmp_e[,ct])
          if(sum(complete.cases(tmp_trait))>0){
            points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                   main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1.5)
            points(dat_TD_d[cn,],dat_evlp_d[cn,],pch=16,
                   main=paste0(groups[g],ObsSpec_now,test_data_now),cex=1,col=col_now[cn])
          }
        }  
      }
      
      
    }
    dev.off()
    
  }
}
