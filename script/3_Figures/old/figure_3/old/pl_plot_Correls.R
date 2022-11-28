
plot_Correls <-function(){

  ##### Master Script 
  if(.Platform$OS.type=="windows"){
    origin ="C:/Profiles/traveller/Documents/TMP/2016_Cluster"
  } else{
    origin =  "/Net/Groups/BGI/work_1/2016_GapFilling"
  }
  res <- read.table(file.path(origin,"runs","META","results_TOTAL.csv"),sep=",",dec=".")
  
  
  #res <- read.table(file.path(origin,"runs","META","results_matrix.csv"),sep=",",dec=".")
  #res <- read.table(file.path(origin,"runs","META","results_TOTAL.csv"),sep=",",dec=".")
  
  correl.cols = c(23,25:33)#c(23,25:which(colnames(res)=="LDMC_LeafArea"),which(colnames(res)=="SLA_ConduitDens"):which(colnames(res)=="SLen_DispULen"))
  Correls <- res[,correl.cols]
  Correls <- Correls[,colSums(!is.na(Correls))!=length(correl.cols)]
  # Correls <- Correls[,-(which(colnames(Correls)=="mean_r2"))]
  colz = rainbow(ncol(Correls))
  
  pe2 <- which(colnames(res)=="GapPercent")
  pe <- which(colnames(res)=="act.gaps")

 
  for(td in 1:4){
   for(ts in 1:2){
    
    TD.choice = TD.choices[td]
    tsub = tsubs[ts]
    
    Correls2 <- as.matrix(res[res[,2]==tsub &
                        res[,3] == TD.choice,
                      which(colnames(res)%in%colnames(Correls))])
    Correls2 <- Correls2[,colSums(!is.na(Correls2))!=0]
    
    if(!is.null(dim(Correls2))){
    Percent <- as.numeric(res[res[,2]==tsub&res[,3]==TD.choice,pe])
  
    pdf(file=file.path(origin,"plots","Correlation",paste0("Correls_",TD.choice,"_",tsub,".pdf")), height = 8, width = 10)
    par(mfrow=c(6,ncol(Correls)/5),mar=c(0,0,2,0))
    
    miny <- 0 #min(na.omit(PCA1,PCA2,PCA3)) *1.2
    maxy <- 1 #max(na.omit(PCA1,PCA2,PCA3)) *1.2
    CorTOT = rep(NA,length(Percent))
    
    for(i in 1:ncol(Correls2)){
      
      Correl_now <- res[res[,2]==tsub &
                          res[,3] == TD.choice,
                        which(colnames(res)==colnames(Correls2)[i])]
      
      mode(Correl_now) <- "numeric"
      
      if(sum(!is.na(Correl_now))>=2){
        
            plot(Percent,Correl_now,col=colz[i],pch=16,cex=.5,yaxt='n',ylim=c(miny,maxy),xaxt="n",
           xlab="% Gaps",ylab="% of explained variance")
            mtext(side=3, colnames(Correls2)[i],cex=.5)
            print("hi")
          CorTOT <- cbind(CorTOT,Correl_now)
      }
    }     
    
    if(sum(colnames(res)=="mean_r2")!=1){
      new.col = rep(NA,nrow(res))
      res <- cbind(res, new.col)
      colnames(res)[ncol(res)] == "mean_r2"
    }
    
    if(!is.null(ncol(CorTOT))){
      res[res[,2]==tsub&res[,3]==TD.choice,which(colnames(res)=="mean_r2")]<- rowMeans(CorTOT[,2:ncol(CorTOT)])
    
      plot(Percent,as.numeric(res[res[,2]==tsub&res[,3]==TD.choice,which(colnames(res)=="mean_r2")]),ylim=c(0,.5), yaxt='n',ylab="mean_r2")
      
      dev.off()  
    }
      
      Percent2 <- as.numeric(res[res[,2]==tsub&res[,3]==TD.choice,pe2])
      
      
      if(!is.null(ncol(CorTOT))){
        
      pdf(file.path(origin,"plots","Gaps",paste0("MeanR2_",TD.choice,"_",tsub,".pdf")), height = 4, width = 4)
      par(mfrow=c(1,1),mar=c(5,6,3,2))
        
      meanC <- rowMeans(CorTOT[,2:ncol(CorTOT)])
      dat <- cbind(Percent2,meanC)
      
      miny = min(meanC[!is.na(meanC)])
      
      boxplot(meanC~Percent2, data=dat,ylab="mean r2",xlab="% gaps", ylim=c(miny,miny+.2),
              main=paste0("Mean r2 ",TD.choice,"_",tsub),cex.axis=0.9)
      dev.off()
      }
    }  
    }
  }




  
  
}
