
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
  
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  

  res_matrix_name="res_20201020"#"res_20201112"
  res_matrix_name= "res_20210303"
  res <- read.table(file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
  res <- as.data.frame(res)
  
  
  TD_choice="Obs_obs_TD"
  test_data="data"
  GapPercent=-1
  ylims=c(-.5,.5)
  nreps=10
  cor_names <- colnames(res)[grep(colnames(res),pattern = "cor")]
  for(TD_choice in TD_choices){
    
    pdf(file=file.path(origin,"figures","figure_3",paste0("figure_3_a_",TD_choice,".pdf")),
        width=10,height=10)
    par(mfrow=c(2,2),mar=c(4,5,3,2))
  
      for(test_data in c("data","data_2")){
      if(test_data=="data_2"){col_now=colz1}
      if(test_data=="data"){col_now=colz}
      i=1
      
      plot(0:10,col="white",ylim=c(ylims),cex.lab=2,ylab="Correlation coeff (pred-observed)",xlab="Gaps filled [%]",xaxt="n")
      axis(side=1,labels = c("obs","0","1","5","10","20","30","40","50","60","70"),at = 1:11,las=2)
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      i=1
      for(i in 1:length(cor_names)){
        res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
                    colnames(res)%in%c("GapPercent",cor_names[i])]
        if(sum(!is.na(res_sub[,2]))>10){
          cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
          cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
          obs_dat=rep(mean(res_sub[which(res_sub$GapPercent==-1)[1:nreps],2],na.rm = TRUE),10)
          dat_now1 <- cbind(obs_dat-obs_dat,
                            res_sub[which(res_sub$GapPercent==0)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==1)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==5)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==10)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==20)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==30)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==40)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==50)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==60)[1:nreps],2]-obs_dat,
                            res_sub[which(res_sub$GapPercent==70)[1:nreps],2]-obs_dat)
          # MEDIAN
          dat_now <- apply(dat_now1,2,FUN=median,na.rm=TRUE)
          x_now=c(1:length(dat_now))
          x_now <- x_now[!is.na(dat_now)]
          dat_now <- dat_now[!is.na(dat_now)]
          lines(x_now,dat_now,
                  ylim=c(ylims),cex.lab=2,
                  col=col_now[5],ylab="",xaxt="n")

        }
      }
      
      plot(0:11,col="white",ylim=c(0,ylims[2]),cex.lab=2,ylab="Range correlation coeff",
             xlab="Gaps filled [%]",xaxt="n")
        axis(side=1,labels = c("obs","0","1","5","10","20","30","40","50","60","70"),at = 1:11,las=2)
        abline(h=seq(-1,1,by=.2),col="darkgray",lty=2)
        abline(h=seq(-1.1,1,by=.2),col="lightgray",lty=2)
        abline(h=0)
        i=1 
        for(i in 1:length(cor_names)){
          res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
                      colnames(res)%in%c("GapPercent",cor_names[i])]
          if(sum(!is.na(res_sub[,2]))>10){
            cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
            cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
            dat_now1 <- cbind(res_sub[which(res_sub$GapPercent==-1)[1:nreps],2],
                              res_sub[which(res_sub$GapPercent==0)[1:nreps],2],
                              res_sub[which(res_sub$GapPercent==1)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==5)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==10)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==20)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==30)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==40)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==50)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==60)[1:nreps],2] ,
                              res_sub[which(res_sub$GapPercent==70)[1:nreps],2] )

            # Range
            dat_now <- rep(NA,ncol(dat_now1))
            for(d in 1:ncol(dat_now1)){
              dat_now[d] <- quantile(x = dat_now1[,d],probs = .90,na.rm = TRUE)-quantile(x = dat_now1[,d],probs = .10,na.rm = TRUE)
            }        
            x_now=c(1:length(dat_now))
            x_now <- x_now[dat_now!=0&!is.na(dat_now)]
            dat_now <- dat_now[dat_now!=0&!is.na(dat_now)]
            lines(x_now,dat_now,
                  ylim=c(ylims),cex.lab=2,
                  col=col_now[5],ylab="",xaxt="n")
          }
          
        }
        
      }
    dev.off()
    
    
  }
  
  
  # against the observed
  TD_choice="Obs_obs_TD"
  test_data="guido"
  GapPercent=-1
  ylims=c(-.5,.5)
  cor_names <- colnames(res)[grep(colnames(res),pattern = "cor")]
  which(cor_names%in%"cor_LeafN_SLA")
  i=34
  for(TD_choice in TD_choices){
    pdf(file=file.path(origin,"figures","figure_3",paste0("figure_3_b_",TD_choice,".pdf")),width=15,height=10)
    par(mfrow=c(5,10),mar=c(3,2,0,1))
    
    
    for(test_data in c("guido","rainfor")){
      if(test_data=="guido"){col_now=colz1}
      if(test_data=="rainfor"){col_now=colz}
      i=1
      for(i in 1:length(cor_names)){
        res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
                    colnames(res)%in%c("GapPercent",cor_names[i])]
        if(sum(!is.na(res_sub[,2]))>10){
          cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
          cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
          obs_dat=rep(mean(res_sub[which(res_sub$GapPercent==-1)[1:10],2],na.rm = TRUE),10)
          print(res_sub[res_sub$GapPercent==-1,2])
          boxplot(cbind(res_sub[which(res_sub$GapPercent==0)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==1)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==5)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==10)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==20)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==30)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==40)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==50)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==60)[1:10],2]-obs_dat,
                        res_sub[which(res_sub$GapPercent==70)[1:10],2]-obs_dat),
                  ylim=c(ylims),cex.lab=2,
                  col=col_now[5],ylab="",xaxt="n")
          text(5,0.4,cor_name_now)
          axis(side=1,labels = c("obs","1","5","10","20","30","40","50","60","70"),at = 1:10,las=2)
          abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
          abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
          abline(h=0)
          #abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col=col_now[8])
        }
      }
    }
    dev.off()
    
    
  }
  
  
  
  # Into the supplementary
  TD_choice="Obs_obs_TD"
  test_data="guido"
  GapPercent=-1
  ylims=c(-.8,.8)
  cor_names <- colnames(res)[grep(colnames(res),pattern = "cor")]
  i=1
  for(TD_choice in TD_choices){
    pdf(file=file.path(origin,"figures","figure_3",paste0("figure_3_c_",TD_choice,".pdf")),width=15,height=10)
    par(mfrow=c(5,10),mar=c(3,2,0,1))
    
    
    for(test_data in c("guido","rainfor")){
      if(test_data=="guido"){col_now=colz1}
      if(test_data=="rainfor"){col_now=colz}
      i=1
      for(i in 1:length(cor_names)){
        res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
                    colnames(res)%in%c("GapPercent",cor_names[i])]
        if(sum(!is.na(res_sub[,2]))>10){
          cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
          cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
          boxplot(cbind(res_sub[res_sub$GapPercent==0,2],
                        res_sub[res_sub$GapPercent==1,2],
                        res_sub[res_sub$GapPercent==5,2],
                        res_sub[res_sub$GapPercent==10,2],
                        res_sub[res_sub$GapPercent==20,2],
                        res_sub[res_sub$GapPercent==30,2],
                        res_sub[res_sub$GapPercent==40,2],
                        res_sub[res_sub$GapPercent==50,2],
                        res_sub[res_sub$GapPercent==60,2],
                        res_sub[res_sub$GapPercent==70,2]),ylim=c(ylims),cex.lab=2,
                  col=col_now[5],ylab="",xaxt="n")
          text(5,0.65,cor_name_now)
          axis(side=1,labels = c("obs","1","5","10","20","30","40","50","60","70"),at = 1:10,las=2)
          abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
          abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
          abline(h=0)
          abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col=col_now[8])
        }
      }
    }
    dev.off()
    
    
  }
  
  
  
  # Load observed correlations
  TD_choice="Obs_obs_TD"
  test_data="guido"
  GapPercent=-1
  ylims=c(-.7,.7)
  cor_names <- colnames(res)[grep(colnames(res),pattern = "cor")]
  i=1
  for(TD_choice in TD_choices){
  pdf(file=file.path(origin,"figures","figure_3",paste0("figure_3_d_",TD_choice,".pdf")))
  par(mfrow=c(5,11),mar=c(0,0,2,0))
  
  
  for(testdata in c("guido","rainfor")){
  for(i in 1:length(cor_names)){
    res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
                colnames(res)%in%c("GapPercent",cor_names[i])]
    if(sum(!is.na(res_sub[,2]))>10){
      cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
      cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
      plot(rnorm(100),col="white",xaxt="n",yaxt="n",xlim=c(-2,2),ylim=c(-2,2))
      text(-1.6,0,test_data,srt=90)
      text(-.9,0,TD_choice,srt=90)
      text(0.1,0,cor_name_now,srt=90)
      boxplot(res_sub[res_sub$GapPercent==0,2],ylim=c(ylims),cex.lab=2,main="observed")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==1,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="1%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==5,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="5%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==10,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="10%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==20,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="20%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==30,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="30%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==40,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="40%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==50,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="50%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==60,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="60%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
      boxplot(res_sub[res_sub$GapPercent==70,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="70%")
      abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
      abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
      abline(h=0)
      abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    }
  }
  }
  dev.off()
  
  
  }
  
    # Load observed correlations
  TD_choice="Obs_obs_TD"
  test_data="guido"
  GapPercent=-1
  ylims=c(-.7,.7)
  cor_names <- colnames(res)[grep(colnames(res),pattern = "cor")]
  i=1
  pdf(file=file.path(origin,"figure","figure_3","figure_3_e.pdf"))
  par(mfrow=c(5,11),mar=c(0,0,2,0))
  
  for(i in 1:length(cor_names)){
  res_sub=res[res$Obs_or_Spec==TD_choice&res$TraitChoice==test_data,
              colnames(res)%in%c("GapPercent",cor_names[i])]
  if(sum(!is.na(res_sub[,2]))>10){
    cor_name_now <- gsub(cor_names[i],pattern = "cor_",replacement = "")
    cor_name_now <- gsub(cor_name_now,pattern = "_",replacement = " ")
    plot(rnorm(100),col="white",xaxt="n",yaxt="n",xlim=c(-2,2),ylim=c(-2,2))
    text(-1.6,0,test_data,srt=90)
    text(-.9,0,TD_choice,srt=90)
    text(0.1,0,cor_name_now,srt=90)
    boxplot(res_sub[res_sub$GapPercent==0,2],ylim=c(ylims),cex.lab=2,main="observed")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==1,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="1%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==5,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="5%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==10,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="10%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==20,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="20%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==30,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="30%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==40,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="40%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==50,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="50%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==60,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="60%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
    boxplot(res_sub[res_sub$GapPercent==70,2],ylim=c(ylims),yaxt="n",cex.lab=2,main="70%")
    abline(h=seq(-1,1,by=.4),col="darkgray",lty=2)
    abline(h=seq(-1.2,1,by=.4),col="lightgray",lty=2)
    abline(h=0)
    abline(h=mean(res_sub[res_sub$GapPercent==0,2],na.rm = TRUE),col="green")
  }
  }
  
  dev.off()
  }
  
  