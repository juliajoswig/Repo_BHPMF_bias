# atm S7: Relative error to trait range for test data 2. Error (absolute) calculated
# as imputed - observed for test data imputed with and without the envelope (back-transformed data 

#------------------------------------------------------------
# define path
#------------------------------------------------------------
setwd("/..")
origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git"
originData = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data"
list.files(file.path(origin,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","helper_scripts","fn_load_functions.R"))
load_functions(origin)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)
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
  GapPercent=50
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  
  
#-------------------------------------------------------------------
# chose trait data   
#-------------------------------------------------------------------
t_choice="data_2"
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}
if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}
if(t_choice=="data"){trait_names=trait_rainfor}
if(t_choice=="data_2"){trait_names=trait_guido}

txs=c("indiv","spec","gen","fam","clad")
tx=0 #(0=indiv,1=spec,2=genus,3=family,4=clades)
for(tx in 0:4){
  tx=tx+1
  par(mfrow=c(4,4))
  {
    
    #-------------------------------------------------------------------
    # load trait data   
    #-------------------------------------------------------------------
#    file.path(originData,"analyes","Point_wise","res.csv")
    colz=colz1
    #load data
  {
      # load Envelope data
      # load TDenvelope
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),"Obs_obs_TD"))
      # load TD data
      # OBSERVED
      TD <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                             "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
      TD_sparse <- as.data.frame(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_80"),
                                                    "Obs_obs_TD","data","traitInfo.csv"),header=TRUE))[,-c(1,2)]
      TD_tax <- as.data.frame(read.table(file = file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_0"),
                                                         "Obs_obs_TD","data","taxInfo.csv"),
                                        sep=",",col.names = c("ID","Species","Genus","Family","Clade")))
      ID_TD <- TD_tax[,1]
      head(TD_tax)
      dim(TD_tax)
      dim(TD)
      dim(TD_sparse)
      plot(TD[,1],TD_sparse[,1])
      abline(0,1)
      
      list.files(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),"Obs_obs_TD","data"))
      
      # predicted 
      TDtd <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                           "Obs_obs","data","traitInfoTD_pred.csv")))[,-c(1,2)]
      TDenv <- as.matrix(read.csv(file.path(originData,"_runs",paste0("Rep_",RepNum),t_choice,paste0("p_","80"),
                                            "Obs_obs_TD","data","traitInfoTD_pred.csv")))[,-c(1,2)]
      head(TDtd)
      head(TDenv)
      
      plot(TD[,1],TDtd[,1])
      abline(0,1)
      plot(TD[,1],TDenv[,1])
      abline(0,1)
      plot(TD_sparse[,1],TD[,1])
      abline(0,1)
    }
    
    trait_names=colnames(TD)
    
    mean.fun <- function(input){
      if(sum(!is.na(input))==0){return(NA)}else{
        return(mean(input[!is.na(input)]))}
    }
    count.fun <- function(input){
      return(sum(!is.na(input)))
    }
    range.fun <- function(input){
      if(sum(!is.na(input))==0){return(NA)}else{
        return(max(input[!is.na(input)])-min(input[!is.na(input)]))
      }
    }
    
    
    
    #_________________________________________________
    # zlog based on TD
    #-------------------------------------------------
    {
      res <- list()
      res$TDtd <- TDtd
      res$TD <- TD
      res$TDenv <- TDenv
      res$TD_sparse <- TD_sparse
      #res$Env_sparse_tx <- Env_sparse_tx
      for(t in 1:length(trait_names)){
        whichcol=which(colnames(TD_tx)==trait_names[t])
        
        TDmn=mean(log(TD[,whichcol]),na.rm = TRUE)
        TDsd=sd(log(TD[,whichcol]),na.rm = TRUE)
        
        # zlog transform on the basis of TD        
        res$TD[,whichcol] <-         (log(res$TD[,whichcol])        -TDmn)/TDsd #zlog TD
        res$TD_sparse[,whichcol] <-  (log(res$TD_sparse[,whichcol]) -TDmn)/TDsd #zlog TD_sparse
        res$TDtd[,whichcol] <-       (log(res$TDtd[,whichcol])      -TDmn)/TDsd #zlog TDtd
        res$TDenv[,whichcol] <-      (log(res$TDenv[,whichcol])     -TDmn)/TDsd #zlog TDenv
        #res$Env_sparse_tx[,whichcol] <- (log(Env_sparse_tx[,whichcol])-TDmn)/TDsd #zlog TD_sparse
      }
      
      par(mfrow=c(1,2))
      plot(res$TD[,1],res$TDtd[,1])
      abline(0,1)
      plot(TD[,1],TDtd[,1])
      abline(0,1)
      
      plot(res$TD[,1],res$TDenv[,1])
      abline(0,1)
      plot(TD[,1],TDenv[,1])
      abline(0,1)
      
      plot(res$TD_sparse[,1],res$TD[,1])
      abline(0,1)
      plot(TD_sparse[,1],TD[,1])
      abline(0,1)
      
    }
    
    
    
    TDtd_zlog <- res$TDtd
    TD_zlog <- res$TD
    TDenv_zlog <- res$TDenv
    TD_sparse_zlog <- res$TD_sparse
    #Env_sparse_now <- res$Env_sparse_tx
    
    
    
    
    
    #_________________________________________________
    # aggregate TD to taxon mean
    #-------------------------------------------------
    TD_tx <- aggregate(TD_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TD_sparse_tx <- aggregate(TD_sparse_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDtd_tx <- aggregate(TDtd_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    TDenv_tx <- aggregate(TDenv_zlog,by=list(TD_tax[,tx]),FUN = mean.fun)
    #Env_sparse_tx <- aggregate(Env_sp,by=list(Env_tax_c[,tx]),FUN = mean.fun)
    
    TD_tx <- TD_tx[,-which(colnames(TD_tx)=="Group.1")]
    TD_sparse_tx <- TD_sparse_tx[,-which(colnames(TD_sparse_tx)=="Group.1")]
    TDtd_tx <- TDtd_tx[,-which(colnames(TDtd_tx)=="Group.1")]
    TDenv_tx <- TDenv_tx[,-which(colnames(TDenv_tx)=="Group.1")]
    #Env_sparse_tx <- Env_sparse_tx[,-which(colnames(Env_sparse_tx)=="Group.1")]
    
    
    
 #------------------------------------------------------------------------
      # residuals of TDtd to TD: 
      bxplTDtd_TD = TDtd_tx - TD_tx
      summary(bxplTDtd_TD)
      #        plot(TDtd_now,TD_now);abline(0,1)
      # residuals of TDtd to TDsparse: 
      bxplTDtd_TDsparse=TDtd_tx- TD_sparse_tx
      summary(bxplTDtd_TDsparse)
      #        plot(TDtd_now,TD_sparse_now);abline(0,1)
      # residuals of TDtd to EXT: 
#      dim(TD_tax$ID%in%Env_tax_c$ID)
#      bxplTDtd_EXT=TDtd_now- Env_sparse_now[TD_tax$ID%in%Env_tax_c$ID,]
#      #        plot(TDtd_now,Env_sparse_now);abline(0,1)
#      # residuals of TDenv to EXT: 
      bxplTDenv_TDsparse=TDenv_tx - TD_sparse_tx
      summary(bxplTDenv_TDsparse)
      #        plot(TDenv_now,Env_sparse_now);abline(0,1)
      # residuals of TDenv to TD: 
      bxplTDenv_TD = TDenv_tx - TD_tx
      summary(bxplTDenv_TD)
      #        plot(TD_now,TDenv_now);abline(0,1)
    
      
    
try(dev.off())
    par(mfrow=c(1,1),mar=c(10,4,1,1))
    pdf(file=file.path(origin,"figures","Figure_2","boxplots",paste0(txs[tx],"_boxplots_",t_choice,".pdf")),width=8,height=2.3,pointsize = 9)
    par(mfrow=c(1,6),mar=c(12,4,1,1))
    t=1

    for(t in 1:ncol(TD)){
      head(bxplTDenv_TD)

      # (bxplTDenv_EXT),  (bxplTDenv_TD))[seq(from = t,to = ncol(TD)*5,by = ncol(TD))]
      dat_plot= cbind(      (bxplTDtd_TD),        (bxplTDtd_TDsparse),      (bxplTDenv_TD),  (bxplTDenv_TDsparse))[seq(from = t,to = ncol(TD)*4,by = ncol(TD))]
      #boxplot(cbind(rnorm(mean=0.01,sd=.1,n = 100),rnorm(mean=10,sd=.1,n = 100),rnorm(mean=20,sd=.1,n = 100)))
      colnames(dat_plot)=  c("IMPobs - OBS",        "IMPobs - OSBsparse",  "IMPobsExt - OBS",  "IMPobsext - OBSsparse")
      dat_plot <- abs(dat_plot)
      
      boxplot(dat_plot,main=trait_names[t],ylim=c(0,
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.95,na.rm=TRUE))),col=colz,las=2,ylab=paste("abs(zlog) distance"))
      rect(xleft = 0,ybottom = min(apply(dat_plot,MARGIN = 2,quantile,prob=0,na.rm=TRUE))-1,xright = 6,ytop = median(dat_plot[,1],na.rm = TRUE),col=colz2[1])
      rect(xleft = 0,ybottom = median(dat_plot[,1],na.rm = TRUE),xright = 6,ytop = max(apply(dat_plot,MARGIN = 2,quantile,prob=1,na.rm=TRUE))+1,col=colz2[3])
      abline(h = median(dat_plot[,1],na.rm = TRUE),col=colz2[8],lwd=1)
#      abline(v = 3.5,col="white",lwd=1.5)
      boxplot(dat_plot,main=trait_names[t],ylim=c(min(apply(dat_plot,MARGIN = 2,quantile,prob=.01,na.rm=TRUE)),
                                                  max(apply(dat_plot,MARGIN = 2,quantile,prob=.99,na.rm=TRUE))),
              col=colz,las=2,ylab=paste("zlog(",units[t],")"),
              add=TRUE,cex.lab=1.5)
    }   
    dev.off()
    
  }
  
}



pdf(file=file.path(origin,"figures","Figure_2",paste0(txs[tx],"_xy_",t_choice,".pdf")),width=8,height=6,pointsize = 9)

plot(unlist(TD_tx),unlist(TDtd_tx), pch=16, col = "blue",
     xlab="OBS",ylab="IMPobs", cex.lab=2,
     xlim=c(min(c(unlist(TD_tx))),max(c(unlist(TD_tx)))),
     ylim=c(min(c(unlist(TD_tx))),max(c(unlist(TD_tx)))))
abline(0,1)
cor_now <- cor(unlist(TD_tx),unlist(TDtd_tx))
text(0,2,paste0("Pearson cor = ",round(cor_now,digits = 2)),cex=2)
plot(unlist(TD_tx),unlist(TDenv_tx),pch=16, col = "blue",
     xlab="OBS",ylab="IMPobsExt",cex.lab=2,
     xlim=c(min(c(unlist(TD_tx))),max(c(unlist(TD_tx)))),
     ylim=c(min(c(unlist(TD_tx))),max(c(unlist(TD_tx)))))
abline(0,1)
cor_now <- cor(unlist(TD_tx),unlist(TDenv_tx))
text(0,2,paste0("Pearson cor = ",round(cor_now,digits = 2)),cex=2)

dev.off()
