

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
gappercents= c("1","5","10","20","30","40","50","60")
repnums=1:3

GapPercent=50
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))

RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
colz=c("#b2182b","#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")


  res <- res[,colSums(!is.na(res))!=0]
  trait_names=as.vector(unique(res$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  missingness = unique(as.vector(res$missingness))
  missingness <- missingness[!is.na(missingness)]
  summary(res$value_obs)
  m=1
  t=2
  w=1

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
    print(head(bxpl_sp))
    

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
                       cbind(res$Sil_spec_obs[ix_trait&!duplicated(res$Species[ix_trait])],
                             res$Sil_spec_pred[ix_trait&!duplicated(res$Species[ix_trait])]))
      bxpl_gen_op=cbind(bxpl_gen_op,
                        cbind(res$Sil_gen_obs[ix_trait&!duplicated(res$Genus[ix_trait])],
                              res$Sil_gen_pred[ix_trait&!duplicated(res$Genus[ix_trait])]))
      bxpl_fam_op=cbind(bxpl_fam_op,
                        cbind(res$Sil_fam_obs[ix_trait&!duplicated(res$Family[ix_trait])],
                              res$Sil_fam_pred[ix_trait&!duplicated(res$Family[ix_trait])]))
      bxpl_clad_op=cbind(bxpl_clad_op,
                         cbind(res$Sil_clad_obs[ix_trait&!duplicated(res$Clade[ix_trait])],
                               res$Sil_clad_pred[ix_trait&!duplicated(res$Clade[ix_trait])]))
      bxpl_GF_op=cbind(bxpl_GF_op,
                       cbind(res$Sil_GF_obs[ix_trait&!duplicated(res$GF[ix_trait])],
                             res$Sil_GF_pred[ix_trait&!duplicated(res$GF[ix_trait])]))
      bxpl_PFT_op=cbind(bxpl_PFT_op,
                        cbind(res$Sil_PFT_obs[ix_trait&!duplicated(res$PFT[ix_trait])],
                              res$Sil_PFT_pred[ix_trait&!duplicated(res$PFT[ix_trait])]))
      
    
    
    bxpl_sp <- bxpl_sp[,colSums(!is.na(bxpl_sp))!=0]
    bxpl_gen <- bxpl_gen[,colSums(!is.na(bxpl_gen))!=0]
    bxpl_fam <- bxpl_fam[,colSums(!is.na(bxpl_fam))!=0]
    bxpl_clad <- bxpl_clad[,colSums(!is.na(bxpl_clad))!=0]
    bxpl_GF <- bxpl_GF[,colSums(!is.na(bxpl_GF))!=0]
    bxpl_PFT <- bxpl_PFT[,colSums(!is.na(bxpl_PFT))!=0]

    colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
           "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
    bxpl_sp_op <- bxpl_sp_op[bxpl_sp_op[,2]!=0,colSums(!is.na(bxpl_sp_op))!=0]
    bxpl_gen_op <- bxpl_gen_op[bxpl_gen_op[,2]!=0,colSums(!is.na(bxpl_gen_op))!=0]
    bxpl_fam_op <- bxpl_fam_op[bxpl_fam_op[,2]!=0,colSums(!is.na(bxpl_fam_op))!=0]
    bxpl_clad_op <- bxpl_clad_op[bxpl_clad_op[,2]!=0,colSums(!is.na(bxpl_clad_op))!=0]
    bxpl_GF_op <- bxpl_GF_op[bxpl_GF_op[,2]!=0,colSums(!is.na(bxpl_GF_op))!=0]
    bxpl_PFT_op <- bxpl_PFT_op[bxpl_PFT_op[,2]!=0,colSums(!is.na(bxpl_PFT_op))!=0]
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Sil",missingness[m],".pdf")),width=15,height=15)
    
    par(mar=c(10,0,4,0),mfrow=c(1,7))
    {
      ymax=1
      ymin=-1
      plot(1:10,col="white",ylim=c(0,ymax),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -20,cex.axis=2)
      axis(side = 2,at = .5,line = -14,labels = "Silhouette Index",tick = FALSE,cex.axis=7,cex.lab=5)

      names_now=c("Species","Genus","Family","Clades","GF","PFT")
      bxpl=bxpl_sp_op
      t=1
      for(t in 1:6){
        if(t==1){bxpl=bxpl_sp_op}
        if(t==2){bxpl=bxpl_gen_op}
        if(t==3){bxpl=bxpl_fam_op}
        if(t==4){bxpl=bxpl_clad_op}
        if(t==5){bxpl=bxpl_GF_op}
        if(t==6){bxpl=bxpl_PFT_op}
        dat_now=bxpl
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz_solid,main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
        axis(1,at = 1:ncol(bxpl),labels = c("obs","pred"),las=2,cex.axis=5)
        yBot_now=median(dat_now[,1],na.rm = TRUE)
        i=1
        rect(xleft = .5,ybottom = ymin,xright = 2.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = .5,ybottom = yBot_now,xright = 2.5,ytop = ymax,col=colz[10],border = NA)
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz_solid,main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n",add=TRUE)
        abline(v=seq(from = .5,to = 12.5,by = 2),lwd=4,col="white")

      }
 
    }
    
    dev.off()
    
  
    
    #--------------------------------------
    # zlog
    #--------------------------------------
    
    
    
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
    print(head(bxpl_sp))
    
    
    ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
    ix_trait[is.na(ix_trait)] <- FALSE
    print(sum(ix_trait))
    
    bxpl_sp=cbind(bxpl_sp,
                  res$Sil_spec_pred_zlog[ix_trait&!duplicated(res$Species[ix_trait])]-
                    res$Sil_spec_obs_zlog[ix_trait&!duplicated(res$Species[ix_trait])])
    bxpl_gen=cbind(bxpl_gen,
                   res$Sil_gen_pred_zlog[ix_trait&!duplicated(res$Genus[ix_trait])]-
                     res$Sil_gen_obs_zlog[ix_trait&!duplicated(res$Genus[ix_trait])])
    bxpl_fam=cbind(bxpl_fam,
                   res$Sil_fam_pred_zlog[ix_trait&!duplicated(res$Family[ix_trait])]-
                     res$Sil_fam_obs_zlog[ix_trait&!duplicated(res$Family[ix_trait])])
    bxpl_clad=cbind(bxpl_clad,
                    res$Sil_clad_pred_zlog[ix_trait&!duplicated(res$Clade[ix_trait])] -
                      res$Sil_clad_obs_zlog[ix_trait&!duplicated(res$Clade[ix_trait])])
    bxpl_GF=cbind(bxpl_GF,
                  res$Sil_GF_pred_zlog[ix_trait&!duplicated(res$GF[ix_trait])]-
                    res$Sil_GF_obs_zlog[ix_trait&!duplicated(res$GF[ix_trait])])
    bxpl_PFT=cbind(bxpl_PFT,
                   res$Sil_PFT_pred_zlog[ix_trait&!duplicated(res$PFT[ix_trait])]-
                     res$Sil_PFT_obs_zlog[ix_trait&!duplicated(res$PFT[ix_trait])])
    
    bxpl_sp_op=cbind(bxpl_sp_op,
                     cbind(res$Sil_spec_obs_zlog[ix_trait&!duplicated(res$Species[ix_trait])],
                           res$Sil_spec_pred_zlog[ix_trait&!duplicated(res$Species[ix_trait])]))
    bxpl_gen_op=cbind(bxpl_gen_op,
                      cbind(res$Sil_gen_obs_zlog[ix_trait&!duplicated(res$Genus[ix_trait])],
                            res$Sil_gen_pred_zlog[ix_trait&!duplicated(res$Genus[ix_trait])]))
    bxpl_fam_op=cbind(bxpl_fam_op,
                      cbind(res$Sil_fam_obs_zlog[ix_trait&!duplicated(res$Family[ix_trait])],
                            res$Sil_fam_pred_zlog[ix_trait&!duplicated(res$Family[ix_trait])]))
    bxpl_clad_op=cbind(bxpl_clad_op,
                       cbind(res$Sil_clad_obs_zlog[ix_trait&!duplicated(res$Clade[ix_trait])],
                             res$Sil_clad_pred_zlog[ix_trait&!duplicated(res$Clade[ix_trait])]))
    bxpl_GF_op=cbind(bxpl_GF_op,
                     cbind(res$Sil_GF_obs_zlog[ix_trait&!duplicated(res$GF[ix_trait])],
                           res$Sil_GF_pred_zlog[ix_trait&!duplicated(res$GF[ix_trait])]))
    bxpl_PFT_op=cbind(bxpl_PFT_op,
                      cbind(res$Sil_PFT_obs_zlog[ix_trait&!duplicated(res$PFT[ix_trait])],
                            res$Sil_PFT_pred_zlog[ix_trait&!duplicated(res$PFT[ix_trait])]))
    
    
    
    bxpl_sp <- bxpl_sp[,colSums(!is.na(bxpl_sp))!=0]
    bxpl_gen <- bxpl_gen[,colSums(!is.na(bxpl_gen))!=0]
    bxpl_fam <- bxpl_fam[,colSums(!is.na(bxpl_fam))!=0]
    bxpl_clad <- bxpl_clad[,colSums(!is.na(bxpl_clad))!=0]
    bxpl_GF <- bxpl_GF[,colSums(!is.na(bxpl_GF))!=0]
    bxpl_PFT <- bxpl_PFT[,colSums(!is.na(bxpl_PFT))!=0]
    
    colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
           "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
    bxpl_sp_op <- bxpl_sp_op[,colSums(!is.na(bxpl_sp_op))!=0]
    bxpl_gen_op <- bxpl_gen_op[,colSums(!is.na(bxpl_gen_op))!=0]
    bxpl_fam_op <- bxpl_fam_op[,colSums(!is.na(bxpl_fam_op))!=0]
    bxpl_clad_op <- bxpl_clad_op[,colSums(!is.na(bxpl_clad_op))!=0]
    bxpl_GF_op <- bxpl_GF_op[,colSums(!is.na(bxpl_GF_op))!=0]
    bxpl_PFT_op <- bxpl_PFT_op[,colSums(!is.na(bxpl_PFT_op))!=0]
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_2_Sil_zlog_",missingness[m],".pdf")),width=15,height=15)
    
    par(mar=c(10,0,4,0),mfrow=c(1,7))
    {
      ymax=1
      ymin=-1
      plot(1:10,col="white",ylim=c(0,ymax),yaxt="n",xaxt="n",ylab="",xlab="",frame=FALSE)
      axis(side = 2,line = -20,cex.axis=2)
      axis(side = 2,at = .5,line = -14,labels = "Silhouette Index",tick = FALSE,cex.axis=7,cex.lab=5)
      
      names_now=c("Species","Genus","Family","Clades","GF","PFT")
      bxpl=bxpl_sp_op
      t=1
      for(t in 1:6){
        if(t==1){bxpl=bxpl_sp_op}
        if(t==2){bxpl=bxpl_gen_op}
        if(t==3){bxpl=bxpl_fam_op}
        if(t==4){bxpl=bxpl_clad_op}
        if(t==5){bxpl=bxpl_GF_op}
        if(t==6){bxpl=bxpl_PFT_op}
        dat_now=bxpl
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz[c(1,3)],main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
        axis(1,at = 1:ncol(bxpl),labels = c("obs","pred"),las=2,cex.axis=5)
        yBot_now=median(dat_now[,1],na.rm = TRUE)
        i=1
        rect(xleft = .5,ybottom = ymin,xright = 2.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = .5,ybottom = yBot_now,xright = 2.5,ytop = ymax,col=colz[10],border = NA)
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz[c(1,3)],main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n",add=TRUE)
        abline(v=seq(from = .5,to = 12.5,by = 2),lwd=4,col="white")
        
      }
      
    }
    
    dev.off()
    
    