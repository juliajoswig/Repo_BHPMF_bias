

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
# load processed data   
#-------------------------------------------------------------------
  colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
  colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255),"#b2182b",rgb(178/255,24/255,43/255,alpha = .8),"#b97880","#b2182b")
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  colz=c("#b2182b","#b2182b","#ef8a62","#ef8a62","#fddbc7","#fddbc7","#f7f7f7","#f7f7f7",
         "#d1e5f0","#d1e5f0","#67a9cf","#67a9cf","#2166ac","#2166ac")
  
  
  RepNum=1
  t_choice="data_2"
  ObsOrTD="Obs_obs_TD"
  Percent=80
  res_sing <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  RepNum=1
  t_choice="data_2"
  ObsOrTD="Obs_obs"
  Percent=80
  res_singENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
  trait_names=as.vector(unique(res_singENV$trait))
  trait_names <- trait_names[!is.na(trait_names)]
  
  head(res_singENV$Sil_spec_obs_zlog)
  head(res_sing$Sil_spec_obs_zlog)
  
  #-------------------------------------------------------------------

  t=2
  w=1
  #--------------
  res_now=res_sing
  #--------------
  par(mfrow=c(1,1),mar=c(10,6,2,2))
    t=1
    ix_trait=res_now$trait==trait_names[t] 
    ix_trait[is.na(ix_trait)] <- FALSE
    sum(ix_trait)
    

    bxpl_sp_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$Species)[ix_trait]),ncol=1)
    bxpl_gen_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$Genus)[ix_trait]),ncol=1)
    bxpl_fam_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$Family)[ix_trait]),ncol=1)
    bxpl_clad_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$Clade)[ix_trait]),ncol=1)
    bxpl_GF_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$GF)[ix_trait]),ncol=1)
    bxpl_PFT_op=matrix(NA,nrow=sum(ix_trait&!duplicated(res_now$PFT)[ix_trait]),ncol=1)
    
    ix_trait=res_now$trait==trait_names[t] 
    ix_trait[is.na(ix_trait)] <- FALSE
    print(sum(ix_trait))
     bxpl_sp_op=cbind(bxpl_sp_op,
                       cbind(res_now$Sil_spec_obs[ix_trait&!duplicated(res_now$Species[ix_trait])],
                             res_now$Sil_spec_pred[ix_trait&!duplicated(res_now$Species[ix_trait])],
                             res_singENV$Sil_spec_pred[ix_trait&!duplicated(res_singENV$Species[ix_trait])]))
      bxpl_gen_op=cbind(bxpl_gen_op,
                        cbind(res_now$Sil_gen_obs[ix_trait&!duplicated(res_now$Genus[ix_trait])],
                              res_now$Sil_gen_pred[ix_trait&!duplicated(res_now$Genus[ix_trait])],
                        res_singENV$Sil_gen_pred[ix_trait&!duplicated(res_singENV$Genus[ix_trait])]))
      bxpl_fam_op=cbind(bxpl_fam_op,
                        cbind(res_now$Sil_fam_obs[ix_trait&!duplicated(res_now$Family[ix_trait])],
                              res_now$Sil_fam_pred[ix_trait&!duplicated(res_now$Family[ix_trait])],
                        res_singENV$Sil_fam_pred[ix_trait&!duplicated(res_singENV$Family[ix_trait])]))
      bxpl_clad_op=cbind(bxpl_clad_op,
                         cbind(res_now$Sil_clad_obs[ix_trait&!duplicated(res_now$Clade[ix_trait])],
                               res_now$Sil_clad_pred[ix_trait&!duplicated(res_now$Clade[ix_trait])],
                         res_singENV$Sil_clad_pred[ix_trait&!duplicated(res_singENV$Clade[ix_trait])]))
      bxpl_GF_op=cbind(bxpl_GF_op,
                       cbind(res_now$Sil_GF_obs[ix_trait&!duplicated(res_now$GF[ix_trait])],
                             res_now$Sil_GF_pred[ix_trait&!duplicated(res_now$GF[ix_trait])],
                       res_singENV$Sil_GF_pred[ix_trait&!duplicated(res_singENV$GF[ix_trait])]))
      bxpl_PFT_op=cbind(bxpl_PFT_op,
                        cbind(res_now$Sil_PFT_obs[ix_trait&!duplicated(res_now$PFT[ix_trait])],
                              res_now$Sil_PFT_pred[ix_trait&!duplicated(res_now$PFT[ix_trait])],
                        res_singENV$Sil_PFT_pred[ix_trait&!duplicated(res_singENV$PFT[ix_trait])]))
    
    
    bxpl_sp_op <- bxpl_sp_op[bxpl_sp_op[,2]!=0,colSums(!is.na(bxpl_sp_op))!=0]
    bxpl_gen_op <- bxpl_gen_op[bxpl_gen_op[,2]!=0,colSums(!is.na(bxpl_gen_op))!=0]
    bxpl_fam_op <- bxpl_fam_op[bxpl_fam_op[,2]!=0,colSums(!is.na(bxpl_fam_op))!=0]
    bxpl_clad_op <- bxpl_clad_op[bxpl_clad_op[,2]!=0,colSums(!is.na(bxpl_clad_op))!=0]
    bxpl_GF_op <- bxpl_GF_op[bxpl_GF_op[,2]!=0,colSums(!is.na(bxpl_GF_op))!=0]
    bxpl_PFT_op <- bxpl_PFT_op[bxpl_PFT_op[,2]!=0,colSums(!is.na(bxpl_PFT_op))!=0]
    
    
    pdf(file=file.path(origin,"_2021","figures","Figure_2",paste0("Figure_S2_Sil_ENV.pdf")),width=17,height=17)
    
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
        dat_now=bxpl[rowSums(bxpl,na.rm = TRUE)!=0,colSums(!is.na(bxpl))!=0]
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz_solid,main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n")
        axis(1,at = 1:ncol(dat_now),labels = c("obs","pred","env"),las=2,cex.axis=5)
        yBot_now=median(dat_now[,1],na.rm = TRUE)
        i=1
        rect(xleft = .5,ybottom = ymin,xright = 3.5,ytop = yBot_now[i],col=colz[6],border = NA)
        rect(xleft = .5,ybottom = yBot_now,xright = 3.5,ytop = ymax,col=colz[10],border = NA) 
        boxplot(dat_now,ylim=c(ymin,ymax),ylab="Sil pred - obs",col=colz_solid,main=names_now[t],cex.main=5,
                las=2,cex.axis=1.5,cex.lab=2,lwd=2,yaxt="n",frame=FALSE,xaxt="n",add=TRUE)
        abline(v=seq(from = .5,to = 12.5,by = 3),lwd=4,col="white")

      }
 
    }
    
    dev.off()
    
  