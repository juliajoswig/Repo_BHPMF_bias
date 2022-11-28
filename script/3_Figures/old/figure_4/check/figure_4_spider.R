

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
t_choice="data"
ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))

#-------------------------------------------------------------------
# for correl: dist_lm
# for tax: dist_
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

colzTRAITS <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .7),rgb(103/255,169/255,207/255,alpha = .7))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255))

lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]&res$missingness==missingness[m]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  vals_op[[t]] <- cbind(res$value_pred_zlog[ix_trait],res$value_obs_zlog[ix_trait])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(res$dist_spec_obs_zlog[ix_trait],
                          res$dist_spec_pred_zlog[ix_trait],
                          res$dist_gen_obs_zlog[ix_trait],
                          res$dist_gen_pred_zlog[ix_trait],
                          res$dist_fam_obs_zlog[ix_trait],
                          res$dist_fam_pred_zlog[ix_trait],
                          res$dist_clad_obs_zlog[ix_trait],
                          res$dist_clad_pred_zlog[ix_trait],
                          res$dist_GF_obs_zlog[ix_trait],
                          res$dist_GF_pred_zlog[ix_trait],
                          res$dist_PFT_obs_zlog[ix_trait],
                          res$dist_PFT_pred_zlog[ix_trait])
  dim(group_op[[t]])
  colnames(group_op[[t]]) <- c("Species","Species","Genus","Genus","Family","Family",
                               "Clade","Clade","Growth form","Growth form","PFT","PFT")
  
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
  lm_op[[t]]<- cbind(obs,pred)
  names(lm_op)[t] <- trait_names[t]
  
  colz[[t]] <- rep(colzTRAITS[match(unique(gsub(colnames(bxpl_tot[[t]])[1:5],pattern = "_zlog_pred",replacement = "")),trait_names)],2)[c(1,1,2,2,3,3,4,4,5,5)]
}







t=1
  #install.packages("fmsb")
  library(fmsb)
  pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_4_Spider_zlog.pdf"),width = 6,height=5)
  par(mar=c(2,2,2,2),mfrow=c(1,1))
  for(t in 1:6){
    data_now <- c(apply(group_op[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                   apply(lm_op[[t]][,1:5],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    
    data_now <- c(apply(group_op[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                   apply(lm_op[[t]][,6:10],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    data_plot <- abs(rbind(data_now_obs,data_now_pred))
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep(1,ncol(data_plot)) , data_plot)
    radarchart(data_chart,pcol = colz_solid,axistype = 2,plty = 1,plwd = 4,pty = 32,axislabcol = "gray",
               title = ,vlcex = 1.8,calcex = 1.5)
    text(0,0,trait_names[t])
  }
    dev.off()  
  
    
    
    