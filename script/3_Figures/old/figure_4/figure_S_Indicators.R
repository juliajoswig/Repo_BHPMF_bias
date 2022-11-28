

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
trait_names <- unique(res$trait)
trait_names <- trait_names[!is.na(trait_names)]


RepNum=1
t_choice="data"
ObsOrTD="Obs_obs"
Percent=80
resENV <- read.csv(file=file.path(origin,"_2021","data","analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021.csv")))
trait_names <- unique(res$trait)
trait_names <- trait_names[!is.na(trait_names)]


#-------------------------------------------------------------------
# for correl: dist_lm
# for tax: dist_
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

colzTRAITS <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .5),rgb(103/255,169/255,207/255,alpha = .5))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255),
             rgb(202/255,0/255,32/255),rgb(5/255,113/255,176/255))


lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
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
}
lm_TD <- lm_op 
group_TD <- group_op
vals_TD <- vals_op

lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]&res$dist_spec_obs_zlog
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(res$dist_spec_obs_zlog[ix_trait],decreasing = TRUE)[1:25]
  ix_trait[ix_trait][which_now]
  
  vals_op[[t]] <- cbind(res$value_pred_zlog[ix_trait][which_now],res$value_obs_zlog[ix_trait][which_now])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(res$dist_spec_obs_zlog[ix_trait][which_now],
                         res$dist_spec_pred_zlog[ix_trait][which_now],
                         res$dist_gen_obs_zlog[ix_trait][which_now],
                         res$dist_gen_pred_zlog[ix_trait][which_now],
                         res$dist_fam_obs_zlog[ix_trait][which_now],
                         res$dist_fam_pred_zlog[ix_trait][which_now],
                         res$dist_clad_obs_zlog[ix_trait][which_now],
                         res$dist_clad_pred_zlog[ix_trait][which_now],
                         res$dist_GF_obs_zlog[ix_trait][which_now],
                         res$dist_GF_pred_zlog[ix_trait][which_now],
                         res$dist_PFT_obs_zlog[ix_trait][which_now],
                         res$dist_PFT_pred_zlog[ix_trait][which_now])
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
}
lm_diss<- lm_op 
group_diss <- group_op
vals_diss <- vals_op


lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]&res$dist_spec_obs_zlog
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(res$dist_spec_obs_zlog[ix_trait],decreasing = FALSE)[1:25]
  ix_trait[ix_trait][which_now]
  
  vals_op[[t]] <- cbind(res$value_pred_zlog[ix_trait][which_now],res$value_obs_zlog[ix_trait][which_now])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(res$dist_spec_obs_zlog[ix_trait][which_now],
                         res$dist_spec_pred_zlog[ix_trait][which_now],
                         res$dist_gen_obs_zlog[ix_trait][which_now],
                         res$dist_gen_pred_zlog[ix_trait][which_now],
                         res$dist_fam_obs_zlog[ix_trait][which_now],
                         res$dist_fam_pred_zlog[ix_trait][which_now],
                         res$dist_clad_obs_zlog[ix_trait][which_now],
                         res$dist_clad_pred_zlog[ix_trait][which_now],
                         res$dist_GF_obs_zlog[ix_trait][which_now],
                         res$dist_GF_pred_zlog[ix_trait][which_now],
                         res$dist_PFT_obs_zlog[ix_trait][which_now],
                         res$dist_PFT_pred_zlog[ix_trait][which_now])
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
 }
lm_sim<- lm_op 
group_sim <- group_op
vals_sim <- vals_op


#------------------------------------------
# ENvelope

lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=2
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  vals_op[[t]] <- cbind(resENV$value_pred_zlog[ix_trait],resENV$value_obs_zlog[ix_trait])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(resENV$dist_spec_obs_zlog[ix_trait],
                         resENV$dist_spec_pred_zlog[ix_trait],
                         resENV$dist_gen_obs_zlog[ix_trait],
                         resENV$dist_gen_pred_zlog[ix_trait],
                         resENV$dist_fam_obs_zlog[ix_trait],
                         resENV$dist_fam_pred_zlog[ix_trait],
                         resENV$dist_clad_obs_zlog[ix_trait],
                         resENV$dist_clad_pred_zlog[ix_trait],
                         resENV$dist_GF_obs_zlog[ix_trait],
                         resENV$dist_GF_pred_zlog[ix_trait],
                         resENV$dist_PFT_obs_zlog[ix_trait],
                         resENV$dist_PFT_pred_zlog[ix_trait])
  dim(group_op[[t]])
  colnames(group_op[[t]]) <- c("Species","Species","Genus","Genus","Family","Family",
                               "Clade","Clade","Growth form","Growth form","PFT","PFT")
  
  colz_now = grep(x = colnames(resENV),pattern = paste0(trait_names[t],"_from_"))
  colz_now = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "zlog")]
  colz_obs = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "obs")]
  colz_pred = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "pred")]
  colnames(resENV)[colz_obs]
  colnames(resENV)[colz_pred]
  colnames(resENV)[colz_now]
  
  obs  =    resENV$value_obs_zlog - resENV[ix_trait,colz_obs]
  pred =    resENV$value_pred_zlog - resENV[ix_trait,colz_pred]
  
  colnames(obs) <- gsub(colnames(obs),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(obs) <- gsub(colnames(obs),pattern = "_lm_obs",replacement = "")
  colnames(obs) <- paste0(colnames(obs),"_obs")
  colnames(pred) <- gsub(colnames(pred),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(pred) <- gsub(colnames(pred),pattern = "_lm_pred",replacement = "")
  colnames(pred) <- paste0(colnames(pred),"_pred")
  lm_op[[t]]<- cbind(obs,pred)
  names(lm_op)[t] <- trait_names[t]
}
lm_ENV<- lm_op 
group_ENV <- group_op
vals_ENV <- vals_op


lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=3
for(t in 1:length(trait_names)){
  ix_trait=resENV$trait==trait_names[t]&resENV$dist_spec_obs_zlog
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(resENV$dist_spec_obs_zlog[ix_trait],decreasing = TRUE)[1:25]
  ix_trait[ix_trait][which_now]

  vals_op[[t]] <- cbind(resENV$value_pred_zlog[ix_trait][which_now],resENV$value_obs_zlog[ix_trait][which_now])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(resENV$dist_spec_obs_zlog[ix_trait][which_now],
                         resENV$dist_spec_pred_zlog[ix_trait][which_now],
                         resENV$dist_gen_obs_zlog[ix_trait][which_now],
                         resENV$dist_gen_pred_zlog[ix_trait][which_now],
                         resENV$dist_fam_obs_zlog[ix_trait][which_now],
                         resENV$dist_fam_pred_zlog[ix_trait][which_now],
                         resENV$dist_clad_obs_zlog[ix_trait][which_now],
                         resENV$dist_clad_pred_zlog[ix_trait][which_now],
                         resENV$dist_GF_obs_zlog[ix_trait][which_now],
                         resENV$dist_GF_pred_zlog[ix_trait][which_now],
                         resENV$dist_PFT_obs_zlog[ix_trait][which_now],
                         resENV$dist_PFT_pred_zlog[ix_trait][which_now])
  dim(group_op[[t]])
  colnames(group_op[[t]]) <- c("Species","Species","Genus","Genus","Family","Family",
                               "Clade","Clade","Growth form","Growth form","PFT","PFT")
  
  colz_now = grep(x = colnames(resENV),pattern = paste0(trait_names[t],"_from_"))
  colz_now = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "zlog")]
  colz_obs = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "obs")]
  colz_pred = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "pred")]
  colnames(resENV)[colz_obs]
  colnames(resENV)[colz_pred]
  colnames(resENV)[colz_now]
  
  obs  =    resENV$value_obs_zlog - resENV[ix_trait,colz_obs]
  pred =    resENV$value_pred_zlog - resENV[ix_trait,colz_pred]
  
  colnames(obs) <- gsub(colnames(obs),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(obs) <- gsub(colnames(obs),pattern = "_lm_obs",replacement = "")
  colnames(obs) <- paste0(colnames(obs),"_obs")
  colnames(pred) <- gsub(colnames(pred),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(pred) <- gsub(colnames(pred),pattern = "_lm_pred",replacement = "")
  colnames(pred) <- paste0(colnames(pred),"_pred")
  lm_op[[t]]<- cbind(obs,pred)
  names(lm_op)[t] <- trait_names[t]
}
lm_DissENV<- lm_op 
group_DissENV <- group_op
vals_DissENV <- vals_op


lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=resENV$trait==trait_names[t]&resENV$dist_spec_obs_zlog
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(resENV$dist_spec_obs_zlog[ix_trait],decreasing = FALSE)[1:25]
  ix_trait[ix_trait][which_now]
  
  vals_op[[t]] <- cbind(resENV$value_pred_zlog[ix_trait][which_now],resENV$value_obs_zlog[ix_trait][which_now])
  colnames(vals_op[[t]]) <- c("obs","pred")
  
  group_op[[t]] <- cbind(resENV$dist_spec_obs_zlog[ix_trait][which_now],
                         resENV$dist_spec_pred_zlog[ix_trait][which_now],
                         resENV$dist_gen_obs_zlog[ix_trait][which_now],
                         resENV$dist_gen_pred_zlog[ix_trait][which_now],
                         resENV$dist_fam_obs_zlog[ix_trait][which_now],
                         resENV$dist_fam_pred_zlog[ix_trait][which_now],
                         resENV$dist_clad_obs_zlog[ix_trait][which_now],
                         resENV$dist_clad_pred_zlog[ix_trait][which_now],
                         resENV$dist_GF_obs_zlog[ix_trait][which_now],
                         resENV$dist_GF_pred_zlog[ix_trait][which_now],
                         resENV$dist_PFT_obs_zlog[ix_trait][which_now],
                         resENV$dist_PFT_pred_zlog[ix_trait][which_now])
  dim(group_op[[t]])
  colnames(group_op[[t]]) <- c("Species","Species","Genus","Genus","Family","Family",
                               "Clade","Clade","Growth form","Growth form","PFT","PFT")
  
  colz_now = grep(x = colnames(resENV),pattern = paste0(trait_names[t],"_from_"))
  colz_now = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "zlog")]
  colz_obs = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "obs")]
  colz_pred = colz_now[grep(x = colnames(resENV)[colz_now],pattern = "pred")]
  colnames(resENV)[colz_obs]
  colnames(resENV)[colz_pred]
  colnames(resENV)[colz_now]
  
  obs  =    resENV$value_obs_zlog - resENV[ix_trait,colz_obs]
  pred =    resENV$value_pred_zlog - resENV[ix_trait,colz_pred]
  
  colnames(obs) <- gsub(colnames(obs),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(obs) <- gsub(colnames(obs),pattern = "_lm_obs",replacement = "")
  colnames(obs) <- paste0(colnames(obs),"_obs")
  colnames(pred) <- gsub(colnames(pred),pattern = paste0(trait_names[t],"_from_"),replacement = "")
  colnames(pred) <- gsub(colnames(pred),pattern = "_lm_pred",replacement = "")
  colnames(pred) <- paste0(colnames(pred),"_pred")
  lm_op[[t]]<- cbind(obs,pred)
  names(lm_op)[t] <- trait_names[t]
}
lm_simENV<- lm_op 
group_simENV <- group_op
vals_simENV <- vals_op


par(mfrow=c(1,1))
  #install.packages("fmsb")
  library(fmsb)
pdf(file=file.path(origin,"_2021","figures","Figure_4","Figure_S4_3.pdf"),width = 11,height=8)
par(mar=c(2,2,4,2),xpd = TRUE)
layout(mat = matrix(c(1,1,1,2,
                      1,1,1,3),nrow=2,byrow = TRUE))

  t=2
  for(t in 1:5){
    data_now <- c(apply(group_TD[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_TD[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_ENV[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_ENV[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_TD[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_TD[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_ENV[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_ENV[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
    #data_plot <- abs(rbind(data_now_obs,data_now_pred))
    #data_plot[data_plot>1] <- 1
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep(1,ncol(data_plot)) , data_plot)
    radarchart(data_chart,pcol = colz_solid,axistype = 2,plty = c(1,1,2,2),plwd = 8,pty = 16,axislabcol = "gray",
               pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Total ",trait_names[t]),vlcex = 2,calcex = 1.5,cex.main = 5,cglcol = "gray")
    text(0,0,trait_names[t])
    #--------------------------------------------------------------------
    #--------------------------------------------------------------------
    
    data_now <- c(apply(group_diss[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_diss[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_diss[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_diss[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_DissENV[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_DissENV[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    colnames(group_DissENV[[t]])
    data_now <- c(apply(group_DissENV[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_DissENV[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
#    data_plot <- abs(rbind(data_now_obs,data_now_pred))
    #data_plot[data_plot>1] <- 1
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep(1,ncol(data_plot)) , data_plot)
    radarchart(data_chart,pcol = colz_solid,axistype = 2,plty = c(1,1,2,2),plwd = 8,pty = 16,axislabcol = "gray",
               pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Diverse"),vlcex = 1,calcex = 1.5,cex.main = 2,cglcol = "gray")
    text(0,0,trait_names[t])
    
    #--------------------------------------------------------------------
    
    head(group_sim[[t]])
    data_now <- c(apply(group_sim[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_sim[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_simENV[[t]][,seq(2,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,5:8],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(group_simENV[[t]][,seq(1,12,2)],MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,1:4],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
    #data_plot <- abs(rbind(data_now_obs,data_now_pred))
    #data_plot[data_plot>1] <- 1
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep(1,ncol(data_plot)) , data_plot)
    radarchart(data_chart,pcol = colz_solid,axistype = 2,plty = c(1,1,2,2),plwd = 8,pty = 16,axislabcol = "gray",
               pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Homogenous"),vlcex = 2,calcex = 1.5,cex.main = 2,cglcol = "gray")
    text(0,0,trait_names[t])
    
  }
    dev.off()  
  
    
    
    