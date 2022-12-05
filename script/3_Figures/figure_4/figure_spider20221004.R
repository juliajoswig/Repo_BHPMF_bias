# Julia Joswig 
# 20221003, last change 202212


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
  
  
  t_choice="data"
  RepNum=1
  
#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(originData,"analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

colzTRAITS <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#d1e5f0", "#67a9cf", "#2166ac")
colz_alpha=c(rgb(239/255,138/255,98/255,alpha = .5),rgb(103/255,169/255,207/255,alpha = .5))
colz_solid=c(rgb(103/255,169/255,207/255),rgb(239/255,138/255,98/255),
             rgb(202/255,0/255,32/255),rgb(5/255,113/255,176/255))

ObsOrTD="Obs_obs_TD"
Percent=80

list.files(file.path(originData,"analyses","Point_wise",
                     RepNum,t_choice,ObsOrTD,Percent))

res <- read.csv(file=file.path(originData,"analyses","Point_wise",RepNum,
                               t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))
res_matrix_name="res_20221203"
res <- read.csv(file=file.path(originData,"analyses","TOTAL",paste0(res_matrix_name,".csv")))
trait_names <- unique(res$trait)
trait_names <- trait_names[!is.na(trait_names)]


ObsOrTD="Obs_obs"
Percent=80
list.files(file.path(originData,"analyses","Point_wise",
                     RepNum,t_choice,ObsOrTD,Percent))

resENV <- read.csv(file=file.path(originData,"analyses","Point_wise",
                                  RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))

#-------------------------------------------------------------------
# for correl: dist_lm
# for tax: dist_
#-------------------------------------------------------------------
{
lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_NO1 = res$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
  ix_trait[is.na(ix_trait)] <- FALSE
  sum(ix_trait)
  
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
lm_TD<- lm_op 
group_TD <- group_op
vals_TD <- vals_op

lm_op <- list()
group_op <- list()
vals_op <- list()
colz <- list()
t=1
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_NO1 = res$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(abs(res$dist_spec_obs_zlog[ix_trait]),decreasing = TRUE)[1:25]
  abs(res$dist_spec_obs_zlog[ix_trait])[which_now]
  
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
  ix_trait=res$trait==trait_names[t]
  ix_NO1 = res$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
  ix_trait[is.na(ix_trait)] <- FALSE
  which_now <- order(abs(res$dist_spec_obs_zlog[ix_trait]),decreasing = FALSE)[1:25]
  abs(res$dist_spec_obs_zlog[ix_trait])[which_now]
  
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
t=1
for(t in 1:length(trait_names)){
  ix_trait=resENV$trait==trait_names[t]
  ix_NO1 = resENV$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
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
t=1
for(t in 1:length(trait_names)){
  ix_trait=resENV$trait==trait_names[t]
  ix_NO1 = resENV$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
  ix_trait[is.na(ix_trait)] <- FALSE  
  which_now <- order(abs(resENV$dist_spec_obs_zlog[ix_trait]),decreasing = TRUE)[1:25]
  abs(res$dist_spec_obs_zlog[ix_trait])[which_now]
  
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
  ix_trait=resENV$trait==trait_names[t]
  ix_NO1 = resENV$nb_spec>1
  ix_trait <- ix_trait&ix_NO1
  ix_trait[is.na(ix_trait)] <- FALSE  
  which_now <- order(abs(resENV$dist_spec_obs_zlog[ix_trait]),decreasing = FALSE)[1:25]
  abs(res$dist_spec_obs_zlog[ix_trait])[which_now]
  
  
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
}

t=2
seqs=seq(2,12,2)
if(t_choice=="data"){columns_now1=1:5;columns_now=6:10}
if(t_choice=="data_2"){columns_now1=1:4;columns_now=5:8}


#install.packages("fmsb")
library(fmsb)
{
  if(t_choice=="data"){
  pdf(file=file.path(originData,"figures","Figure_6","Figure_6_Spider_r.pdf"),width = 8,height=11)
  }
  if(t_choice=="data_2"){
    pdf(file=file.path(originData,"figures","Figure_6","Figure_6_SpiderTD2_r.pdf"),width = 8,height=11)}
  layout(mat = matrix(c(1,1,
                        1,1,
                        2,3),nrow=3,byrow = TRUE))
  border=1
  t=1
  for(t in 1:length(trait_names)){
    data_now <- c(apply(abs(group_TD[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_TD[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))

    data_now <- c(apply(abs(group_ENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_ENV[[t]][,1:length(trait_names)],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_TD[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_TD[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_ENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_ENV[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))[,c(1:6)]
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 2,axistype = 2,plty = c(1,1,1,1),plwd = c(20,8,8,8),pty = 16,axislabcol = "gray",
               #  pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Total ",trait_names[t]),vlcex = 3,calcex = 1.5,cex.main = 5,cglcol = "gray")
    #text(0,0,trait_names[t])
    
    #--------------------------------------------------------------------
    #--------------------------------------------------------------------
    
    data_now <- c(apply(abs(group_diss[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_diss[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_diss[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_diss[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_DissENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_DissENV[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    colnames(group_DissENV[[t]])
    data_now <- c(apply(abs(group_DissENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_DissENV[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))[,c(1:6)]
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 1,axistype = 2,plty = c(1,1,1,1),plwd = c(8,5,5,5),pty = 16,axislabcol = "gray",
               #pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Outliers"),vlcex = 1.4,calcex = 1.5,cex.main = 2,cglcol = "gray")
    #text(0,0,trait_names[t])
    
    #--------------------------------------------------------------------
    
    
    data_now <- c(apply(abs(group_sim[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_sim[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_simENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,columns_now],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_simENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_sim[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))[,c(1:6)]
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 1,axistype = 2,plty = c(1,1,1,1),plwd = c(8,5,5,5),pty = 16,axislabcol = "gray",
               #pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Aligned"),vlcex = 1.4,calcex = 1.5,cex.main = 2,cglcol = "gray")
    #text(0,0,trait_names[t])
    
  }
  dev.off()  
}




#par(mar=c(2,2,4,2),xpd = TRUE,mfrow=c(1,1))


{
  if(t_choice=="data"){
    pdf(file=file.path(originData,"figures","Figure_4","Figure_SpiderTraitWise_r.pdf"),width = 8,height=11)
  }
  if(t_choice=="data_2"){
    pdf(file=file.path(originData,"figures","Figure_4","Figure_S_SpiderTD2TraitWise_r.pdf"),width = 8,height=11)}
  par(mar=c(2,2,4,2),xpd = TRUE)
  layout(mat = matrix(c(1,1,
                        1,1,
                        2,3),nrow=3,byrow = TRUE))
  
  border=1
  t=1
  for(t in 1:length(trait_names)){
    head(group_TD[[1]][,1:5])
    colnames(group_TD[[t]])
    data_now <- c(apply(abs(group_TD[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_ENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_ENV[[t]][,1:length(trait_names)],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_TD[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_ENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 2,axistype = 2,plty = c(1,1,1,1),plwd = c(20,8,8,8),pty = 16,axislabcol = "gray",
               #  pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Total ",trait_names[t]),vlcex = 3,calcex = 1.5,cex.main = 5,cglcol = "gray")
    #text(0,0,trait_names[t])
    #--------------------------------------------------------------------
    #--------------------------------------------------------------------
    data_now <- c(apply(abs(group_diss[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_diss[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_DissENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    colnames(group_DissENV[[t]])
    data_now <- c(apply(abs(group_DissENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE),
                  apply(lm_DissENV[[t]][,columns_now1],MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 1,axistype = 2,plty = c(1,1,1,1),plwd = c(8,5,5,5),pty = 16,axislabcol = "gray",
               #pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Outliers"),vlcex = 1.4,calcex = 1.5,cex.main = 2,cglcol = "gray")
    #text(0,0,trait_names[t])
    
    #--------------------------------------------------------------------
    
    
    data_now <- c(apply(abs(group_sim[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obs <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_sim[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_pred <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_simENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_pred",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_predENV <- as.data.frame(t(data_now))
    
    data_now <- c(apply(abs(group_simENV[[t]][,seqs]),MARGIN = 2,FUN = median,na.rm=TRUE))
    names(data_now) <-gsub(names(data_now),pattern = "_obs",replacement = "") 
    names(data_now) <-gsub(names(data_now),pattern = "_zlog",replacement = "") 
    data_now_obsENV <- as.data.frame(t(data_now))
    
    data_plot <- abs(rbind(data_now_obs,data_now_pred,data_now_predENV))
    data_plot[data_plot>border] <- border
    data_chart <- rbind(rep(0,ncol(data_plot)) , rep((border-.1),ncol(data_plot)) , data_plot)
    par(mar=c(0,0,4,0),xpd = TRUE)
    radarchart(data_chart,pcol = colz_solid,cglty = 1,cglwd = 1,axistype = 2,plty = c(1,1,1,1),plwd = c(8,5,5,5),pty = 16,axislabcol = "gray",
               #pfcol = c(colz_alpha[2],NA,NA),
               title = paste0("Aligned"),vlcex = 1.4,calcex = 1.5,cex.main = 2,cglcol = "gray")
    #text(0,0,trait_names[t])
    
  }
  dev.off()  
}



data_chart <- matrix(0,ncol=6,nrow=4)
data_chart <- rbind(data_chart,rep(1,6))

#data_charto <- data_chart
#data_chart <- data_chart[1:3,]
data_chart[3,] <- c(3,3,3,3,3,3)
data_chart <- data.frame(data_chart)
# ---------------------------------------------------------------------
# explanatory figure
# ---------------------------------------------------------------------
getwd()
list.files(file.path(origin,"figures"))

pdf(file = file.path(origin,"figures","Figure_5","Figure_SpiderExplain_r.pdf"), width = 18,height=6)

  par(mar=c(0,0,0,0),xpd = TRUE,mfrow=c(1,2))
  
  radarchart(data_chart, pcol = "white", cglty = 1, 
             cglwd = 2.5, axistype = 3,
             plty = c(1,1,1,1), plwd = c(.01,8,8,8),pty = 16,
             axislabcol = "black",
             #  pfcol = c(colz_alpha[2],NA,NA),
             vlcex =1.5, calcex = 1.5,cglcol = "gray",
             caxislabels = c(">=1 little clustering",
                             "0.75",
                             "0.5",
                             "0.25",
                             "0 strong clustering"),
             vlabels = c("species",
                          "genera",
                          "families",
                          "clades",
                          "GF",
                             "PFT"
))
  
  plot(x=1:10,y=1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",ylab="",xlab="")
  legend(0.5,10,legend = c("OBS - reference line for:","IMPobs","IMPobsExt"),
         col=colz_solid[1:3],lwd=8,lty=1,bty="n",cex=3)
  text(2.6,5,labels = "x: individual trait value",cex=2)
   
  text(2.75,4,labels = "Little clustering = centre",cex=2)
  text(3.4,3,labels = "Strong Clustering = circle border",cex=2)
#  text(4.75,2,labels = "lm(trait_n): residual of predicted target trait value ",cex=2)
#  text(2.5,1,labels = "from lm, input trait_n",cex=2)
  
  dev.off()

  
  pdf(file = file.path(origin,"figures","Figure_5","Figure_SpiderExplainTD2_r.pdf"), width = 18, height=6)
  par(mar=c(0,0,0,0),xpd = TRUE,mfrow=c(1,2))
  
  radarchart(data_chart, pcol = "white", cglty = 1, 
             cglwd = 2.5, axistype = 3,
             plty = c(1,1,1,1), plwd = c(.01,8,8,8),pty = 16,
             axislabcol = "black",
             #  pfcol = c(colz_alpha[2],NA,NA),
             vlcex =1.5, calcex = 1.5,cglcol = "gray",
             caxislabels = c(">=1 little clustering",
                             "0.75",
                             "0.5",
                             "0.25",
                             "0 strong clustering"),
             vlabels = c("species",
                         "genera",
                         "families",
                         "clades",
                         "GF",
                         "PFT"
             ))
  
  plot(x=1:10,y=1:10,col="white",frame=FALSE,xaxt="n",yaxt="n",ylab="",xlab="")
  legend(0.5,10,legend = c("OBS2 - reference line for:","IMP2obs","IMP2obsExt"),
         col=colz_solid[1:3],lwd=8,lty=1,bty="n",cex=3)
  text(2.6,5,labels = "x: individual trait value",cex=2)
  
  text(2.75,4,labels = "Little clustering = centre",cex=2)
  text(3.4,3,labels = "Strong Clustering = circle border",cex=2)
  
  dev.off()
  
  