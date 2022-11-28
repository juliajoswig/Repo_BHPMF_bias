setwd("/..")
origin = "Volumes/bgi/work_1/2016_GapFilling"
#origin="2016_GapFilling"
res_matrix_name="res_20201020"
#res_o=res

res <- read.table(file.path(origin,"runs","META",paste0(res_matrix_name,".csv")),sep=",",dec=".")
colnames(res)
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
# Silhouette increases with gap-size
#--------------------------------------
sel_ObSpe1="Obs_obs_TD"
pattern_now="RMSE"
#chose guido
ix1=res[,colnames(res)=="TraitChoice"]=="guido"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent1 <- res[ix1,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse <- res[ix1,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,colSums(!is.na(dat_rmse))!=0]
colnames(dat_rmse)

dat_rmse <- as.matrix(dat_rmse)
mode(dat_rmse)="numeric"
colnames(dat_rmse)

ix=res$GapPercent==0
boxplot(dat_rmse,las=2)
#dat_plot <- aggregate(dat_rmse,by=list(Percent1),FUN=mean,na.rm = TRUE)
#dat_plot <- dat_plot[complete.cases(dat_plot),]

#chose rainfor
ix2=res[,colnames(res)=="TraitChoice"]=="rainfor"&res[,colnames(res)=="Obs_or_Spec"]==sel_ObSpe1
Percent2 <- res[ix2,colnames(res)=="GapPercent"]
mode(Percent1) <- "numeric"

dat_rmse2 <- res[ix2,grep(colnames(res),pattern = pattern_now)]
dat_rmse <- dat_rmse[,colSums(!is.na(dat_rmse))!=0]
colnames(dat_rmse2)

dat_rmse2 <- as.matrix(dat_rmse2)
mode(dat_rmse2)="numeric"
colnames(dat_rmse2)

boxplot(dat_rmse2,las=2)
#dat_plot2 <- aggregate(dat_rmse2,by=list(Percent1),FUN=mean,na.rm = TRUE)
#dat_plot2 <- dat_plot2[complete.cases(dat_plot2),]

