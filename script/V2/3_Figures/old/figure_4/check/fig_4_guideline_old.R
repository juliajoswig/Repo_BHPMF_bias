

# for 0% gaps and for 70% gaps:

# load the RMSE 
# calculate per cluster an average
# load the Silhouette data per cluster
# load the mean(sd) data per cluster within silhouettes I think
# load the coefficience of variance data per cluster
# load the number of observations somehow
GapPercent1=70
RepNum=1
TD_choice="Obs_obs_TD"
trait_sub="guido"
# load the RMSE 
rmse_0_general <- read.csv(file=file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all.csv")))
rmse_0_rawclust <- read.csv(file= file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent1,RepNum,paste0("RMSE_all_PlusCluster.csv")))

# calculate per cluster an average
rmse_0_spec <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,2]),FUN=mean,na.rm=TRUE)
rmse_0_gen <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,3]),FUN=mean,na.rm=TRUE)
rmse_0_fam <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,4]),FUN=mean,na.rm=TRUE)
rmse_0_pg <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,5]),FUN=mean,na.rm=TRUE)
rmse_0_gf <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,6]),FUN=mean,na.rm=TRUE)
rmse_0_pft <- aggregate(rmse_0_rawclust[,8:ncol(rmse_0_rawclust)],by=list(rmse_0_rawclust[,7]),FUN=mean,na.rm=TRUE)

# load the Silhouette data per cluster
# load the mean(sd) data per cluster within silhouettes I think
# load the number of observations somehow
sil_0_rawclust <- read.csv(file=file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")),row.names = NULL)
sil_0_spec <- sil_0_rawclust[sil_0_rawclust[,2]=="Species",3:ncol(sil_0_rawclust)]
sil_0_spec <- sil_0_spec[complete.cases(sil_0_spec),]
sil_0_gen <- sil_0_rawclust[sil_0_rawclust[,2]=="Genus",3:ncol(sil_0_rawclust)]
sil_0_gen <- sil_0_gen[complete.cases(sil_0_gen),]
sil_0_fam <- sil_0_rawclust[sil_0_rawclust[,2]=="Family",3:ncol(sil_0_rawclust)]
sil_0_fam <- sil_0_fam[complete.cases(sil_0_fam),]
sil_0_pg <- sil_0_rawclust[sil_0_rawclust[,2]=="PG",3:ncol(sil_0_rawclust)]
sil_0_pg <- sil_0_pg[complete.cases(sil_0_pg),]
sil_0_gf <- sil_0_rawclust[sil_0_rawclust[,2]=="GF",3:ncol(sil_0_rawclust)]
sil_0_gf <- sil_0_gf[complete.cases(sil_0_gf),]
sil_0_pft <- sil_0_rawclust[sil_0_rawclust[,2]=="PFT",3:ncol(sil_0_rawclust)]
sil_0_pft <- sil_0_pft[complete.cases(sil_0_pft),]

# load the coefficience of variance data per cluster
cova_0_general <- read.csv(file=file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all.csv")))
cova_0_rawclust <- read.csv(file=file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")))
cova_0_spec <- cova_0_rawclust[cova_0_rawclust[,1]=="Species",2:ncol(cova_0_rawclust)]
cova_0_spec <- cova_0_spec[complete.cases(cova_0_spec),]
cova_0_gen <- cova_0_rawclust[cova_0_rawclust[,1]=="Genus",2:ncol(cova_0_rawclust)]
cova_0_gen <- cova_0_gen[complete.cases(cova_0_gen),]
cova_0_fam <- cova_0_rawclust[cova_0_rawclust[,1]=="Family",2:ncol(cova_0_rawclust)]
cova_0_fam <- cova_0_fam[complete.cases(cova_0_fam),]
cova_0_pg <- cova_0_rawclust[cova_0_rawclust[,1]=="PG",2:ncol(cova_0_rawclust)]
cova_0_pg <- cova_0_pg[complete.cases(cova_0_pg),]
cova_0_gf <- cova_0_rawclust[cova_0_rawclust[,1]=="GF",2:ncol(cova_0_rawclust)]
cova_0_gf <- cova_0_gf[complete.cases(cova_0_gf),]
cova_0_pft <- cova_0_rawclust[cova_0_rawclust[,1]=="PFT",2:ncol(cova_0_rawclust)]
cova_0_pft <- cova_0_pft[complete.cases(cova_0_pft),]


# For Species:
names(sil_0_spec)[1] <- "cluster"
names(rmse_0_spec)[1] <- "cluster"
names(cova_0_spec)[1] <- "cluster"
# merge
tot_spec <- merge(sil_0_spec,cova_0_spec,by = "cluster",all.y = TRUE,all.x = TRUE)
tot_spec <- merge(tot_spec,rmse_0_spec,by = "cluster",all.y = TRUE,all.x = TRUE)
colnames(tot_spec) <- gsub(colnames(tot_spec),pattern = ".x",replacement = "_cova")
colnames(tot_spec) <- gsub(colnames(tot_spec),pattern = ".y",replacement = "_rmse")
tot_spec$Cluster_size_log <- log(tot_spec$Cluster_size)
tot_spec[,grep(names(tot_spec),pattern = "rmse")] <- tot_spec[,grep(names(tot_spec),pattern = "rmse")]^(-.5)
tot_spec$rmse_mean <- rowMeans(tot_spec[,grep(names(tot_spec),pattern = "rmse")])
tot_spec$cova_mean <- rowMeans(tot_spec[,grep(names(tot_spec),pattern = "cova")])

plot(tot_spec)
require(FactoMineR)
tot_spec_m <- as.matrix(tot_spec[,2:ncol(tot_spec)])
mode(tot_spec_m) <- "numeric"
PCA(tot_spec_m[complete.cases(tot_spec_m),])
tot_spec_m3 <- tot_spec_m[,colnames(tot_spec_m)%in%c("Cluster_size","Silhouette_Ind_cova","std","rmse_mean")]
PCA(tot_spec_m3[tot_spec_m3[,1]>1,])
PCA(tot_spec_m3)
PCA_3a <- PCA(tot_spec_m3[tot_spec_m3[,1]>1,],graph = FALSE)
PCA_3b <- PCA(tot_spec_m3,graph = FALSE)
tot_spec_mrest <- tot_spec_m[,!colnames(tot_spec_m)%in%c("Cluster_size","Silhouette_Ind_cova","std")]
PCA(tot_spec_mrest[complete.cases(tot_spec_mrest),])


pdf(file=file.path(origin,"figures","PCA_70.pdf"),widt=10,height=6)
par(mfrow=c(1,2))
  plot(PCA_3a$ind$coord[,1:2],pch=16,main="Including clusters n=1",xlim = c(-4,4),ylim = c(-4,4),col="gray",
       xlab=paste0(round(PCA_3a$eig[1,2],digits=1),"%"),ylab=paste0(round(PCA_3a$eig[2,2],digits=1),"%")) 
  points(0,0,cex=1.2,col="blue",pch=15)
  i=1
  for(i in 1:nrow(PCA_3a$var$coord)){
    lines(x=c(0,PCA_3a$var$coord[i,1]*2),y=c(0,PCA_3a$var$coord[i,2]*2),col="blue")
  }
  text(PCA_3a$var$coord[,1:2]*2,labels = c("Cluster size","Silhouette","SD","Error"),col="blue",cex=1.2)
  
  plot(PCA_3b$ind$coord[,1:2],pch=16,main="Excluding clusters n=1",xlim = c(-4,4),ylim = c(-4,4),col="gray",
       xlab=paste0(round(PCA_3b$eig[1,2],digits=1),"%"),ylab=paste0(round(PCA_3b$eig[2,2],digits=1),"%"))
  points(0,0,cex=1.2,col="blue",pch=15)
  i=1
  for(i in 1:nrow(PCA_3b$var$coord)){
    lines(x=c(0,PCA_3b$var$coord[i,1]*2),y=c(0,PCA_3b$var$coord[i,2]*2),col="blue")
  }
  text(PCA_3b$var$coord[,1:2]*2,labels = c("Cluster size","Silhouette","SD","Error"),col="blue",cex=1.2)
  
dev.off()

t=1
par(mfrow=c(2,2))
  plot(tot_spec$Cluster_size,tot_spec$Silhouette_Ind_cova,ylab="Silhouette index",xlab="cluster size",pch=16)
  plot(tot_spec$Cluster_size,tot_spec$rmse_mean,xlab="cluster size",ylab="error",pch=16,log="y")
  #  plot(tot_spec$Cluster_size,tot_spec$std,ylab="std",xlab="cluster size",pch=16)
  plot(tot_spec$std,tot_spec$Silhouette_Ind_cova,xlab="std",ylab="Silhouette index",pch=16,log="x")
  plot(tot_spec$std,tot_spec$rmse_mean,xlab="std",ylab="error",pch=16,log="y")
  
  par(mfrow=c(2,3))
  for(t in 2:length(grep(colnames(tot_spec),pattern = "cova"))){
    plot(tot_spec$Silhouette_Ind_cova,tot_spec[,grep(colnames(tot_spec),pattern = "cova")[t]],
         main=colnames(tot_spec)[grep(colnames(tot_spec),pattern = "cova")][t],
         xlab="Silhouette index",ylab="error",pch=16,ylim = c(-5,5))
  }

  par(mfrow=c(2,3))
  for(t in 1:5){
    plot(tot_spec$Silhouette_Ind_cova,tot_spec[,grep(colnames(tot_spec),pattern = "rmse")[t]],
         main=colnames(tot_spec)[grep(colnames(tot_spec),pattern = "rmse")[t]],
         xlab="Silhouette index",ylab="error",pch=16,log="y")
  }
  
  par(mfrow=c(2,3))
  for(t in 1:5){
    plot(tot_spec[,grep(colnames(tot_spec),pattern = "cova")[t]]+1000,tot_spec[,grep(colnames(tot_spec),pattern = "rmse")[t]],
         main=colnames(tot_spec)[grep(colnames(tot_spec),pattern = "rmse")[t]],
         xlab="cova",ylab="error",pch=16,log="xy")
  }
  
plot(rmse_0_spec[,2:ncol(rmse_0_spec)],cova_0_spec)
plot(rmse_0_gen[,2:ncol(rmse_0_gen)],cova_0_gen)
plot(rmse_0_fam[,2:ncol(rmse_0_fam)],cova_0_fam)
plot(rmse_0_fam[,2:ncol(rmse_0_fam)],cova_0_fam)
plot(rmse_0_gf[,2:ncol(rmse_0_gf)],cova_0_gf)
plot(sil_0_spec[,2:ncol(sil_0_spec)],rmse_0_spec[,2:ncol(rmse_0_spec)])
plot(sil_0_gen[,2:ncol(sil_0_spec)],rmse_0_gen)
plot(sil_0_fam[,2:ncol(sil_0_spec)],rmse_0_fam)
plot(sil_0_pg[,2:ncol(sil_0_spec)],rmse_0_pg)
plot(sil_0_gf[,2:ncol(sil_0_spec)],rmse_0_gf)
plot(sil_0_pft[,2:ncol(sil_0_spec)],rmse_0_pft)
