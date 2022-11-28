
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
Version_now="V3"
list.files(file.path(origin,"script","analysis",Version_now))
#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions()

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
TD_choices =  out$TD_choices
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido =out$trait_guido
trait_rainfor =out$trait_rainfor

source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_add_col_to_res.R"))
source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_logscale.R"))
#source(file.path(origin,"script","analysis",Version_now,"index_functions","fn_RMSE_all.R"))
#calculate_RMSE_function<- function(origin,repnums,whichDataSet,obsspec,gappercents,tsubs,TD_choices){
  
  #-------------------------------------------------------------------
  # Loop through all runs and calculate the indices
  #-------------------------------------------------------------------
  wds=2
  wtd=1
  RepNum=2
  GapPercent=0
  #gappercents=-1
  repnums=2:10
  #whichDataSet=1
  for(RepNum in repnums){
    print("------------------------------")
    print(RepNum)
    print("------------------------------")
    for(wds in whichDataSet){
      for(wtd in obsspec){

        trait_sub = tsubs[wds]
        TD_choice = TD_choices[wtd]
         
        for(GapPercent in gappercents){
          
          # define dat path
          if(GapPercent==-1){GapPercent1=0}else{GapPercent1=GapPercent}
            dat_path_org1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","TestData_org.RData")
            dat_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","HPMF_filled.csv")
            gap_ix_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","Negap.ix.csv")
            std_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","std.csv")
            
          if((file.exists(dat_path)&file.exists(dat_path_org1))|
             (GapPercent==-1&file.exists(dat_path_org1))){

            #--------------------------------------------------------------
            # load data
            #--------------------------------------------------------------
            # load predicted data
            traitInfo_hpmf <- as.matrix(read.table(file = dat_path, sep=",", dec="."))
            # load original test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org1) 
            # load function information
            GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
            PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")
            
            if(GapPercent==-1){
              # fit colnames
              traitInfo_observed <- TestData_tot[,colnames(TestData_tot)%in%colnames(traitInfo_hpmf)]
              NEW_gap_ix <- as.matrix(read.table(file = gap_ix_path,sep=",", dec="."))
              
              # ----------------------------------------------------------------------------------------------------------
              # Rescale the original data
              # ----------------------------------------------------------------------------------------------------------
              # load scaling factors
              load(file.path(origin,"runs",paste0("Rep",RepNum),paste0("Meta_",trait_sub),"scaling_factors.RData"))
                
                if(TD_choice =="Spec_spec"){    
                  for(i in 1:ncol(traitInfo_observed)){
                    traitInfo_observed[,i] <- (log(traitInfo_observed[,i]) - scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_hpmf)[i]),4])/# changed traitInfo to traitInfo_hpmf
                      scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_observed)[i]),3]
                  }
                }else{
                  for(i in 1:ncol(traitInfo_observed)){
                    traitInfo_observed[,i] <- (log(traitInfo_observed[,i]) - scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_hpmf)[i]),2])/# changed traitInfo to traitInfo_hpmf
                      scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_observed)[i]),1]
                  }
                }
            }
            
            if(GapPercent!=-1){traitInfo <- traitInfo_hpmf}
            if(GapPercent==-1){traitInfo <- traitInfo_observed}
            
            taxInfo <- TestData_tot[,colnames(TestData_tot)%in%c("IDs","AccSpeciesName","Genus","Family","PhylogeneticGroup")]
  
            # Growth Forms
            match_tax <- match(x = taxInfo$AccSpeciesName,table = GF$AccSpeciesName)
            GF_now <-  GF$GF[match_tax]
            
            # PlantFunctionalTypes
            match_tax <- match(x = taxInfo$AccSpeciesName,table = PFT$AccSpeciesName)
            PFT_now <-  PFT$PFT[match_tax]
            
            #----------------------------------------------------------------
            
            require(cluster)
            all_out=rep(NA,4)
            data_now=traitInfo
            
            # Function
            {
            print("Function silhouette")
              
            fac_noNA <- rep("is_NA",length(GF_now))
            fac_noNA[!is.na(GF_now)] <- as.vector(GF_now[!is.na(GF_now)])
            fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
            fac_now <- fac_now[fac_now!=0]
              
            cats <- factor(fac_noNA)
            ranks <- rank(-table(cats), ties.method="first")
            DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
            taxnow <- DF[,2]
            
            try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
            try(si_sry <- summary(si2))

            try(silhouette_out <- cbind(rep("GF",length(si_sry$clus.sizes)),
                                        names(fac_now),
                                        si_sry$clus.sizes,
                                        si_sry$clus.avg.widths))
            colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
            
            all_out=rbind(all_out,silhouette_out)
            all_out <- as.matrix(all_out)
            rm("silhouette_out")
            
            # PFT
            fac_noNA <- rep("is_NA",length(PFT_now))
            fac_noNA[!is.na(PFT_now)] <- as.vector(PFT_now[!is.na(PFT_now)])
            fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
            fac_now <- fac_now[fac_now!=0]
            
            cats <- factor(fac_noNA)
            ranks <- rank(-table(cats), ties.method="first")
            DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
            taxnow <- DF[,2]

            try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
            try(si_sry <- summary(si2))
            try(silhouette_out <- cbind(rep("PFT",length(si_sry$clus.sizes)),
                                        names(fac_now),
                                        si_sry$clus.sizes,
                                        si_sry$clus.avg.widths))
            colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
            
            all_out=rbind(all_out,silhouette_out)
            all_out <- as.matrix(all_out)
            
            rm("silhouette_out")
          }

            # Taxonomy
            {
              #Species
              if(TD_choice=="Obs_obs_TD"|TD_choice=="Obs_obs"){
                if(length(unique(taxInfo[,2]))!=nrow(taxInfo)){
                  print("Species silhouette")
                  
                  fac_noNA <- rep("is_NA",length(taxInfo[,2]))
                  fac_noNA[!is.na(taxInfo[,2])] <- as.vector(taxInfo[,2][!is.na(taxInfo[,2])])
                  fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
                  fac_now <- fac_now[fac_now!=0]
                  
                  cats <- factor(fac_noNA)
                  ranks <- rank(-table(cats), ties.method="first")
                  DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
                  taxnow <- DF[,2]
                  head(data_now)
                  
                  try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
                  try(si_sry <- summary(si2))
                  try(silhouette_out <- cbind(rep("Species",length(si_sry$clus.sizes)),
                                              names(fac_now),
                                              si_sry$clus.sizes,
                                              si_sry$clus.avg.widths))
                  colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
                  
                  all_out=rbind(all_out,silhouette_out)
                  all_out <- as.matrix(all_out)
                  rm("silhouette_out")
                }
              }
              
              
              # Genera
              if(length(unique(taxInfo[,3]))!=nrow(taxInfo)){
                print("Genus silhouette")
                
                fac_noNA <- rep("is_NA",length(taxInfo[,3]))
                fac_noNA[!is.na(taxInfo[,3])] <- as.vector(taxInfo[,3][!is.na(taxInfo[,3])])
                fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
                fac_now <- fac_now[fac_now!=0]
                
                cats <- factor(fac_noNA)
                ranks <- rank(-table(cats), ties.method="first")
                DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
                taxnow <- DF[,2]
                
                try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
                try(si_sry <- summary(si2))
                try(silhouette_out <- cbind(rep("Genus",length(si_sry$clus.sizes)),
                                            names(fac_now),
                                            si_sry$clus.sizes,
                                            si_sry$clus.avg.widths))
                colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
                
                all_out=rbind(all_out,silhouette_out)
                all_out <- as.matrix(all_out)
                rm("silhouette_out")
              }
              
              
              # Family
              if(length(unique(taxInfo[,4]))!=nrow(taxInfo)){
                print("Family silhouette")
                
                fac_noNA <- rep("is_NA",length(taxInfo[,4]))
                fac_noNA[!is.na(taxInfo[,4])] <- as.vector(taxInfo[,4][!is.na(taxInfo[,4])])
                fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
                fac_now <- fac_now[fac_now!=0]
                
                cats <- factor(fac_noNA)
                ranks <- rank(-table(cats), ties.method="first")
                DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
                taxnow <- DF[,2]
                
                
                try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
                try(si_sry <- summary(si2))
                try(silhouette_out <- cbind(rep("Family",length(si_sry$clus.sizes)),
                                            names(fac_now),
                                            si_sry$clus.sizes,
                                            si_sry$clus.avg.widths))
                colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
                
                all_out=rbind(all_out,silhouette_out)
                all_out <- as.matrix(all_out)
                
                rm("silhouette_out")
              }
              
              
              # PhylogGroup
              if(length(unique(taxInfo[,5]))!=nrow(taxInfo)){
                print("PG silhouette")
                
                fac_noNA <- rep("is_NA",length(taxInfo[,5]))
                fac_noNA[!is.na(taxInfo[,5])] <- as.vector(taxInfo[,5][!is.na(taxInfo[,5])])
                fac_now <- table(fac_noNA)[order(table(fac_noNA),decreasing = TRUE)]
                fac_now <- fac_now[fac_now!=0]
                
                cats <- factor(fac_noNA)
                ranks <- rank(-table(cats), ties.method="first")
                DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
                taxnow <- DF[,2]
                
                try(si2 <- silhouette(taxnow, dist(data_now,method = "canberra")))
                try(si_sry <- summary(si2))
                try(silhouette_out <- cbind(rep("PG",length(si_sry$clus.sizes)),
                                            names(fac_now),
                                            si_sry$clus.sizes,
                                            si_sry$clus.avg.widths))
                colnames(silhouette_out) <- c("Level","Cluster","Cluster_size","Silhouette_Index")
                
                all_out=rbind(all_out,silhouette_out)
                all_out <- as.matrix(all_out)
                rm("silhouette_out")
                
              }
            }
            #-------------------------------------------------
            # create folder
            if(GapPercent==-1){GapPercent1="org"}else{GapPercent1=GapPercent}
            if(!file.exists(file.path(origin,"data_output","Sil",trait_sub))){
              dir.create(file.path(origin,"data_output","Sil",trait_sub))}
            if(!file.exists(file.path(origin,"data_output","Sil",trait_sub,TD_choice))){
              dir.create(file.path(origin,"data_output","Sil",trait_sub,TD_choice))}
            if(!file.exists(file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1))){
              dir.create(file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1))}
            if(!file.exists(file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1,RepNum))){
              dir.create(file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1,RepNum))}
            
            write.table(all_out,file= 
                          file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")),
                        sep = ",",row.names = TRUE,col.names = TRUE)
            
            
            print(all_out)
            rm("all_out")
            
            
      }
      }
      }
      }
    }


  
  