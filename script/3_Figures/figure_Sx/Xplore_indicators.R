
#------------------------------------------------------------
# define path
#------------------------------------------------------------
is.it.on.cluster=TRUE
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
require(cluster)

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
  
  whichDataSet=1:2
  wds=2
  wtd=3
  repnums=sample(1:30,30)
  GapPercent=10
  RepNum=1
  for(RepNum in repnums){
    print("------------------------------")
    print(RepNum)
    print("------------------------------")
    for(wds in whichDataSet){
      for(wtd in obsspec){

        trait_sub = tsubs[wds]
        TD_choice = TD_choices[wtd]
        
         print(paste0(trait_sub," ",TD_choice))
         print("------------------------------")
         for(GapPercent in gappercents){
         
          { 
          # define dat path
          if(GapPercent==-1){GapPercent1=0}else{GapPercent1=GapPercent}
            dat_path_tot <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","traitInfo.csv")
            dat_path_org <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),TD_choice,"data","TestData_org.RData")
            dat_path_org1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","TestData_org.RData")
            dat_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","mean.csv")
            list.files(file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data"))
            dat_path1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","HPMF_filled.csv")
            gap_ix_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","Negap.ix.csv")
            tax_path=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","taxInfo.csv")
            trait_path_org <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),TD_choice,"data", "traitInfo_TD.csv")
            trait_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data", "traitInfo_TD.csv")
            print(trait_path)
          if((file.exists(dat_path)&file.exists(dat_path_org1)&file.exists(trait_path))|
             (GapPercent==-1&file.exists(dat_path_org1)&file.exists(trait_path))){

            #--------------------------------------------------------------
            # load data
            #--------------------------------------------------------------
            # load predicted data
            traitInfo_hpmf_TD <- as.matrix(read.table(file = dat_path1, sep=",", dec="."))
            traitInfo_hpmf <- as.matrix(read.csv(file = dat_path, sep="\t", dec="."))
            traitInfo_input <- as.matrix(read.csv(file = dat_path_tot, sep=",", dec="."))
            
            # load original tax data
            tax_path_hpmf <- as.data.frame(read.csv(file = tax_path, sep=",", dec="."))
            
            # load original 0% test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org)
            head(traitInfo_observed)
            
            #Observed data
            traitInfo_observed <- TestData_tot[,colnames(TestData_tot)%in%colnames(traitInfo_hpmf)]
            #needs to get scaled!
            # Rescale the original data
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
            
            if(GapPercent!=-1){traitInfo <- traitInfo_hpmf}
            if(GapPercent==-1){traitInfo <- traitInfo_observed}
            
            # load original test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org1) 
            # load function information
            GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
            PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")
            # load gapIX
            trait_tot <- as.matrix(read.table(file = dat_path_tot, sep=",", dec="."))
            org_tot <- as.matrix(read.table(file = trait_path_org, sep=",", dec="."))
            test_tot <- as.matrix(read.table(file = trait_path, sep=",", dec="."))
            gap_org <- is.na(org_tot)
            gap_now <- is.na(test_tot)
            gap_ix=gap_org!=gap_now
            gap_ix <- gap_ix[,colnames(gap_ix)%in%colnames(traitInfo_hpmf)]
            out_list <- list()
            
            
          }
            
            #----------------------------------------------------------------
            # ok we want to get 
            # 1. CV obs data set 1/2 ENVELOPE
            # 2. deviation CV pred-obs (data set_envelope) Irgendwo schon gemacht!
            # 3. CV obs data set 1/2
            covaTD <- read.csv(file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")))
            #SilTD <- read.csv(file.path(origin,"data_output","Sil",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")),row.names = NULL)
            #RMSETD <- read.csv(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent1,RepNum,paste0("RMSE.csv")),row.names = NULL)
            list.files(file.path(origin,"data_output"))
            {           
            all_out=rep(NA,6)
            data_now=traitInfo_scaled
            data_now2=traitGappy_scaled
            
            
            new_mean_fun=function(x){
              if(sum(!is.na(x))>0){
                out=mean(x[!is.na(x)])      
              }else{out=NA}
            }
            new_sd_fun=function(x){
              if(sum(!is.na(x))>0){
                out=sd(x[!is.na(x)])      
              }else{out=NA}
            }
            
            group_names=c("GF","PFT","Species","Genus","Family","PG")
            fctrs=list(GF_now=GF_now,
                         PFT_now=PFT_now,
                       taxInfo_spec=taxInfo[,2],
                       taxInfo_gen=taxInfo[,3],
                       taxInfo_fam=taxInfo[,4],
                       taxInfo_pg=taxInfo[,5])

            g=1
            for(g in 1:length(group_names)){
            print(paste0(group_names[g]," CoVa"))
            
            cats <- factor(fctrs[[g]])
            ranks <- rank(-table(cats), ties.method="first")
            DF <- data.frame(category=cats, rank=ranks[as.character(cats)])
            taxnow <- DF[,2]
            
            mean_tax <- aggregate(data_now,by=list(fctrs[[g]]),FUN=new_mean_fun)
            sd_tax <- aggregate(data_now,by=list(fctrs[[g]]),FUN=new_sd_fun)
            mean_tax2 <- aggregate(data_now2,by=list(fctrs[[g]]),FUN=new_mean_fun)
            sd_tax2 <- aggregate(data_now2,by=list(fctrs[[g]]),FUN=new_sd_fun)
            ix=!is.na(sd_tax)
            gaps <- aggregate(gap_ix,by=list(fctrs[[g]]),FUN=sum)
            
            cova_tmp <- sd_tax[,2:ncol(sd_tax)]/mean_tax[,2:ncol(mean_tax)]
            cova_tmp <- cbind(rep(group_names[g],nrow(cova_tmp)),mean_tax$Group.1,cova_tmp)
            colnames(cova_tmp)[1:2] <- c("Group","Cluster")
            out_list[[g]] <- group_names[g]
            out_list[[g]]$pred <- list(cova_tmp=cova_tmp)
            
            cova_tmp <- sd_tax2[,2:ncol(sd_tax2)]/mean_tax2[,2:ncol(mean_tax2)]
            cova_tmp <- cbind(rep(group_names[g],nrow(cova_tmp)),mean_tax2$Group.1,cova_tmp)
            colnames(cova_tmp)[1:2] <- c("Group","Cluster")
            out_list[[g]]$obs <- list(cova_tmp=cova_tmp)

            out_list[[g]]$gaps <- list(gaps=gaps)
            
            for(clmns in 2:ncol(sd_tax)){
              try(cv_sry <- summary(sd_tax[ix[,clmns],clmns]/mean_tax[ix[,clmns],clmns]))
              try(cova_out <- t(as.matrix(cv_sry)))
              all_out=rbind(all_out,cova_out)
              all_out <- as.matrix(all_out)
              rownames(all_out)[nrow(all_out)] <- paste0(group_names[g],"_",colnames(sd_tax)[clmns])
            }
            }
            
            names(out_list) <-  group_names
            
            
            for(g in 1:length(group_names)){
              colnames(out_list[[g]]$obs$cova_tmp) <- paste0(colnames(out_list[[g]]$obs$cova_tmp),"_obs")
            }
            
            if(trait_sub=="rainfor"){out_matrix=rep(NA,20)}
            if(trait_sub=="guido"){out_matrix=rep(NA,17)}
            
              for(g in 1:length(group_names)){
                out_matrix <- rbind(out_matrix,
                                cbind(out_list[[g]]$pred$cova_tmp, 
                                out_list[[g]]$obs$cova_tmp[,3:ncol(out_list[[g]]$obs$cova_tmp)],
                                out_list[[g]]$gaps$gaps[,2:ncol(out_list[[g]]$gaps$gaps)]))
                }
            
            print(GapPercent)
            print(summary(out_matrix[,c(3,8,13)]))
            print(summary(out_matrix[,c(3,9,15)]))
            }
            #----------------------------------------------------------------
            # create folder
            if(GapPercent==-1){GapPercent1="org"}else{GapPercent1=GapPercent}
            if(!file.exists(file.path(origin,"data_output","CoVa"))){
              dir.create(file.path(origin,"data_output","CoVa"))}
            if(!file.exists(file.path(origin,"data_output","CoVa",trait_sub))){
              dir.create(file.path(origin,"data_output","CoVa",trait_sub))}
            if(!file.exists(file.path(origin,"data_output","CoVa",trait_sub,TD_choice))){
              dir.create(file.path(origin,"data_output","CoVa",trait_sub,TD_choice))}
            if(!file.exists(file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1))){
              dir.create(file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1))}
            if(!file.exists(file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum))){
              dir.create(file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum))}
            
            write.table(all_out,file= 
                          file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all.csv")),
                        sep = ",",row.names = TRUE,col.names = TRUE)
            write.table(out_matrix,file= 
                          file.path(origin,"data_output","CoVa",trait_sub,TD_choice,GapPercent1,RepNum,paste0("all_clusters.csv")),
                        sep = ",",row.names = TRUE,col.names = TRUE)
           
            rm("all_out")
          }
            
      }
      }
      }
      }
    


  
  