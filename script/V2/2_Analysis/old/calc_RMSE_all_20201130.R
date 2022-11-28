

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
load_functions(origin = origin,Version_now = Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
  tsubs <- out$tsubs
  ObsOrTDs = out$ObsOrTDs
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
source(file.path(origin,"_2021","script",Version_now,"helper_scripts","fn_add_col_to_res.R"))
#source(file.path(origin,"_2021","script",Version_now,"index_functions","fn_RMSE_all.R"))




#calculate_RMSE_function<- function(origin,repnums,whichDataSet,obsspec,gappercents,tsubs,TD_choices){
  
  #-------------------------------------------------------------------
  # Loop through all runs and calculate the indices
  #-------------------------------------------------------------------

  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  
  GapPercent=50
  RepNum=1
  td=1
  p=2
  TDno=1
  t_choice <- t_choices[td]
  Percent = gappercents[p] 
  ObsOrTD <- TDnos[TDno]

for(RepNum in 3){
  for(td in 1){
    t_choice <- t_choices[td]
    for(p in length(gappercents):2){
      Percent = gappercents[p] 
      for(TDno in 1){
        ObsOrTD <- TDnos[TDno]
        
        print(paste(RepNum,ObsOrTD,t_choice,Percent))
          
          if(GapPercent!=(-1)){
          # define dat path
          dat_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","HPMF_filled.csv")
          std_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","std.csv")
          dat_path_org1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),TD_choice,"data","TestData_org.RData")
          gap_ix_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",trait_sub),paste0(TD_choice),"data","Negap.ix.csv")
          
          if(file.exists(dat_path)&GapPercent!=(-1)&(file.exists(dat_path_org1))){

            #--------------------------------------------------------------
            # load data
            #--------------------------------------------------------------
            # load data
            #try(rm("data_pred"),silent = TRUE)
            try(rm("traitInfo_hpmf"),silent = TRUE)
            
            # load predicted data
            traitInfo_hpmf <-as.matrix(read.table(file = dat_path, sep=",", dec="."))
            # load test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org1) 
            # fit colnames
            traitInfo_observed <- TestData_tot[,colnames(TestData_tot)%in%colnames(traitInfo_hpmf)]
            NEW_gap_ix <- as.matrix(read.table(file = gap_ix_path,sep=",", dec="."))
            # ----------------------------------------------------------------------------------------------------------
            # Rescale
            # ----------------------------------------------------------------------------------------------------------
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
                
            # ----------------------------------------------------------------------------------------------------------
            # RMSE
            # ----------------------------------------------------------------------------------------------------------
            
            if(GapPercent!=(-1)){
              Filled_tot <- traitInfo_hpmf
              Observed_tot <- traitInfo_observed

              
            if(sum(!is.na(Filled_tot))!=0){

              rmse_input <- cbind(Filled_tot,Observed_tot)
              rmse_tot_tot <- rmse(as.vector(Filled_tot),as.vector(unlist((Observed_tot))))
            
              if(GapPercent!=0){
                Filled_tot <- traitInfo_hpmf[NEW_gap_ix]
              Observed_tot <- traitInfo_observed[NEW_gap_ix]
              
              rmse_input <- cbind(Filled_tot,Observed_tot)
              rmse_tot_gap <- rmse(Observed_tot,Filled_tot)
              }else{rmse_tot_gap=NA}
            
            #trait-wise
            rmse_trait_gap <- list()
            rmse_trait_tot <- list()
            tr=1
            for(tr in 1:ncol(traitInfo_hpmf)){
              
              if(GapPercent!=0){
                Filled <- as.vector(traitInfo_hpmf[as.vector(NEW_gap_ix[,tr]),tr])
                Observed <- as.vector(unlist(traitInfo_observed[as.vector(NEW_gap_ix[,tr]),tr]))
                if(sum(!is.na(Filled))!=0){
                  try(rm("rmse_now"))
                  try(rmse_now <- rmse(Observed,Filled))
                  rmse_trait_gap[[tr]] <- rmse_now
                }else{rmse_trait_gap[[tr]] <- NA}
              }else{rmse_trait_gap=rep(NA,ncol(traitInfo_hpmf))}
              
              Filled <- as.vector(traitInfo_hpmf[,tr])
              Observed <- as.vector(unlist(traitInfo_observed[,tr]))
              if(sum(!is.na(Filled))!=0){
                try(rm("rmse_now"))
                plot(as.vector(Filled),as.vector(unlist((Observed))))
                abline(0,1)
                try(rmse_now <- rmse(Observed,Filled))
                rmse_trait_tot[[tr]] <- rmse_now
              }else{rmse_trait_tot[[tr]] <- NA}
            }
            
            rmse_man <- matrix(NA,ncol=ncol(traitInfo_hpmf),nrow=4)
            rmse_man[1,] <- unlist(rmse_trait_tot)
            rmse_man[2,] <- unlist(rmse_trait_gap)
            rownames(rmse_man) <- c("rmse_trait_tot","rmse_trait_gap","rmse_tot_gap","rmse_tot_tot")
            rmse_man[3,1] <- unlist(rmse_tot_gap)
            rmse_man[4,1] <- unlist(rmse_tot_tot)
            colnames(rmse_man) <- colnames(traitInfo_hpmf)
            
            # create folder
            if(!file.exists(file.path(origin,"data_output","RMSE",trait_sub))){
                dir.create(file.path(origin,"data_output","RMSE",trait_sub))}
            if(!file.exists(file.path(origin,"data_output","RMSE",trait_sub,TD_choice))){
                dir.create(file.path(origin,"data_output","RMSE",trait_sub,TD_choice))}
            if(!file.exists(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent))){
                dir.create(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent))}
            if(!file.exists(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent,RepNum))){
                dir.create(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent,RepNum))}
              
            print(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent,RepNum))
            print(rmse_man)
            
            write.csv(rmse_man,
                      file=file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent,RepNum,paste0("RMSE.csv")))
            print(file.path(origin,"data_output","RMSE",trait_sub,TD_choice,GapPercent,RepNum))
            }     
      }
      }
      }
      }
    }

  }
  }
  
  
  