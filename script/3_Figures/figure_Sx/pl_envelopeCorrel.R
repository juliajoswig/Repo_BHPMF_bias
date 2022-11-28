
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
            taxInfo_hpmf <- as.data.frame(read.csv(file = tax_path, sep=",", dec="."))
            
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
            
            # load function information
            GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
            PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")
            
            # Growth Forms
            match_tax <- match(x = taxInfo_hpmf$AccSpeciesName,table = GF$AccSpeciesName)
            GF_now <-  GF$GF[match_tax]
            
            # PlantFunctionalTypes
            match_tax <- match(x = taxInfo_hpmf$AccSpeciesName,table = PFT$AccSpeciesName)
            PFT_now <-  PFT$PFT[match_tax]
            
            # load original test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org1) 
            taxInfo_TD <- TestData_tot[,1:5]
            head(taxInfo_TD)
            
            
    


  
  