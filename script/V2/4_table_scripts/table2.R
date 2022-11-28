
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
colz1 =out$colz1
colz2 =out$colz2



#----------------------------------------------------------------
# get nb of unique individuals first columns
#----------------------------------------------------------------
trait_sub="guido"


  RepNum=1
  GapPercent1=0
  GapPercent=0
  TD_choice="Obs_obs_TD"
  
  
  #load tax and test data
  
    dat_path_tot <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","traitInfo.csv")
    dat_path_org <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),TD_choice,"data","TestData_org.RData")
    dat_path_org1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","TestData_org.RData")
    dat_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","HPMF_filled.csv")
    gap_ix_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","Negap.ix.csv")
    tax_path=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","taxInfo.csv")
    trait_path_org <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),TD_choice,"data", "traitInfo_TD.csv")
    trait_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data", "traitInfo_TD.csv")
    print(trait_path)
  
    #--------------------------------------------------------------
    # load data
    #--------------------------------------------------------------
    # load predicted data
    traitInfo_hpmf <- as.matrix(read.table(file = dat_path, sep=",", dec="."))
    # load original 0% test data "TestData_tot" containing tax, info and ID data
    load(dat_path_org) 
    #Observed data
    traitInfo_observed <- TestData_tot[,colnames(TestData_tot)%in%colnames(traitInfo_hpmf)]
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
    
    
    if(GapPercent!=-1){traitInfo_TD <- traitInfo_hpmf}
    if(GapPercent==-1){traitInfo_TD <- traitInfo_observed}
    
    traitGappy_TD <- traitInfo_observed
    traitGappy_TD[gap_ix] <- NA
    
    print("scale")
    traitInfo_TD_scaled <- (traitInfo_TD-apply(traitInfo_TD,2,mean))/apply(traitInfo_TD,2,sd)
    traitGappy_TD_scaled <- (traitGappy_TD-apply(traitGappy_TD,2,mean,na.rm=TRUE))/apply(traitGappy_TD,2,sd,na.rm=TRUE)
    
    taxInfo_TD <- TestData_tot[,colnames(TestData_tot)%in%c("IDs","AccSpeciesName","Genus","Family","PhylogeneticGroup")]
    
    # Growth Forms
    match_tax <- match(x = taxInfo_TD$AccSpeciesName,table = GF$AccSpeciesName)
    GF_now <-  GF$GF[match_tax]
    taxInfo_TD$GF <- GF_now
    
    # PlantFunctionalTypes
    match_tax <- match(x = taxInfo_TD$AccSpeciesName,table = PFT$AccSpeciesName)
    PFT_now <-  PFT$PFT[match_tax]
    taxInfo_TD$PFT <- PFT_now
    
    
  
    new_count <- function(input){
      return(sum(!is.na(input)))
    }
  #----------------------------------------------------------------
  # get nb of matching individuals second
  #----------------------------------------------------------------
  table_final <- matrix(NA,ncol = 6,nrow=100)
  colnames(table_final) <- c("trait","Nb taxa data set 1","Matching taxa","Additional individuals","Matching indiv(nb)","Nb indiv data set")
  trait_names <- trait_guido
  list.files(file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),TD_choice,"data"))
  #----------------------------------------------------------------
  # load
  #----------------------------------------------------------------
#  TD_choice = "Obs_obs_TD"
#  traitInfo_TD <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
#                                       TD_choice,"data","traitInfo.csv"))
#  taxInfo_TD <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
#                                       TD_choice,"data","taxInfo.csv"))
#  taxInfo_TD <- cbind(taxInfo_TD,GF_now,PFT_now)
  TD_choice = "Obs_obs"
  traitInfo <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                       TD_choice,"data","traitInfo.csv"))
  taxInfo <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                     TD_choice,"data","taxInfo.csv"))
  taxInfo$PFT <- PFT$PFT
  taxInfo$GF <- GF$GF
  taxInfo <- cbind(taxInfo,GF,PFT)
  head(taxInfo)
  dim(traitInfo)
  dim(traitInfo_TD)
  #----------------------------------------------------------------
  # nb of ind/spec/gen/fam/PG/
  #----------------------------------------------------------------
  #Species
require("xtable")
  if(!exists("res")){
    res <- matrix(NA,ncol=4,nrow=20)
    colnames(res) <- c("","Data set 1", "Data set 2", "Envelope data")  }


res[1,1] <-"trait number"
res[2,1] <-"trait names" 
res[3,1] <-"Individuals" 
res[4,1] <-"Species" 
res[5,1] <-"Genera" 
res[6,1] <-"Families" 
res[7,1] <-"Phylogenetic groups" 
res[8,1] <-"" 
res[9,1] <-"Gowth forms" 
res[10,1] <-"tree" 
res[11,1] <-"shrub" 
res[12,1] <-"herb" 
res[13,1] <-"graminoid" 
res[14,1] <-"fern" 
res[15,1] <-"none" 
res[16,1] <-"other" 
res[17,1] <-"PFTs"  

if(trait_sub=="guido"){cn=2}
if(trait_sub=="rainfor"){cn=3}

res[1,cn] <- ncol(traitInfo_TD)
res[2,cn] <- paste(colnames(traitInfo_TD),collapse=" ")
res[3,cn] <- length(unique(taxInfo_TD$IDs))
res[4,cn] <-length(unique(taxInfo_TD$AccSpeciesName))
res[5,cn] <- length(unique(taxInfo_TD$Genus))
res[6,cn] <- length(unique(taxInfo_TD$Family))
res[7,cn] <- length(unique(taxInfo_TD$PhylogeneticGroup))
res[9,cn] <-length(unique(taxInfo_TD$GF))
res[10,cn] <-sum(taxInfo_TD$GF=="tree")
res[11,cn] <-sum(taxInfo_TD$GF=="shrub") 
res[12,cn] <- sum(taxInfo_TD$GF=="herb")
res[13,cn] <- sum(taxInfo_TD$GF=="graminoid")
res[14,cn] <- sum(taxInfo_TD$GF=="fern")
res[15,cn] <- sum(taxInfo_TD$GF=="none")
res[16,cn] <- nrow(taxInfo_TD)-sum(as.numeric(res[10:15,cn]))
res[17,cn] <-length(unique(taxInfo_TD$PFT))
  
cn=4  
  res[1,cn] <- ncol(traitInfo)
  res[2,cn] <- paste(colnames(traitInfo),collapse=" ")
  res[3,cn] <- length(unique(taxInfo$IDs))
  res[4,cn] <-length(unique(taxInfo$AccSpeciesName))
  res[5,cn] <- length(unique(taxInfo$Genus))
  res[6,cn] <- length(unique(taxInfo$Family))
  res[7,cn] <- length(unique(taxInfo$PhylogeneticGroup))
  res[8,cn] <- sum(!is.na(traitInfo))/(18*nrow(traitInfo))*100
  res[9,cn] <-length(unique(taxInfo$GF))
  res[10,cn] <-sum(taxInfo$GF=="tree")
  res[11,cn] <-sum(taxInfo$GF=="shrub") 
  res[12,cn] <- sum(taxInfo$GF=="herb")
  res[13,cn] <- sum(taxInfo$GF=="graminoid")
  res[14,cn] <- sum(taxInfo$GF=="fern")
  res[15,cn] <- sum(taxInfo$GF=="none")
  res[16,cn] <- nrow(taxInfo)-sum(as.numeric(res[10:15,cn]))
  res[17,cn] <-length(unique(taxInfo$PFT))
  
  
length(unique(taxInfo_TD$IDs))  
  length(unique(taxInfo_TD$IDs))  
  length(unique(taxInfo_TD$IDs))  
  print(xtable(res),include.rownames = FALSE)
  
  