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
t_choices <- out$tsubs
TDnos = out$TD_choices
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
gappercents=c(1,5,10,20,30,40,50,60,70,80)
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")


#----------------------------------------------------------------
# get nb of unique individuals first columns
#----------------------------------------------------------------
t_choice="data_2"

{
  RepNum=1
  GapPercent1=0
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
    
    
    if(GapPercent!=-1){traitInfo <- traitInfo_hpmf}
    if(GapPercent==-1){traitInfo <- traitInfo_observed}
    
    traitGappy <- traitInfo_observed
    traitGappy[gap_ix] <- NA
    
    print("scale")
    traitInfo_scaled <- (traitInfo-apply(traitInfo,2,mean))/apply(traitInfo,2,sd)
    traitGappy_scaled <- (traitGappy-apply(traitGappy,2,mean,na.rm=TRUE))/apply(traitGappy,2,sd,na.rm=TRUE)
    
    taxInfo <- TestData_tot[,colnames(TestData_tot)%in%c("IDs","AccSpeciesName","Genus","Family","PhylogeneticGroup")]
    
    # Growth Forms
    match_tax <- match(x = taxInfo$AccSpeciesName,table = GF$AccSpeciesName)
    GF_now <-  GF$GF[match_tax]
    
    # PlantFunctionalTypes
    match_tax <- match(x = taxInfo$AccSpeciesName,table = PFT$AccSpeciesName)
    PFT_now <-  PFT$PFT[match_tax]
    
    
    
  
    new_count <- function(input){
      return(sum(!is.na(input)))
    }
    trait_now=trait_rainfor
    trait_now=trait_guido
    #Species  
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(TestData_tot$AccSpeciesName),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)
    #Genus  
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(TestData_tot$Genus),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
                   round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)  
    #Family  
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(TestData_tot$Family),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
                   round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)  
    #PG 
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(TestData_tot$PhylogeneticGroup),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
                   round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)  
    #GF  
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(GF_now),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
                   round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)  
    #PFT  
    TestData_nb <- aggregate(TestData_tot[,-c(1:5)],by=list(PFT_now),FUN=new_count)
    nb_traits <- TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]
    colSums(nb_traits!=0)  
    cbind(round(colMeans(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)]),digits=2),
          round(apply(TestData_nb[,which(colnames(TestData_nb)%in%trait_now)],2,FUN = sd),digits = 2))
    cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
          round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2))
    colMeans(cbind(round(colMeans(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]),digits=2),
                   round(apply(TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]],2,FUN = sd),digits = 2)))
    add_traits <- TestData_nb[,which(!colnames(TestData_nb)%in%trait_now)[-1]]
    colSums(add_traits!=0)  
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
  TD_choice = "Obs_obs_TD"
  traitInfo_TD <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                       TD_choice,"data","traitInfo.csv"))
  taxInfo_TD <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                       TD_choice,"data","taxInfo.csv"))
  taxInfo_TD <- cbind(taxInfo_TD,GF_now,PFT_now)
  TD_choice = "Obs_obs"
  traitInfo <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                       TD_choice,"data","traitInfo.csv"))
  taxInfo <- read.csv(file=file.path(origin,"runs",paste0("Rep",RepNum), paste0("p","0","_",trait_sub),
                                     TD_choice,"data","taxInfo.csv"))
  taxInfo <- cbind(taxInfo,GF,PFT)
  head(taxInfo)
  dim(traitInfo)
  dim(traitInfo_TD)
  #----------------------------------------------------------------
  # matching taxa in subset and envelope
  #----------------------------------------------------------------
  #Species
  g=1
  for(g in 1:4){
  groups=c("Species","Genus","Family","PhylogeneticGroups")
  group_now=groups[g]
  if(length(which(table_final[,1]%in%group_now))!=1){
    table_final[which(is.na(table_final[,1]))[1],1] <- group_now
  }

  tI_s <- traitInfo[which(taxInfo[,(g+1)]%in%taxInfo_TD[,(g+1)]),]
  xI_s <- taxInfo[which(taxInfo[,(g+1)]%in%taxInfo_TD[,(g+1)]),]
  dim(tI_s)
  tI_match <- tI_s[which(xI_s$IDs%in%taxInfo_TD$IDs),]
  xI_match <- xI_s[which(xI_s$IDs%in%taxInfo_TD$IDs),]
  tI_additionals <- tI_s[-which(xI_s$IDs%in%taxInfo_TD$IDs),]
  xI_additionals <- xI_s[-which(xI_s$IDs%in%taxInfo_TD$IDs),]
  dim(xI_additionals)
  
  table_final[table_final[,1]==group_now,3] <- length(unique(xI_additionals[,(g+1)]))
  table_final[table_final[,1]==group_now,2] <- length(unique(as.vector(taxInfo_TD[,(g+1)])))
  table_final[table_final[,1]==group_now,6] <- length(unique(as.vector(taxInfo_TD$IDs)))
  table_final[table_final[,1]==group_now,4] <- length(unique(as.vector(xI_additionals$IDs)))
  table_final[table_final[,1]==group_now,5] <- length(unique(as.vector(xI_match$IDs)))


  i=1
  for(i in 1:ncol(tI_additionals)){
    trait_now <- colnames(tI_additionals)[i]
    if(length(which(table_final[,1]%in%paste0(trait_now,"_",group_now)))!=1){
    table_final[which(is.na(table_final[,1]))[1],1] <- paste0(trait_now,"_",group_now)
    }
    whichrow = which(table_final[,1]%in%paste0(trait_now,"_",group_now))
    table_final[whichrow,4]  <- colSums(!is.na(tI_additionals))[i]
    table_final[whichrow,3]  <- length(unique(xI_additionals[,(g+1)][!is.na(tI_additionals[,i])]))
    table_final[whichrow,5] <- colSums(!is.na(tI_match))[i]
    if(sum(colnames(traitInfo_TD)%in%trait_now)==1){
      table_final[whichrow,6] <- length(as.vector(traitInfo_TD[,colnames(traitInfo_TD)%in%trait_now]))
      table_final[whichrow,2] <- length(unique(as.vector(taxInfo_TD[,(g+1)])))
    }else{table_final[whichrow,2] <- 0}
  }
  }
  table_final
  #install.packages("xtable")
  require(xtable)
  colnames(table_final)
  col1 <- paste0(table_final[,c(colnames(table_final)=="Nb taxa data set")]," - ",
                 table_final[,c(colnames(table_final)=="Matching taxa")])
  col2 <- paste0(table_final[,c(colnames(table_final)=="Nb indiv data set")]," - ",
                 table_final[,c(colnames(table_final)=="Matching indiv(nb)")])
  col3 <- table_final[,c(colnames(table_final)=="Additional individuals")]
  tablex<- cbind(table_final[,1],col1,col2,col3)

  if(!exists("tab_compare")){  
  tab_compare <- matrix(NA,ncol=13,nrow=5)
  colnames(tab_compare) <- c("Group","Data set 1 (tax nb)","Data set 1 (tax %)","Data set 2 (tax nb)","Data set 2 (tax %)",
                             "Data set 1 (ind nb)","Data set 1 (ind %)","Data set 2 (ind nb)","Data set 2 (ind %)",
                             "Additionall values 1 (nb)","Additionall values 1(%)","Additionall values 2 (nb)","Additionall values 2 (%)")
  }
  if(trait_sub=="rainfor"){
    tab_compare[1,1] <- "Species"
    tab_compare[1,4] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[1,5] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Species",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)
    tab_compare[1,8] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[1,9] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[1,12] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[1,13] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Species",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
    
    tab_compare[2,1] <- "Genus"
    tab_compare[2,4] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[2,5] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)
    tab_compare[2,8] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[2,9] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[2,12] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[2,13] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Genus",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    
    #-------------------------------------
    tab_compare[3,1] <- "Family"
    tab_compare[3,4] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[3,5] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Family",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)
    tab_compare[3,8] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[3,9] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Family",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[3,12] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[3,13] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Family",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
    
    tab_compare[4,1] <- "PhyloGroups"
    tab_compare[4,4] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[4,5] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"PhylogeneticGroups",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)    
    tab_compare[4,8] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[4,9] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"PhylogeneticGroups",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[4,12] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[4,13] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"PhylogeneticGroups",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
  }
  if(trait_sub=="guido"){
    tab_compare[1,2] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[1,3] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Species",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)    
    tab_compare[1,6] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[1,7] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[1,10] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[1,11] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Species")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Species",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
    
    tab_compare[2,2] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[2,3] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)    
    tab_compare[2,6] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[2,7] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[2,10] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[2,11] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Genus")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Genus",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
    
    tab_compare[3,2] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[3,3] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Family",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)    
    tab_compare[3,6] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[3,7] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[3,10] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[3,11] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="Family")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"Family",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
    #-------------------------------------
    
    tab_compare[4,2] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE),digits=0)
    tab_compare[4,3] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching taxa")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"PhylogeneticGroups",c(colnames(table_final)=="Matching taxa")])*100,
                              digits=1)    
    tab_compare[4,6] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE),digits=0)
    tab_compare[4,7] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups"),
                                                          c(colnames(table_final)=="Matching indiv(nb)")]),na.rm = TRUE)/
                                as.numeric(table_final[table_final[,1]%in%"Genus",c(colnames(table_final)=="Matching indiv(nb)")])*100,
                              digits=1)
    tab_compare[4,10] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups")[-1],
                                                          c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE),digits=0)
    tab_compare[4,11] <- round(mean(as.numeric(table_final[grep(table_final[,1],pattern="PhylogeneticGroups")[-1],
                                                           c(colnames(table_final)=="Additional individuals")]),na.rm = TRUE)/
                                 as.numeric(table_final[table_final[,1]%in%"PhylogeneticGroups",colnames(table_final)=="Nb indiv data set"])
                               ,digits=1)
  #-------------------------------------
  }
  col1 <- paste0(tab_compare[,2]," - ", tab_compare[,3])
  col2 <- paste0(tab_compare[,4]," - ", tab_compare[,5])
  col3 <- paste0(tab_compare[,6]," - ", tab_compare[,7])
  col4 <- paste0(tab_compare[,8]," - ", tab_compare[,9])
  col5 <- paste0(tab_compare[,10]," - ", tab_compare[,11])
  col5 <- paste0(tab_compare[,12]," - ", tab_compare[,13])
  tab_compare
  # average 
  xtable(table_final,include.rownames=FALSE)
  print(xtable(tab_compare), include.rownames=FALSE)
  print(xtable(tablex), include.rownames=FALSE)
  length(unique(xI_additionals$AccSpeciesName))
  
  table_final[table_final[,1]=="Species",2] <- length(unique(xI_additionals$AccSpeciesName))
  table_final[table_final[,1]=="Genus",2] <- length(unique(xI_additionals$Genus))
  table_final[table_final[,1]=="Family",2] <- length(unique(xI_additionals$Family))
  table_final[table_final[,1]=="PG",2] <- length(unique(xI_additionals$PhylogeneticGroup))
  table_final[table_final[,1]=="GF",2] <- length(unique(xI_additionals$GF))
  table_final[table_final[,1]=="PFT",2] <- length(unique(xI_additionals$PFT))
  length(unique(xI_additionals$Genus))
  length(unique(xI_additionals$Family))
  length(unique(xI_additionals$PhylogeneticGroup))
  length(unique(xI_additionals$GF))
  length(unique(xI_additionals$PFT))
  
  dim(tI_additionals)
  colSums(!is.na(tI_additionals))
  table_final  
  