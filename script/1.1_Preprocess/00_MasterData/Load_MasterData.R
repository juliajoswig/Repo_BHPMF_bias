

is.it.on.cluster=FALSE
if(is.it.on.cluster){
  setwd("/..")
  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("work_1","2016_GapFilling")}
if(!is.it.on.cluster){
  setwd("/..")
  origin = "Volumes/bgi-1/work_1/2016_GapFilling"
  origin = "Volumes/bgi/work_1/2016_GapFilling"
}
Version_now="V1"


convert_to_point <- function(input){
  return(as.numeric(gsub(",", ".", gsub("\\.", "", input))))
}
trait_rainfor =  c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
trait_guido =  c("SLA","PlantHeight","SeedMass","LDMC","LeafArea")#c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")

list.files(file.path(origin,"_2021","data","MasterData","TRY_HPMF_2013"))
list.files(file.path(origin,"_2021","data","MasterData","TRY_HPMF_2013","Data"))
#install.packages("readxl")
require("readxl")
#------------------------------------------
# load  
#------------------------------------------

list.files(file.path(origin,"_2021","data","MasterData","TRY_HPMF_2013","Data"))
MasterData2013cleaned <- read_excel(path = file.path(origin,"_2021","data","MasterData","TRY_HPMF_2013","Data","TRY_pmf_DataRelease_cleaned_2013_06_19.xlsx"))

colnames(MasterData2013cleaned)
MasterData2013_c <- MasterData2013cleaned[,colnames(MasterData2013cleaned)%in%c(
  "ObservationID","Lat","Lon","Alt","PlantDevState",                                      
  "Exposition","AccSpeciesID","AccSpeciesName","Type", "Genus","Family", "PhylogeneticGroup","PlantGrowthForm" ,
  "Succulent","Climbing" ,                                               
 "LeafType","LeafPhenology", "PhotosyntheticPathway","Woodiness","LeafCompoundness","KoeppenGeigerCode")]
 MasterData2013_t <- MasterData2013cleaned[, colnames(MasterData2013cleaned)%in%
                                             c("SLA","PlantHeight", "SeedMass","LDMC","StemSpecificDensity","LeafArea","LeafN","LeafP","Stem conduit density (vessels and tracheids)", 
                                             "Seed number per reproducton unit","Wood vessel element length",             
 "Leaf nitrogen (N) content per area","Leaf fresh mass","Leaf nitrogen/phosphorus (N/P) ratio",
 "Leaf carbon (C) content per dry mass","Leaf delta 15N","Seed length","Dispersal unit length")]
MasterData2013 <- read.csv(file = file.path(origin,"_2021","data","MasterData","TRY_HPMF_2013","Data","TRY_pmf_DataRelease2012_08_25.csv"))
mtch <- match(MasterData2013$ObservationID,MasterData2013_c$ObservationID)
MasterData <- cbind(MasterData2013_c[mtch,],MasterData2013)
summary(MasterData)

#------------------------------------------
# Rename
#------------------------------------------
  colnames(MasterData)[colnames(MasterData)=="StemSpecificDensity"] <- "SSD"
  colnames(MasterData)[colnames(MasterData)=="N"] <- "LeafN"
  colnames(MasterData)[colnames(MasterData)=="P"] <- "LeafP"
  colnames(MasterData)[colnames(MasterData)=="NperArea"] <- "LeafNArea"
#  colnames(MasterData)[colnames(MasterData)=="Leaf nitrogen (N) content per area"] <- "LeafNArea"
  colnames(MasterData)[colnames(MasterData)=="PlantHeight"] <- "PlantHeight"
  colnames(MasterData)[colnames(MasterData)=="N_Pratio"] <- "LeafNPratio"
  colnames(MasterData)[colnames(MasterData)=="Area"] <- "LeafArea"
#  colnames(MasterData)[colnames(MasterData)=="Leaf nitrogen/phosphorus (N/P) ratio"] <- "LeafNPratio"
  #  colnames(MasterData)[colnames(MasterData)=="Leaf carbon (C) content per dry mass"] <- "LeafC"
  colnames(MasterData)[colnames(MasterData)=="CperDryMass"] <- "LeafC"
  colnames(MasterData)[colnames(MasterData)=="Seed length"] <- "SeLen"
  colnames(MasterData)[colnames(MasterData)=="FreshMass"] <- "LeafFMass"
#  colnames(MasterData)[colnames(MasterData)=="Leaf fresh mass"] <- "LeafFMass"
  colnames(MasterData)[colnames(MasterData)=="DUL"] <- "DispUlen"
#  colnames(MasterData)[colnames(MasterData)=="Dispersal unit length"] <- "DispUlen"
#  colnames(MasterData)[colnames(MasterData)=="Leaf delta 15N"] <- "Leafdelta15N"
  colnames(MasterData)[colnames(MasterData)=="delta15N"] <- "Leafdelta15N"
#  colnames(MasterData)[colnames(MasterData)=="Seed number per reproducton unit"] <- "SeNbU"
  colnames(MasterData)[colnames(MasterData)=="SeedNr_reproductonUnit"] <- "SeNbU"
  colnames(MasterData)[colnames(MasterData)=="StemConduitDensity"] <- "ConduitDens"
#  colnames(MasterData)[colnames(MasterData)=="Stem conduit density (vessels and tracheids)"] <- "ConduitDens"
  #  colnames(MasterData)[colnames(MasterData)=="SeedNr_reproductonUnit"] <- "SeNbU"
#  colnames(MasterData)[colnames(MasterData)=="Wood vessel element length"] <- "VesLen"
  colnames(MasterData)[colnames(MasterData)=="WVEL"] <- "VesLen"
  colnames(MasterData)[colnames(MasterData)=="Area"] <- "LeafArea"
  colnames(MasterData)
  summary(as.data.frame(MasterData2013_t))
  
  #------------------------------------------
  # correct
  #------------------------------------------
  MasterData <- as.data.frame(MasterData)
  MasterData$AccSpeciesName[which(MasterData$PlantHeight==160)]
  MasterData[which(MasterData$PlantHeight==160),colnames(MasterData)=="PlantHeight"] <- as.numeric(1.6)
  MasterData$AccSpeciesName[which(MasterData$PlantHeight==150)]
  MasterData[which(MasterData$PlantHeight==150),colnames(MasterData)=="PlantHeight"] <- as.numeric(0.015)
  MasterData$AccSpeciesName[which(MasterData$PlantHeight==120)]
  
  #------------------------------------------
  # kick out zeros
  #------------------------------------------
  # cut down if tax missing
  MasterTax <-  MasterData[,colnames(MasterData)%in%
                              c("ObservationID","AccSpeciesName","Genus",
                                "Family","PhylogeneticGroup")]
  cc_tax <- rowSums(!is.na(MasterTax))!=0
  MasterData <- MasterData[cc_tax,]

  MasterTrait <-  MasterData[,colnames(MasterData)%in%
                             c("SLA","PlantHeight",
                                "SeedMass","LDMC","SSD","LeafArea",             
                                "LeafN","LeafP","ConduitDens","SeNbU" ,               
                                "VesLen","LeafNArea","LeafFMass","LeafNPratio" ,         
                                "LeafC","Leafdelta15N","SeedLength","DispUlen")]
  cc_trait <- rowSums(!is.na(MasterTrait))!=0
  MasterData <- MasterData[cc_trait,]
  
  guido <- MasterData[,colnames(MasterData)%in%c("ObservationID","Lat","Lon","Alt","PlantDevState",                                      
                                                 "Exposition","AccSpeciesID","AccSpeciesName","Type", "Genus","Family", "PhylogeneticGroup","PlantGrowthForm" ,
                                                 "Succulent","Climbing" ,                                               
                                                 "LeafType","LeafPhenology", "PhotosyntheticPathway","Woodiness","LeafCompoundness","KoeppenGeigerCode",trait_guido)]
  dim(guido[complete.cases(guido[,colnames(guido)%in%trait_guido]),])
  guido <- guido[complete.cases(guido[,colnames(guido)%in%trait_guido]),]
  colnames(guido)
  head(guido)
  rainfor <- MasterData[,colnames(MasterData)%in%c("ObservationID","Lat","Lon","Alt","PlantDevState",                                      
                                                   "Exposition","AccSpeciesID","AccSpeciesName","Type", "Genus","Family", "PhylogeneticGroup","PlantGrowthForm" ,
                                                   "Succulent","Climbing" ,                                               
                                                   "LeafType","LeafPhenology", "PhotosyntheticPathway","Woodiness","LeafCompoundness","KoeppenGeigerCode",trait_rainfor)]
  rainfor <- rainfor[complete.cases(rainfor[,colnames(rainfor)%in%trait_rainfor]),]
  dim(rainfor)
  head(rainfor)
  colnames(rainfor)

  #------------------------------------------
  # write
  #------------------------------------------
  write.csv(x = rainfor,file=file.path(origin,"_2021","data","MasterData","data","Obs_obs_TD","MasterData.csv"))
  write.csv(x = MasterData,file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData.csv"))

  write.csv(x = guido,file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs_TD","MasterData.csv"))
  write.csv(x = MasterData,file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs","MasterData.csv"))

# same but old Version: Without Taxonomic information



# LOAD TAXA !
