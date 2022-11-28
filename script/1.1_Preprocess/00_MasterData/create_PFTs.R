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
#Version_now="V1"
gappercents=c(0,1,5,10,20,30,40,50,60,70,80,90,100)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
trait_rainfor =  c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
trait_guido =  c("SLA","PlantHeight","SeedMass","LDMC","LeafArea")#c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")

data_rainforTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs_TD","MasterData.csv"))
data_rainforTD <- data_rainforTD[,-grep(colnames(data_rainforTD),pattern = "X")]
colnames(data_rainforTD)
data_rainfor <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData.csv"))
data_rainfor <- data_rainfor[,-grep(colnames(data_rainfor),pattern = "X")]
colnames(data_rainfor)

data_guidoTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs_TD","MasterData.csv"))
data_guidoTD <- data_guidoTD[,-grep(colnames(data_guidoTD),pattern = "X")]
colnames(data_guidoTD)
data_guido <- read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs","MasterData.csv"))
data_guido <- data_guido[,-grep(colnames(data_guido),pattern = "X")]
colnames(data_guido)

taxcols=colnames(data_guido)%in%c("ObservationID","AccSpeciesName","Genus","Family","PhylogeneticGroup")
funcols=colnames(data_guido)%in%c("PlantGrowthForm")
#Do PFT
PFT_guido=rep(NA,nrow(data_guido))
PFT_guidoTD=rep(NA,nrow(data_guidoTD))
PFT_rainfor=rep(NA,nrow(data_rainfor))
PFT_rainforTD=rep(NA,nrow(data_rainforTD))


GF=as.vector(unique(data_guido$PlantGrowthForm))
PW=c("C3","C4","CAM")

i=1
j=1
for(i in 1:2){
for(j in 1:2){
  t_choice=t_choices[i]
  TDno=TDnos[j]
  
  if(t_choice=="data"&TDno=="Obs_obs_TD"){
    PFT_now=PFT_rainforTD
    data_now=data_rainforTD
  }
  if(t_choice=="data"&TDno=="Obs_obs"){
    PFT_now=PFT_rainfor
    data_now=data_rainfor
  }
  if(t_choice=="data_2"&TDno=="Obs_obs_TD"){
    PFT_now=PFT_guidoTD
    data_now=data_guidoTD
  }
  if(t_choice=="data_2"&TDno=="Obs_obs"){
    PFT_now=PFT_guido
    data_now=data_guido
  }
  print(dim(data_now))
  print(length(PFT_now))
  g=1
p=1
for(g in 1:length(GF)){
  for(p in 1:length(PW)){
    PFT_now[data_now$PlantGrowthForm==GF[g]&
                data_now$LeafPhenology=="evergreen"&
                data_now$LeafType=="needleleaved"&is.na(PFT_now)] <- paste0("NeedleleafEvergreen",GF[g])
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafPhenology=="evergreen"&
                  data_now$LeafType=="broadleaved"&is.na(PFT_now)] <- paste0("BroadleafEvergreen",GF[g])
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafPhenology=="deciduous"&
                  data_now$LeafType=="broadleaved"&is.na(PFT_now)] <- paste0("NeedleleafDeciduous",GF[g])
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafPhenology=="deciduous"&
                  data_now$LeafType=="needleleaved"&is.na(PFT_now)] <- paste0("NeedleleleafDeciduous",GF[g])
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafPhenology=="deciduous"&
                  data_now$LeafType=="needleleaved"&is.na(PFT_now)] <- paste0("NeedleleleafDeciduous",GF[g])
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$Climbing=="climber"&
                  is.na(PFT_now)]     <- paste0("climber",GF[g]) 
  
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafType=="needleleaved"&
                  data_now$PhotosyntheticPathway==PW[p]&
                  is.na(PFT_now)] <-paste0("needleleaved",PW[p],GF[g]) 
  PFT_now[data_now$PlantGrowthForm==GF[g]&
                  data_now$LeafType=="broadleaved"&
                  data_now$PhotosyntheticPathway==PW[p]&
                  is.na(PFT_now)] <-paste0("broadleaved",PW[p],GF[g]) 
}
}

print(dim(data_now))
print(length(PFT_now))
# Rest lower rank
for(g in 1:length(GF)){
  for(p in 1:length(PW)){
    PFT_now[data_now$PlantGrowthForm==GF[g]&
                    data_now$PhotosyntheticPathway==PW[p]&
                    is.na(PFT_now)] <-paste0(PW[p],GF[g]) 
    PFT_now[data_now$PlantGrowthForm==GF[g]&
                    data_now$PhotosyntheticPathway==PW[p]&
                    is.na(PFT_now)] <-paste0(PW[p],GF[g]) 
    PFT_now[data_now$PlantGrowthForm==GF[g]&
                data_now$LeafType=="needleleaved"&is.na(PFT_now)] <- paste0("Needleleaf",GF[g])
    PFT_now[data_now$PlantGrowthForm==GF[g]&
                data_now$LeafType=="broadleaved"&is.na(PFT_now)] <- paste0("Broadleaf",GF[g])
}
}
print(dim(data_now))
print(length(PFT_now))

for(g in 1:length(GF)){
  #lowest rank == GF
  PFT_now[data_now$PlantGrowthForm==GF[g]&is.na(PFT_now)] <- paste0(GF[g])
}

if(t_choice=="data"&TDno=="Obs_obs_TD"){PFT_rainforTD=PFT_now;data_rainforTD <- cbind(PFT_rainforTD,data_rainforTD)}
if(t_choice=="data"&TDno=="Obs_obs"){PFT_rainfor=PFT_now;data_rainfor <- cbind(PFT_rainfor,data_rainfor)}
if(t_choice=="data_2"&TDno=="Obs_obs_TD"){PFT_guidoTD=PFT_now;data_guidoTD <- cbind(PFT_guidoTD,data_guidoTD)}
if(t_choice=="data_2"&TDno=="Obs_obs"){PFT_guido=PFT_now;data_guido <- cbind(PFT_guido,data_guido)}
  
}
}

data_rainforTD <- data_rainforTD[,c(2,1,3:ncol(data_rainforTD))]
data_rainfor <- data_rainfor[,c(2,1,3:ncol(data_rainfor))]
data_guidoTD <- data_guidoTD[,c(2,1,3:ncol(data_guidoTD))]
data_guido <- data_guido[,c(2,1,3:ncol(data_guido))]

colnames(data_guido)[2] <- "PFT"
colnames(data_guidoTD)[2] <- "PFT"
colnames(data_rainfor)[2] <- "PFT"
colnames(data_rainforTD)[2] <- "PFT"


write.csv(data_rainforTD,file=file.path(origin,"_2021","data","MasterData","data","Obs_obs_TD","MasterData_inclPFT.csv"))
write.csv(data_rainfor,file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData_inclPFT.csv"))

write.csv(data_guidoTD,file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs_TD","MasterData_inclPFT.csv"))
write.csv(data_guido,file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs","MasterData_inclPFT.csv"))

