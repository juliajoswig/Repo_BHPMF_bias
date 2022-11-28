is.it.on.cluster=FALSE
is.it.on.cluster=TRUE
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
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
trait_rainfor =  c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
trait_guido =  c("SLA","PlantHeight","SeedMass","LDMC","LeafArea")#c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")

MasterData <- as.data.frame(read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData_inclPFT.csv")))
  colnames(MasterData[,25:ncol(MasterData)])
  summary(MasterData[,25:ncol(MasterData)])
  hist(MasterData[,colnames(MasterData)%in%"Leafdelta15N"])
  MasterData_log <- cbind(log(MasterData[,c(25:39)]),
                          MasterData[,40],
                          log(MasterData[,c(41:ncol(MasterData))]))
  colnames(MasterData_log)[16] <- "Leafdelta15N"
  colSums(is.na(MasterData_log)==is.na(MasterData[,25:ncol(MasterData)]))==nrow(MasterData)
  MasterData_mean <- apply(MasterData_log,2,mean,na.rm=TRUE)
  MasterData_sd <- apply(MasterData_log,2,sd,na.rm=TRUE)
  MeanAndSd <- rbind(MasterData_mean,MasterData_sd)
  save(MeanAndSd,file=file.path(origin,"_2021","data","helper_data","scaling","factors.Rdata"))
  

  

td=1
TDno=1
RepNum=1
rn=1
p=8
t_choice <- t_choices[td]
Percent = gappercents[p]
gappercents=0
ObsOrTD <- TDnos[TDno]

#list.files(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum)))
           

for(RepNum in 1:3){
  
  print(RepNum)
  for(td in 1:2){
    t_choice <- t_choices[td]
    p=1
    for(p in 1:length(gappercents)){
      Percent = gappercents[p]
     
        for(TDno in 1:2){
          ObsOrTD <- TDnos[TDno]
          traitInfo <- read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                             paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
          
          traitInfo_c <- traitInfo[,-grep(colnames(traitInfo),pattern = "ObservationID")]
          traitInfo_c <- traitInfo_c[,-grep(colnames(traitInfo_c),pattern = "X")]
          traitInfo_log <- log(traitInfo_c)

          zlog_trans = rbind(apply(traitInfo_log,2,mean,na.rm=TRUE),
                       apply(traitInfo_log,2,sd,na.rm=TRUE))
          rownames(zlog_trans) <- c("mean","sd")
          meanM <- matrix(NA,ncol=ncol(traitInfo_log),nrow=nrow(traitInfo_log))
          sdM <- meanM
          for(i in 1:nrow(traitInfo_log)){
            meanM[i,] <- zlog_trans[1,]
            sdM[i,] <- zlog_trans[2,]
          }
          zlog_trans <- list()
          zlog_trans$meanM <- meanM
          zlog_trans$sdM <- sdM
          traitInfo_zlog1 <- (traitInfo_log - meanM) /sdM
          traitInfo_zlog <- cbind(traitInfo[,colnames(traitInfo)== "ObservationID"],traitInfo_zlog1)
          colnames(traitInfo_zlog)[1]="ObservationID"
          
          print("*********")
          print(paste0(sum(rowSums(!is.na(traitInfo_zlog[,-1]))==0)," missing values per row"))
          
          unlink(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                           paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.txt"),recursive = TRUE)
          unlink(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                           paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"),recursive = TRUE)
          write.table(traitInfo_zlog,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                                   paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"), sep="," ,dec = ".")
#          read.csv(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","traitInfo_zlog.csv"))[,-1]
          write.table(traitInfo_log,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                                    paste0("p_",Percent),ObsOrTD,"data","traitInfo_log.csv"), sep="," ,dec = ".")
          save(zlog_trans,file = file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                           paste0("p_",Percent),ObsOrTD,"data","zlog_transform.RData"))
          
          print(file.path(paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD)) 
        }
      }
    }
  }



























