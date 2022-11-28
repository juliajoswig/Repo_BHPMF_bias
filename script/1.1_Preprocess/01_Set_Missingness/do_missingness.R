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

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
trait_rainfor =  c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
trait_guido =  c("SLA","PlantHeight","SeedMass","LDMC","LeafArea")#c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")

data_rainforTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs_TD","MasterData_inclPFT.csv"))
colnames(data_rainforTD)
data_rainfor <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData_inclPFT.csv"))
colnames(data_rainfor)

data_guidoTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs_TD","MasterData_inclPFT.csv"))
colnames(data_guidoTD)
data_guido <- read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs","MasterData_inclPFT.csv"))
colnames(data_guido)



data_rainforTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs_TD","MasterData_inclPFT.csv"))
data_rainforTD <- data_rainforTD[,-grep(colnames(data_rainforTD),pattern = "X")]
data_rainforTD <- data_rainforTD[,-grep(colnames(data_rainforTD),pattern = ".1")]
data_rainfor <- read.csv(file=file.path(origin,"_2021","data","MasterData","data","Obs_obs","MasterData_inclPFT.csv"))
data_rainfor <- data_rainfor[,-grep(colnames(data_rainfor),pattern = "X")]
data_rainfor <- data_rainfor[,-grep(colnames(data_rainfor),pattern = ".1")]

data_guidoTD <- read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs_TD","MasterData_inclPFT.csv"))
data_guidoTD <- data_guidoTD[,-grep(colnames(data_guidoTD),pattern = ".1")]
data_guidoTD <- data_guidoTD[,-grep(colnames(data_guidoTD),pattern = "X")]
data_guido <-   read.csv(file=file.path(origin,"_2021","data","MasterData","data_2","Obs_obs","MasterData_inclPFT.csv"))
data_guido <- data_guido[,-grep(colnames(data_guido),pattern = ".1")]
data_guido <- data_guido[,-grep(colnames(data_guido),pattern = "X")]

taxcols=colnames(data_guido)%in%c("ObservationID","AccSpeciesName","Genus","Family","PhylogeneticGroup")
funcols=colnames(data_guido)%in%c("ObservationID","PlantGrowthForm","PFT")
dim(data_guido)
dim(data_guidoTD)

colnames(data_rainforTD)
td=1
TDno=1
RepNum=1
rn=1
p=5
t_choice <- t_choices[td]
Percent = gappercents[p]
ObsOrTD <- TDnos[TDno]

for(RepNum in repnums:1){
  
  print(RepNum)
  for(td in 1:2){
    t_choice <- t_choices[td]
    print(t_choice)
    if(t_choice=="data"){
      fun_envelope=data_rainfor[,funcols]
      fun_now=data_rainforTD[,funcols]
      tax_envelope=data_rainfor[,taxcols]
      tax_now=data_rainforTD[,taxcols]
      data_now=data_rainforTD[,c(1,23:ncol(data_rainforTD))]
      data_envelope=data_rainfor[,c(1,23:ncol(data_rainfor))]}
    if(t_choice=="data_2"){
      fun_envelope=data_guido[,funcols]
      fun_now=data_guidoTD[,funcols]
      tax_envelope=data_guido[,taxcols]
      tax_now=data_guidoTD[,taxcols]
      data_now = data_guidoTD[,c(1,23:ncol(data_guidoTD))]
      data_envelope = data_guido[,c(1,23:ncol(data_guido))]
      }

    # add according to % 
    # select entries
    for(p in 1:length(gappercents)){
      
      #Missingness select base entries
      # retain 1 per row
      # TRUE = set NA i.e. FALSE = value
      ix_NA <- matrix(FALSE,ncol=(ncol(data_now)-1), nrow=nrow(data_now))
      for(i in 1:nrow(data_now)){ix_NA[i,sample(1:(ncol(data_now)-1),1)] <- TRUE}
      keep_base <- sum(ix_NA)/(sum(ix_NA)+sum(!ix_NA))
      entries_base <- sum(!ix_NA)/(sum(ix_NA)+sum(!ix_NA))
      ix_NA2=ix_NA
      print(sum(rowSums(!ix_NA)==0))
      
      Percent = gappercents[p]
      total_vals=(ncol(data_now)-1)*nrow(data_now)
      nb_total_entries = round(total_vals*((100-Percent)/100),digits = 0)
      nb_entries_add <- nb_total_entries-sum(ix_NA)
      if(nb_entries_add<=0){nb_entries_add=1}
      length(as.vector(ix_NA2[!ix_NA]))
      
      ix_NA2[!ix_NA][sample(1:(total_vals-sum(ix_NA)),nb_entries_add)] <- TRUE
      
      print("-------------------------")  
      print(round(sum(ix_NA2)/(sum(ix_NA2)+sum(!ix_NA2)),digits = 2))
      print((100-Percent)/100)
      print(sum(rowSums(!ix_NA2)==0))
      
      if(round(sum(ix_NA2)/(sum(ix_NA2)+sum(!ix_NA2)),digits = 2)==(100-Percent)/100){
        
        TDno=2
        for(TDno in 1:2){
          ObsOrTD <- TDnos[TDno]
          
          ix_NA3 <- cbind(rep(TRUE,nrow=nrow(data_now)),ix_NA2)
          data_now[!ix_NA3] <- NA
          if(ObsOrTD=="Obs_obs"){
            mtch <- match(data_now$ObservationID,data_envelope$ObservationID)
            data_envelope$ObservationID[mtch]==data_now$ObservationID
            if(t_choice=="data"){trait_subset_names=trait_rainfor}
            if(t_choice=="data_2"){trait_subset_names=trait_guido}
            mtch_traits=match(trait_subset_names,colnames(data_envelope))
            data_envelope[mtch,mtch_traits] <- data_now[,-grep(colnames(data_now),pattern = "ObservationID")]
            colnames(data_envelope[mtch,match(trait_subset_names,colnames(data_envelope))] )
            data_save=data_envelope
          }
          if(ObsOrTD=="Obs_obs_TD"){data_save=data_now}
          dim(data_save)
          
          #Save
          if(!dir.exists(file.path(origin,"_2021","data","_runs"))){
            dir.create(file.path(origin,"_2021","data","_runs"))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum)))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum)))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent)))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent)))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data"))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","tmp"))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","tmp"))}
          if(!dir.exists(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"R"))){
            dir.create(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"R"))}
          
          write.csv(data_save,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                             paste0("p_",Percent),ObsOrTD,"data","traitInfo.csv"))
          if(ObsOrTD=="Obs_obs"){
            write.table(tax_envelope,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,
                                                      paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"), sep="," ,dec = ".", col.names = F, row.names = F)
            write.csv(fun_envelope,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","funInfo.csv"))
          }
          if(ObsOrTD=="Obs_obs_TD"){
            write.table(tax_now,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","taxInfo.csv"),
                        sep="," ,dec = ".", col.names = F, row.names = F)
            write.csv(fun_now,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data","funInfo.csv"))
          }
          print(file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD)) 
          #rm("data_save")
        }
      }
    }
    rm("data_now")
  }
}



























