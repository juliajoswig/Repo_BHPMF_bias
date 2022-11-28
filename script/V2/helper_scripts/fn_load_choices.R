choices <- function(){
  
  out <- list()
  out$repnums = 1:30
  out$gappercents = c(-1,0,1,5,10,20,30,40,50,60,70,80,100)
  out$whichDataSet = 1:2
  out$ObsSpec = 1:4
  out$obsSpec = 1:4
  out$preparation = "no"

  out$tsubs = c("data","data_2")
  out$TD_choices =  c("Obs_obs_TD","Obs_obs")
  out$ObsOrTDs <- c("Obs_obs_TD","Obs_obs")
  out$trait_rainfor =  c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
  out$trait_guido =  c("SLA","PlantHeight","SeedMass","LDMC","LeafArea")#c("SLA","PlantHeight","SSD","LeafN","LeafP","LeafNArea")
  out$colz2<- c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#8c2d04")
  out$colz1<- c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")
  
  # Load functional data  
  list.files(file.path(origin,"_2021","data","helper_data","functional_information"))
  GF <- read.table(file.path(origin,"_2021","data","helper_data","functional_information","GF_Obs.csv"),sep=",",dec=".")
  PFT <- read.table(file.path(origin,"_2021","data","helper_data","functional_information","PFT_Obs.csv"),sep=",",dec=".")

  #------------------------------
  #------------------------------

  out$data <- read.csv(file.path(origin,"_2021","data","helper_data","ObservationIDs","data","traitInfo.csv"))[,-1]
  out$data_2 <- read.csv(file.path(origin,"_2021","data","helper_data","ObservationIDs","data_2","traitInfo.csv"))[,-1]
  
  out$new.mean.fun <- function(x){return(mean(x,na.rm = TRUE))}
  out$new.sd.fun <- function(x){return(sd(x,na.rm = TRUE))}
  out$new.count.fun <- function(x){return(sum(!is.na(x)))}

  out$add_col_to_res <- function(new.col.names,input){
    for(icn in 1:length(new.col.names)){
      new_col=rep(NA,nrow(input))
      if(sum(colnames(input)%in%new.col.names[icn])==0){
        input$new_col <-  new_col
        colnames(input)[ncol(input)] = new.col.names[icn]
        print(paste(new.col.names[icn],"added"))
      }
    }
    return(input)
  }
  
  return(out)
}