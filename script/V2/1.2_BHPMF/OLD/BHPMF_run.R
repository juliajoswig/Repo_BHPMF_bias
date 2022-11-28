require(BHPMF)

is.it.on.cluster=TRUE
#  setwd("/..")
#  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("Net","Groups","BGI","work_1","2016_GapFilling")
  origin2=file.path("work_1","2016_GapFilling","_2021","data","_runs")
  Version_now="V2"
  
  # Packages
#  require(doParallel)
#  require(foreach) # parallelisierte Vorschleife
#  require(plyr) # evtl.

  
  
  #  TuneBhpmf(X, hierarchy.info, num.folds=5, num.samples=sampleNUM, burn=burnNUM, tmp.dir="../data/fold1/Tunning")
  
  print("packages loaded")

  # necessary infos
  sampleNUM= 1000
  gapNUM = 20
  burnNUM= 200
  gappercents=c(1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3


# start parallel
#  cl <- makeCluster(2)
#  registerDoParallel(cl)


td=1
TDno=1
RepNum=3
rn=1
p=3
t_choice <- t_choices[td]
Percent = gappercents[p] 
ObsOrTD <- TDnos[TDno]
template=TRUE

#for(RepNum in 1){
  for(td in 1){
    t_choice <- t_choices[td]
    for(p in 2:length(gappercents)){
    Percent = gappercents[p] 
    for(TDno in 1){
      ObsOrTD <- TDnos[TDno]
      
      print(paste(RepNum,ObsOrTD,t_choice,Percent))
 
      if(template){
        setwd("/..")
#        setwd(file.path("Volumes","BGI"))
        setwd(file.path("Net","Groups","BGI"))
        setwd(file.path(origin2,"template"))
        XT <- as.matrix(read.table("data/trait_info.csv", sep=","))
        hierarchy.infoT <- read.table("data/hierarchy_info.csv", sep=",")
      }
      
        setwd("/..")
        setwd(file.path("Net","Groups","BGI"))
#        setwd(file.path("Volumes","BGI"))
        setwd(file.path(origin2,paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))
        X <- as.matrix(read.table("data/traitInfo_zlog.txt", sep=","))[,-1]
#        colnames(X) <- paste0("V",1:ncol(X))
        rownames(X) <- NULL
        hierarchy.info <- read.table("data/taxInfo.csv", sep=",")
        rownames(hierarchy.info) <- NULL
      
        unlink("data/tmp", recursive=TRUE)
        dir.create("data/tmp")
        
      XT[1:1140,1:5]<-X
      XT<-XT[1:1140,1:5]
      mode(hierarchy.info[,1:5])
      hierarchy.infoT[1:1140,1] <- hierarchy.info[,1]
      str(hierarchy.infoT[,2])
      hierarchy.infoT[1:1140,2] <- as.factor(hierarchy.info[,2])
      hierarchy.infoT[1:1140,3] <- hierarchy.info[,3]
      hierarchy.infoT[1:1140,4] <- hierarchy.info[,4]
      hierarchy.infoT[1:1140,5] <- hierarchy.info[,5]
      
      GapFilling(XT, hierarchy.infoT, 
                 num.samples=sampleNUM, 
                 burn=burnNUM,
                 gaps=gapNUM, 
                 num.folds=5,
                 num.latent.feats=15,
                 tmp.dir="data/tmp", 
                 mean.gap.filled.output.path="data/mean.csv",
                 std.gap.filled.output.path="data/std.csv")
      
      GapFilling(X, hierarchy.info, 
                 num.samples=sampleNUM, 
                 burn=burnNUM,
                 gaps=gapNUM, 
                 num.folds=5,
                 num.latent.feats=15,
                 tmp.dir="data/tmp", 
                 mean.gap.filled.output.path="data/mean.csv",
                 std.gap.filled.output.path="data/std.csv")
if(1==1){      
      rmse_out <- CalculateCvRmse(X, hierarchy.info, 
                                  num.samples=sampleNUM, burn=burnNUM, 
                                  gaps=gapNUM, num.folds=5, num.latent.feats=15,
                                  tmp.dir="data/tmp")
      save(rmse_out,file=file.path("data/rmse_out.RData"))
}
      
    print("done")
    
    }
  }
  }

#}

