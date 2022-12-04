require(BHPMF)

is.it.on.cluster=FALSE
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
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
# gappercents=c(5,10)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3


# start parallel
#  cl <- makeCluster(2)
#  registerDoParallel(cl)


td=1
TDno=2
RepNum=1
rn=1
p=3
t_choice <- t_choices[td]
Percent = gappercents[p] 
ObsOrTD <- TDnos[TDno]
template=FALSE
gappercents=c(1)
Percent=50

for(RepNum in 2){#sample(1:3,1)
  for(td in 1){
    t_choice <- t_choices[td]
    for(p in 2){ #sample(1:length(gappercents),2) 
      Percent = gappercents[p] 
      for(TDno in 2){
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
        setwd(file.path(origin2,paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))
        
        X <- as.matrix(read.csv("data/traitInfo_zlog.csv")[,-1])
        #        X <- as.matrix(read.table("data/traitInfo.csv", sep=","))[-1,-1]
        rownames(X) <- NULL
        hierarchy.info <- read.table("data/taxInfo.csv", sep=",")
        rownames(hierarchy.info) <- NULL
        print(paste0(sum(rowSums(!is.na(X))==0)," missing row values"))
        if(sum(rowSums(!is.na(X))==0)==0){
          
          
          #if(!file.exists("data/mean.csv")){# do run only if it does not exist yet...
          try(unlink("data/tmp", recursive=TRUE))
          dir.create("data/tmp")
          GapFilling(X, hierarchy.info, 
                     num.samples=sampleNUM, 
                     burn=burnNUM,
                     gaps=gapNUM, 
                     num.folds=5,
                     num.latent.feats=15,
                     tmp.dir="data/tmp", 
                     mean.gap.filled.output.path="data/mean.csv",
                     std.gap.filled.output.path="data/std.csv")
          
          if(TDno==1){      
            try(unlink("data/tmpRMSE", recursive=TRUE))
            dir.create("data/tmpRMSE")
            rmse_out <- CalculateCvRmse(X, hierarchy.info, 
                                        num.samples=sampleNUM, burn=burnNUM, 
                                        gaps=gapNUM, num.folds=5, num.latent.feats=15,
                                        tmp.dir="data/tmpRMSE")
            save(rmse_out,file=file.path("data/rmse_out.RData"))
          }
          
          print("done")
          #}
        }
      }
    }
  }
  
}
# 2 running data TD 
print("done :)")
