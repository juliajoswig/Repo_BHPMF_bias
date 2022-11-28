

RunTheHPMF <- function(origin2,Percent,t_choice,ObsOrTD,RepNum){
  
#  setwd(file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent,"_",t_choice),TD_choice,"R"))
  setwd(file.path(origin2,paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))
  require(BHPMF)
  
#  X <- as.matrix(read.table("../data/traitInfo.csv", sep=","))
#  hierarchy.info <- read.table("../data/taxInfo.csv", sep=",")
  X <- as.matrix(read.csv(file = "data/traitInfo.csv"))
  X <- X[,-1]
  hierarchy.info <- as.matrix(read.csv(file = "data/taxInfo.csv"))
  hierarchy.info <- hierarchy.info[,-1]

  print(paste0(paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD))
  
  GapFilling(X, hierarchy.info, 
             num.samples=sampleNUM, 
             burn=burnNUM,
             gaps=gapNUM, 
             num.folds=5,
             num.latent.feats=15,
             tmp.dir="../data/tmp", 
             mean.gap.filled.output.path="../data/mean.csv",
             std.gap.filled.output.path="../data/std.csv")
  
  rmse_out <- CalculateCvRmse(X, hierarchy.info, num.samples=sampleNUM, burn=burnNUM, gaps=gapNUM, num.folds=5, num.latent.feats=15, tmp.dir="../data")
  save(rmse_out,file=file.path(origin,"_2021","data","_runs",paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"data",
                               "rmse_out.RData"))
}


#  TuneBhpmf(X, hierarchy.info, num.folds=5, num.samples=sampleNUM, burn=burnNUM, tmp.dir="../data/fold1/Tunning")

