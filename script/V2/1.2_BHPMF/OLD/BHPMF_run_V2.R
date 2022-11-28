is.it.on.cluster=TRUE
#  setwd("/..")
#  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("Net","Groups","BGI","work_1","2016_GapFilling")
  origin2=file.path("work_1","2016_GapFilling","_2021","data","_runs")
  Version_now="V2"
  
  # Packages
# require(BHPMF)
# read all required R scripts (?)

    
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
RepNum=1
rn=1
p=3
t_choice <- t_choices[td]
Percent = gappercents[p] 
ObsOrTD <- TDnos[TDno]
#for(RepNum in 1){
  for(td in 1){
    t_choice <- t_choices[td]
    for(p in 2:length(gappercents)){
    Percent = gappercents[p] 
    for(TDno in 1){
      ObsOrTD <- TDnos[TDno]
      
      print(paste(RepNum,ObsOrTD,t_choice,Percent))
 
        setwd("/..")
        setwd(file.path("Net","Groups","BGI"))
#        setwd(file.path("Volumes","BGI"))
        setwd(file.path(origin2,paste0("Rep_",RepNum),t_choice,paste0("p_",Percent),ObsOrTD,"R"))
        X <- as.matrix(read.csv("../data/traitInfo.csv", sep=","))[,-1]
#        colnames(X) <- paste0("V",1:ncol(X))
        rownames(X) <- NULL
        hierarchy.info <- read.table("data/taxInfo.csv", sep=",")
        rownames(hierarchy.info) <- NULL



      ###################
      # this script has to be run from the R directory of the dataset
      
      
      ########### this does not work: ################
      # set env variable OMP_NUM_THREADS=$max_threads
      Sys.setenv(OMP_NUM_THREADS = 4)
      ###############################################
      
      ########## this neither: #####################
      #args <- commandArgs(trailingOnly = T)
      #tmp_dir <- args[1]
      #message("tmp_dir: ", tmp_dir)
      #############################################
      
      
      message(Sys.time(), ": Loading data")
      source("transform_back.R")
      source("gap_filling.R")
#      X    <- as.matrix(read.table("../data/traitInfo.txt", header = F, sep = "\t"))[,-1]
      X    <- as.matrix(read.csv("../data/traitInfo.csv", header = F, sep = ","))[-1,-c(1:2)]
#      hier <-           read.table("../data/phyAllInfo.txt",sep = "\t")
      hier <-           read.csv("../data/taxInfo.csv",sep = ",")
      
      try(unlink("../data/output"))# Julia new
      if(!file.exists("../data/output")) dir.create("../data/output")
      
      if(typeof(X) != "double")       stop("X must be of type double")
      if(class(hier) != "data.frame") stop("hier must be of class data.frame")
      if(nrow(X) != nrow(hier))       stop("hier and X must have the same number of rows")
      
      #apply(hier,2, function(x) sum(is.na(x)))
      
      
      back_trans_pars <- list()
      # transform data and 
      # TODO!!! save data for backtransformation
      rm_col <- c()
      for(i in 1:ncol(X)){
        x <- X[,i]
        min_x <- min(x,na.rm = T)
        x <- x - min_x + 1  # make this optional?
        logx <- log(x)
        mlogx <- mean(logx, na.rm = T)
        slogx <- sd(logx, na.rm = T)
        x <- (logx - mlogx)/slogx
        back_trans_pars[[i]] <- list(min_x = min_x,
                                     mlogx = mlogx,
                                     slogx = slogx)
        X[,i] <- x
      }
      
      write.csv(X,file="../data/traitInfo_zlog.csv")#Julia add
      save(back_trans_pars, file = "../data/output/back_trans_pars.RData")
      
      # message("Intermediate results are saved in: ", tmp_dir)
      # 
      # message(Sys.time(), ": Starting PreprocessCv()")
      
      # PreprocessCv(X, 
      #              hier, 
      #              num.folds = 5, 
      #              tmp.dir   = tmp_dir)
      
      message(Sys.time(), ": Starting GapFilling()")
      
      pdf("../data/output/plot.pdf")
      
      GapFilling(X, 
                 hier, 
                 num.samples                 = 10, #1000, 
                 burn                        = 5, #100, 
                 gaps                        = 2,
                 num.latent.feats            = 10, 
                 tuning                      = FALSE, 
                 num.folds                   = 5, 
                 tmp.dir                     = "../data",
                 mean.gap.filled.output.path = "../data/output/mean_gap_filled.txt", 
                 std.gap.filled.output.path  = "../data/output/std_gap_filled.txt",
                 rmse.plot.test.data         = TRUE)
      
      dev.off()
      
      X_ <- read.table("../data/output/mean_gap_filled.txt", header = T, sep = "\t")
      
      X_ <- transform_back(X_, back_trans_pars)
      
      write.table(X_, file = "../data/output/mean_gap_filled_back_trans.txt", 
                  sep = "\t", row.names = F,
                  quote = F)
      
      # message(Sys.time(), ": Starting CalculateCvRmse()")
      
      # CalculateCvRmse(X, 
      #                hier, 
      #                num.folds        = 10, 
      #                num.samples      = 1000,
      #                burn             = 100, 
      #                gaps             = 2,
      #                num.latent.feats = 15,
      #                tuning           = FALSE, 
      #                tmp.dir          = tmp_dir)
      
      
      
      message(Sys.time(), ": DONE")
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
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

