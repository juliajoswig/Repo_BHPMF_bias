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
X    <- as.matrix(read.table("../data/traitInfo.txt", header = F, sep = "\t"))[,-1]
hier <-           read.table("../data/phyAllInfo.txt",sep = "\t")

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




