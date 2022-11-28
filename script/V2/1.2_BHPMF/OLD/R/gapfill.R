###################
# this script has to be run from the R directory of the dataset


message(Sys.time(), ": Loading data")
#source("transform_back.R")
source("gap_filling.R")
X    <- as.matrix(read.table("../data/traitInfo.txt", header = F, sep = "\t"))[,-1]
hier <-           read.table("../data/phyAllInfo.txt",sep = "\t")

if(!file.exists("../data/output")) dir.create("../data/output")

if(typeof(X) != "double")       stop("X must be of type double")
if(class(hier) != "data.frame") stop("hier must be of class data.frame")
if(nrow(X) != nrow(hier))       stop("hier and X must have the same number of rows")

#apply(hier,2, function(x) sum(is.na(x)))


#back_trans_pars <- list()
# transform data and
# TODO!!! save data for backtransformation
#rm_col <- c()

### some possible change:
# this part should not be done here, instead it should be done directly after the thinning
# thinning should write a file with unprocessed in its name
# preprocessing should make the transformation and write the actual file
for(i in 1:ncol(X)){
  x <- X[,i]
  minx <- min(x, na.rm = T)
  x_rng <- quantile(x,probs = c(.1,.9), na.rm = T)
  x[ x < x_rng[1] | x > x_rng[2] ] <- NA # make it robust, ie remove tails

  # if there are no negative values take logarithm first and then z-transform
  # this should be the case for all variables, except Leaf_nitrogen_phosphorus_N_P_ratio
  if(minx > 0){
    logx <- log(x)
    mlogx <- mean(logx, na.rm = T)
    slogx <- sd(logx, na.rm = T)
    #back_trans_funs[[i]] <- (function(s,m) function(x) exp(x*s + m))(slogx,mlogx)
    X[,i] <- (log(X[,i]) - mlogx) / slogx
    # if there are negative values do only z-transformation
    # this should be the case only for Leaf_nitrogen_phosphorus_N_P_ratio
  } else {
    m <- mean(x, na.rm = T)
    s <- sd(x, na.rm = T)
    #back_trans_funs[[i]] <- (function(s,j) function(x) x*s + m)(s,m)
    X[,i] <- (X[,i] - m) / s
  }
}

#save(back_trans_funs, file = "../data/output/back_trans_pars.RData")

tmp_dir = "../data"
#message("Intermediate results are saved in: ", tmp_dir)



message(Sys.time(), ": Starting GapFilling()")

pdf("../data/output/plot.pdf")

GapFilling(X,
           hier,
           num.samples                 = 3000,
           burn                        = 500,
           gaps                        = 20,
           num.latent.feats            = 10,
           tuning                      = FALSE,
           num.folds                   = 5,
           tmp.dir                     = "../data",
           mean.gap.filled.output.path = "../data/output/mean_gap_filled.txt",
           std.gap.filled.output.path  = "../data/output/std_gap_filled.txt",
           rmse.plot.test.data         = TRUE)

dev.off()

#X_ <- read.table("../data/output/mean_gap_filled.txt", header = T, sep = "\t")

#X_ <- transform_back(X_, back_trans_pars)

#write.table(X_, file = "../data/output/mean_gap_filled_back_trans.txt",
#            sep = "\t", row.names = F,
#            quote = F)

message(Sys.time(), ": DONE")
