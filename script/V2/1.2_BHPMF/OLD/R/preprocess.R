###################
# this script has to be run from the R directory of the dataset


message(Sys.time(), ": Loading data")
source("preprocess_cv.R")
X    <- as.matrix(read.table("../data/traitInfo.txt", header = F, sep = "\t"))[,-1]
hier <-           read.table("../data/phyAllInfo.txt",sep = "\t")

if(!file.exists("../data/output")) dir.create("../data/output")

if  (typeof(X) != "double")     stop("X must be of type double")
if(class(hier) != "data.frame") stop("hier must be of class data.frame")
if    (nrow(X) != nrow(hier))   stop("hier and X must have the same number of rows")


back_trans_funs <- list()
# transform data:
# make a robust logz transforation
rm_col <- c()
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
    back_trans_funs[[i]] <- (function(s,m) function(x) exp(x*s + m))(slogx,mlogx)
    X[,i] <- (log(X[,i]) - mlogx) / slogx
    # if there are negative values do only z-transformation
    # this should be the case only for Leaf_nitrogen_phosphorus_N_P_ratio
  } else {
    m <- mean(x, na.rm = T)
    s <- sd(x, na.rm = T)
    back_trans_funs[[i]] <- (function(s,j) function(x) x*s + m)(s,m)
    X[,i] <- (X[,i] - m) / s
  }
}

save(back_trans_funs, file = "../data/output/back_trans_funs.RData")
write.csv(X, file = "../data/output/TraitInfo.csv")

tmp_dir = "../data"
message("Intermediate results are saved in: ", tmp_dir)

message(Sys.time(), ": Starting PreprocessCv()")

PreprocessCv(X,
             hier,
             num.folds = 5,
             tmp.dir   = tmp_dir)


message(Sys.time(), ": DONE")
