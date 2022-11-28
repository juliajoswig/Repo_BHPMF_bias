GapFilling <- function(X, hierarchy.info, prediction.level,
        used.num.hierarchy.levels, num.samples=1000, burn=100, gaps=2,
        num.latent.feats=10, tuning=FALSE, num.folds=5, tmp.dir,
		mean.gap.filled.output.path, std.gap.filled.output.path,
		rmse.plot.test.data=TRUE) {
# This function calculate the average RMSE in cross validation.
  #
  # Args:
  #   X: A matrix containing the missing values, rows are the observations,
  #      and columns are the features. Missing values are filled by NA.
  #   hierarchy.info: A matrix containing the hierarchical info
  #      in the format of taxa-table.
  #   prediction.level: The level at which gaps are filled.
  #      e.g., =4, filling gaps at leaves with a hierarchy informaiton of 3 levels.
  #		 The default value is at the observation level.
  #   used.num.hiearchy.levels: Number of hierarchy levels that is used for gap filling.
  #      e.g., =2 means using only the first and second level of hierarchy.
  #      It should between 1 and total number of hierarchy.
  #		 The default value is total number of hierarchy level.
  #   num.samples: Total number of generated samples at each fold using gibbs sampling.
  #      It is not the effective number of samples.
  #		 The default value is 1000.
  #   burn: Number of initial sampled parameters discarded. The default value is 100.
  #   gaps: Gap between sampled parameters kept. The default value is 2.
  #	  num.latent.feats: Size of latent vectors in BHPMF. Set it if tuning is False.
  #	  	 The default value is 10.
  #	  tuning: If set true, first tune BHPMF to choose the best value of num.latent.feats.
  #	  	 The default value is False.
  #   num.folds: Number of cross validation (CV) folds used in tuning. The default value is 10.
  #	  tmp.dir: A temporary directory used to save preprocessing files.
  #	  	 If not provided, a tmp directory in /R/tmp will be created.
  #		 If provided, each time calling a function from this package, will use the saved
  #		 preprocessing files saved in this directory, it helps to avoid running preprocessing
  #		 functions for the same input.
  #		 WARNING: This directory should be empty for the first time call on a new dataset.
  #	  mean.gap.filled.output.path: A file path for saving the mean value of gap filled data.
  #	  	 In the same format as the input X.
  #	  	 It contain the prediction for both missing and observed values.
  #	  std.gap.filled.output.path: A file path for saving the std value of gap filled data.
  #	  	 In the same format as the input X.
  #	  rmse.plot.test.data: If TRUE, provide the rmse vs std plot for test data.
  #
  # Returns:
  #   The average RMSE in cross validation

    source("file_exists.R")
    source("preprocess_cv.R")
	source("plot_rmse_std.R")
	num.hierarchy.levels <- ncol(hierarchy.info)
	num.cols <- ncol(X)

    # Check missing arguments
    if (missing(X)) {
	    stop("Missing X!")
    }
    if (missing(hierarchy.info)) {
	    stop("Missing hierarchy.info!")
    }
    if (missing(mean.gap.filled.output.path)) {
       stop("A file path for mean.gap.filled.output.path should be provided.")
    }
	if (missing(std.gap.filled.output.path)) {
	   stop("A file path for std.gap.filled.output.path should be provided.")
	}
    if (missing(used.num.hierarchy.levels)) {
        used.num.hierarchy.levels <- num.hierarchy.levels-1
    }
    if (missing(prediction.level)) {
       prediction.level <- num.hierarchy.levels
    }
	preprocess.flag <- FALSE
    if (missing(tmp.dir)) {
       tmp.dir <- paste(tempdir(), "/BHPMFAuthorsTmp", sep = "")
       if (file.exists(tmp.dir)) {
          unlink(tmp.dir, recursive = TRUE, force = TRUE)
       }
       dir.create(tmp.dir)
       preprocess.flag <- TRUE
    } else if (!file.exists(tmp.dir)) {
       stop("tmp.dir: ", tmp.dir, "  does not exist")
    }

    if (!preprocess.flag) {   # tmp directory is provided by user
       if (!CheckPreprocessFilesExist(tmp.dir, 1, used.num.hierarchy.levels, prediction.level)) { # return TRUE if files exists
          preprocess.flag <- TRUE
       }
    }

	cat("preprocess.flag: ",preprocess.flag, "\n")
	num.folds.cv1 <- 2
	if (tuning) {
	   num.folds.cv2 <- num.folds
	} else {
	   num.folds.cv2 <- 0
	}
    if (preprocess.flag) {
       PreprocessCv(X, hierarchy.info, num.folds.cv1, tmp.dir)
       tmp.tune.dir = paste(tmp.dir, "/fold1/Tunning", sep="")
	   cat(tmp.tune.dir)
       dir.create(tmp.tune.dir)
    }

	load(paste(tmp.dir, "/processed_hierarchy_info.Rda", sep=""))

	#Tune the num.latent.feats parameter
	#By choosing the best parameter that minimize the
	#CV RMSE over training data in the first fold
	#Assumption: the best parameter is the same for other folds
    if (tuning) {
       cat("tuning: ")
       source("tune_BHPMF.R")
       # find the X matrix for fold 1
       file.name <- paste(tmp.dir, '/fold1/Ytrain', prediction.level, ".txt", sep = "")
       Y.tune <- as.matrix(read.table(file.name, sep="\t", header=F))
       # X.tune <- sparseMatrix (Y.tune[,1], Y.tune[,2], x=Y.tune[,3])
       X.tune <- matrix(data=NA, nrow=nrow(X), ncol=ncol(X))
       X.tune[cbind(Y.tune[, 1], Y.tune[, 2])] <- Y.tune[, 3]
       rm(Y.tune)

       tmp.tune.dir = paste(tmp.dir, "/fold1/Tunning", sep="")

       out <- TuneBhpmf(X.tune, hierarchy.info, prediction.level, used.num.hierarchy.levels,
                    num.folds, num.samples, burn, gaps, tmp.tune.dir)
       rm(X.tune)
       num.latent.feats <- out$BestNumLatentFeats
    }

	source("utillity.R")
	dyn.load("../src/HPMF.so")

	save.file.flag <- 0
	out.whole.flag <- 1
	fold <- 1
	opt <- 2
	tmp.env <- new.env()
    args <- list("NumSamples" = as.integer(num.samples),
                 "InputDir" = tmp.dir,
                 "DatasetId" = as.integer(fold),
                 "Gaps" = as.integer(gaps),
                 "Burn" = as.integer(burn),
                 "SaveFileFlag" = as.integer(save.file.flag),
                 "OutWholeFlag" = as.integer(out.whole.flag),
                 "NumTraits" = as.integer(num.cols),
                 "NumFeats" = as.integer(num.latent.feats),
                 "NumHierarchyLevel" = as.integer(num.hierarchy.levels-1),
                 "PredictLevel" = as.integer(prediction.level),
                 "UsedNumHierarchyLevel" = as.integer(used.num.hierarchy.levels),
                 "Opt" = as.integer(opt),
                 "Env" = tmp.env
                 )

	out <- .Call("DemoHPMF", args, mean.gap.filled.output.path, std.gap.filled.output.path, num.nodes.per.level)
	
	gap.filled.dat = as.matrix(read.table(mean.gap.filled.output.path, sep="\t"))
	gap.filled.dat = gap.filled.dat[,1:ncol(X)]
	colnames(gap.filled.dat) = colnames(X)
	write.table(gap.filled.dat, file = mean.gap.filled.output.path, sep="\t", col.names = T, row.names = F)

	std.filled.dat = as.matrix(read.table(std.gap.filled.output.path, sep="\t"))
	std.filled.dat = std.filled.dat[,1:ncol(X)]
	colnames(std.filled.dat) = colnames(X)
	write.table(std.filled.dat, file = std.gap.filled.output.path, sep="\t", col.names = T, row.names = F)

	if (rmse.plot.test.data) {

	    file.name <- paste(tmp.dir, "/fold1/Ytest", prediction.level, ".txt", sep = "")
        test.data <- as.matrix(read.table(file.name))
		row.idx <- test.data[,1]
		col.idx <- test.data[,2]
		nrows <- nrow(X)
	    test.idx <- (col.idx - 1) * nrows + row.idx

		test.mean <- gap.filled.dat[test.idx]
		test.res <- test.mean - test.data[, 3]
        test.rmse <- sqrt(mean(test.res^2));
		cat("RMSE for the test data: ", test.rmse)

		test.std <- std.filled.dat[test.idx]
	   	PlotRmseVsStd(test.res, test.std)
	}	
}