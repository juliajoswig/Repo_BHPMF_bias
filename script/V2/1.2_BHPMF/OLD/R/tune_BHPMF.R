TuneBhpmf <- function(X, hierarchy.info, prediction.level, used.num.hierarchy.levels,
		  	 		num.folds=10, num.samples=1000, burn=100, gaps=2, tmp.dir) {

# Tune the BHPMF parameter, number of latent factors.
  #
  # Args:
  #   X: A matrix containing the missing values, rows are the observations,
  #      and columns are the features. Missing values are filled by NA.
  #   hierarchy.info: A matrix containing the hierarchical info
  #      in the format of taxa-table.
  #   num.folds: Number of cross validation (CV) folds. The default value is 10.
  #   prediction.level: The level at which gaps are filled.
  #      e.g., =4, filling gaps at leaves with a hierarchy informaiton of 3 levels.
  #		  The default value is at the observation level.
  #   used.num.hierarchy.levels: Number of hierarchy levels that is used for gap filling.
  #      e.g., =2 means using only the first and second level of hierarchy.
  #		  It should between 1 and total number of hierarchy.
  #      The default value is total number of hierarchy level.
  #   num.samples: Total number of generated samples at each fold using gibbs sampling.
  #      It is not the effective number of samples.
  #		  The default value is 1000.
  #   burn: Number of initial sampled parameters discarded. The	default	value is 100.
  #   gaps: Gap between sampled parameters kept. The default value is 2.
  #   tmp.dir: A temporary directory used to save preprocessing	files.
  #      If not provided, a tmp directory in /R/tmp will be created.
  #		 If provided, each time calling	a function from	this package, will use the saved
  #		 preprocessing files saved in this directory, it helps to avoid running preprocessing
  #		 functions for the same	input.
  #		 WARNING: This directory should	be empty for the first time call on a new dataset.
  #
  # Returns:
  #   The tuned number of latent factors paramter for BHMPF.


    source("file_exists.R")
    source("preprocess_cv.R")
    source("utillity.R")
    dyn.load("../src/HPMF.so")

    num.hierarchy.levels <- ncol(hierarchy.info)
    num.cols <- ncol(X)

    # Check missing arguments
    if (missing(X)) {
        stop("Missing X!")
    }
    if (missing(hierarchy.info)) {
        stop("Missing hierarchy.info!")
    }
    if (missing(used.num.hierarchy.levels)) {
        used.num.hierarchy.levels <- num.hierarchy.levels - 1
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
       if (!CheckPreprocessFilesExist(tmp.dir, num.folds, used.num.hierarchy.levels, prediction.level)) { # return TRUE if files exists
          preprocess.flag <- TRUE
       }
    }

	if (preprocess.flag) {
	   cat("preprocessing: \n")
	   PreprocessCv(X, hierarchy.info, num.folds, tmp.dir)
	}

	load(paste(tmp.dir, "/processed_hierarchy_info.Rda", sep=""))
	tmp.env <- new.env()
	min_rmse <- Inf
	save_file_flag <- 0
	opt <- 2
	out.whole.flag <- 0
	num_latent_feats_set <- seq(5, num.cols+6, by=5);
	rmse_vec <- matrix(0, length(num_latent_feats_set), num.folds)

	for (ind in 1 : length(num_latent_feats_set)) {
		for (fold in 1 : num.folds) {
			num_latent_feats = num_latent_feats_set[ind];
			args <- list(
				  "NumSamples" = as.integer(num.samples),
				  "InputDir" = tmp.dir,
				  "DatasetId" = as.integer(fold),
				  "Gaps" = as.integer(gaps),
				  "Burn" = as.integer(burn),
                  "SaveFileFlag" = as.integer(save_file_flag),
                  "OutWholeFlag" = as.integer(out.whole.flag),
				  "NumTraits" = as.integer(num.cols),
				  "NumFeats" = as.integer(num_latent_feats),
				  "NumHierarchyLevel" = as.integer(num.hierarchy.levels-1),
				  "PredictLevel" = as.integer(prediction.level),
				  "UsedNumHierarchyLevel" = as.integer(used.num.hierarchy.levels),
				  "Opt" = as.integer(opt),
				  "Env" = tmp.env
				  )

			 out <- .Call("DemoHPMF", args, "", "", num.nodes.per.level);
			 cat(out$RMSE);
			 cat("\n")
			 rmse_vec[ind, fold] = out$RMSE;	 
		}

		avg_rmse = mean(rmse_vec[ind,]);
		if (avg_rmse < min_rmse) {
		   min_rmse = avg_rmse;
		   best_num_latent_feats = num_latent_feats_set[ind]
		}
	}
	
	result <- list("MinRMSE" = as.numeric(min_rmse),
		   	  	   "BestNumLatentFeats" = as.integer(best_num_latent_feats)
				   )
	return(result);
}