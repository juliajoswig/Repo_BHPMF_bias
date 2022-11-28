readFileNumeric <- function(args)
{
	X = as.matrix(read.table(args[[1]], sep = args[[2]], colClasses = "numeric")); 
	retList <- list("matrix" = t(X),
                 "ncols" = as.integer(ncol(X)),
                 "nrows" = as.integer(nrow(X)))
	return(retList);
}

readFileInteger <- function(args)
{
	X = as.matrix(read.table(args[[1]], sep = args[[2]], colClasses = "integer")); 
	retList <- list("matrix" = t(X),
                 "ncols" = as.integer(ncol(X)),
                 "nrows" = as.integer(nrow(X)))
	return(retList);
}

# CheckPreprocessFilesExist <- function(tmp.dir) {
#	files.exist.flag = FALSE
#	if (file.exists())
#
#}