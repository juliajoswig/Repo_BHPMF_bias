# transforms the matrix back after gapfilling

transform_back <- function(X, pars){
  if(length(pars) != ncol(X)) stop("length(pars) must be equal to number of cols X")
  for(i in 1:ncol(X)){
    X[,i] <- exp( X[,1] * pars[[i]]$slogx + pars[[i]]$mlogx ) + min_x - 1
  }
  X
}

# back_trans_pars <- list()
# # transform data and save data for backtransformation
# for(i in 1:ncol(X)){
#   x <- X[,i]
#   min_x <- min(x,na.rm = T)
#   x <- x - min_x + 1  # make this optional?
#   logx <- log(x)
#   mlogx <- mean(logx, na.rm = T)
#   slogx <- sd(logx, na.rm = T)
#   x <- (logx - mlogx)/slogx
#   back_trans_pars[[i]] <- list(min_x = min_x,
#                                mlogx = mlogx,
#                                slogx = slogx)
#   X[,i] <- x
# }