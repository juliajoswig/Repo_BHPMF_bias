

add_col_to_res <- function(new.col.names,res){
  res = as.matrix(res)

  for(i in 1:length(new.col.names)){
    new.col=rep(NA,nrow(res))
    if(sum(colnames(res)%in%new.col.names[i])<1){
      res <-  cbind(res,new.col)  
      colnames(res)[ncol(res)] = new.col.names[i]
      print(paste(new.col.names[i],"added"))
    }
  }
  
  res = as.matrix(res)
  return(res)
}