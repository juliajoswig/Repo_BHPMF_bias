create_res <- function(res_matrix_name, tsubs,TD_choices,repnums = out$repnums,gappercents,
                       whichDataSet,ObsSpec){
  
  res <- matrix(NA,ncol=4,nrow=10000)
  colnames(res)[1:4] = c("Repetition","Obs_or_Spec","TraitChoice","GapPercent")        
  
  RepNum=1
  n=0
  ts=2
  for(RepNum in repnums){
    for(ts in whichDataSet){
      for(td in ObsSpec){
        for(GapPercent in gappercents){
          
          
          trait_sub = tsubs[ts]
          TD_choice = TD_choices[td]
          
          n=n+1
          res[n,1:4] <- c(RepNum,TD_choice,trait_sub,GapPercent)  
        }
      }
    }
  }    
  
  res <- as.matrix(res[!is.na(res[,4]),])
  dim(res)
  head(res)
  write.table(res,file=file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
}