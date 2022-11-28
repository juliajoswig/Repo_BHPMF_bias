
#------------------------------------------------------------
# define path
#------------------------------------------------------------
is.it.on.cluster=FALSE
if(is.it.on.cluster){
  setwd("/..")
  setwd(file.path("Net","Groups","BGI"))
  origin=file.path("work_1","2016_GapFilling")}
if(!is.it.on.cluster){
  setwd("/..")
  origin = "Volumes/bgi/work_1/2016_GapFilling"
}
Version_now="V2"
list.files(file.path(origin,"_2021","script",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin,Version_now)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
t_choices <- out$tsubs
TDnos = out$TD_choices
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido = out$trait_guido
trait_rainfor = out$trait_rainfor
colz1 = out$colz1
colz2 = out$colz2
new.mean.fun = out$new.mean.fun
new.sd.fun = out$new.sd.fun
add_col_to_res <- out$add_col_to_res
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)



#------------------------------------------------------------
# create results list
#------------------------------------------------------------
res_matrix_name="res_20210303"#"res_20201112"
#create_res(res_matrix_name = res_matrix_name,tsubs = t_choices,TD_choices = TDnos,repnums = 1:3,gappercents = gappercents,whichDataSet =whichDataSet,ObsSpec =  1:2)
#-------------------------------------------------------------------
# Loop through all runs and load the indices
#-----------------------------------------------------------------
  res <- read.table(file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
  RepNum=1
  GapPercent=5
  ts=1
  td=1
  repnums=1:3
  
  Indeces=c("RMSE","Sil","PCA","Corr")
  #Indeces=c("RMSE")
  
  GapPercent=50
  RepNum=1
  t_choice="data"
  ObsOrTD="Obs_obs_TD"
  gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
  t_choices=c("data","data_2")
  TDnos=c("Obs_obs_TD","Obs_obs")
  repnums=3
  td=1
  tc=1
  p=3
  
  for(RepNum in 1:3){
    for(tc in 1:2){
      for(p in 1:length(gappercents)){
        for(td in 1:2){
          
          
        
        print("#####################")
        Percent = gappercents[p] 
        t_choice =t_choices[tc]
        ObsOrTD = TDnos[td]
        print(paste(RepNum,ObsOrTD,t_choice,Percent))
        print("#####################")
        
        res[,4] <- as.numeric(res[,4])
        
        
        n <- which( as.numeric(res$Repetition)== RepNum & 
                      res[,2]== ObsOrTD & 
                      res[,3]== t_choice & 
                      as.numeric(res[,4])== Percent)
        I=4
        for(I in 1:length(Indeces)){
       
        Index_now = Indeces[I]
        
        #--------------------------------------------------------------
        # check if data available
        #--------------------------------------------------------------
        # PCA
        if(Index_now == "PCA"){
          pca_path = file.path(origin,"_2021","data","analyses","PCA",t_choice,ObsOrTD,Percent,RepNum,paste0("PCA.csv"))
          res <- add_col_to_res("PCA1",res)
          res <- add_col_to_res("PCA2",res)
          res <- add_col_to_res("PCA3",res)
          
          if(file.exists(pca_path)){
            print("PCA")
            print(Percent)
            pca_out <- read.csv(pca_path)
            whichcol = which(colnames(res)%in%c("PCA1","PCA2","PCA3"))
            try(res[n,whichcol] <-  as.numeric(as.vector(pca_out[1,])))
            try(rm(pca_out))
          }
        } 
        
        # sil   
        # set path
        if(Index_now=="Sil"){
            print("Silhouette Index")
            if(Percent==-1){Percent1="org"}else{Percent1=Percent}
            sil_path=file.path(origin,"_2021","data","analyses","Silhouette",t_choice,ObsOrTD,Percent,RepNum,"Silhouette.RData")
            if(file.exists(sil_path)){
              load(sil_path)

              if(ObsOrTD!="Spec_spec"&ObsOrTD!="Spec_spec_TD"){
                cl=1
                res <- add_col_to_res(paste0("Sil_Spec_md"),res) 
                res <- add_col_to_res(paste0("Sil_Spec_mdNO1"),res) 
                
                  res[n, which(colnames(res)==paste0("Sil_Spec_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
                  res[n, which(colnames(res)==paste0("Sil_Spec_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
                  
                print("SPEC")
            }
              cl=2
              res <- add_col_to_res(paste0("Sil_Gen_md"),res) 
              res <- add_col_to_res(paste0("Sil_Gen_mdNO1"),res) 
              
              res[n, which(colnames(res)==paste0("Sil_Gen_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
              res[n, which(colnames(res)==paste0("Sil_Gen_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
              
              cl=3
              res <- add_col_to_res(paste0("Sil_Fam_md"),res) 
              res <- add_col_to_res(paste0("Sil_Fam_mdNO1"),res) 
              
              res[n, which(colnames(res)==paste0("Sil_Fam_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
              res[n, which(colnames(res)==paste0("Sil_Fam_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
              
              cl=4
              res <- add_col_to_res(paste0("Sil_Clad_md"),res) 
              res <- add_col_to_res(paste0("Sil_Clad_mdNO1"),res) 
              
              res[n, which(colnames(res)==paste0("Sil_Clad_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
              res[n, which(colnames(res)==paste0("Sil_Clad_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
              
              cl=5
              res <- add_col_to_res(paste0("Sil_GF_md"),res) 
              res <- add_col_to_res(paste0("Sil_GF_mdNO1"),res) 
              
              res[n, which(colnames(res)==paste0("Sil_GF_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
              res[n, which(colnames(res)==paste0("Sil_GF_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
              
              cl=6
              res <- add_col_to_res(paste0("Sil_PFT_md"),res) 
              res <- add_col_to_res(paste0("Sil_PFT_mdNO1"),res) 
              
              res[n, which(colnames(res)==paste0("Sil_PFT_md" ))] <- median(as.numeric(sil_l[[cl]][,4]))
              res[n, which(colnames(res)==paste0("Sil_PFT_mdNO1" ))] <-median(as.numeric(sil_l[[cl]][as.numeric(sil_l[[cl]][,3])>1,4]))
              try(rm("sil_l"))
              
              
            }
          }
          


      # RMSE
      if(Index_now == "RMSE"){
        print("RMSE")            

        path_rmse <-  file.path(origin,"_2021","data","analyses",
                                      "RMSE",t_choice,ObsOrTD,Percent,RepNum,"RMSE.csv")
        
        if(file.exists(path_rmse)){
          all_out <- read.csv(file.path(origin,"_2021","data","analyses",
                             "RMSE",t_choice,ObsOrTD,Percent,RepNum,"RMSE.csv"),header = TRUE,row.names = NULL)
          
          for(i in 1:ncol(all_out)){
            res <- add_col_to_res(paste0("RMSE_",colnames(all_out)[i]),res)
            res <- add_col_to_res(paste0("RMSE_gap_",colnames(all_out)[i]),res)
            res <- add_col_to_res(paste0("RMSE_zlog_",colnames(all_out)[i]),res)
            res <- add_col_to_res(paste0("RMSE_gap_zlog_",colnames(all_out)[i]),res)
          }
        
        # add entries
            whichcol = which(colnames(res)==paste0("RMSE_Total"))
            res[n,whichcol] <-  all_out[which(all_out[,1]=="no_trans"),which(colnames(all_out)=="Total")]
            whichcol = which(colnames(res)==paste0("RMSE_gap_Total"))
            res[n,whichcol] <-  all_out[which(all_out[,1]=="no_transNA"),which(colnames(all_out)=="Total")]
            
            whichcol = which(colnames(res)==paste0("RMSE_zlog_Total"))
            res[n,whichcol] <-  all_out[which(all_out[,1]=="zlog"),which(colnames(all_out)=="Total")]
            whichcol = which(colnames(res)==paste0("RMSE_gap_zlog_Total"))
            res[n,whichcol] <-  all_out[which(all_out[,1]=="zlogNA"),which(colnames(all_out)=="Total")]
            
            i=2
            for(i in 2:ncol(all_out)){
              whichcol = which(colnames(res)==paste0("RMSE_",colnames(all_out)[i]))
              res[n,whichcol] <-  all_out[which(all_out[,1]=="no_trans"),which(colnames(all_out)==colnames(all_out)[i])]
              whichcol = which(colnames(res)==paste0("RMSE_gap_",colnames(all_out)[i]))
              res[n,whichcol] <-  all_out[which(all_out[,1]=="no_transNA"),which(colnames(all_out)==colnames(all_out)[i])]
              
              whichcol = which(colnames(res)==paste0("RMSE_zlog_",colnames(all_out)[i]))
              res[n,whichcol] <-  all_out[which(all_out[,1]=="zlog"),which(colnames(all_out)==colnames(all_out)[i])]
              whichcol = which(colnames(res)==paste0("RMSE_gap_zlog_",colnames(all_out)[i]))
              res[n,whichcol] <-  all_out[which(all_out[,1]=="zlogNA"),which(colnames(all_out)==colnames(all_out)[i])]
              }
            }
        } 

      
        # Corr
      if(Index_now == "Corr"){
        cor_path=       file.path(origin,"_2021","data","analyses","Correl",t_choice,ObsOrTD,Percent,RepNum,paste0("Correl.csv"))
        cor_zlog_path = file.path(origin,"_2021","data","analyses","Correl",t_choice,ObsOrTD,Percent,RepNum,paste0("Correl_zlog.csv"))
        
        if(file.exists(cor_path)){
            cor_dat <- read.csv(cor_path)
            i=1
            j=2
            for(i in 1:ncol(cor_dat)){
              for(j in 1:ncol(cor_dat)){
                if(i!=j){
              res <-              add_col_to_res(paste0("cor_",colnames(cor_dat)[i],"_",colnames(cor_dat)[j]),res)
              whichcols = which(colnames(res)%in%paste0("cor_",colnames(cor_dat)[i],"_",colnames(cor_dat)[j]))
              res[n,whichcols] <-  cor_dat[i,j]
              }}}
            rm("cor_dat")
            print("Correlations")
          }

        if(file.exists(cor_zlog_path)){
          cor_dat <- read.csv(cor_zlog_path)
          i=1
          j=2
          for(i in 1:ncol(cor_dat)){
            for(j in 1:ncol(cor_dat)){
              if(i!=j){
                res <-              add_col_to_res(paste0("corzlog_",colnames(cor_dat)[i],"_",colnames(cor_dat)[j]),res)
                whichcols = which(colnames(res)%in%paste0("corzlog_",colnames(cor_dat)[i],"_",colnames(cor_dat)[j]))
                res[n,whichcols] <-  cor_dat[i,j]
                if(Percent==-1){print(res[n,c(4,whichcols)])}              
              }}}
          rm("cor_dat")
        }
        
        
        } 
        } 
        
        }
        }
        }

    }
      
  
  res=as.data.frame(res)

  write.table(res,file=file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")),sep=",",dec=".")
  
  