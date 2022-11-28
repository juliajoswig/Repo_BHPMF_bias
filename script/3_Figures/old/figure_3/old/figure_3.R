figuree_correl_scatter <- function(){
  #install.packages("psych")
  library(psych)
  
  colz1=c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
  colz=c("#fff7ec","#fee8c8","#fdd49e","#fdbb84","#fc8d59","#ef6548","#d7301f","#b30000","#7f0000")
  
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
  Version_now="V3"
  list.files(file.path(origin,"script","analysis",Version_now))
  
  #------------------------------------------------------------
  # load some functions
  #------------------------------------------------------------
  source(file.path(origin,"script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
  load_functions()
  
  #------------------------------------------------------------
  # define data set approaches/choices
  #------------------------------------------------------------
  out <- choices()
  tsubs <- out$tsubs
  TD_choices =  out$TD_choices
  repnums = out$repnums
  gappercents = out$gappercents
  whichDataSet = out$whichDataSet
  ObsSpec = out$ObsSpec
  obsspec = ObsSpec
  preparation = out$preparation
  trait_guido =out$trait_guido
  trait_rainfor =out$trait_rainfor
  colz1 =out$colz1
  colz2 =out$colz2
  Indeces <- c("rmse_man")

  
  
  
  # Load observed correlations
  TD_choice="Obs_obs_TD"
  GapPercent=-1
  RepNum=2 
  wds=1
  wtd=1
  whichDataSet=1:2
  obsspec=1
  repnums=c(sample(1:30,3))
  gappercents=c(-1,0,20,30,40,50,60,70)
  for(RepNum in repnums){
    for(wds in whichDataSet){
      for(wtd in obsspec){
        trait_sub = tsubs[wds]
        TD_choice = TD_choices[wtd]
       for(GapPercent in gappercents){
          
          # define dat path
          if(GapPercent==-1){GapPercent1=0}else{GapPercent1=GapPercent}
          dat_path_org1 <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","TestData_org.RData")
          dat_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),TD_choice,"data","HPMF_filled.csv")
          gap_ix_path <- file.path(origin,"runs",paste0("Rep",RepNum), paste0("p",GapPercent1,"_",trait_sub),paste0(TD_choice),"data","Negap.ix.csv")
          
          if((file.exists(dat_path)&file.exists(dat_path_org1))|
             (GapPercent==-1&file.exists(dat_path_org1))){
            
            #--------------------------------------------------------------
            # load data
            #--------------------------------------------------------------
            # load predicted data
            traitInfo_hpmf <- as.matrix(read.table(file = dat_path, sep=",", dec="."))
            # load original test data "TestData_tot" containing tax, info and ID data
            load(dat_path_org1) 
            # load function information
            GF <- read.table(file.path(origin,"runs","META","GF_Obs.csv"),sep=",",dec=".")
            PFT <- read.table(file.path(origin,"runs","META","PFT_Obs.csv"),sep=",",dec=".")
            
            if(GapPercent==-1){
              # fit colnames
              traitInfo_observed <- TestData_tot[,colnames(TestData_tot)%in%colnames(traitInfo_hpmf)]
              NEW_gap_ix <- as.matrix(read.table(file = gap_ix_path,sep=",", dec="."))
              
              # ----------------------------------------------------------------------------------------------------------
              # Rescale the original data
              # ----------------------------------------------------------------------------------------------------------
              # load scaling factors
              load(file.path(origin,"runs",paste0("Rep",RepNum),paste0("Meta_",trait_sub),"scaling_factors.RData"))
              
              if(TD_choice =="Spec_spec"){    
                for(i in 1:ncol(traitInfo_observed)){
                  traitInfo_observed[,i] <- (log(traitInfo_observed[,i]) - scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_hpmf)[i]),4])/# changed traitInfo to traitInfo_hpmf
                    scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_observed)[i]),3]
                }
              }else{
                for(i in 1:ncol(traitInfo_observed)){
                  traitInfo_observed[,i] <- (log(traitInfo_observed[,i]) - scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_hpmf)[i]),2])/# changed traitInfo to traitInfo_hpmf
                    scaling.factors[which(rownames(scaling.factors)==colnames(traitInfo_observed)[i]),1]
                }
              }
            }
            
            if(GapPercent!=-1){traitInfo <- traitInfo_hpmf}
            if(GapPercent==-1){traitInfo <- traitInfo_observed}
    
            pdf(file=file.path(origin,"plots","figure_3",paste0("Pairs",GapPercent,trait_sub,TD_choice,".pdf")),
                width = 8,height=8)
            par(mfrow=c(1,1))
            
            pairs.panels(traitInfo,pch=16,cex=.4, cex.cor = 1,
                         method = "pearson", # correlation method
                         hist.col = colz1[8],ylim=c(-4,4),
                         density = TRUE,  # show density plots
                         ellipses = TRUE # show correlation ellipses
            )
            dev.off()
            
          }
       }
        }
    }
  }

  }
  