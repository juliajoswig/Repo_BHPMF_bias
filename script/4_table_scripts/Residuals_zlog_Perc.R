# small error for SSD or LDMC, correct.

#------------------------------------------------------------
# define path
#------------------------------------------------------------
setwd("/..")
origin = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_git"
originData = "Volumes/Data_JJoswig/BGC/projects_BGC/2016_GapFilling/Repo_data"
list.files(file.path(origin,"script"))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"script","helper_scripts","fn_load_functions.R"))
load_functions(origin)

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices(originData)
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
GapPercent=50
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3


#----------------------------------------
# chose data
#----------------------------------------
t_choice="data"
if(t_choice=="data"){units <- c("mm2 mg-1","m","mm2 mg-1","mg g-1","mg g-1","g m-2")}
if(t_choice=="data_2"){units <- c("mm2 mg-1","m","","","")}

if(t_choice=="data"){RepNum=1}
if(t_choice=="data_2"){RepNum=2}



#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(originData,"analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")

ObsOrTD="Obs_obs_TD"
Percent=80
res <- read.csv(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))


ObsOrTD="Obs_obs"
Percent=80
resENV <- read.csv(file.path(originData,"analyses","Point_wise",RepNum,t_choice,ObsOrTD,Percent,paste0("res_2021_08.csv")))

colz=c("#ef8a62","#fddbc7","#f7f7f7","#d1e5f0","#67a9cf","#2166ac")
colzENV="#b2182b"

colz=c("#f4a582","#fddbc7","#d1e5f0","#92c5de","#4393c3","#2166ac","#b2182b")
res <- res[,colSums(!is.na(res))!=0]
trait_names=as.vector(unique(res$trait))
trait_names <- trait_names[!is.na(trait_names)]
missingness = unique(as.vector(res$missingness))
missingness <- missingness[!is.na(missingness)]
m=1
t=2
w=1




# get some numbers for the paper:
table_now <- matrix(NA,ncol=8,nrow=length(trait_names))
table_2 <- matrix(NA,ncol=8,nrow=length(trait_names))
ix_trait=res$trait==trait_names[t]
ix_trait[is.na(ix_trait)] <- FALSE
t=3
for(t in 1:length(trait_names)){
  ix_trait=res$trait==trait_names[t]
  ix_trait[is.na(ix_trait)] <- FALSE
  
  bxpl=res$value_pred_zlog[ix_trait]- res$value_obs_zlog[ix_trait]
  bxplENV=resENV$value_pred_zlog[ix_trait]- resENV$value_obs_zlog[ix_trait]
  
  table_now[t,1] <- quantile(abs(bxpl),probs = .25)
  table_now[t,2] <- quantile(abs(bxpl),probs = .5)
  table_now[t,3] <- quantile(abs(bxpl),probs = .75)
  table_now[t,4] <- max(abs(bxpl))
  table_now[t,5] <- quantile(abs(bxplENV),probs = .25)
  table_now[t,6] <- quantile(abs(bxplENV),probs = .5)
  table_now[t,7] <- quantile(abs(bxplENV),probs = .75)
  table_now[t,8] <- max(abs(bxplENV))
  
  print(quantile(abs(bxpl),probs = 1)<c(dist(range(abs(res$value_obs)[ix_trait]))))
  print(quantile(abs(bxpl),probs = 1)<c(dist(range(abs(res$value_obs)[ix_trait]))))
  table_2[t,1] <- quantile(abs(bxpl),probs = .25)/c(dist(range(res$value_obs[ix_trait])))*100
  table_2[t,2] <- quantile(abs(bxpl),probs = .5)/c(dist(range(res$value_obs[ix_trait])))*100
  table_2[t,3] <- quantile(abs(bxpl),probs = .75)/c(dist(range(res$value_obs[ix_trait])))*100
  table_2[t,4] <- quantile(abs(bxpl),probs = 1)/c(dist(range(res$value_obs[ix_trait])))*100
  table_2[t,5] <- quantile(abs(bxplENV),probs = .25)/c(dist(range(resENV$value_obs[ix_trait])))*100
  table_2[t,6] <- quantile(abs(bxplENV),probs = .5)/c(dist(range(resENV$value_obs[ix_trait])))*100
  table_2[t,7] <- quantile(abs(bxplENV),probs = .75)/c(dist(range(resENV$value_obs[ix_trait])))*100
  table_2[t,8] <- quantile(abs(bxplENV),probs = 1)/c(dist(range(resENV$value_obs[ix_trait])))*100
}   

print(table_now[,c(4,8)])
print(table_2[,c(4,8)])

#table_now <- table_now[c(-4,8)]
#table_2 <- table_2[c(-4,8)]
if(t_choice=="data"){
colnames(table_now) <- c("25th quantile TD","Median TD","75th quantile TD","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
rownames(table_now) <- trait_names
colnames(table_2) <- c("25th quantile TD","Median TD","75th quantile TD","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
}
if(t_choice=="data_2"){
  colnames(table_now) <- c("25th quantile TD2","Median TD2","75th quantile TD2","Max TD2",
                           "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
  rownames(table_now) <- trait_names
  colnames(table_2) <- c("25th quantile TD2","Median TD2","75th quantile TD2","Max TD",
                         "25th quantile envelope","Median  envelope","75th quantile  envelope","Max envelope")
}
rownames(table_2) <- trait_names
require(xtable)

#Absolute residuals of imputed - observed for test data imputed with and without the envelope (back-transformed data).
xtable(table_now)

#Residuals relative to trait range. Residuals (absolute) calculated as  imputed - observed for test data imputed with and without the envelope (back-transformed data).
xtable(table_2)

