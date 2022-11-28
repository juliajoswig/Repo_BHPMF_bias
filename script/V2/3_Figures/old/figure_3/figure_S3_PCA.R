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
Version_now="V1"
list.files(file.path(origin,"script","analysis",Version_now))

#------------------------------------------------------------
# load some functions
#------------------------------------------------------------
source(file.path(origin,"_2021","script","analysis",Version_now,"helper_scripts","fn_load_functions.R"))
load_functions(origin = origin,Version_now = "V1")

#------------------------------------------------------------
# define data set approaches/choices
#------------------------------------------------------------
out <- choices()
tsubs <- out$tsubs
ObsOrTDs =  out$ObsOrTDs
repnums = out$repnums
gappercents = out$gappercents
whichDataSet = out$whichDataSet
ObsSpec = out$ObsSpec
obsspec = ObsSpec
preparation = out$preparation
trait_guido =out$trait_guido
trait_rainfor =out$trait_rainfor
colz1 = out$colz1
colz2 = out$colz2
gappercents= c("1","5","10","20","30","40","50","60")
repnums=1:3

GapPercent=50
RepNum=1

gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3

#-------------------------------------------------------------------
# load trait data   
#-------------------------------------------------------------------
file.path(origin,"_2021","data","analyes","Point_wise","res.csv")
colorset1 <- c("#b2e2e2","#66c2a4","#2ca25f","#006d2c")
colz=c("#b2182b","#92c5de","#2166ac")
aph=.5
colz_add_alpha=c(rgb(178, 24, 43,maxColorValue = 255,alpha = aph),
                 rgb(214, 96, 77,alpha = aph,maxColorValue = 255),
                 rgb(244, 165, 130,alpha = aph,maxColorValue = 255),
                 rgb(253, 219, 199,alpha = aph,maxColorValue = 255),
                 rgb(146, 197, 222,alpha = aph,maxColorValue = 255),
                 rgb(67, 147, 195,alpha = aph,maxColorValue = 255),
                 rgb(33, 102, 172,alpha = aph,maxColorValue = 255))

colz_add_alpha <- rep(rgb(red = 1,green = 1,blue = 1,alpha = .1),6)


GapPercent=50
RepNum=1
t_choice="data"
ObsOrTD="Obs_obs_TD"
gappercents=c(0,1,5,10,20,30,40,50,60,70,80)
t_choices=c("data","data_2")
TDnos=c("Obs_obs_TD","Obs_obs")
repnums=3
td=1
p=3
res_matrix_name="res_20210303"
res <- read.csv(file=file.path(origin,"_2021","data","analyses","TOTAL",paste0(res_matrix_name,".csv")))

res <- res[,colSums(!is.na(res))!=0]
res_now <- res[res$TraitChoice=="data_2",]
  trait_names=as.vector(unique(res_now$trait))
  trait_names <- trait_names[!is.na(trait_names)]
m=1
t=2
w=1


  t=1
  pdf(file=file.path(origin,"_2021","Figures","Figure_3","Figure_S3_PCA.pdf"),width=6,height=8)
  par(mfrow=c(1,1),mar=c(6,8,2,2))
  
    plot(res_now$GapPercent,res_now$PCA1,ylim=c(0,1),ylab="Variance Explained",xlab="",
         main="",col=colz[1],pch=16,cex.axis=3,cex.lab=2.2,cex=1,yaxt="n",xaxt="n")
    axis(side = 1,line = 2,at = seq(0,80,40),labels = c("","Missingness [%]",""),tick = FALSE,cex.axis=2)
    axis(side = 1,at = seq(0,80,20),labels = seq(0,80,20),tick = TRUE,cex.axis=2)
    
    axis(side = 2,at = seq(0,1,.2),labels = seq(0,1,.2),tick = TRUE,cex.axis=2)
#    axis(side = 2,line = 1,at = seq(0,.8,.4),labels = c("","Variance Explained",""),tick = FALSE,cex.axis=2)
    abline(v = seq(0,80,20),col="gray",lty=2)
    abline(h = seq(0,.80,.20),col="gray",lty=2)
    
    lines(aggregate(x = res_now$PCA1,list(res_now$GapPercent),FUN=median,na.rm=TRUE),ylim=c(0,1),
         main="",col=colz[1],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=2,lwd=4)
    points(res_now$GapPercent,res_now$PCA2,ylim=c(0,1),
         main="",col=colz[2],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=1)
    lines(aggregate(x = res_now$PCA2,list(res_now$GapPercent),FUN=median,na.rm=TRUE),ylim=c(0,1),
          main="",col=colz[2],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=2,lwd=4)
    points(res_now$GapPercent,res_now$PCA3,ylim=c(0,1),
         main="",col=colz[3],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=1)
    lines(aggregate(x = res_now$PCA3,list(res_now$GapPercent),FUN=median,na.rm=TRUE),ylim=c(0,1),
          main="",col=colz[3],pch=16,cex.main=4,cex.axis=3,cex.lab=2.2,cex=2,lwd=4)
    
    legend(x = 40,y = .98,pch = 16,col =colz[1:3],legend = c("Axis 1","Axis 2","Axis 3"),cex=2,border = FALSE,bty = "n",fill = "white")
  dev.off()
  
  