
load("../data/ProcessedData_pmf_org.Rda")
hierarchy.info <- read.table("../data/phyinfoAllMod", sep="\t")
source("cv_rmse.R")
source("gap_filling.R")
#out <- CalculateCvRmse(X, hierarchy.info, used.num.hierarchy.levels=2, num.latent.feats=15, tmp.dir="../data", num.folds = 2, num.samples=100, burn=10, gaps=2, tuning=FALSE)
pdf("../data/output/plot%03d.pdf")
GapFilling(X, hierarchy.info, 
           num.latent.feats=15, 
           mean.gap.filled.output.path="../data/output/mean_gap_filled.txt", 
           std.gap.filled.output.path="../data/output/std_gap_filled.txt", 
           num.samples=500, 
           burn=100, 
           gaps=2, 
           tmp.dir="../data")
dev.off()

