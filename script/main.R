

## created by Julia S. Joswig
## email: juliajoswigjj@gmail.com
## For the manuscript 
## Gap-filling of plant trait data with BHPMF induces taxonomic patterns
## Joswig, Kattge, Kramer, Mahecha, Rüger, Schaepman, Schrodt, Schuman 
## date 20221128 


# Content
# 1. Create the data sets
# 2. Set gaps (3 Repetitions)
# 3. zlog transform 
# 4. BHPMF impute
# 5. Back transform
# 6. Analyse & plot


# ------------------------------------------------------------------------------
# 1. Create the data sets
# ------------------------------------------------------------------------------

# 1.1
# run Repo_git/script/1.1_Preprocess/00_MasterData/Load_MasterData.R
# requires "readxl"
# loads "TRY_pmf_DataRelease_cleaned_2013_06_19.xlsx"
# writes "data/MasterData/data/Obs_obs/MasterData.csv"  # ExTD
# writes "data/MasterData/data/Obs_obs_TD/MasterData.csv" # TD
# writes "data/MasterData/data_2/Obs_obs_TD/MasterData.csv"# TD2
# writes "data/MasterData/data_2/Obs_obs/MasterData.csv"# ExTD2

# 1.2
# run Repo_git/script/1.1_Preprocess/00_MasterData/create_PFTs.R
# reads "data/MasterData/data/Obs_obs/MasterData.csv"  # ExTD
# reads "data/MasterData/data/Obs_obs_TD/MasterData.csv" # TD
# reads "data/MasterData/data_2/Obs_obs_TD/MasterData.csv"# TD2
# reads "data/MasterData/data_2/Obs_obs/MasterData.csv"# ExTD2
# adds PFT
# creates "data/MasterData/data/Obs_obs_TD/MasterData_inclPFT.csv"
# creates "data/MasterData/data/Obs_obs/MasterData_inclPFT.csv"
# creates "data/MasterData/data_2/Obs_obs_TD/MasterData_inclPFT.csv"
# creates "data/MasterData/data_2/Obs_obs/MasterData_inclPFT.csv"

# ------------------------------------------------------------------------------
# 2. Set gaps (3 Repetitions)
# ------------------------------------------------------------------------------
# 2.1
# run Repo_git/script/1.1_Preprocess/01_Set_Missingness/do_missingness.R
# loads "data/MasterData/data/Obs_obs_TD/MasterData_inclPFT.csv"
# loads "data/MasterData/data/Obs_obs/MasterData_inclPFT.csv"
# loads "data/MasterData/data_2/Obs_obs_TD/MasterData_inclPFT.csv"
# loads "data/MasterData/data_2/Obs_obs/MasterData_inclPFT.csv"

# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/taxInfo.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/funInfo.csv"  
# * "data" and "data_2" refers to TD and TD2
# ** missingness levels from 0 to 80%
# *** Obs_obs refers to the ExTD (ExTD2) and Obs_obs_TD refers to TD (TD2)


# ------------------------------------------------------------------------------
# 3. zlog transform
# ------------------------------------------------------------------------------
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo_log.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo_zlog.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/zlog_transform.RData"


# ------------------------------------------------------------------------------
# 4.  BHPMF imputes
# ------------------------------------------------------------------------------
# requires "BHPMF"
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/taxInfo.csv"
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/traitInfo_zlog.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/mean.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/sd.csv"
# * "data" and "data_2" refers to TD and TD2
# ** missingness levels from 0 to 80%
# *** Obs_obs refers to the ExTD (ExTD2) and Obs_obs_TD refers to TD (TD2)


# ------------------------------------------------------------------------------
# 5.  Back transform (crucial!)
# ------------------------------------------------------------------------------
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo.csv"
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/traitInfo_zlog.csv"
# writes "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/mean.csv"
# loads "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/taxInfo.csv"

# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo_pred.csv"
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfo_pred_zlog.csv"
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/TDzlog_transform.RData“
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfoTD_pred_REzlog.csv“
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfoTD_obs.csv“
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfoTD_obs_zlog.csv“
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfoTD_pred_zlog.csv“
# write "data/_runs/Rep_<RepNum>/data[or data_2]*/p_<0-80>**/Obs_obs[or Obs_obs_TD]***/data/traitInfoTD_pred.csv“
# * "data" and "data_2" refers to TD and TD2
# ** missingness levels from 0 to 80%
# *** Obs_obs refers to the ExTD (ExTD2) and Obs_obs_TD refers to TD (TD2)

# ------------------------------------------------------------------------------
# 6.Analyse 
# ------------------------------------------------------------------------------




# ------------------------------------------------------------------------------
# 6.Plot figures and tables
# ------------------------------------------------------------------------------

# Figure 1
# Figure 2
# Figure 3
# Figure 4
# Figure 5
# Figure 6

