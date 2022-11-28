

## created by Julia S. Joswig
## email: juliajoswigjj@gmail.com
## For the manuscript 
## Gap-filling of plant trait data with BHPMF induces taxonomic patterns
## Joswig, Kattge, Kramer, Mahecha, Rüger, Schaepman, Schrodt, Schuman 
## date 20221128 


# Content
# 1. Create the data sets
# 2. Set gaps (3 Repetitions)
# 3. BHPMF impute
# 4. Analyse & plot


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
# run Repo_git/script/1.1_Preprocess/xxx.R
# writes "data/MasterData/data/Obs_obs_TD/MasterData_inclPFT.csv"
# writes "data/MasterData/data_2/Obs_obs_TD/MasterData_inclPFT.csv"
# requires "xx"
# loads "xxx"
