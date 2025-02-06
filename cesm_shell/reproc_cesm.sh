######################################
# reproc_cesm
# Sometimes you process CESM data, and then realize you need to run
# something else (often selyear or selmonth) later on. This script 
# makes that process easy! Designed for LENS2.
# Output files are named manually in the CUST_NAME variable. Watch out 
# for this!
# 
# Globals:
#   None
# Written by Daniel Hueholt
# Graduate Research Assistant at Colorado State University
#######################################

### Input variables
# IN_PATH="/Users/dhueholt/Documents/gddt_data/LENS2/monthly_SST/1935_1965/" #Path to data
IN_PATH="/Volumes/Polycrystal/Data/LENS2/merge_PSL/AMJJAS/" #Path to data
# IN_PATH="/barnes-scratch/DATA/CESM2-LE/raw_data/monthly/PSL/" # ASHA
IN_TOKEN="*.nc" #Token to match data files, e.g. "*.001.*.nc"
# OUT_PATH="/Volumes/Polycrystal/Data/LENS2/monthly_SST/merge/1935_1965/" #Path to save full dataset files
# OUT_PATH="/barnes-scratch/DATA/CESM2-LE/processed_data/monthly/PSL/" #Path to save full dataset files
OUT_PATH="/Volumes/Polycrystal/Data/LENS2/merge_PSL/ann/"

### Script
IN_CARD="$IN_PATH$IN_TOKEN"
PATH_LENGTH=${#IN_PATH}
echo $PATH_LENGTH

for f in $IN_CARD; do
  ACTIVE_FNAME=${f:$PATH_LENGTH}
  NO_EXT=$(echo $ACTIVE_FNAME | cut -d'.' -f1)
  NO_TIME=$(echo $NO_EXT | rev | cut -d'_' -f2- | rev)
  CUST_NAME=${OUT_PATH}${NO_TIME}_1850-2100_AMJJAS.nc # WATCH OUT
  echo $CUST_NAME
  # cdo -L -selyear,1935/1965 $f $CUST_NAME
  cdo -L -yearmonmean $f $CUST_NAME
done