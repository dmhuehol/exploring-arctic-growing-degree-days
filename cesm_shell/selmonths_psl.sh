#######################################
# selmonths_psl
# Select months
# Output files are named:
# [original_name]_AMJJAS.nc
# 
# Globals:
#   None
# Written by Daniel Hueholt
# Graduate Research Assistant at Colorado State University
#######################################

### Input variables
IN_PATH="/Users/dhueholt/Documents/gddt_data/LENS2/merge_PSL/" #Path to data
IN_TOKEN="*.nc" #Token to match data files, e.g. "*.001.*.nc"
SEL_MONTHS=4,5,6,7,8,9
OUT_APPENDIX="_AMJJAS.nc"
OUT_PATH="/Users/dhueholt/Documents/gddt_data/LENS2/merge_PSL/AMJJAS/" #Path to save full dataset files

### Script
IN_CARD="$IN_PATH$IN_TOKEN"
PATH_LENGTH=${#IN_PATH}
echo $PATH_LENGTH

### Bits for filename
for f in $IN_CARD; do
  ACTIVE_FNAME=${f:$PATH_LENGTH} 
  ACTIVE_NONC=$(echo $ACTIVE_FNAME | cut -d'.' -f1)
  OUT_NAME=$ACTIVE_NONC$OUT_APPENDIX
  # echo $OUT_NAME
  OUT_CARD=$OUT_PATH$OUT_NAME
  # echo $OUT_CARD
  cdo -selmonth,$SEL_MONTHS $f $OUT_CARD
done