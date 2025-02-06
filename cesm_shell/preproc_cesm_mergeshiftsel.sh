######################################
# preproc_cesm_mergeshiftsel
# Merge files, shift times, and select variables and months, and do
# annual average from monthly data for CESM netCDF files. Designed for 
# LENS2. This is different enough from daily temperature processing that 
# I broke out the scripts for now.
# Output files are named:
# LE2_macro-micro_var_1850-2100_loc.nc
# See stackoverflow.com/a/13648438 for explanation of ENS_UNQ parsing.
# 
# Globals:
#   None
# Written by Daniel Hueholt
# Graduate Research Assistant at Colorado State University
#######################################

### Input variables
#  Local path to data
# IN_PATH="/Volumes/Polycrystal/Data/LENS2/monthly_SST/"
#  COE-HPC path to data
IN_PATH="/barnes-engr-scratch1/DATA/CMIP6/raw_data/CMIP/piControl/day_1/TREFHT/"
IN_TOKEN="*.nc" #Token to match data files, e.g. "*.001.*.nc"
SEL_VAR=TREFHT
#  Local path to save processed data
# OUT_PATH="/Volumes/Polycrystal/Data/LENS2/monthly_SST/merge/"
#  COE-HPC path to save processed data
OUT_PATH="/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/day_1/gdd_base5/"

### Script
IN_CARD="$IN_PATH$IN_TOKEN"
PATH_LENGTH=${#IN_PATH}
echo $PATH_LENGTH

### Bits for filename
for f in $IN_CARD; do
  ACTIVE_FNAME=${f:$PATH_LENGTH} 
  EXP_ENSMAC=$(echo $ACTIVE_FNAME | cut -d'.' -f5)
  ENSMIC=$(echo $ACTIVE_FNAME | cut -d'.' -f6)
  MATCH_ENS="${EXP_ENSMAC}.${ENSMIC}"
  ACTIVE_EXPENS+=( $MATCH_ENS )
  VARNAME=$(echo $ACTIVE_FNAME | cut -d'.' -f9)
done

### Merge files, select variables, mask ocean, reduce region, and convert unit
ENS_UNQ="$(tr ' ' '\n' <<< "${ACTIVE_EXPENS[@]}" | sort -u | tr '\n' ' ')"
echo $ENS_UNQ
RUN_TIME=()
for me in ${ENS_UNQ[@]}; do
    MATCH_MERGE="${IN_PATH}*${me}*"
    OUT_ENS=${me//-/_}
    OUT_ENS=${OUT_ENS//./-}
    MERGE_NAME=${OUT_ENS}_${VARNAME}_2015-2024.nc # WATCH OUT
    OUT_NAME=${OUT_PATH}$MERGE_NAME
    # cdo -L -yearmonmean -selmonth,4/9 -selvar,$SEL_VAR -shifttime,'-1days' -mergetime $MATCH_MERGE $OUT_NAME
    cdo -L -yearmonmean -sellonlatbox,0,360,50,90 -selvar,$SEL_VAR -shifttime,'-1days' -mergetime $MATCH_MERGE $OUT_NAME
done

### Remove temporary files
# rm ${TEMP_PATH}*.nc
