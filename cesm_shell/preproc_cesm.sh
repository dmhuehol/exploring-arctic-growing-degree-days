#######################################
# preproc_cesm
# Merge files, select variables, apply landmask, reduce region, and 
# convert unit for CESM TREFHT netCDF files. Designed for LENS2.
# Individual processing steps can be enabled/disabled by commenting
# relevant CDO line(s) without breaking automatic filename construction.
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
# IN_PATH="/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/" #Path to data
IN_PATH="/barnes-engr-scratch1/DATA/CMIP6/raw_data/CMIP/piControl/day_1/TREFHT/" # ASHA
IN_TOKEN="*.nc" #Token to match data files, e.g. "*.001.*.nc"
# MASK_F="/Users/dhueholt/Documents/gddt_data/mask/mask_cdo.nc"
MASK_F="/scratch/dhueholt/mask/mask_cdo.nc" # ASHA
SEL_VAR=TREFHT
# TEMP_PATH="/Volumes/Polycrystal/Data/LENS2/temp/" #Location for intermediate files
TEMP_PATH="/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/temp/" # ASHA
# FULLD_PATH="/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/fulld/" #Path to save full dataset files
FULLD_PATH="/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/temp/fulld/" # ASHA
# RED_PATH="/Volumes/Polycrystal/Data/LENS2/daily_TREFHT/refined/" #Path to save reduced dataset files
RED_PATH="/barnes-engr-scratch1/DATA/CMIP6/processed/CMIP/piControl/temp/refined/" # ASHA

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
    MERGE_NAME=${OUT_ENS}_${VARNAME}_1850-2100.nc
    TEMP_NAME=${TEMP_PATH}$MERGE_NAME
    cdo -mergetime $MATCH_MERGE $TEMP_NAME
    # FULLD_NAME=$FULLD_PATH$MERGE_NAME
    # cdo -L -div -selvar,$SEL_VAR $TEMP_NAME $MASK_F $FULLD_NAME
    ARC_NAME=${MERGE_NAME//.nc/_arc.nc}
    ANTAR_NAME=${MERGE_NAME//.nc/_antar.nc}
    OUT_ARC=$RED_PATH$ARC_NAME
    OUT_ANTAR=$RED_PATH$ANTAR_NAME
    cdo -L -setattribute,$SEL_VAR@units="degC" -subc,273.15 -selmonth,3,4,5,6,7,8,9,10 -sellonlatbox,0,360,50,90 $FULLD_NAME $OUT_ARC
    # cdo -L -setattribute,$SEL_VAR@units="degC" -subc,273.15 -sellonlatbox,0,360,-90,-60 $FULLD_NAME $OUT_ANTAR
done

### Remove temporary files
# rm ${TEMP_PATH}*.nc