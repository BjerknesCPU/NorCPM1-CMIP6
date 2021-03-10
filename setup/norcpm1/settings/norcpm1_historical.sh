# EXPERIMENT DEFAULT SETTINGS 
# DO NOT CHANGE, USE VARNAME=VALUE ARGUMENT WHEN CALLING SCRIPT TO OVERRIDE DEFAULTS 

# experiment settings
: ${CASE_PREFIX:=`basename $1 .sh`} # case prefix, not including _YYYYMMDD_XX suffix ; extracted from name of this file
: ${REST_PREFIX:=N1850AERCNOC_f19_g16_CMIP6forcings_15} 
: ${REST_PATH_REMOTE:=}
: ${REST_PATH_LOCAL:=$INPUTDATA/ccsm4_init}
: ${START_YEARS:=1850} # multiple start dates only for prediction
: ${START_MONTHS:=01} # multiple start dates only for prediction
: ${START_DAYS:=01} # multiple start dates only for prediction
: ${REF_YEARS:=0051} # multiple reference dates only for RUN_TYPE=hybrid
: ${REF_MONTHS:=01} # multiple reference dates only for RUN_TYPE=hybrid
: ${REF_DAYS:=01} # multiple reference dates only for RUN_TYPE=hybrid
: ${REF_ENSEMBLE:=1} # set to 1 if ensemble of perturbed intial conditions with same start date, only for RUN_TYPE=hybrid
: ${RUN_TYPE:=hybrid} # use "branch" if unspecified 
: ${ENSSIZE:=2} # number of members 
: ${STOP_OPTION:=nmonths} # units for run length specification STOP_N 
: ${STOP_N:=1} # run continuesly for this length 
: ${RESTART:=59} # restart this many times (set to 5 years here) 
: ${WALLTIME:='96:00:00'}  

# general settings 
: ${VERSION:=`basename $SETUPROOT`}
: ${CASESROOT:=$SETUPROOT/../../cases/$VERSION}
: ${CCSMROOT:=$SETUPROOT/../../model/$VERSION}
: ${COMPSET:=N20TREXTAERCNOCCMIP6}
: ${PECOUNT:=L} # T=32, S=64, M=96, L=128, X1=502
: ${RES:=f19_g16}
: ${ACCOUNT:=nn9039k}
: ${ASK_BEFORE_REMOVE:=0} # 1=will ask before removing existing cases 
: ${MEMBERTAG:=mem} # leave empty or set to 'mem' 
: ${MAX_PARALLEL_STARCHIVE:=30} 
: ${VERBOSE:=1} # set -vx option in all scripts

# derived settings
: ${START_YEAR1:=`echo $START_YEARS | cut -d" " -f1`}
: ${START_MONTH1:=`echo $START_MONTHS | cut -d" " -f1`}
: ${START_DAY1:=`echo $START_DAYS | cut -d" " -f1`}
: ${REF_YEAR1:=`echo $REF_YEARS | cut -d" " -f1`}
: ${REF_MONTH1:=`echo $REF_MONTHS | cut -d" " -f1`}
: ${REF_DAY1:=`echo $REF_DAYS | cut -d" " -f1`}
: ${SCRIPTSROOT:=$CCSMROOT/scripts}
