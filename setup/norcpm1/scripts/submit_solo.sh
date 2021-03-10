#!/bin/sh -evx

SETUPROOT=`readlink -f \`dirname $0\``
. $SETUPROOT/scripts/source_settings.sh $*

echo PROPAGATE INDIVIDUAL ENSEMBLE MEMBERS AS SEPARATE JOBS

echo + IF CONTINUE_RUN NOT SET: CHECK IF SIMULATION SHOULD BE CONTINUED
if [ -z $CONTINUE_RUN ]
then 
  ENSEMBLE_PREFIX1=${CASE_PREFIX}_${START_YEAR1}${START_MONTH1}${START_DAY1}
  CASE1=${ENSEMBLE_PREFIX1}_${MEMBERTAG}01
  if [ $CASE1 == `head -1 $WORK/noresm/$ENSEMBLE_PREFIX1/$CASE1/run/rpointer.atm | cut -d. -f1` ] 
  then 
    CONTINUE_RUN=TRUE
  else
    CONTINUE_RUN=FALSE
  fi
fi
echo ++ CONTINUE_RUN set to $CONTINUE_RUN

echo + BEGIN LOOP OVER START DATES 
for START_YEAR in $START_YEARS
do
for START_MONTH in $START_MONTHS
do
for START_DAY in $START_DAYS
do

  echo +++ SET CONTINUE_RUN, STOP_OPTION AND STOP_N
  ENSEMBLE_PREFIX=${CASE_PREFIX}_${START_YEAR}${START_MONTH}${START_DAY} 
  : ${MEMBERS:=`seq -w 01 $ENSSIZE`}
  for MEMBER in $MEMBERS 
  do 
    CASE=${ENSEMBLE_PREFIX}_${MEMBERTAG}$MEMBER
    cd $CASESROOT/${ENSEMBLE_PREFIX}/${CASE}
    ./xmlchange -file env_run.xml -id STOP_OPTION -val $STOP_OPTION 
    ./xmlchange -file env_run.xml -id STOP_N -val $STOP_N 
    ./xmlchange -file env_run.xml -id CONTINUE_RUN -val $CONTINUE_RUN 
    sed -i "s/SBATCH --time=.*/SBATCH --time=${WALLTIME}/" ${CASE}.${MACH}.run
    ./${CASE}.${MACH}.submit
  done

done ; done ; done 
echo + END LOOP OVER START DATES 
