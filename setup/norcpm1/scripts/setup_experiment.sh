#!/bin/sh -e

SETUPROOT=`readlink -f \`dirname $0\`` 
. $SETUPROOT/scripts/source_settings.sh $*

echo + BEGIN LOOP OVER START DATES AND MEMBERS 
for START_YEAR in $START_YEARS
do
for START_MONTH in $START_MONTHS
do
for START_DAY in $START_DAYS
do
for MEMBER in `seq -w 01 $ENSSIZE`
do

  echo ++ PICK REFERENCE DATE
  COUNT=0
  for YEAR in $REF_YEARS
  do
  for MONTH in $REF_MONTHS
  do
  for DAY in $REF_DAYS
  do
    COUNT=`expr $COUNT + 1`
    if [ $COUNT -eq 1 ] 
    then 
      REF_YEAR=$YEAR
      REF_MONTH=$MONTH
      REF_DAY=$DAY
      REF_YEAR1=$YEAR
      REF_MONTH1=$MONTH
      REF_DAY1=$DAY
    fi  
    if [ $COUNT -eq $MEMBER ] 
    then 
      REF_YEAR=$YEAR
      REF_MONTH=$MONTH
      REF_DAY=$DAY
    fi
  done ; done ; done  

  ENSEMBLE_PREFIX=${CASE_PREFIX}_${START_YEAR}${START_MONTH}${START_DAY}
  CASE=${ENSEMBLE_PREFIX}_${MEMBERTAG}${MEMBER}
  ENSEMBLE_PREFIX1=${CASE_PREFIX}_${START_YEAR1}${START_MONTH1}${START_DAY1}
  CASE1=${ENSEMBLE_PREFIX1}_${MEMBERTAG}01
  if [[ $SKIP_CASE1 && $CASE == $CASE1 ]]
  then
    echo ++ SKIP CASE $CASE
    continue 
  else 
    echo ++ PREPARE CASE $CASE 
  fi 
 
  echo ++ REMOVE OLD CASE IF NEEDED
  CASEROOT=$CASESROOT/$ENSEMBLE_PREFIX/$CASE
  EXEROOT=$WORK/noresm/$ENSEMBLE_PREFIX/$CASE
  DOUT_S_ROOT=$WORK/archive/$ENSEMBLE_PREFIX/$CASE
  for ITEM in $CASEROOT $EXEROOT $DOUT_S_ROOT
  do
    if [ -e $ITEM ]
    then
      if [ $ASK_BEFORE_REMOVE -eq 1 ] 
      then 
        echo "remove existing $ITEM? (y/n)"
        if [ `read line ; echo $line` == "y" ]
        then
          rm -rf $ITEM
        fi
      else
        rm -rf $ITEM
      fi
    fi
  done

  if [ $CASE == $CASE1 ] 
  then 
    echo +++ PREPARE MEMBER 1 FROM SCRATCH

    echo +++ CREATE CASES DIRECTORY 
    mkdir -p $CASESROOT/$ENSEMBLE_PREFIX $WORK/noresm/$ENSEMBLE_PREFIX $WORK/archive/$ENSEMBLE_PREFIX

    echo +++ CREATE MEMBER 1 CASE
    $SCRIPTSROOT/create_newcase -case $CASEROOT -compset $COMPSET -res $RES -mach $MACH -pecount $PECOUNT

    echo +++ SET INITIALISATION 
    cd $CASEROOT
    ./xmlchange -file env_build.xml -id EXEROOT -val $EXEROOT
    ./xmlchange -file env_run.xml -id DOUT_S_ROOT -val $DOUT_S_ROOT
    if [[ $RUN_TYPE && $RUN_TYPE == "hybrid" ]]
    then
      REST_CASE=${REST_PREFIX}
      if (( $REF_ENSEMBLE )) 
      then 
        REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${REF_YEAR}-${REF_MONTH}-${REF_DAY}_${MEMBERTAG}${MEMBER}
      else
        REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${REF_YEAR}-${REF_MONTH}-${REF_DAY}
      fi
      if [ ! -e $REST_PATH ]
      then
        echo cannot locate restart data in $REST_PATH  
        exit
      fi
      ./xmlchange -file env_conf.xml -id RUN_TYPE -val hybrid
      ./xmlchange -file env_conf.xml -id RUN_REFDATE -val ${REF_YEAR}-${REF_MONTH}-${REF_DAY}
      ./xmlchange -file env_conf.xml -id RUN_STARTDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
    else
      REST_CASE=${REST_PREFIX}${MEMBER}
      REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${START_YEAR}-${START_MONTH}-${START_DAY}
      if [ ! -e $REST_PATH ]
      then
        echo cannot locate restart data in $REST_PATH  
        exit
      fi
      ./xmlchange -file env_conf.xml -id BRNCH_RETAIN_CASENAME -val TRUE
      ./xmlchange -file env_conf.xml -id RUN_TYPE -val branch
      ./xmlchange -file env_conf.xml -id RUN_REFDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
      ./xmlchange -file env_conf.xml -id RUN_STARTDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
    fi
    ./xmlchange -file env_conf.xml -id RUN_REFCASE -val ${REST_CASE}
    ./xmlchange -file env_conf.xml -id GET_REFCASE -val FALSE

    echo +++ CONFIGURE MEMBER 1 CASE 
    ./configure -case

    echo +++ DEACTIVATE MICOM RESTART COMPRESSION
    sed -i s/" RSTCMP   =".*/" RSTCMP   = 0"/ Buildconf/micom.buildnml.csh

    echo +++ BUILD MEMBER 1 CASE 
    set +e 
    ./${CASE}.${MACH}.build
    set -e 

  else
    echo +++ CLONE MEMBER 1 
    $SCRIPTSROOT/create_clone -clone $CASESROOT/$ENSEMBLE_PREFIX1/$CASE1 -case $CASEROOT

    echo +++ LINK BUILD OBJECTS AND EXECUTABLE FROM MEMBER 1
    mkdir -p $EXEROOT/run
    cd $EXEROOT
    for ITEM in atm cpl ccsm csm_share glc ice ocn pio lib
    do
      ln -s  $WORK/noresm/$ENSEMBLE_PREFIX1/$CASE1/$ITEM .
    done
    cd run
    ln -s $WORK/noresm/$ENSEMBLE_PREFIX1/$CASE1/run/ccsm.exe . 

    echo +++ CONFIGURE CASE 
    cd $CASEROOT
    ./xmlchange -file env_build.xml -id EXEROOT -val $EXEROOT
    ./xmlchange -file env_run.xml -id DOUT_S_ROOT -val $DOUT_S_ROOT
    ./xmlchange -file env_conf.xml -id RUN_REFCASE -val ${REST_PREFIX}${MEMBER}
    ./xmlchange -file env_conf.xml -id GET_REFCASE -val FALSE
    if [[ $RUN_TYPE && $RUN_TYPE == "hybrid" ]]
    then
      ./xmlchange -file env_conf.xml -id RUN_TYPE -val hybrid
      ./xmlchange -file env_conf.xml -id RUN_REFDATE -val ${REF_YEAR}-${REF_MONTH}-${REF_DAY}
      ./xmlchange -file env_conf.xml -id RUN_STARTDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
    else
      ./xmlchange -file env_conf.xml -id BRNCH_RETAIN_CASENAME -val TRUE
      ./xmlchange -file env_conf.xml -id RUN_TYPE -val branch
      ./xmlchange -file env_conf.xml -id RUN_REFDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
      ./xmlchange -file env_conf.xml -id RUN_STARTDATE -val ${START_YEAR}-${START_MONTH}-${START_DAY}
    fi 
    ./configure -case

    echo +++ MODIFY RESTART PATH IN NAMELISTS
    if [[ $RUN_TYPE && $RUN_TYPE == "hybrid" ]]
    then
      REST_CASE=${REST_PREFIX}
      if (( $REF_ENSEMBLE ))
      then
        REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${REF_YEAR}-${REF_MONTH}-${REF_DAY}_${MEMBERTAG}${MEMBER}
      else
        REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${REF_YEAR}-${REF_MONTH}-${REF_DAY}
      fi
      if [ ! -e $REST_PATH ]
      then
        echo cannot locate restart data in $REST_PATH  
        exit
      fi
      sed -i "s%${REST_PREFIX}01.cice.r.${REF_YEAR1}-${REF_MONTH1}-${REF_DAY1}%${REST_PREFIX}${REF_MEMBER}.cice.r.${REF_YEAR}-${REF_MONTH}-${REF_DAY}%" Buildconf/cice.buildnml.csh
      sed -i "s%${REST_PREFIX}01.clm2.r.${REF_YEAR1}-${REF_MONTH1}-${REF_DAY1}%${REST_PREFIX}${REF_MEMBER}.clm2.r.${REF_YEAR}-${REF_MONTH}-${REF_DAY}%" Buildconf/clm.buildnml.csh
      #Fanf: typically ifile does not have the same date than restart file
      ifile=$(basename `find $REST_PATH/ -name '*cam2.i*'`)
      sed -i s/"ncdata".*/"ncdata = ${ifile} "/g Buildconf/cam.input_data_list
      sed -i s/"ncdata".*/"ncdata = '${ifile}'"/g Buildconf/cam.buildnml.csh
      sed -i s/"ncdata".*/"ncdata = '${ifile}'"/g Buildconf/camconf/ccsm_namelist
    else
      REST_CASE=${REST_PREFIX}${MEMBER}
      REST_PATH=$REST_PATH_LOCAL/${REST_CASE}/${START_YEAR}-${START_MONTH}-${START_DAY}
      if [ ! -e $REST_PATH ]
      then
        echo cannot locate restart data in $REST_PATH  
        exit
      fi
      sed -i "s%${ENSEMBLE_PREFIX1}/${CASE1}/run/${REST_PREFIX}01.cam2.r.${START_YEAR1}-${START_MONTH1}-${START_DAY1}%${ENSEMBLE_PREFIX}/${CASE}/run/${REST_PREFIX}${MEMBER}.cam2.r.${START_YEAR}-${START_MONTH}-${START_DAY}%" Buildconf/cam.buildnml.csh 
      sed -i "s%${REST_PREFIX}01.cice.r.${START_YEAR1}-${START_MONTH1}-${START_DAY1}%${REST_PREFIX}${MEMBER}.cice.r.${START_YEAR}-${START_MONTH}-${START_DAY}%" Buildconf/cice.buildnml.csh
      sed -i "s%${REST_PREFIX}01.clm2.r.${START_YEAR1}-${START_MONTH1}-${START_DAY1}%${REST_PREFIX}${MEMBER}.clm2.r.${START_YEAR}-${START_MONTH}-${START_DAY}%" Buildconf/clm.buildnml.csh
    fi 
    sed -i "s/start_ymd      =.*/start_ymd      = ${START_YEAR}${START_MONTH}${START_DAY}/" Buildconf/cpl.buildnml.csh 

    echo +++ DUMMY BUILD
    ./xmlchange -file env_build.xml -id BUILD_COMPLETE -val TRUE 
    sed -i '/source $CASETOOLS\/ccsm_buildexe/d' ${CASE}.${MACH}.build
    ./${CASE}.${MACH}.build 

  fi 
  echo ++ FINISHED PREPARING CASE

  echo ++ STAGE RESTART DATA 
  cd $EXEROOT/run 
  ln -sf $REST_PATH/*nc . 
  cp -f $REST_PATH/rpointer* .  

  echo ++ CREATE SUBDIRS FOR TIMING AND SHORT TERM ARCHIVING
  mkdir -p timing/checkpoints $WORK/archive/$ENSEMBLE_PREFIX

done ; done ; done ; done 
echo + END LOOP OVER START DATES AND MEMBERS 

echo + BUILD ASSIMILATION CODE IF NEEDED
cd $SETUPROOT
source scripts/build_assimcode.sh $*
