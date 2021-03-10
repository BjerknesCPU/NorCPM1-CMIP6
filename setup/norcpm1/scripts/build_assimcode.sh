#!/bin/sh -e

SETUPROOT=`readlink -f \`dirname $0\``
. $SETUPROOT/scripts/source_settings.sh $*

if [[ $ENKF_VERSION && $ENKF_CODE ]]
then
  echo + build EnKF
  mkdir -p $ANALYSISROOT/bld/EnKF/TMP
  cd $ANALYSISROOT/bld/EnKF
  cp -f $ENKFROOT/shared/* . 
  cp -f $ENKF_CODE/* . 
  make clean
  make  
fi

if [[ $ENKF_VERSION && $PREP_OBS_CODE ]]
then
  echo + build prep_obs
  mkdir -p $ANALYSISROOT/bld/prep_obs/TMP
  cd $ANALYSISROOT/bld/prep_obs
  cp -f $ENKFROOT/shared/* . 
  cp -f $PREP_OBS_CODE/* . 
  make clean
  make
fi

if [[ $ENKF_VERSION && $ENSAVE_FIXENKF_CODE ]]
then
  echo + build ensave and fixenkf
  mkdir -p $ANALYSISROOT/bld/ensave_fixenkf/TMP
  cd $ANALYSISROOT/bld/ensave_fixenkf
  cp -f $ENKFROOT/shared/* . 
  cp -f $ENSAVE_FIXENKF_CODE/* . 
  make clean
  make
fi

if [[ $ENKF_VERSION && $MICOM_INIT_CODE ]]
then
  echo + build micom_init
  mkdir -p $ANALYSISROOT/bld/micom_init
  cd $ANALYSISROOT/bld/micom_init
  cp -f $MICOM_INIT_CODE/* . 
  make clean
  make
fi
