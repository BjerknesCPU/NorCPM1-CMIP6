#!/bin/sh -e

if [[ ! $1 || $1 == "-h" || $1 == "--help"  ]]
then 
  echo "USAGE: ./`basename $0` <path to settings file> [var1=value var2=value]" 
  echo 
  echo "EXAMPLE: ./`basename $0` settings/predictiontest.sh MACH=betzy" 
  exit 
fi 

# OVERRIDE DEFAULTS
for ARG in $* 
do
 if [ `echo $ARG | grep =` ] 
 then
   declare `echo $ARG | cut -d= -f1`=`echo $ARG | cut -d= -f2`
 fi
done 

# SET MACHINE SPECIFIC SETTINGS
. $SETUPROOT/settings/setmach.sh

# READING SETTINGS FROM FILE 
if [ ! $SETUPROOT/settings/`basename $1` ]
then
  echo cannot read settings file $SETUPROOT/settings/`basename $1` 
  exit
fi
. $SETUPROOT/settings/`basename $1`

# SET VERBOSITY
if (( $VERBOSE ))
then
  set -vx
  echo set logging verbose
fi
