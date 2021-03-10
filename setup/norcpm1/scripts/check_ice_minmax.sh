#!/bin/sh -e

SETUPROOT=`readlink -f \`dirname $0\``
. $SETUPROOT/scripts/source_settings.sh $*

echo CHECK VALUE RANGES IN RESTARTS OF SEA ICE AREA AND VOLUME

cd $ANALYSISROOT
rm -f aice_minmax.txt aicen_minmax.txt vice_minmax.txt vicen_minmax.txt
for MEM in `seq -w 01 $ENSSIZE`
do
  CASE=${CASE_PREFIX}_${START_YEAR1}${START_MONTH1}${START_DAY1}_${MEMBERTAG}${MEM}
  REST=`head -1 ../$CASE/run/rpointer.ice`

  ncwa -O -N -a ncat -v aicen ../$CASE/run/$REST aice.nc 
  ncwa -O -y min aice.nc aicemin.nc
  AICEMIN=`ncdump aicemin.nc | grep "aicen =" | awk '{print $3}'`
  ncwa -O -y max aice.nc aicemax.nc
  AICEMAX=`ncdump aicemax.nc | grep "aicen =" | awk '{print $3}'`
  echo $MEM $AICEMIN $AICEMAX >> aice_minmax.txt

  ncwa -O -y min -v aicen ../$CASE/run/$REST aicenmin.nc
  AICENMIN=`ncdump aicenmin.nc | grep "aicen =" | awk '{print $3}'`
  ncwa -O -y max -v aicen ../$CASE/run/$REST aicenmax.nc
  AICENMAX=`ncdump aicenmax.nc | grep "aicen =" | awk '{print $3}'`
  echo $MEM $AICENMIN $AICENMAX >> aicen_minmax.txt

  ncwa -O -N -a ncat -v vicen ../$CASE/run/$REST vice.nc 
  ncwa -O -y min vice.nc vicemin.nc
  VICEMIN=`ncdump vicemin.nc | grep "vicen =" | awk '{print $3}'`
  ncwa -O -y max vice.nc vicemax.nc
  VICEMAX=`ncdump vicemax.nc | grep "vicen =" | awk '{print $3}'`
  echo $MEM $VICEMIN $VICEMAX >> vice_minmax.txt

  ncwa -O -y min -v vicen ../$CASE/run/$REST vicenmin.nc
  VICENMIN=`ncdump vicenmin.nc | grep "vicen =" | awk '{print $3}'`
  ncwa -O -y max -v vicen ../$CASE/run/$REST vicenmax.nc
  VICENMAX=`ncdump vicenmax.nc | grep "vicen =" | awk '{print $3}'`
  echo $MEM $VICENMIN $VICENMAX >> vicen_minmax.txt
done 
rm -f aice.nc aicemin.nc aicemax.nc aicenmin.nc aicenmax.nc
rm -f vice.nc vicemin.nc vicemax.nc vicenmin.nc vicenmax.nc

echo 
echo RANGE OF ICE AREA SUMMED OVER CATAGORIES
cat aice_minmax.txt

echo 
echo RANGE OF ICE AREA FROM INDIVIDUAL CATAGORIES
cat aicen_minmax.txt

echo 
echo RANGE OF ICE VOLUME SUMMED OVER CATAGORIES
cat vice_minmax.txt

echo 
echo RANGE OF ICE VOLUME FROM INDIVIDUAL CATAGORIES
cat vicen_minmax.txt
