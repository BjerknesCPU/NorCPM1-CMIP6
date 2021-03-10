#!/bin/bash
set -e

sn=`basename $0`

#-----------------------------------------------------------------------
# Create remote archive root directory.
#-----------------------------------------------------------------------

ssh $CCSMUSER@$DOUT_L_MSHOST "mkdir -p $DOUT_L_MSROOT"

#-----------------------------------------------------------------------
# Archive source, object, library, and executable files related to the
# latest build.
#-----------------------------------------------------------------------

cd $EXEROOT

if ! ls *.ccsm.exe.* > /dev/null 2>&1; then
  echo "$sn: No build to archive."
else
  lastexe=`ls -1 *.ccsm.exe.* | tail -1`
  lid=${lastexe##*.}
  tarfile=build.$lid.tar
  if ssh $CCSMUSER@$DOUT_L_MSHOST "ls $DOUT_L_MSROOT/$tarfile.gz" > /dev/null 2>&1; then
    echo "$sn: Build $lid already archived."
  else
    tar cf $tarfile atm ccsm cpl csm_share glc ice lib lnd mct ocn pio $lastexe
    gzip -f $tarfile
    scp -p -oNoneSwitch=yes -oNoneEnabled=yes $tarfile.gz $CCSMUSER@$DOUT_L_MSHOST:$DOUT_L_MSROOT
    remote_cksum=`ssh $CCSMUSER@$DOUT_L_MSHOST "cksum $DOUT_L_MSROOT/$tarfile.gz"`
    local_cksum=`cksum $tarfile.gz`
    if [[ ${remote_cksum% *} == ${local_cksum% *} ]]; then
      echo "$sn: Build $lid archived."
      rm -f $tarfile.gz
    else
      echo "$sn: Archiving of build $lid failed!"
      exit 1
    fi
  fi
fi

#-----------------------------------------------------------------------
# Archive model output
#-----------------------------------------------------------------------

echo "$sn: Archive model output and restart files..."
rsync -a --progress --rsh='ssh -oNoneSwitch=yes -oNoneEnabled=yes' $DOUT_S_ROOT/ $CCSMUSER@$DOUT_L_MSHOST:$DOUT_L_MSROOT

echo "$sn: Archiving completed."
exit 0
