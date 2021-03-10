#! /bin/csh -f

# -------------------------------------------------------------------------
# Check for successful run
# -------------------------------------------------------------------------

set sdate = `date +"%Y-%m-%d %H:%M:%S"`

cd $RUNDIR
set CplLogFile = `ls -1t cpl.log* | head -1` 
if ($CplLogFile == "") then
  echo "Model did not complete - no cpl.log file present - exiting"
  exit -1
endif
grep 'SUCCESSFUL TERMINATION' $CplLogFile  || echo "Model did not complete - see $RUNDIR/$CplLogFile" && echo "run FAILED $sdate" >>& $CASEROOT/CaseStatus && exit -1

echo "run SUCCESSFUL $sdate" >>& $CASEROOT/CaseStatus

# -------------------------------------------------------------------------
# Update env variables in case user changed them during run
# -------------------------------------------------------------------------

cd $CASEROOT
source ./Tools/ccsm_getenv

# -------------------------------------------------------------------------
# Save model output stdout and stderr 
# -------------------------------------------------------------------------

cd $EXEROOT
gzip */*.$LID
if ($LOGDIR != "") then
  if (! -d $LOGDIR/bld) mkdir -p $LOGDIR/bld || echo " problem in creating $LOGDIR/bld"
  cp -p */*build.$LID.* $LOGDIR/bld  
  cp -p */*log.$LID.*   $LOGDIR      
endif

# -------------------------------------------------------------------------
# Perform short term archiving of output
# -------------------------------------------------------------------------

if ($DOUT_S == 'TRUE') then
  echo "Archiving ccsm output to $DOUT_S_ROOT"
  echo "Calling the short-term archiving script st_archive.sh"
  cd $RUNDIR; $CASETOOLS/st_archive.sh
endif

# -------------------------------------------------------------------------
# Submit longer term archiver if appropriate
# -------------------------------------------------------------------------

cd $CASEROOT
if ($DOUT_L_MS == 'TRUE' && $DOUT_S == 'TRUE') then
  echo "Long term archiving ccsm output using the script $CASE.$MACH.l_archive"
  set num = 0
  if ($LBQUERY == "TRUE") then
     set num = `$BATCHQUERY | grep $CASE.l_archive | wc -l`
  endif
  if ($LBSUBMIT == "TRUE" && $num < 1) then
cat > templar <<EOF
   $BATCHSUBMIT ./$CASE.$MACH.l_archive
EOF
    source templar
    if ($status != 0) then
      echo "ccsm_postrun error: problem sourcing templar " 
    endif
    rm templar
  endif 
endif

# -------------------------------------------------------------------------
# Resubmit another run script or restart model
# -------------------------------------------------------------------------

cd $CASEROOT
if ($RESUBMITNOW == "TRUE") then
  if ($RESUBMIT > 0) then
    @ RESUBMIT = $RESUBMIT - 1
    echo RESUBMIT is now $RESUBMIT

    #tcraig: reset CONTINUE_RUN on RESUBMIT if NOT doing timing runs
    #use COMP_RUN_BARRIERS as surrogate for timing run logical
    if (${COMP_RUN_BARRIERS} == "FALSE") then
       ./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
    endif
    ./xmlchange -file env_run.xml -id RESUBMIT     -val $RESUBMIT

    if ($LBSUBMIT == "TRUE") then
cat > tempres <<EOF
   $BATCHSUBMIT ./$CASE.$MACH.run
EOF
      source tempres
      if ($status != 0) then
        echo "ccsm_postrun error: problem sourcing tempres " 
      endif
      rm tempres
    endif 
  endif
else
  if ($RESTART > 0) then
    ./xmlchange -file env_run.xml -id CONTINUE_RUN -val TRUE
    echo "Restarts model without resubmit"
  endif
endif

if ($CHECK_TIMING == 'TRUE') then
  cd $CASEROOT
  if !(-d timing) mkdir timing
  $CASETOOLS/perf_summary.pl $RUNDIR/timing/ccsm_timing >& timing/ccsm_timing_summary.$LID
  $CASETOOLS/getTiming.csh -lid $LID -mach $MACH
  gzip timing/ccsm_timing_summary.$LID
endif

if ($SAVE_TIMING == 'TRUE') then
  cd $RUNDIR
  mv timing timing.$LID
  cd $CASEROOT
endif


