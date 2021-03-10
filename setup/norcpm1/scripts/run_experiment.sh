#!/bin/sh -evx

. $SETUPROOT/scripts/source_settings.sh $*

echo + PREPARE ANALYSIS DIRECTORY IF NEEDED 
if [ $ENKF_VERSION ] 
then
  cd $ANALYSISROOT
  ln -sf $OCNGRIDFILE grid.nc 
  if [ $COMPENSATE_ICE_FRESHWATER -eq 1 ] 
  then
    touch compensate_ice_freshwater
  fi   
fi

echo + BEGIN LOOP OVER START DATES 
for START_YEAR in $START_YEARS
do
for START_MONTH in $START_MONTHS
do
for START_DAY in $START_DAYS
do

  echo ++ BEGIN RESTART LOOP 
  for RESTART_COUNT in `seq 0 $RESTART`
  do 

    echo +++ CHECK IF SIMULATION SHOULD BE CONTINUED
    ENSEMBLE_PREFIX=${CASE_PREFIX}_${START_YEAR}${START_MONTH}${START_DAY}
    CASE1=${ENSEMBLE_PREFIX}_${MEMBERTAG}01
    if [ $CASE1 == `head -1 $WORK/noresm/$ENSEMBLE_PREFIX/$CASE1/run/rpointer.atm | cut -d. -f1` ]
    then
      CONTINUE_RUN=TRUE
    else
      CONTINUE_RUN=FALSE
    fi
    echo ++++ CONTINUE_RUN set to $CONTINUE_RUN


    echo +++ SET CONTINUE_RUN, STOP_OPTION AND STOP_N
    for MEMBER in `seq -w 01 $ENSSIZE`
    do 
      CASE=${ENSEMBLE_PREFIX}_${MEMBERTAG}$MEMBER
      cd $CASESROOT/${ENSEMBLE_PREFIX}/${CASE}
      ./xmlchange -file env_run.xml -id STOP_OPTION -val $STOP_OPTION 
      ./xmlchange -file env_run.xml -id STOP_N -val $STOP_N 
      ./xmlchange -file env_run.xml -id CONTINUE_RUN -val $CONTINUE_RUN 
      cd $WORK/noresm/${ENSEMBLE_PREFIX}/${CASE}/run/
      sed -i "s/stop_option    =.*/stop_option    ='${STOP_OPTION}'/" drv_in 
      sed -i "s/restart_option =.*/restart_option ='${STOP_OPTION}'/" drv_in 
      sed -i "s/stop_n         =.*/stop_n         =${STOP_N}/" drv_in 
      sed -i "s/restart_n      =.*/restart_n      =${STOP_N}/" drv_in 
      if [ $CONTINUE_RUN == "TRUE" ] 
      then 
        sed -i "s/start_type    =.*/start_type    = 'continue'/" drv_in 
      else
      if [[ $RUN_TYPE && $RUN_TYPE == "hybrid" ]] 
        then 
          sed -i "s/start_type    =.*/start_type    = 'startup'/" drv_in 
        else 
          sed -i "s/start_type    =.*/start_type    = 'branch'/" drv_in 
        fi
      fi 
    done
    yr=`head -1 rpointer.atm | cut -d. -f4 | cut -c1-4`
    mm=`head -1 rpointer.atm | cut -d. -f4 | cut -c6-7`
    y=$((10#$yr)) # strip leading zeros
    m=$((10#$mm)) # strip leading zeros
    echo ++++ year $yr month $mm

    if [[ $yr -eq $START_YEAR1 && $mm -eq $START_MONTH1 ]]
    then 
      SKIPASSIM=1
    fi 
    if [[ $ENKF_VERSION && $SKIPASSIM -eq 0 ]]
    then
      echo +++ PERFORM ASSIMILATION UPDATE `date`
      cd $ANALYSISROOT

      #every second month we reduced rfactor by 1
      nmonth=$(( (y - $START_YEAR1) * 12 + m - $START_MONTH1 )) 
      RFACTOR=$(( $RFACTOR_START - $nmonth / 2 ))
      if [ $RFACTOR -lt 1 ] 
      then 
        RFACTOR=1
      fi  
      echo ++++ SET RFACTOR TO $RFACTOR

      echo ++++ LINK FORECASTS
      for MEMBER in `seq -w 01 $ENSSIZE`
      do
        CASE=${ENSEMBLE_PREFIX}_${MEMBERTAG}$MEMBER
        ln -sf $WORK/noresm/${ENSEMBLE_PREFIX}/${CASE}/run/${CASE}.micom.r.${yr}-${mm}-15-00000.nc forecast0${MEMBER}.nc
        ln -sf $WORK/noresm/${ENSEMBLE_PREFIX}/${CASE}/run/${CASE}.cice.r.${yr}-${mm}-15-00000.nc forecast_ice0${MEMBER}.nc
        ncks -O -v aicen forecast_ice0${MEMBER}.nc aiceold0${MEMBER}.nc
        ncks -O -v vicen forecast_ice0${MEMBER}.nc viceold0${MEMBER}.nc
      done

      ENKF_CNT=0                           # Counter of EnKF sequential call
      echo ++++ PREPARE OBSERVATIONS AND DO SEQUENTIAL/CONCURRENT ASSIMILATION
      for iobs in ${!OBSLIST[*]}
      do
        OBSTYPE=${OBSLIST[$iobs]}
        PRODUCER=${PRODUCERLIST[$iobs]}
        REF_PERIOD=${REF_PERIODLIST[$iobs]}
        COMB_ASSIM=${COMBINE_ASSIM[$iobs]}    #sequential/joint observation assim 
        ln -sf ${INPUTDATA}/obs/${OBSTYPE}/${PRODUCER}/${yr}_${mm}.nc . || { echo "${INPUTDATA}/obs/${OBSTYPE}/${PRODUCER}/${yr}_${mm}.nc missing, we quit" ; exit 1 ; }
        ln -sf ${INPUTDATA}/enkf/${RES}/CMIP6/Free-average${mm}-${REF_PERIOD}.nc mean_mod.nc || { echo "Error ${INPUTDATA}/enkf/${RES}/Free-average${mm}-${REF_PERIOD}.nc missing, we quit" ; exit 1 ; }
        if [ -f ${INPUTDATA}/enkf/${RES}/${PRODUCER}/${RES}_${OBSTYPE}_obs_unc_anom.nc ]
        then
          ln -sf ${INPUTDATA}/enkf/${RES}/${PRODUCER}/${RES}_${OBSTYPE}_obs_unc_anom.nc  obs_unc_${OBSTYPE}.nc
        fi
        ln -sf ${INPUTDATA}/obs/${OBSTYPE}/${PRODUCER}/${OBSTYPE}_avg_${mm}-${REF_PERIOD}.nc mean_obs.nc || { echo "Error ${INPUTDATA}/obs/${OBSTYPE}/${PRODUCER}/${OBSTYPE}_avg_${mm}-${REF_PERIOD}.nc missing, we quit" ; exit 1 ; }
        cat ${INPUTDATA}/enkf/infile.data.${OBSTYPE}.${PRODUCER} | sed  "s/yyyy/${yr}/" | sed  "s/mm/${mm}/" > infile.data
        time mpirun -n 1 ./prep_obs
        mv observations.uf observations.uf_${OBSTYPE}.${PRODUCER}
        if (( $COMB_ASSIM ))
        then
          let ENKF_CNT=ENKF_CNT+1
          cat observations.uf_* > observations.uf
          rm -f observations.uf_*
          cp -f ${INPUTDATA}/enkf/analysisfields_V${ENKF_VERSION}_${ENKF_CNT}.in analysisfields.in
          cat ${INPUTDATA}/enkf/enkf.prm_V${ENKF_VERSION}_${ENKF_CNT} | sed  "s/XXX/${RFACTOR}/" > enkf.prm
          sed -i s/"enssize =".*/"enssize = "${ENSSIZE}/g enkf.prm

          if (( $ENSAVE ))
          then 
            echo +++++ COMPUTE PRE-ASSIMILATION ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
            time ./ensave forecast $ENSSIZE 
            mv forecast_avg.nc forecast_avg_${ENKF_CNT}.nc
            if [ $ENKF_VERSION -eq 2 ]
            then 
              time ./ensave_ice forecast_ice $ENSSIZE 
              mv forecast_ice_avg.nc forecast_ice_avg_${ENKF_CNT}.nc
            fi 
          fi

          echo +++++ CALL ENKF
          time mpirun -n $ENKF_NTASKS ./EnKF enkf.prm
          mv enkf_diag.nc enkf_diag_${ENKF_CNT}.nc      
          mv tmpX5.uf tmpX5_${ENKF_CNT}.uf

          if (( $ENSAVE ))
          then
            echo +++++ COMPUTE POST-ASSIMILATION ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
            time mpirun -n $ENSSIZE ./ensave forecast $ENSSIZE 
            mv forecast_avg.nc analysis_avg_${ENKF_CNT}.nc
            if [ $ENKF_VERSION -eq 2 ]
            then
              time ./ensave_ice forecast_ice $ENSSIZE 
              mv forecast_ice_avg.nc analysis_ice_avg_${ENKF_CNT}.nc
            fi
          fi

          echo 'Finished with EnKF; call number :' $ENKF_CNT
          date
        fi
      done  #OBS list

      echo ++++ ARCHIVE ASSIMILATION FILES
      mkdir -p $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm}  || { echo "Could not create $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm} dir" ; exit 1 ; }
      mv enkf_diag_*.nc observations-*.nc tmpX5_*.uf $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm}
      if (( $ENSAVE ))
      then
        mv analysis_*avg_*.nc forecast_*avg_*.nc $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm}
      fi 

      if [ $ENKF_VERSION -eq 2 ]
      then
        echo ++++ FIXENKF_ICE - POST-ASSIMILATION CORRECTION OF SEA ICE STATE
        for MEMBER in `seq -w 01 $ENSSIZE`
        do
          time ./fixenkf_cice $MEMBER & # process members in parallel on first node
        done
        wait # wait until all members are finished
      fi 

      echo ++++ MICOM_INIT - POST-ASSIMILATION CORRECTION OF OCEAN STATE    
      time mpirun -n $(( $MICOM_INIT_NTASKS_PER_MEMBER * $ENSSIZE )) ./micom_init $ENSSIZE

      if (( $ENSAVE ))
      then
	echo ++++ COMPUTE FINAL ENSEMBLE MEANS FOR OCEAN AND SEA ICE 
	time mpirun -n $ENSSIZE ./ensave forecast $ENSSIZE
	mv forecast_avg.nc $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm}/fix_analysis_avg.nc
	if [ $ENKF_VERSION -eq 2 ]
	then
	  time ./ensave_ice forecast_ice $ENSSIZE
	  mv forecast_ice_avg.nc $WORK/noresm/${ENSEMBLE_PREFIX}/RESULT/${yr}_${mm}/fix_analysis_ice_avg.nc
	fi
      fi

      echo ++++ FINISH ASSIM POST-PROCESSING `date`
      rm -f forecast???.nc forecast_ice???.nc aiceold???.nc viceold???.nc
      rm -f observations.uf enkf.prm* infile.data*

      echo +++ FINISHED ASSIMILATION UPDATE
    fi 
    SKIPASSIM=0

    echo +++ LAUNCH FIRST MEMBER - WILL RUN THE ENTIRE ENSEMBLE
    CASE=${ENSEMBLE_PREFIX}_${MEMBERTAG}01 
    cd $CASESROOT/${ENSEMBLE_PREFIX}/${CASE} 
    ./${CASE}.${MACH}.run

    echo +++ SHORT TERM ARCHIVING OF REMAINING MEMBERS 
    N_PARALLEL_STARCHIVE=0 
    for MEMBER in `seq -w 02 $ENSSIZE`
    do 
      export MACH
      export CASE=${CASE_PREFIX}_${START_YEAR}${START_MONTH}${START_DAY}_${MEMBERTAG}$MEMBER
      export RUNDIR=$WORK/noresm/${ENSEMBLE_PREFIX}/${CASE}/run
      export DOUT_S_ROOT=$WORK/archive/${ENSEMBLE_PREFIX}/${CASE}
      cd ${RUNDIR}
      $CASESROOT/${ENSEMBLE_PREFIX}/${CASE}/Tools/st_archive.sh &
      N_PARALLEL_STARCHIVE=`expr ${N_PARALLEL_STARCHIVE} + 1` 
      if [ ${N_PARALLEL_STARCHIVE} -eq $MAX_PARALLEL_STARCHIVE ] 
      then 
        N_PARALLEL_STARCHIVE=0
        wait 
      fi 
    done 
    wait 

    CONTINUE_RUN=TRUE 
  done  
  echo ++ END RESTART LOOP 

done ; done ; done 
echo + END LOOP OVER START DATES 
