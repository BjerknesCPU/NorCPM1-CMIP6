      SUBROUTINE AUFR_BGC(kpie,kpje,kpke,pddpo,kplyear,kplmon     &
     &                    ,kplday,kpldtoce,omask                  &
     &                    ,rstfnm_ocn,path,path_len)

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/aufr_bgc.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!****************************************************************
!
!**** *AUFR_BGC* - reads marine bgc restart data.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!     - extra SBR for reading bgc data from the restart file.
!     S.Legutke,        *MPI-MaD, HH*    15.08.01
!     - netCDF version (with cond.comp. PNETCDF)
!     - no use of chemc values from netCDF restart
!
!     Patrick Wetzel,    *MPI-Met, HH*    16.04.02
!     - read chemcm(i,j,7,12) from netCDF restart
!
!     Filippa Fransner & Ingo Bethke    *GFI, Bergen Uni*  19.04.15       
!     - addition of accumulated diagnostics
!     PROBLEM: omegaC at certain grid points
!
!     Purpose
!     -------
!     Read restart data to continue an interrupted integration.
!
!     Method
!     -------
!     The bgc data are read from an extra file, other than the ocean data.
!     The time stamp of the bgc restart file (idate) is specified from the
!     ocean time stamp through the SBR parameter list of AUFW_BGC. The only 
!     time control variable proper to the bgc is the time step number 
!     (idate(5)). It can differ from that of the ocean (idate(4)) by the 
!     difference of the offsets of restart files.
!
!**   Interface.
!     ----------
!
!     *CALL*       *AUFR_BGC(kpie,kpje,kpke,pddpo
!                            ,kplyear,kplmon,kplday,kpldtoce)*
!
!     *COMMON*     *MO_PARAM1_BGC* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of grid cell (3rd dimension) [m].
!     *INTEGER* *kplyear*   - year  in ocean restart date
!     *INTEGER* *kplmon*  - month in ocean restart date
!     *INTEGER* *kplday*    - day   in ocean restart date
!     *INTEGER* *kpldtoce*  - step  in ocean restart date
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
!ik IK introduced ocetra array elements for phyto, grazer, poc (=det), calciu 
!ik array indices are: iphy, izoo, idet, icalc
!iktodo IK introduced new variable opal (index iopal)
!ik nocetra is the number of all BGC element (size of ocetra(,,,l))
!ik nosedi is the number of all elements interacting with the sediment

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc
! diagnostic accumulation - start
      use mo_bgcmean, ncdims_nctools => ncdims
! diagnostic accumulation - end

      implicit none
      
      INTEGER  kpie,kpje,kpke
      INTEGER  kplyear,kplmon,kplday,kpldtoce
      REAL pddpo(kpie,kpje,kpke)
      REAL chemcm_t(kpie,kpje,8)
      REAL omask(kpie,kpje)
      INTEGER  i,j,k,l,kmon

      INTEGER idate(5)

      INTEGER :: restyear            !  year of restart file
      INTEGER :: restmonth           !  month of restart file
      INTEGER :: restday             !  day of restart file
      INTEGER :: restdtoce           !  time step number from bgc ocean file

      character rstfnm*80
      character*(*) rstfnm_ocn,path
      integer path_len
! diagnostic accumulation - start
      integer :: n
      real :: rtmp 
      character(len=2) :: c2
! diagnostic accumulation - end
#define PNETCDF 
#ifdef PNETCDF                
      INCLUDE 'netcdf.inc'
      INTEGER ncid,ncstat,ncvarid

!
! Open netCDF data file
!
      IF(mnproc==1) THEN

#ifdef CCSMCOUPLED
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'.micom.r.')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'.micom.rbgc.'//rstfnm_ocn(i+9:)
#else
        i=1
        do while (rstfnm_ocn(i:i+8).ne.'_restphy_')
          i=i+1
          if (i+8.gt.len(rstfnm_ocn)) then
            write (io_stdo_bgc,*)                                    &
     &        'Could not generate restart file name!'
            call xchalt('(aufr_bgc)')
            stop '(aufr_bgc)'
          endif
        enddo
        rstfnm=rstfnm_ocn(1:i-1)//'_rest_b_'//rstfnm_ocn(i+9:)
#endif

        ncstat = NF_OPEN(path(1:path_len)//rstfnm,NF_NOWRITE, ncid)
        IF ( ncstat .NE. NF_NOERR ) THEN
             CALL xchalt('(AUFR: Problem with netCDF1)')
                    stop '(AUFR: Problem with netCDF1)'
        ENDIF

!
! Read restart data : date
!

        ncstat = NF_GET_ATT_INT(ncid, NF_GLOBAL,'date', idate)
        IF ( ncstat .NE. NF_NOERR ) THEN
          CALL xchalt('(AUFR: Problem reading date of restart file)')
                 stop '(AUFR: Problem reading date of restart file)'
        ENDIF
        restyear  = idate(1)
        restmonth = idate(2)
        restday   = idate(3)
        restdtoce = idate(4)
        ldtbgc = idate(5)
        WRITE(io_stdo_bgc,*) ' '
        WRITE(io_stdo_bgc,*) 'Date of bgc restart file : '
        WRITE(io_stdo_bgc,*) ' year  = ',restyear
        WRITE(io_stdo_bgc,*) ' month = ',restmonth
        WRITE(io_stdo_bgc,*) ' day   = ',restday
        WRITE(io_stdo_bgc,*) ' dtoce = ',restdtoce
        WRITE(io_stdo_bgc,*) ' dtbgc = ',ldtbgc
        WRITE(io_stdo_bgc,*) ' '
      ENDIF

!
! Compare with date read from ocean restart file
!
! As the ocean is already in its first step, its counter has 
! gone up one step already for the year and month. The ocean day
! counter is still at its restart date. Therefore:

!      restmonth = restmonth + 1
!
!      IF (restmonth .GT. 12) THEN
!         restmonth=1
!         restyear=restyear+1
!      ENDIF

! Memorize ocean start date :   
      bgcstartyear  = kplyear
      bgcstartmonth = kplmon
      bgcstartday   = kplday

      IF ( kplyear  .NE. restyear  ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart years in oce/bgc are not the same : '  &
     &   ,kplyear,'/',restyear,' !!!'
         ENDIF
      ENDIF

      IF ( kplmon .NE. restmonth ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart months in oce/bgc are not the same : '   &
     &   ,kplmon,'/',restmonth,' !!!'
         ENDIF
!         STOP 'Stop : restart months in oce/bgc are not the same.'
      ENDIF

      IF ( kplday   .NE. restday   ) THEN
         IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)                                     &
     &   'WARNING: restart days in oce/bgc are not the same : '   &
     &   ,kplday,'/',restday,' !!!'
         ENDIF
!         STOP 'Stop : restart days in oce/bgc are not the same.'
      ENDIF 

!      IF ( kpldtoce .NE. ldtbgc   ) THEN
!         WRITE(io_stdo_bgc,*)                                       & 
!     &   'WARNING: restart step no.  in oce/bgc are not the same : '&
!     &   ,kpldtoce,'/',ldtbgc,' !!!'
!      ENDIF

!
! Read restart data : ocean aquateous tracer
!                
      CALL read_netcdf_var(ncid,'sco212',ocetra(1,1,1,isco212),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco212),kpke,'sco212')
! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sco213',ocetra(1,1,1,isco213),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco213),kpke,'sco213')
      CALL read_netcdf_var(ncid,'sco214',ocetra(1,1,1,isco214),kpke,0)
!      call chksumbgc(ocetra(1,1,1,isco214),kpke,'sco214')
#endif
      CALL read_netcdf_var(ncid,'alkali',ocetra(1,1,1,ialkali),kpke,0)
!      call chksumbgc(ocetra(1,1,1,ialkali),kpke,'alkali')
      CALL read_netcdf_var(ncid,'phosph',ocetra(1,1,1,iphosph),kpke,0)
      CALL read_netcdf_var(ncid,'oxygen',ocetra(1,1,1,ioxygen),kpke,0)
      CALL read_netcdf_var(ncid,'gasnit',ocetra(1,1,1,igasnit),kpke,0)
      CALL read_netcdf_var(ncid,'ano3',ocetra(1,1,1,iano3),kpke,0)
      CALL read_netcdf_var(ncid,'silica',ocetra(1,1,1,isilica),kpke,0)
      CALL read_netcdf_var(ncid,'doc',ocetra(1,1,1,idoc),kpke,0)
      CALL read_netcdf_var(ncid,'phyto',ocetra(1,1,1,iphy),kpke,0)
      CALL read_netcdf_var(ncid,'grazer',ocetra(1,1,1,izoo),kpke,0)
      CALL read_netcdf_var(ncid,'poc',ocetra(1,1,1,idet),kpke,0)
! comment out for first run (no 13c in restart file)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'poc13',ocetra(1,1,1,idet13),kpke,0)
      CALL read_netcdf_var(ncid,'poc14',ocetra(1,1,1,idet14),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'calciu',ocetra(1,1,1,icalc),kpke,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'calciu13',ocetra(1,1,1,icalc13),kpke,0)
      CALL read_netcdf_var(ncid,'calciu14',ocetra(1,1,1,icalc14),kpke,0)
#endif
      CALL read_netcdf_var(ncid,'opal',ocetra(1,1,1,iopal),kpke,0)
      CALL read_netcdf_var(ncid,'n2o',ocetra(1,1,1,ian2o),kpke,0)
      CALL read_netcdf_var(ncid,'dms',ocetra(1,1,1,idms),kpke,0)
      CALL read_netcdf_var(ncid,'fdust',ocetra(1,1,1,ifdust),kpke,0)
      CALL read_netcdf_var(ncid,'iron',ocetra(1,1,1,iiron),kpke,0)

#ifdef AGG
      CALL read_netcdf_var(ncid,'snos',ocetra(1,1,1,inos),kpke,0)
      CALL read_netcdf_var(ncid,'adust',ocetra(1,1,1,iadust),kpke,0)
#endif /*AGG*/

#ifdef ANTC14
      CALL read_netcdf_var(ncid,'antc14',ocetra(1,1,1,iantc14),kpke,0)
#endif
#ifdef PCFC
      CALL read_netcdf_var(ncid,'cfc11',ocetra(1,1,1,icfc11),kpke,0)
      CALL read_netcdf_var(ncid,'cfc12',ocetra(1,1,1,icfc12),kpke,0)
#endif

!
!Check aquateous restart data for topography
! 

!      DO i    =1,kpie
!      DO j    =1,kpje
!      DO k    =1,kpke      
!      DO l    =1,nocetra 
!	IF (omask(i,j) .le. 0.5 ) THEN
!	  IF ( ocetra(i,j,k,l) .NE. rmasko ) THEN
!             WRITE(io_stdo_bgc,*) 'ocetra not properly masked at :'  &
!     &                             ,i,j,k,l,ocetra(i,j,k,l)
!	  ENDIF
!	ELSE
!	  IF ( ocetra(i,j,k,l) .EQ. rmasko ) THEN
!             WRITE(io_stdo_bgc,*) 'land mask values at wet points :'  &
!     &                             ,i,j,k,l,ocetra(i,j,k,l)
!             call STOP_ALL('Stop : Restart file with different topography:     &
!     &	     land mask values at wet points')
!
!	  ENDIF	 
!	ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

!
! Read restart data : other fields
!
     CALL read_netcdf_var(ncid,'hi',hi(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'co3',co3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'akw3',akw3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'akb3',akb3(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'ak13',ak13(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'ak23',ak23(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'aksp',aksp(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satoxy',satoxy(1,1,1),kpke,0)
     CALL read_netcdf_var(ncid,'satn2o',satn2o(1,1),1,0)

!
! Read restart data : chemical constants
!
      DO kmon =1,12
      CALL read_netcdf_var(ncid,'chemcm',chemcm_t,8,kmon)
      DO i    =1,kpie
      DO j    =1,kpje
	 IF (omask(i,j) .lt. 0.5 ) THEN
	    DO k=1,8
	       chemcm(i,j,k,kmon)  =   rmasko
	    ENDDO
	 ELSE
	   chemcm(i,j,1,kmon) = chemcm_t(i,j,1)
	   chemcm(i,j,2,kmon) = chemcm_t(i,j,2)
	   chemcm(i,j,3,kmon) = chemcm_t(i,j,3)
	   chemcm(i,j,4,kmon) = chemcm_t(i,j,4)
	   chemcm(i,j,5,kmon) = chemcm_t(i,j,5)
	   chemcm(i,j,6,kmon) = chemcm_t(i,j,6)
	   chemcm(i,j,7,kmon) = chemcm_t(i,j,7)
	   chemcm(i,j,8,kmon) = chemcm_t(i,j,8)
	 ENDIF
      ENDDO
      ENDDO
      ENDDO
!      call chksumbgc(chemcm,8*12,'chemcm')

!
! Read restart data : sediment variables.
! js: reading of terrigenous sediment was missing until 02.08.2005
      CALL read_netcdf_var(ncid,'ssso12',sedlay(1,1,1,issso12),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'ssso13',sedlay(1,1,1,issso13),ks,0)
      CALL read_netcdf_var(ncid,'ssso14',sedlay(1,1,1,issso14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'sssc12',sedlay(1,1,1,isssc12),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'sssc13',sedlay(1,1,1,isssc13),ks,0)
      CALL read_netcdf_var(ncid,'sssc14',sedlay(1,1,1,isssc14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'ssssil',sedlay(1,1,1,issssil),ks,0)
      CALL read_netcdf_var(ncid,'ssster',sedlay(1,1,1,issster),ks,0)
      CALL read_netcdf_var(ncid,'bur_o12',burial(1,1,issso12),1,0)
      CALL read_netcdf_var(ncid,'bur_c12',burial(1,1,isssc12),1,0)
      CALL read_netcdf_var(ncid,'bur_sil',burial(1,1,issssil),1,0)
      CALL read_netcdf_var(ncid,'bur_clay',burial(1,1,issster),1,0)
      CALL read_netcdf_var(ncid,'sedhpl',sedhpl(1,1,1),ks,0)
      CALL read_netcdf_var(ncid,'powaic',powtra(1,1,1,ipowaic),ks,0)
!aufsetz!
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'powc13',powtra(1,1,1,ipowc13),ks,0)
      CALL read_netcdf_var(ncid,'powc14',powtra(1,1,1,ipowc14),ks,0)
#endif
      CALL read_netcdf_var(ncid,'powaal',powtra(1,1,1,ipowaal),ks,0)
      CALL read_netcdf_var(ncid,'powaph',powtra(1,1,1,ipowaph),ks,0)
      CALL read_netcdf_var(ncid,'powaox',powtra(1,1,1,ipowaox),ks,0)
      CALL read_netcdf_var(ncid,'pown2',powtra(1,1,1,ipown2),ks,0)
      CALL read_netcdf_var(ncid,'powno3',powtra(1,1,1,ipowno3),ks,0)
      CALL read_netcdf_var(ncid,'powasi',powtra(1,1,1,ipowasi),ks,0)

#ifdef DIFFAT 
!
! Read restart data : co2 diffusion
!
      CALL read_netcdf_var(ncid,'atmco2',atm(1,1,iatmco2),1,0)
      CALL read_netcdf_var(ncid,'atmo2',atm(1,1,iatmo2),1,0)
      CALL read_netcdf_var(ncid,'atmn2',atm(1,1,iatmn2),1,0)
#ifdef __c_isotopes
      CALL read_netcdf_var(ncid,'atmc13',atm(1,1,iatmc13),1,0)
      CALL read_netcdf_var(ncid,'atmc14',atm(1,1,iatmc14),1,0)
#endif

#endif

! read partly accumulated diagnostic fields from restart - start
      do n=1,nbgc
        write(c2,'(i2.2)') n
        rtmp=0. 
        if (mnproc.eq.1) then 
          if (ncinqa('nacc_bgc'//c2)) then 
            call ncgeti('nacc_bgc'//c2,nacc_bgc(n))
          else
            nacc_bgc(n)=0 
          endif
          rtmp=nacc_bgc(n)
        endif 
        call xcmaxr(rtmp) 
        nacc_bgc(n)=nint(rtmp)
        if (nacc_bgc(n).ne.0) then
!         if (jsrfphyc(n).ne.0) then 
!           if (ncinqv('srfphyc_acc'//c2)) call read_netcdf_var( & 
!    &       ncid,'srfphyc_acc'//c2,bgcm2d(1,1,jsrfphyc(n)),1,0)
!         endif
          if (jintphyc(n).ne.0) call read_netcdf_var(ncid,              &
     &      'intphyc_acc'//c2,bgcm2d(1:kpie,1:kpje,jintphyc(n)),1,0)
          if (jintphosy(n).ne.0) call read_netcdf_var(ncid,             &
     &      'intpp_acc'//c2,bgcm2d(1:kpie,1:kpje,jintphosy(n)),1,0) 
          if (jsrfphyc(n).ne.0) call read_netcdf_var(ncid,              &
     &      'srfphyc_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfphyc(n)),1,0)
          if (jsrfphosy(n).ne.0) call read_netcdf_var(ncid,             &
     &      'srfpp_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfphosy(n)),1,0)
          if (jsrfdic(n).ne.0) call read_netcdf_var(ncid,               &
     &      'srfdic_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfdic(n)),1,0)
          if (jsrfsilica(n).ne.0) call read_netcdf_var(ncid,            &
     &      'srfsi_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfsilica(n)),1,0)
          if (jsrfalkali(n).ne.0) call read_netcdf_var(ncid,            &
     &      'srfalk_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfalkali(n)),1,0)
          if (jsrfano3(n).ne.0) call read_netcdf_var(ncid,              &
     &      'srfno3_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfano3(n)),1,0)
          if (jsrfoxygen(n).ne.0) call read_netcdf_var(ncid,            &
     &      'srfoxy_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfoxygen(n)),1,0)
          if (jsrfphosph(n).ne.0) call read_netcdf_var(ncid,            &
     &      'srfpo4_acc'//c2,bgcm2d(1:kpie,1:kpje,jsrfphosph(n)),1,0)
          if (jexposi(n).ne.0) call read_netcdf_var(ncid,               &
     &      'exposi_acc'//c2,bgcm2d(1:kpie,1:kpje,jexposi(n)),1,0)
          if (jexpoca(n).ne.0) call read_netcdf_var(ncid,               &
     &      'expoca_acc'//c2,bgcm2d(1:kpie,1:kpje,jexpoca(n)),1,0)
          if (jexport(n).ne.0) call read_netcdf_var(ncid,               &
     &      'export_acc'//c2,bgcm2d(1:kpie,1:kpje,jexport(n)),1,0)
          if (jdms_uv(n).ne.0) call read_netcdf_var(ncid,               &
     &      'dmsuv_acc'//c2,bgcm2d(1:kpie,1:kpje,jdms_uv(n)),1,0)
          if (jdms_bac(n).ne.0) call read_netcdf_var(ncid,              &
     &      'dmsbac_acc'//c2,bgcm2d(1:kpie,1:kpje,jdms_bac(n)),1,0)
          if (jdmsprod(n).ne.0) call read_netcdf_var(ncid,              &
     &      'dmsprod_acc'//c2,bgcm2d(1:kpie,1:kpje,jdmsprod(n)),1,0)
          if (jdms(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'dms_acc'//c2,bgcm2d(1:kpie,1:kpje,jdms(n)),1,0)
          if (jniflux(n).ne.0) call read_netcdf_var(ncid,               &
     &      'niflux_acc'//c2,bgcm2d(1:kpie,1:kpje,jniflux(n)),1,0)
          if (joxflux(n).ne.0) call read_netcdf_var(ncid,               &
     &      'oxflux_acc'//c2,bgcm2d(1:kpie,1:kpje,joxflux(n)),1,0)
          if (jco2fxu(n).ne.0) call read_netcdf_var(ncid,               &
     &      'co2fxu_acc'//c2,bgcm2d(1:kpie,1:kpje,jco2fxu(n)),1,0)
          if (jco2fxd(n).ne.0) call read_netcdf_var(ncid,               &
     &      'co2fxd_acc'//c2,bgcm2d(1:kpie,1:kpje,jco2fxd(n)),1,0)
          if (jdmsflux(n).ne.0) call read_netcdf_var(ncid,              &
     &      'dmsflux_acc'//c2,bgcm2d(1:kpie,1:kpje,jdmsflux(n)),1,0)
          if (jpco2(n).ne.0) call read_netcdf_var(ncid,                 &
     &      'pco2_acc'//c2,bgcm2d(1:kpie,1:kpje,jpco2(n)),1,0)
          if (jkwco2(n).ne.0) call read_netcdf_var(ncid,                &
     &      'kwco2_acc'//c2,bgcm2d(1:kpie,1:kpje,jkwco2(n)),1,0)
          if (jatmco2(n).ne.0) call read_netcdf_var(ncid,               &
     &      'atmco2_acc'//c2,bgcm2d(1:kpie,1:kpje,jatmco2(n)),1,0)
          if (jatmo2(n).ne.0) call read_netcdf_var(ncid,                &
     &      'atmo2_acc'//c2,bgcm2d(1:kpie,1:kpje,jatmo2(n)),1,0)
          if (jatmn2(n).ne.0) call read_netcdf_var(ncid,                &
     &      'atmn2_acc'//c2,bgcm2d(1:kpie,1:kpje,jatmn2(n)),1,0)
! 3D
          if (jlvlphyto(n).ne.0) call read_netcdf_var(ncid,             &
     &      'phyclvl_acc'//c2,                                          &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlphyto(n)),ddm,0)
          if (jlvlgrazer(n).ne.0) call read_netcdf_var(ncid,            &
     &       'zoolvl_acc'//c2,                                          &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlgrazer(n)),ddm,0)
          if (jlvldoc(n).ne.0) call read_netcdf_var(ncid,               &
     &      'doclvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvldoc(n)),ddm,0)
          if (jlvlphosy(n).ne.0) call read_netcdf_var(ncid,             &
     &      'pplvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlphosy(n)),ddm,0)
          if (jlvlphosph(n).ne.0) call read_netcdf_var(ncid,            &
     &      'po4lvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlphosph(n)),ddm,0)
          if (jlvloxygen(n).ne.0) call read_netcdf_var(ncid,            &
     &      'oxylvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvloxygen(n)),ddm,0)
          if (jlvliron(n).ne.0) call read_netcdf_var(ncid,              &
     &      'felvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvliron(n)),ddm,0)
          if (jlvlano3(n).ne.0) call read_netcdf_var(ncid,              &
     &      'no3lvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlano3(n)),ddm,0)
          if (jlvlalkali(n).ne.0) call read_netcdf_var(ncid,            &
     &      'alklvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlalkali(n)),ddm,0)
          if (jlvlsilica(n).ne.0) call read_netcdf_var(ncid,            &
     &      'silvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlsilica(n)),ddm,0)
          if (jlvldic(n).ne.0) call read_netcdf_var(ncid,               &
     &      'diclvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvldic(n)),ddm,0)
          if (jlvlpoc(n).ne.0) call read_netcdf_var(ncid,               &
     &      'poclvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlpoc(n)),ddm,0)
          if (jlvlcalc(n).ne.0) call read_netcdf_var(ncid,              &
     &      'calvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlcalc(n)),ddm,0)
          if (jlvlopal(n).ne.0) call read_netcdf_var(ncid,              &
     &      'oplvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlopal(n)),ddm,0)
          if (jlvlco3(n).ne.0) call read_netcdf_var(ncid,               &
     &      'co3lvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlco3(n)),ddm,0)
          if (jlvlph(n).ne.0) call read_netcdf_var(ncid,                &
     &      'phlvl_acc'//c2,                                            &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlph(n)),ddm,0)
          if (jlvlomegac(n).ne.0) call read_netcdf_var(ncid,            &
     &      'omegaclvl_acc'//c2,                                        &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlomegac(n)),ddm,0)
          if (jlvldic13(n).ne.0) call read_netcdf_var(ncid,             &
     &      'dic13lvl_acc'//c2,                                         &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvldic13(n)),ddm,0)
          if (jlvldic14(n).ne.0) call read_netcdf_var(ncid,             &
     &      'dic14lvl_acc'//c2,                                         &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvldic14(n)),ddm,0)
          if (jlvlnos(n).ne.0) call read_netcdf_var(ncid,               &
     &      'noslvl_acc'//c2,                                           &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlnos(n)),ddm,0)
          if (jlvlcfc11_t(n).ne.0) call read_netcdf_var(ncid,           &
     &      'cfc11lvl_acc'//c2,                                         &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlcfc11_t(n)),ddm,0)
          if (jlvlcfc12_t(n).ne.0) call read_netcdf_var(ncid,           &
     &      'cfc12lvl_acc'//c2,                                         &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlcfc12_t(n)),ddm,0)
          if (jlvlac14_t(n).ne.0) call read_netcdf_var(ncid,            &
     &      'ac14lvl_acc'//c2,                                          &
     &       bgcm3dlvl(1:kpie,1:kpje,:,jlvlac14_t(n)),ddm,0)
! 3D sigma layers
          if (jdp(n).ne.0) call read_netcdf_var(ncid,                   &
     &      'dp_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jdp(n)),kpke,0)
          if (jphyto(n).ne.0) call read_netcdf_var(ncid,                &
     &      'phyc_acc'//c2,                                             &
     &       bgcm3d(1:kpie,1:kpje,:,jphyto(n)),kpke,0)
          if (jgrazer(n).ne.0) call read_netcdf_var(ncid,               &
     &       'zoo_acc'//c2,                                             &
     &       bgcm3d(1:kpie,1:kpje,:,jgrazer(n)),kpke,0)
          if (jdoc(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'doc_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jdoc(n)),kpke,0)
          if (jphosy(n).ne.0) call read_netcdf_var(ncid,                &
     &      'pp_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jphosy(n)),kpke,0)
          if (jphosph(n).ne.0) call read_netcdf_var(ncid,               &
     &      'po4_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jphosph(n)),kpke,0)
          if (joxygen(n).ne.0) call read_netcdf_var(ncid,               &
     &      'oxy_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,joxygen(n)),kpke,0)
          if (jiron(n).ne.0) call read_netcdf_var(ncid,                 &
     &      'fe_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jiron(n)),kpke,0)
          if (jano3(n).ne.0) call read_netcdf_var(ncid,                 &
     &      'no3_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jano3(n)),kpke,0)
          if (jalkali(n).ne.0) call read_netcdf_var(ncid,               &
     &      'alk_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jalkali(n)),kpke,0)
          if (jsilica(n).ne.0) call read_netcdf_var(ncid,               &
     &      'si_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jsilica(n)),kpke,0)
          if (jdic(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'dic_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jdic(n)),kpke,0)
          if (jpoc(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'poc_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jpoc(n)),kpke,0)
          if (jcalc(n).ne.0) call read_netcdf_var(ncid,                 &
     &      'ca_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jcalc(n)),kpke,0)
          if (jopal(n).ne.0) call read_netcdf_var(ncid,                 &
     &      'op_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jopal(n)),kpke,0)
          if (jco3(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'co3_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jco3(n)),kpke,0)
          if (jph(n).ne.0) call read_netcdf_var(ncid,                   &
     &      'ph_acc'//c2,                                               &
     &       bgcm3d(1:kpie,1:kpje,:,jph(n)),kpke,0)
          if (jomegac(n).ne.0) call read_netcdf_var(ncid,               &
     &      'omegac_acc'//c2,                                           &
     &       bgcm3d(1:kpie,1:kpje,:,jomegac(n)),kpke,0)
          if (jdic13(n).ne.0) call read_netcdf_var(ncid,                &
     &      'dic13_acc'//c2,                                            &
     &       bgcm3d(1:kpie,1:kpje,:,jdic13(n)),kpke,0)
          if (jdic14(n).ne.0) call read_netcdf_var(ncid,                &
     &      'dic14_acc'//c2,                                            &
     &       bgcm3d(1:kpie,1:kpje,:,jdic14(n)),kpke,0)
          if (jnos(n).ne.0) call read_netcdf_var(ncid,                  &
     &      'nos_acc'//c2,                                              &
     &       bgcm3d(1:kpie,1:kpje,:,jnos(n)),kpke,0)
          if (jcfc11_t(n).ne.0) call read_netcdf_var(ncid,              &
     &      'cfc11_acc'//c2,                                            &
     &       bgcm3d(1:kpie,1:kpje,:,jcfc11_t(n)),kpke,0)
          if (jcfc12_t(n).ne.0) call read_netcdf_var(ncid,              &
     &      'cfc12_acc'//c2,                                            &
     &       bgcm3d(1:kpie,1:kpje,:,jcfc12_t(n)),kpke,0)
          if (jac14_t(n).ne.0) call read_netcdf_var(ncid,               &
     &      'ac14_acc'//c2,                                             &
     &       bgcm3d(1:kpie,1:kpje,:,jac14_t(n)),kpke,0)
! Sediments
          if (jpowaic(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powaic_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowaic(n)),ks,0)
          if (jpowaal(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powaal_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowaal(n)),ks,0)
          if (jpowaph(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powaph_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowaph(n)),ks,0)
          if (jpowaox(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powaox_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowaox(n)),ks,0)
          if (jpown2(n).ne.0) call read_netcdf_var(ncid,                &
     &      'pown2_acc'//c2,                                            &
     &       bgct_sed(1:kpie,1:kpje,:,jpown2(n)),ks,0)
          if (jpowno3(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powno3_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowno3(n)),ks,0)
          if (jpowasi(n).ne.0) call read_netcdf_var(ncid,               &
     &      'powasi_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jpowasi(n)),ks,0)
          if (jssso12(n).ne.0) call read_netcdf_var(ncid,               &
     &      'ssso12_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jssso12(n)),ks,0)
          if (jssssil(n).ne.0) call read_netcdf_var(ncid,               &
     &      'ssssil_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jssssil(n)),ks,0)
          if (jsssc12(n).ne.0) call read_netcdf_var(ncid,               &
     &      'sssc12_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jsssc12(n)),ks,0)
          if (jssster(n).ne.0) call read_netcdf_var(ncid,               &
     &      'ssster_acc'//c2,                                           &
     &       bgct_sed(1:kpie,1:kpje,:,jssster(n)),ks,0)
        endif
      enddo
! read partly accumulated diagnostic fields from restart - end

      if(mnproc==1) ncstat = NF_CLOSE(ncid)

#else

     ! Attention - this doesn't work any more for MPI
     ! One big FORTRAN READ like here is a very, very bad idea
     ! for parallel programs !!!!!

#ifndef NOMPI
     call stop_all("can't read in bgc restart with MPI and no NETCDF")
! would need to read these in on p_io and scatter if we were actually i
! going to do this - c.f. read_netcdf_var.F90

#else

      OPEN(io_rsti_bgc,FILE='restart_bgc',STATUS='UNKNOWN'             &
     &               ,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')

      READ(io_rsti_bgc)                                                &
     &           (((ocetra(i,j,k,iphosph),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,isilica),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,ioxygen),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iphy   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,izoo   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,idet   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,idoc   ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,icalc  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iano3  ),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,igasnit),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,iopal  ),i=1,kpie),j=1,kpje),k=1,kpke)


      READ(io_rsti_bgc) chemcm,hi,co3,aksp                             &
     &          ,(((ocetra(i,j,k,isco212),i=1,kpie),j=1,kpje),k=1,kpke)&
     &          ,(((ocetra(i,j,k,ialkali),i=1,kpie),j=1,kpje),k=1,kpke)

      READ(io_rsti_bgc) sedlay,sedhpl

      READ(io_rsti_bgc) 
     &           (((powtra(i,j,k,ipowaic),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaal),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaph),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowaox),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowasi),i=1,kpie),j=1,kpje),k=1,ks)  &
     &          ,(((powtra(i,j,k,ipowno3),i=1,kpie),j=1,kpje),k=1,ks)  &
     &           ,(((powtra(i,j,k,ipown2) ,i=1,kpie),j=1,kpje),k=1,ks)

      CLOSE (io_rsti_bgc)
#endif

#endif

!     
!  Masking aqueous sea water tracer.
!
!      DO l=1,nocetra
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(pddpo(i,j,k) .LT. 0.5) THEN
!         ocetra(i,j,k,l)=rmasko
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO
! only for the first run!!      
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!        ocetra(i,j,k,isco214)=ocetra(i,j,k,isco212)*0.75
!      ENDDO
!      ENDDO
!      ENDDO      
!
!  Masking other sea water tracer.
!
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(pddpo(i,j,k) .LT. 0.5) THEN
!ik      phyto(i,j,k)=rmasko
!ik      grazer(i,j,k)=rmasko
!ik      poc(i,j,k)=rmasko
!         hi(i,j,k)=rmasko
!ik      calciu(i,j,k)=rmasko
!         co3(i,j,k)=rmasko
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

!
!  Masking sediment pore water tracer.
!
!      DO  k=1,ks
!      DO  j=1,kpje
!      DO  i=1,kpie 
!      IF(kbo(i,j) .LE. 1) THEN
!         powtra(i,j,k,ipowaic)=rmasks
!         powtra(i,j,k,ipowaal)=rmasks
!         powtra(i,j,k,ipowaph)=rmasks
!         powtra(i,j,k,ipowaox)=rmasks
!         powtra(i,j,k,ipown2)=rmasks
!         powtra(i,j,k,ipowno3)=rmasks
!         powtra(i,j,k,ipowasi)=rmasks
!         sedlay(i,j,k,issso12)=rmasks
!         sedlay(i,j,k,isssc12)=rmasks
!         sedlay(i,j,k,issssil)=rmasks
!         burial(i,j,issso12)=rmasks
!         burial(i,j,isssc12)=rmasks
!         burial(i,j,issssil)=rmasks
!         sedhpl(i,j,k)=rmasks
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

!     
!  Restrict to positive values (until error is found only !!!!)
!

!      DO l=1,nocetra
!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(omask(i,j) .GT. 0.5) THEN
!	 ocetra(i,j,k,l)=MAX(ocetra(i,j,k,l),0.)
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO
!      ENDDO

!      DO k=1,kpke
!      DO j=1,kpje
!      DO i=1,kpie
!      IF(omask(i,j) .GT. 0.5) THEN
!ik         phyto (i,j,k)=MAX(phyto (i,j,k),0.)
!ik         grazer(i,j,k)=MAX(grazer(i,j,k),0.)
!ik         poc   (i,j,k)=MAX(poc   (i,j,k),0.)
!        hi    (i,j,k)=MAX(hi    (i,j,k),1.e-12)
!ik         calciu(i,j,k)=MAX(calciu(i,j,k),0.)
!ka multiply DIC to fix atmospheric CO2
!         ocetra(i,j,k,isco212)=ocetra(i,j,k,isco212)*0.995
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO


!      DO  k=1,ks
!      DO  j=1,kpje
!      DO  i=1,kpie 
!      IF(omask(i,j) .GT. 0.5) THEN
!         powtra(i,j,k,ipowaic)=MAX(powtra(i,j,k,ipowaic),0.)
!#ifdef __c_isotopes
!         powtra(i,j,k,ipowc13)=MAX(powtra(i,j,k,ipowc13),0.)
!         powtra(i,j,k,ipowc14)=MAX(powtra(i,j,k,ipowc14),0.)
!#endif
!         powtra(i,j,k,ipowaal)=MAX(powtra(i,j,k,ipowaal),0.)
!         powtra(i,j,k,ipowaph)=MAX(powtra(i,j,k,ipowaph),0.)
!         powtra(i,j,k,ipowaox)=MAX(powtra(i,j,k,ipowaox),0.)
!         powtra(i,j,k,ipown2) =MAX(powtra(i,j,k,ipown2) ,0.)
!         powtra(i,j,k,ipowno3)=MAX(powtra(i,j,k,ipowno3),0.)
!         powtra(i,j,k,ipowasi)=MAX(powtra(i,j,k,ipowasi),0.)
!         sedlay(i,j,k,issso12)=MAX(sedlay(i,j,k,issso12),0.)
!         sedlay(i,j,k,isssc12)=MAX(sedlay(i,j,k,isssc12),0.)
!#ifdef __c_isotopes
!         sedlay(i,j,k,isssc13)=MAX(sedlay(i,j,k,isssc13),0.)
!         sedlay(i,j,k,isssc14)=MAX(sedlay(i,j,k,isssc14),0.)
!#endif
!         sedlay(i,j,k,issssil)=MAX(sedlay(i,j,k,issssil),0.)
!         sedhpl(i,j,k)        =MAX(sedhpl(i,j,k)    ,1.e-12)
!      ENDIF
!      ENDDO
!      ENDDO
!      ENDDO

      RETURN
      END
