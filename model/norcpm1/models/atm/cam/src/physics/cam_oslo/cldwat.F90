#undef DEBUG
!#define AEROFFL
module cldwat
!----------------------------------------------------------------------- 
! 
! Purpose: Prognostic cloud water data and methods.
! 
! Public interfaces:
!
! inimc -- Initialize constants
! pcond -- Calculate prognostic condensate
!
! cldwat_fice   -- calculate fraction of condensate in ice phase (radiation partitioning)
!
! Author: P. Rasch, with Modifications by Minghua Zhang
! January 2010, modified by J. Kay to add precip fluxes for COSP simulator
! 
!         Modified by Alf Kirkevåg to include cldwat_par version of CAM4-Oslo/NorESM 
! based on CCSM4, September 2011. Note for DIAGNCDNC (and DIRIND version) of CAM4-Oslo:
! This is old code for diagnostic CCN-calculations. It is kept here as a possibility to 
! address at a later stage (e.g. for testing ourposes). In that case, the code needs to 
! be re-structured and tested, since this option has not been used or checked for some time.
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use spmd_utils,    only: masterproc
   use ppgrid,        only: pcols, pver, pverp
   use wv_saturation, only: estblf, hlatv, tmin, hlatf, rgasv, pcf, &
                            cp, epsqs, ttrice
   use abortutils,    only: endrun
   use cam_logfile,   only: iulog

   implicit none

!-----------------------------------------------------------------------
! PUBLIC: Make default data and interfaces private
!-----------------------------------------------------------------------
   private
   save
   public inimc, pcond          ! Public interfaces
   public cldwat_fice          ! Public interfaces
   public cldwat_readnl
   integer, public::  ktop      ! Level above 10 hPa

   real(r8),public ::  icritc               ! threshold for autoconversion of cold ice
   real(r8),public ::  icritw               ! threshold for autoconversion of warm ice
!!$   real(r8),public,parameter::  conke  = 1.e-6_r8    ! tunable constant for evaporation of precip
!!$   real(r8),public,parameter::  conke  =  2.e-6_r8    ! tunable constant for evaporation of precip
   real(r8),public ::  conke                ! tunable constant for evaporation of precip
   real(r8),public ::  r3lcrit              ! critical radius where liq conversion begins

!-----------------------------------------------------------------------
! PRIVATE: Everything else is private to this module
!-----------------------------------------------------------------------
   real(r8), private:: tmax_fice! max temperature for cloud ice formation
   real(r8), private:: tmin_fice! min temperature for cloud ice formation
   real(r8), private:: tmax_fsnow ! max temperature for transition to convective snow
   real(r8), private:: tmin_fsnow ! min temperature for transition to convective snow
   real(r8), private:: rhonot   ! air density at surface
   real(r8), private:: t0       ! Freezing temperature
   real(r8), private:: cldmin   ! assumed minimum cloud amount
   real(r8), private:: small    ! small number compared to unity
   real(r8), private:: c        ! constant for graupel like snow cm**(1-d)/s
   real(r8), private:: d        ! constant for graupel like snow
   real(r8), private:: esi      ! collection efficient for ice by snow
   real(r8), private:: esw      ! collection efficient for water by snow
   real(r8), private:: nos      ! particles snow / cm**4
   real(r8), private:: pi       ! Mathematical constant
   real(r8), private:: gravit   ! Gravitational acceleration at surface
   real(r8), private:: rh2o
   real(r8), private:: prhonos
   real(r8), private:: thrpd    ! numerical three added to d
   real(r8), private:: gam3pd   ! gamma function on (3+d)
   real(r8), private:: gam4pd   ! gamma function on (4+d)
   real(r8), private:: rhoi     ! ice density
   real(r8), private:: rhos     ! snow density
   real(r8), private:: rhow     ! water density
   real(r8), private:: mcon01   ! constants used in cloud microphysics
   real(r8), private:: mcon02   ! constants used in cloud microphysics
   real(r8), private:: mcon03   ! constants used in cloud microphysics
   real(r8), private:: mcon04   ! constants used in cloud microphysics
   real(r8), private:: mcon05   ! constants used in cloud microphysics
   real(r8), private:: mcon06   ! constants used in cloud microphysics
   real(r8), private:: mcon07   ! constants used in cloud microphysics
   real(r8), private:: mcon08   ! constants used in cloud microphysics

   integer, private ::  k1mb    ! index of the eta level near 1 mb

! Parameters used in findmcnew
   real(r8) :: capnsi               ! sea ice cloud particles / cm3
   real(r8) :: capnc                ! cold and oceanic cloud particles / cm3
   real(r8) :: capnw                ! warm continental cloud particles / cm3
   real(r8) :: kconst               ! const for terminal velocity (stokes regime)
   real(r8) :: effc                 ! collection efficiency
   real(r8) :: alpha                ! ratio of 3rd moment radius to 2nd
   real(r8) :: capc                 ! constant for autoconversion
   real(r8) :: convfw               ! constant used for fall velocity calculation
   real(r8) :: cracw                ! constant used for rain accreting water
   real(r8) :: critpr               ! critical precip rate collection efficiency changes
   real(r8) :: ciautb               ! coefficient of autoconversion of ice (1/s)

#ifdef DEBUG
   integer, private,parameter ::  nlook = 1  ! Number of points to examine
   integer, private ::  ilook(nlook)         ! Longitude index to examine
   integer, private ::  latlook(nlook)       ! Latitude index to examine
   integer, private ::  lchnklook(nlook)     ! Chunk index to examine
   integer, private ::  icollook(nlook)      ! Column index to examine
#endif
  ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)

contains
!===============================================================================
  subroutine cldwat_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input




   ! Namelist variables
   real(r8) :: cldwat_icritw  = unset_r8    !   icritw  = threshold for autoconversion of warm ice  
   real(r8) :: cldwat_icritc  = unset_r8    !   icritc  = threshold for autoconversion of cold ice  
   real(r8) :: cldwat_conke   = unset_r8    !   conke   = tunable constant for evaporation of precip
   real(r8) :: cldwat_r3lcrit = unset_r8    !   r3lcrit = critical radius where liq conversion begins

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'cldwat_readnl'

   namelist /cldwat_nl/ cldwat_icritw, cldwat_icritc, cldwat_conke, cldwat_r3lcrit

   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'cldwat_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, cldwat_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)

      ! set local variables
      icritw  = cldwat_icritw 
      icritc  = cldwat_icritc
      conke   = cldwat_conke
      r3lcrit = cldwat_r3lcrit

   end if



#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(icritw,            1, mpir8,  0, mpicom)
   call mpibcast(icritc,            1, mpir8,  0, mpicom)
   call mpibcast(conke,             1, mpir8,  0, mpicom)
   call mpibcast(r3lcrit,           1, mpir8,  0, mpicom)
#endif

end subroutine cldwat_readnl

!================================================================================================
  subroutine cldwat_fice(ncol, t, fice, fsnow)
!
! Compute the fraction of the total cloud water which is in ice phase.
! The fraction depends on temperature only. 
! This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
!  modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
!-----------------------------------------------------------------------
    implicit none

! Arguments
    integer,  intent(in)  :: ncol                 ! number of active columns
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature

    real(r8), intent(out) :: fice(pcols,pver)     ! Fractional ice content within cloud
    real(r8), intent(out) :: fsnow(pcols,pver)    ! Fractional snow content for convection

! Local variables
    integer :: i,k                                   ! loop indexes

!-----------------------------------------------------------------------

! Define fractional amount of cloud that is ice
    do k=1,pver
       do i=1,ncol

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fice) then
             fice(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fice) then
             fice(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fice(i,k) =(tmax_fice - t(i,k)) / (tmax_fice - tmin_fice)
          end if

! snow fraction partitioning

! If warmer than tmax then water phase
          if (t(i,k) > tmax_fsnow) then
             fsnow(i,k) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(i,k) < tmin_fsnow) then
             fsnow(i,k) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fsnow(i,k) =(tmax_fsnow - t(i,k)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do

    return
  end subroutine cldwat_fice

subroutine inimc( tmeltx, rhonotx, gravitx, rh2ox)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! initialize constants for the prognostic condensate
! 
! Author: P. Rasch, April 1997
! 
!-----------------------------------------------------------------------
   use pmgrid,       only: plev, plevp
   use dycore,       only: dycore_is, get_resolution
   use ref_pres,       only: pref_mid
   use phys_control, only: phys_getopts

   integer k
   real(r8), intent(in) :: tmeltx
   real(r8), intent(in) :: rhonotx
   real(r8), intent(in) :: gravitx
   real(r8), intent(in) :: rh2ox


#ifdef UNICOSMP
   real(r8) signgam              ! variable required by cray gamma function
   external gamma
#endif
   character(len=16)    :: microp_scheme 

! Get microphysics option

   call phys_getopts( microp_scheme_out = microp_scheme )

! Set following for all physics packages

   tmax_fice = tmeltx    - 10._r8
!! tmax_fice = tmeltx
!! tmin_fice = tmax_fice - 20._r8
   tmin_fice = tmax_fice - 30._r8
   tmax_fsnow = tmeltx
   tmin_fsnow = tmeltx   - 5._r8

! Set remaining for RK microphysics

   if( microp_scheme .eq. 'RK' ) then
      rhonot = rhonotx             ! air density at surface (gm/cm3)
      gravit = gravitx
      rh2o   = rh2ox
      rhos = .1_r8                 ! assumed snow density (gm/cm3)
      rhow = 1._r8                 ! water density
      rhoi = 1._r8                 ! ice density
      esi = 1.0_r8                 ! collection efficient for ice by snow
      esw = 0.1_r8                 ! collection efficient for water by snow
      t0 = tmeltx                  ! approximate freezing temp
      cldmin = 0.02_r8             ! assumed minimum cloud amount
      small = 1.e-22_r8            ! a small number compared to unity
      c = 152.93_r8                ! constant for graupel like snow cm**(1-d)/s
      d = 0.25_r8                  ! constant for graupel like snow
      nos = 3.e-2_r8               ! particles snow / cm**4
      pi = 4._r8*atan(1.0_r8)
      prhonos = pi*rhos*nos
      thrpd = 3._r8 + d
      if (d==0.25_r8) then
         gam3pd = 2.549256966718531_r8 ! only right for d = 0.25
         gam4pd = 8.285085141835282_r8
      else
#ifdef UNICOSMP
         call gamma(3._r8+d, signgam, gam3pd)
         gam3pd = sign(exp(gam3pd),signgam)
         call gamma(4._r8+d, signgam, gam4pd)
         gam4pd = sign(exp(gam4pd),signgam)
         write(iulog,*) ' d, gamma(3+d), gamma(4+d) =', gam3pd, gam4pd
#else
         write(iulog,*) ' can only use d ne 0.25 on a cray '
         stop
#endif
      endif
      mcon01 = pi*nos*c*gam3pd/4._r8
      mcon02 = 1._r8/(c*gam4pd*sqrt(rhonot)/(6*prhonos**(d/4._r8)))
      mcon03 = -(0.5_r8+d/4._r8)
      mcon04 = 4._r8/(4._r8+d)
      mcon05 = (3+d)/(4+d)
      mcon06 = (3+d)/4._r8
      mcon07 = mcon01*sqrt(rhonot)*mcon02**mcon05/prhonos**mcon06
      mcon08 = -0.5_r8/(4._r8+d)

!  find the level about 1mb, we wont do the microphysics above this level
      k1mb = 1
      do k=1,pver-1
         if (pref_mid(k) < 1.e2_r8 .and. pref_mid(k+1) >= 1.e2_r8) then
            if (1.e2_r8-pref_mid(k) < pref_mid(k+1)-1.e2_r8) then
               k1mb = k
            else
               k1mb = k + 1
            end if
            goto 20
         end if
      end do
      if (masterproc) then
         write(iulog,*)'inimc: model levels bracketing 1 mb not found'
      end if
!     call endrun
      k1mb = 1
20    if( masterproc ) write(iulog,*)'inimc: model level nearest 1 mb is',k1mb,'which is',pref_mid(k1mb),'pascals'

      if( masterproc ) write(iulog,*) 'cloud water initialization by inimc complete '

! Initialize parameters used by findmcnew
      capnw = 400._r8              ! warm continental cloud particles / cm3
      capnc = 150._r8              ! cold and oceanic cloud particles / cm3
!     capnsi = 40._r8              ! sea ice cloud particles density  / cm3
      capnsi = 75._r8              ! sea ice cloud particles density  / cm3

      kconst = 1.18e6_r8           ! const for terminal velocity

!     effc = 1._r8                 ! autoconv collection efficiency following boucher 96
!     effc = .55*0.05_r8           ! autoconv collection efficiency following baker 93
      effc = 0.55_r8               ! autoconv collection efficiency following tripoli and cotton
!     effc = 0._r8    ! turn off warm-cloud autoconv
      alpha = 1.1_r8**4
      capc = pi**(-.333_r8)*kconst*effc *(0.75_r8)**(1.333_r8)*alpha  ! constant for autoconversion
!r3lcrit: critical radius where liq conversion begins (10.0 micron)
#ifdef AERLIFE
#ifdef DIAGNCDNC 
   r3lcrit = 10.0e-6_r8
#else
   r3lcrit = 14.0e-6_r8
!orig   r3lcrit = 10.0e-6_r8        ! bare midlertidig --> lik met
#endif
#endif

! critical precip rate at which we assume the collector drops can change the
! drop size enough to enhance the auto-conversion process (mm/day)
#ifdef AERLIFE
#ifdef DIAGNCDNC
   critpr = 0.5_r8
#else
   critpr = 5.0_r8 
!orig   critpr = 0.5_r8 
#endif
#else
      critpr = 0.5_r8        ! bare midlertidig --> lik met
#endif
      convfw = 1.94_r8*2.13_r8*sqrt(rhow*1000._r8*9.81_r8*2.7e-4_r8)

! liquid microphysics
!     cracw = 6_r8                 ! beheng
      cracw = .884_r8*sqrt(9.81_r8/(rhow*1000._r8*2.7e-4_r8)) ! tripoli and cotton

! ice microphysics
      ciautb = 5.e-4_r8

      if ( masterproc ) then
         write(iulog,*)'tuning parameters cldwat: icritw',icritw,'icritc',icritc,'conke',conke,'r3lcrit',r3lcrit
         write(iulog,*)'tuning parameters cldwat: capnw',capnw,'capnc',capnc,'capnsi',capnsi,'kconst',kconst
         write(iulog,*)'tuning parameters cldwat: effc',effc,'alpha',alpha,'capc',capc
         write(iulog,*)'tuning parameters cldwat: critpr',critpr,'convfw',convfw,'cracw',cracw,'ciautb',ciautb
      endif
   endif

   return
end subroutine inimc

#ifdef DIRIND
subroutine pcond (lchnk   ,ncol    , &
                  tn      ,ttend   ,qn      ,qtend   ,omega   , &
                  cwat    ,p       ,pdel    ,cldn    ,fice    , fsnow, &
                  cme     ,prodprec,prodsnow,evapprec,evapsnow,evapheat, prfzheat, &     
                  meltheat,precip  ,snowab  ,deltat  ,fwaut   , &
                  fsaut   ,fracw   ,fsacw   ,fsaci   ,lctend  , &
                  rhdfda  ,rhu00   ,seaicef, zi      ,ice2pr, liq2pr, &
                  liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio, &
                  lctendx, qm, precc, cwatx, sigw1, cldo, cnlxin, dlf, &        ! in
                  liq2prx, relca, ntotal, cdnc, cxsind, n1, n2, n3, n4, n5, &   ! out
                  n6, n7, n8, n9, n10, n11, n12, n13, n14, nrmodes, relh, &     ! out
                  cmex, selfx, lossx, nucrat, nrainx, supersat, evapx, freez )  ! out
#else
subroutine pcond (lchnk   ,ncol    , &
                  tn      ,ttend   ,qn      ,qtend   ,omega   , &
                  cwat    ,p       ,pdel    ,cldn    ,fice    , fsnow, &
                  cme     ,prodprec,prodsnow,evapprec,evapsnow,evapheat, prfzheat, &     
                  meltheat,precip  ,snowab  ,deltat  ,fwaut   , &
                  fsaut   ,fracw   ,fsacw   ,fsaci   ,lctend  , &
                  rhdfda  ,rhu00   ,seaicef, zi      ,ice2pr, liq2pr, &
                  liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio)
#endif

!----------------------------------------------------------------------- 
! 
! Purpose: 
! The public interface to the cloud water parameterization
! returns tendencies to water vapor, temperature and cloud water variables
! 
! For basic method 
!  See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
!  model climate using diagnosed and 
!  predicted condensate parameterizations, 1998, J. Clim., 11,
!  pp1587---1614.
! 
! For important modifications to improve the method of determining
! condensation/evaporation see Zhang et al (2001, in preparation)
!
! Authors: M. Zhang, W. Lin, P. Rasch and J.E. Kristjansson
!          B. A. Boville (latent heat of fusion)
!-----------------------------------------------------------------------
   use wv_saturation, only: vqsatd
   use cam_control_mod, only: nlvdry
   use constituents, only: pcnst, cnst_get_ind
   use comsrf,       only: landm
!
!---------------------------------------------------------------------
!
! Input Arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: fice(pcols,pver)     ! fraction of cwat that is ice
   real(r8), intent(in) :: fsnow(pcols,pver)    ! fraction of rain that freezes to snow
   real(r8), intent(in) :: cldn(pcols,pver)     ! new value of cloud fraction    (fraction)
   real(r8), intent(in) :: cwat(pcols,pver)     ! cloud water (kg/kg)
   real(r8), intent(in) :: omega(pcols,pver)    ! vert pressure vel (Pa/s)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure          (K)
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure thickness (Pa)
   real(r8), intent(in) :: qn(pcols,pver)       ! new water vapor    (kg/kg)
   real(r8), intent(in) :: qtend(pcols,pver)    ! mixing ratio tend  (kg/kg/s)
   real(r8), intent(in) :: tn(pcols,pver)       ! new temperature    (K)
   real(r8), intent(in) :: ttend(pcols,pver)    ! temp tendencies    (K/s)
   real(r8), intent(in) :: deltat               ! time step to advance solution over
   real(r8), intent(in) :: lctend(pcols,pver)   ! cloud liquid water tendencies   ====wlin
   real(r8), intent(in) :: rhdfda(pcols,pver)   ! dG(a)/da, rh=G(a), when rh>u00  ====wlin
   real(r8), intent(in) :: rhu00 (pcols,pver)   ! Rhlim for cloud                 ====wlin
   real(r8), intent(in) :: seaicef(pcols)       ! sea ice fraction  (fraction)
   real(r8), intent(in) :: zi(pcols,pverp)      ! layer interfaces (m)
   real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
#ifdef DIRIND
   real(r8), intent(in) :: qm(pcols,pver,pcnst) ! Common aerosol (only!) tracers for indirect and direct calculations
   real(r8), intent(in) :: precc(pcols)          ! convective-scale precipitation rate
   real(r8), intent(in) :: cwatx(pcols,pver)     ! cloud water (kg/kg)
   real(r8) :: cnlx(pcols,pver)     ! cloud droplet number concentration (cm-3)
   real(r8), intent(in) :: cnlxin(pcols,pver)     ! cloud droplet number concentration (cm-3)
   real(r8), intent(in) :: cldo(pcols,pver)     ! TS old value of cloud fraction    (fraction)
   real(r8), intent(in) :: lctendx(pcols,pver)  ! offline cloud liquid water tendencies   ====wlin
!corinna for modified detrainment from convective clouds
   real(r8), intent(in) :: dlf(pcols,pver)       ! detrained water from convective clouds
   real(r8), intent(in) :: sigw1(pcols,pver)     ! Standard deviation, vertical velocity distribution
#endif
!
! Output Arguments
!
   real(r8), intent(out) :: cme     (pcols,pver) ! rate of cond-evap of condensate (1/s)
   real(r8), intent(out) :: prodprec(pcols,pver) ! rate of conversion of condensate to precip (1/s)
   real(r8), intent(out) :: evapprec(pcols,pver) ! rate of evaporation of falling precip (1/s)
   real(r8), intent(out) :: evapsnow(pcols,pver) ! rate of evaporation of falling snow (1/s)
   real(r8), intent(out) :: evapheat(pcols,pver) ! heating rate due to evaporation of precip (W/kg)
   real(r8), intent(out) :: prfzheat(pcols,pver) ! heating rate due to freezing of precip (W/kg)
   real(r8), intent(out) :: meltheat(pcols,pver) ! heating rate due to snow melt (W/kg)
   real(r8), intent(out) :: precip(pcols)        ! rate of precipitation (kg / (m**2 * s))
   real(r8), intent(out) :: snowab(pcols)        ! rate of snow (kg / (m**2 * s))
   real(r8), intent(out) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
   real(r8), intent(out) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
   real(r8), intent(out) :: liq2snow(pcols,pver) ! rate of conversion of liquid to snow
   real(r8), intent(out) :: rkflxprc(pcols,pverp)   ! grid-box mean RK flux_large_scale_cloud_rain+snow at interfaces (kg m^-2 s^-1)
   real(r8), intent(out) :: rkflxsnw(pcols,pverp)   ! grid-box mean RK flux_large_scale_cloud_snow at interfaces (kg m^-2 s^-1)
! intent(out)s here for pcond to pass to stratiform.F90 to be addflded/outflded
   real(r8), intent(out) :: pracwo(pcols,pver)      ! accretion of cloud water by rain (1/s)
   real(r8), intent(out) :: psacwo(pcols,pver)      ! accretion of cloud water by snow (1/s)
   real(r8), intent(out) :: psacio(pcols,pver)      ! accretion of cloud ice by snow (1/s)

   real(r8) nice2pr     ! rate of conversion of ice to snow
   real(r8) nliq2pr     ! rate of conversion of liquid to precip
   real(r8) nliq2snow   ! rate of conversion of liquid to snow
   real(r8) prodsnow(pcols,pver) ! rate of production of snow

#ifdef DIRIND
   real(r8), intent(out) :: ntotal(pcols,pver)   ! total aerosol number concentration (3D)
   real(r8), intent(out) :: relca(pcols,pver)    ! computed effective radius of liquid clouds 
   real(r8), intent(out) :: cdnc(pcols,pver)     ! cloud droplet number concentration
   real(r8), intent(out) :: cxsind(pcols,pver)   ! total excess mass concentration wrt. internal mixing
   real(r8), intent(out) :: n1(pcols,pver)
   real(r8), intent(out) :: n2(pcols,pver)
   real(r8), intent(out) :: n3(pcols,pver)
   real(r8), intent(out) :: n4(pcols,pver)
   real(r8), intent(out) :: n5(pcols,pver)
   real(r8), intent(out) :: n6(pcols,pver)
   real(r8), intent(out) :: n7(pcols,pver)
   real(r8), intent(out) :: n8(pcols,pver)
   real(r8), intent(out) :: n9(pcols,pver)
   real(r8), intent(out) :: n10(pcols,pver)
   real(r8), intent(out) :: n11(pcols,pver)
   real(r8), intent(out) :: n12(pcols,pver)
   real(r8), intent(out) :: n13(pcols,pver)
   real(r8), intent(out) :: n14(pcols,pver)
   real(r8), intent(out) :: nrmodes(pcols,pver,pcnst) ! number concentration in each mode 
   real(r8), intent(out) :: liq2prx(pcols,pver)  ! rate of conversion of liquid to precip
   real(r8), intent(out) :: cmex(pcols,pver)     ! rate of cond-evap within the cloud
   real(r8), intent(out) :: selfx(pcols,pver)    ! rate of selfcollection of cloud droplets (cm-3/s)
   real(r8), intent(out) :: freez(pcols,pver)    ! rate of freezing of cloud droplets (simple) (cm-3/s)
   real(r8), intent(out) :: evapx(pcols,pver) 	 ! Evaporation of cloud droplets (cm-3/s)
   real(r8), intent(out) :: nucrat(pcols,pver)   ! Nucleation rate of cloud droplets (cm-3/s)
   real(r8), intent(out) :: lossx(pcols,pver)    ! Loss rate for cloud droplets (not including precip.) (cm-3/s)
   real(r8), intent(out) :: nrainx(pcols,pver)   ! Loss rate for cloud droplets due to precip. (cm-3/s)
   real(r8), intent(out) :: relh(pcols,pver)     ! relative humidity
   real(r8), intent(out) :: supersat(pcols,pver)
#endif

!
! Local workspace
!
#ifdef DIRIND
!  off-line variables
   real(r8) rhocgs		 ! Air density in cgs units   
   real(r8) capna(pcols)         ! Cloud droplet number concentration
   real(r8) rer3fac
   real(r8) volrad(pcols,pver)   ! Mean volume radius (um)
   real(r8) fracwx(pcols,pver)   ! relative importance of collection of liquid by rain
   real(r8) fsacix(pcols,pver)   ! relative importance of collection of ice by snow
   real(r8) fsacwx(pcols,pver)   ! relative importance of collection of liquid by snow
   real(r8) fsautx(pcols,pver)   ! relative importance of ice auto conversion
   real(r8) fwautx(pcols,pver)   ! relative importance of warm cloud autoconversion
   real(r8) icolx(pcols,pver)    
   real(r8) sigw(pcols)	         ! Same as sigw1, but for one vertical level
   real(r8) nuclx(pcols,pver)    ! Number of droplets nucleated over one timestep (cm-3)
   real(r8) icnlx(pcols)         ! In-cloud droplet number concentration (cm-3)
   real(r8) moistrho(pcols, pver)! Moist air density
   real(r8) ncoefx(pcols)        ! Conversion time scale for cloud droplets to turn into rain 
   real(r8) evaprx(pcols,pver)   ! rate of evaporation of falling precipitation (1/s)
   real(r8) prainx(pcols,pver)   ! rate of conversion of condensate to precipitation (1/s)
   real(r8) nliq2prx             ! rate of conversion of liquid to precip
   real(r8) coefx(pcols)         ! conversion time scale for condensate to rain
   real(r8) cwmx(pcols)          ! cwat mixing ratio at midpoint
   real(r8) cwnx(pcols)          ! cwat mixing ratio at end
   real(r8) ncwnx(pcols)         ! Cloud droplet number conc. at end
   real(r8) icwcx(pcols)         ! in-cloud water content (kg/kg)
   real(r8) precabx(pcols)       ! rate of precipitation (kg / (m**2 * s))
   real(r8) prprovx(pcols)       ! provisional value of precip at btm of layer 
   real(r8) prtmpx               ! work variable
   real(r8) polx                 ! work variable
   real(r8) cdtx                 ! work variable
   real(r8) npolx                ! work variable
   real(r8) ncdtx                ! work variable
   real(r8) cmec1x (pcols)       !c1    of new C - E scheme formulation
   real(r8) cmec2x (pcols)       !c2    of new C - E scheme formulation
   real(r8) cmec3x (pcols)       !c3    of new C - E scheme formulation
   real(r8) cmec4x (pcols)       !c4    of new C - E scheme formulation
   real(r8) csigmax(pcols)       !sigma of new C - E scheme formulation
   real(r8) nummode(pcols,pcnst)
   real(r8) rcwnx(pcols,2,pver)
   real(r8) epsbeta
   real(r8) :: dropmass, reffl   ! mass of a fixed-radius droplet (for detrainment)
#endif

   real(r8) :: precab(pcols)        ! rate of precipitation (kg / (m**2 * s))
   integer i                 ! work variable
   integer iter              ! #iterations for precipitation calculation
   integer k                 ! work variable
   integer l                 ! work variable
#ifdef DIRIND
   integer m                 ! work variable
#endif

   real(r8) cldm(pcols)          ! mean cloud fraction over the time step
   real(r8) cldmax(pcols)        ! max cloud fraction above
   real(r8) coef(pcols)          ! conversion time scale for condensate to rain
   real(r8) cwm(pcols)           ! cwat mixing ratio at midpoint of time step
   real(r8) cwn(pcols)           ! cwat mixing ratio at end
   real(r8) denom                ! work variable
   real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
   real(r8) es(pcols)            ! sat. vapor pressure
   real(r8) fracw(pcols,pver)    ! relative importance of collection of liquid by rain
   real(r8) fsaci(pcols,pver)    ! relative importance of collection of ice by snow
   real(r8) fsacw(pcols,pver)    ! relative importance of collection of liquid by snow
   real(r8) fsaut(pcols,pver)    ! relative importance of ice auto conversion
   real(r8) fwaut(pcols,pver)    ! relative importance of warm cloud autoconversion
   real(r8) gamma(pcols)         ! d qs / dT
   real(r8) icwc(pcols)          ! in-cloud water content (kg/kg)
   real(r8) mincld               ! a small cloud fraction to avoid / zero
   real(r8) omeps                ! 1 minus epsilon
   real(r8),parameter ::omsm=0.99999_r8                 ! a number just less than unity (for rounding)
   real(r8) prprov(pcols)        ! provisional value of precip at btm of layer
   real(r8) prtmp                ! work variable
   real(r8) q(pcols,pver)        ! mixing ratio before time step ignoring condensate
   real(r8) qs(pcols)            ! spec. hum. of water vapor
   real(r8) qsn, esn             ! work variable
   real(r8) qsp(pcols,pver)      ! sat pt mixing ratio
   real(r8) qtl(pcols)           ! tendency which would saturate the grid box in deltat
   real(r8) qtmp, ttmp           ! work variable
   real(r8) relhum1(pcols)        ! relative humidity
   real(r8) relhum(pcols)        ! relative humidity
!!$   real(r8) tc                   ! crit temp of transition to ice
   real(r8) t(pcols,pver)        ! temp before time step ignoring condensate
   real(r8) tsp(pcols,pver)      ! sat pt temperature
   real(r8) pol                  ! work variable
   real(r8) cdt                  ! work variable
   real(r8) wtthick              ! work variable

! Extra local work space for cloud scheme modification       

   real(r8) cpohl                !Cp/Hlatv
   real(r8) hlocp                !Hlatv/Cp
   real(r8) dto2                 !0.5*deltat (delta=2.0*dt)
   real(r8) calpha(pcols)        !alpha of new C - E scheme formulation
   real(r8) cbeta (pcols)        !beta  of new C - E scheme formulation
   real(r8) cbetah(pcols)        !beta_hat at saturation portion 
   real(r8) cgamma(pcols)        !gamma of new C - E scheme formulation
   real(r8) cgamah(pcols)        !gamma_hat at saturation portion
   real(r8) rcgama(pcols)        !gamma/gamma_hat
   real(r8) csigma(pcols)        !sigma of new C - E scheme formulation
   real(r8) cmec1 (pcols)        !c1    of new C - E scheme formulation
   real(r8) cmec2 (pcols)        !c2    of new C - E scheme formulation
   real(r8) cmec3 (pcols)        !c3    of new C - E scheme formulation
   real(r8) cmec4 (pcols)        !c4    of new C - E scheme formulation
   real(r8) cmeres(pcols)        !residual cond of over-sat after cme and evapprec
   real(r8) ctmp                 !a scalar representation of cmeres
   real(r8) clrh2o               ! Ratio of latvap to water vapor gas const
   real(r8) ice(pcols,pver)    ! ice mixing ratio
   real(r8) liq(pcols,pver)    ! liquid mixing ratio
   real(r8) rcwn(pcols,2,pver), rliq(pcols,2,pver), rice(pcols,2,pver)
   real(r8) cwnsave(pcols,2,pver), cmesave(pcols,2,pver)
   real(r8) prodprecsave(pcols,2,pver)
   logical error_found
!
!------------------------------------------------------------
!
   clrh2o = hlatv/rh2o   ! Ratio of latvap to water vapor gas const
   omeps = 1.0_r8 - epsqs
#ifdef PERGRO
   mincld = 1.e-4_r8
   iter = 1   ! number of times to iterate the precipitation calculation
#else
   mincld = 1.e-4_r8
   iter = 2
#endif
!   omsm = 0.99999_r8
   cpohl = cp/hlatv
   hlocp = hlatv/cp
   dto2=0.5_r8*deltat
!
! Constant for computing rate of evaporation of precipitation:
!
!!$   conke = 1.e-5_r8
!!$   conke = 1.e-6_r8
!
! initialize a few single level fields
!
   do i = 1,ncol
      precip(i) = 0.0_r8
      precab(i) = 0.0_r8
      snowab(i) = 0.0_r8
      cldmax(i) = 0.0_r8
#ifdef DIRIND
      precabx(i) = 0.0_r8
      icnlx(i) = 0.0_r8
!This new initialization of standard variable below 
!does not change the meteorology (checked)
      cmec1(i) = 0.0_r8 !corinna: used for progcdnc
#endif
   end do
!
! initialize multi-level fields 
!
   do k = 1,pver
      do i = 1,ncol
         q(i,k) = qn(i,k) 
         t(i,k) = tn(i,k)
!         q(i,k)=qn(i,k)-qtend(i,k)*deltat
!         t(i,k)=tn(i,k)-ttend(i,k)*deltat
    end do
   end do
   cme     (:ncol,:) = 0._r8
   evapprec(:ncol,:) = 0._r8
   prodprec(:ncol,:) = 0._r8
   evapsnow(:ncol,:) = 0._r8
   prodsnow(:ncol,:) = 0._r8
   evapheat(:ncol,:) = 0._r8
   meltheat(:ncol,:) = 0._r8
   prfzheat(:ncol,:) = 0._r8
   ice2pr(:ncol,:)   = 0._r8
   liq2pr(:ncol,:)   = 0._r8
   liq2snow(:ncol,:) = 0._r8
   fwaut(:ncol,:) = 0._r8
   fsaut(:ncol,:) = 0._r8
   fracw(:ncol,:) = 0._r8
   fsacw(:ncol,:) = 0._r8
   fsaci(:ncol,:) = 0._r8
   rkflxprc(:ncol,:) = 0._r8
   rkflxsnw(:ncol,:) = 0._r8

   pracwo(:ncol,:) = 0._r8
   psacwo(:ncol,:) = 0._r8
   psacio(:ncol,:) = 0._r8

#ifdef DIRIND
   cmex    (:ncol,:) = 0._r8
   prainx  (:ncol,:) = 0._r8
   evapx   (:ncol,:) = 0._r8
   evaprx  (:ncol,:) = 0._r8
   lossx   (:ncol,:) = 0._r8
   nrainx  (:ncol,:) = 0._r8
   selfx   (:ncol,:) = 0._r8
   freez   (:ncol,:) = 0._r8
   supersat(:ncol,:) = 0._r8
   nucrat  (:ncol,:) = 0._r8
   nuclx   (:ncol,:) = 0._r8
   liq2prx (:ncol,:) = 0._r8
   fwautx(:ncol,:) = 0._r8
   fsautx(:ncol,:) = 0._r8
   fracwx(:ncol,:) = 0._r8
   fsacwx(:ncol,:) = 0._r8
   fsacix(:ncol,:) = 0._r8
#endif


!
! find the wet bulb temp and saturation value
! for the provisional t and q without condensation
!
   call findsp (lchnk, ncol, qn, tn, p, tsp, qsp)
   do 800 k = k1mb,pver
      call vqsatd (t(1,k), p(1,k), es, qs, gamma, ncol)
      do i = 1,ncol
         relhum(i) = q(i,k)/qs(i)
!
         cldm(i) = max(cldn(i,k),mincld)
!
! the max cloud fraction above this level
!
         cldmax(i) = max(cldmax(i), cldm(i))

! define the coefficients for C - E calculation

         calpha(i) = 1.0_r8/qs(i)
         cbeta (i) = q(i,k)/qs(i)**2*gamma(i)*cpohl
         cbetah(i) = 1.0_r8/qs(i)*gamma(i)*cpohl
         cgamma(i) = calpha(i)+hlatv*cbeta(i)/cp
         cgamah(i) = calpha(i)+hlatv*cbetah(i)/cp
         rcgama(i) = cgamma(i)/cgamah(i)

#ifdef DIRIND
	 moistrho(i,k) = p(i,k)/(287._r8*(1._r8+0.6_r8*q(i,k))*t(i,k))
!        Converting from #/kg to #/cm3 and calculating in-cloud value
 	 cnlx(i,k) = cnlxin(i,k)*moistrho(i,k)*1.e-6_r8
#endif
         if(cldm(i) > mincld) then
            icwc(i) = max(0._r8,cwat(i,k)/cldm(i))
#ifdef DIRIND
	    icwcx(i) = max(0._r8,cwatx(i,k)/cldm(i))
	    icnlx(i) = max(0._r8,(cnlx(i,k)/cldm(i)))
#endif
         else
            icwc(i) = 0.0_r8
#ifdef DIRIND
            icwcx(i) = 0.0_r8
            icnlx(i) = 0.0_r8
#endif
         endif

!PJR the above logic give zero icwc with nonzero cwat, dont like it!
!PJR generates problems with csigma
!PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
!         icwc(i) = max(1.e-8_r8,cwat(i,k)/cldm(i))

!
! initial guess of evaporation, will be updated within iteration
!
         evapprec(i,k) = conke*(1._r8 - cldm(i))*sqrt(precab(i)) &
                        *(1._r8 - min(relhum(i),1._r8))
#ifdef DIRIND
         evaprx(i,k) = conke*(1._r8 - cldm(i))*sqrt(precabx(i)) &
                        *(1._r8 - min(relhum(i),1._r8))
#endif

!
! zero cmeres before iteration for each level
!
         cmeres(i)=0.0_r8

      end do
      do i = 1,ncol
!
! fractions of ice at this level
!
!!$         tc = t(i,k) - t0
!!$         fice(i,k) = max(0._r8,min(-tc*0.05,1.0_r8))
!
! calculate the cooling due to a phase change of the rainwater
! from above
!
         if (t(i,k) >= t0) then
            meltheat(i,k) =  -hlatf * snowab(i) * gravit/pdel(i,k)
            snowab(i) = 0._r8
         else
            meltheat(i,k) = 0._r8
         endif
      end do

!
! calculate cme and formation of precip. 
!
! The cloud microphysics is highly nonlinear and coupled with cme
! Both rain processes and cme are calculated iteratively.
! 
      do 100 l = 1,iter

        do i = 1,ncol

!
! calculation of cme has 4 scenarios
! ==================================
!
           if(relhum(i) > rhu00(i,k)) then
    
           ! 1. whole grid saturation
           ! ========================
              if(relhum(i) >= 0.999_r8 .or. cldm(i) >= 0.999_r8 ) then
                 cme(i,k)=(calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)
#ifdef DIRIND
		 cmex(i,k)=(calpha(i)*qtend(i,k)-cbetah(i)*ttend(i,k))/cgamah(i)  ! = cme
!ts+ No evaporation of cloud droplets in this case, obviously:
	         evapx(i,k) = 0._r8
#endif

           ! 2. fractional saturation
           ! ========================
              else
                 if (rhdfda(i,k) .eq. 0._r8 .and. icwc(i).eq.0._r8) then
                    write (iulog,*) ' cldwat.F90:  empty rh cloud ', i, k, lchnk
                    write (iulog,*) ' relhum, iter ', relhum(i), l, rhu00(i,k), cldm(i), cldn(i,k)
                    call endrun ()
                 endif
                  csigma(i) = 1.0_r8/(rhdfda(i,k)+cgamma(i)*icwc(i))
                  cmec1(i) = (1.0_r8-cldm(i))*csigma(i)*rhdfda(i,k)
                  cmec2(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0_r8-rcgama(i)*cldm(i))*   &
                             csigma(i)*calpha(i)*icwc(i)
                  cmec3(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
                           (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigma(i)*icwc(i)
                  cmec4(i) = csigma(i)*cgamma(i)*icwc(i)

                  ! Q=C-E=-C1*Al + C2*Aq - C3* At + C4*Er
  
                  cme(i,k) = -cmec1(i)*lctend(i,k) + cmec2(i)*qtend(i,k)  &
                             -cmec3(i)*ttend(i,k) + cmec4(i)*evapprec(i,k)
#ifdef DIRIND
		  csigmax(i) = 1.0_r8/(rhdfda(i,k)+cgamma(i)*icwcx(i))
		  cmec1x(i) = (1.0_r8-cldm(i))*csigmax(i)*rhdfda(i,k)
 		  cmec2x(i) = cldm(i)*calpha(i)/cgamah(i)+(1.0_r8-rcgama(i)*cldm(i))*   &
                             csigmax(i)*calpha(i)*icwcx(i)
 		  cmec3x(i) = cldm(i)*cbetah(i)/cgamah(i) +  &
                           (cbeta(i)-rcgama(i)*cldm(i)*cbetah(i))*csigmax(i)*icwcx(i)
		  cmec4x(i) = csigmax(i)*cgamma(i)*icwcx(i)
 		  !cmex(i,k) = -cmec1x(i)*lctend(i,k) + cmec2x(i)*qtend(i,k)  & 
 		  cmex(i,k) = -cmec1x(i)*lctendx(i,k) + cmec2x(i)*qtend(i,k)  & !corinna: use LWX tendency
                             -cmec3x(i)*ttend(i,k) + cmec4x(i)*evaprx(i,k)
!ts+ Here we evaporate droplets if the cloud is shrinking, otherwise not.
                  if(cldn(i,k) < cldo(i,k)) then
                     evapx(i,k) =( (cldn(i,k)-cldo(i,k))/cldo(i,k))*cnlx(i,k)/deltat
                  end if
#endif
               endif

           ! 3. when rh < rhu00, evaporate existing cloud water
           ! ================================================== 
           else if(cwat(i,k) > 0.0_r8)then
              ! liquid water should be evaporated but not to exceed 
              ! saturation point. if qn > qsp, not to evaporate cwat
              cme(i,k)=-min(max(0._r8,qsp(i,k)-qn(i,k)),cwat(i,k))/deltat 
#ifdef DIRIND
              cmex(i,k)=-min(max(0._r8,qsp(i,k)-qn(i,k)),cwatx(i,k))/deltat
!ts+	Here we evaporate droplets in proportion to the cloud water evaporated
#ifdef AEROFFL
        !corinna: changed to homogeneous mixing assumption, evaporate only if cloud disappears
	!if(cwatx(i,k) > 0._r8 .and. cnlx(i,k) > 1.e-10_r8) then
	!evapx(i,k) = cmex(i,k)*(cnlx(i,k)/cwatx(i,k))
        if(cwatx(i,k) > 0._r8 .and. cmex(i,k)==-cwatx(i,k)/deltat .and. cnlx(i,k) > 1.e-10_r8) then
          evapx(i,k) = -cnlx(i,k)/deltat !remove droplets if all water evaporates
        endif
#else
	!if(cwat(i,k) > 0._r8 .and. cnlx(i,k) > 1.e-10_r8) then
	!evapx(i,k) = cme(i,k)*(cnlx(i,k)/cwat(i,k))
        if(cwat(i,k) > 0._r8 .and. cme(i,k)==-cwat(i,k)/deltat .and. cnlx(i,k) > 1.e-10_r8) then
          evapx(i,k) = -cnlx(i,k)/deltat !remove droplets if all water evaporates
        endif
#endif
#endif

           ! 4. no condensation nor evaporation
           ! ==================================
           else
              cme(i,k)=0.0_r8
#ifdef DIRIND
	      cmex(i,k)=0.0_r8
	      evapx(i,k) = 0.0_r8
#endif
           endif

#ifdef DIRIND
! 	Make sure all droplets are frozen if temp. < -40:
	if(cnlx(i,k) > 0._r8 .and. t(i,k) < 233._r8) then
	  freez(i,k) = - cnlx(i,k)/deltat
	end if
#endif
  
        end do    !end loop for cme update

! Because of the finite time step, 
! place a bound here not to exceed wet bulb point
! and not to evaporate more than available water
!
         do i = 1, ncol
            qtmp = qn(i,k) - cme(i,k)*deltat

! possibilities to have qtmp > qsp
!
!   1. if qn > qs(tn), it condenses; 
!      if after applying cme,  qtmp > qsp,  more condensation is applied. 
!      
!   2. if qn < qs, evaporation should not exceed qsp,
    
            if(qtmp > qsp(i,k)) then
              cme(i,k) = cme(i,k) + (qtmp-qsp(i,k))/deltat
#ifdef DIRIND
	      cmex(i,k) = cmex(i,k) + (qtmp-qsp(i,k))/deltat
#endif
            endif

!
! if net evaporation, it should not exceed available cwat
!
            if(cme(i,k) < -cwat(i,k)/deltat)  &
               cme(i,k) = -cwat(i,k)/deltat
#ifdef DIRIND
            if(cmex(i,k) < -cwatx(i,k)/deltat)  &
               cmex(i,k) = -cwatx(i,k)/deltat
#endif
!
! addition of residual condensation from previous step of iteration
!
            cme(i,k) = cme(i,k) + cmeres(i)
#ifdef DIRIND
            cmex(i,k) = cmex(i,k) + cmeres(i) 
#endif

         end do

         !      limit cme for roundoff errors
         do i = 1, ncol
            cme(i,k) = cme(i,k)*omsm
#ifdef DIRIND
            cmex(i,k) = cmex(i,k)*omsm          
#endif
         end do

         do i = 1,ncol
!
! as a safe limit, condensation should not reduce grid mean rh below rhu00
! 
           if(cme(i,k) > 0.0_r8 .and. relhum(i) > rhu00(i,k) )  &
              cme(i,k) = min(cme(i,k), (qn(i,k)-qs(i)*rhu00(i,k))/deltat)
#ifdef DIRIND
           if(cmex(i,k) > 0.0_r8 .and. relhum(i) > rhu00(i,k) )  &    
              cmex(i,k) = min(cmex(i,k), (qn(i,k)-qs(i)*rhu00(i,k))/deltat)
#endif
!
! initial guess for cwm (mean cloud water over time step) if 1st iteration
!
           if(l < 2) then
             cwm(i) = max(cwat(i,k)+cme(i,k)*dto2,  0._r8)
#ifdef DIRIND
	     cwmx(i) = max(cwatx(i,k)+cmex(i,k)*dto2,  0._r8)
#endif
           endif
#ifdef DIRIND
	   sigw(i) = sigw1(i,k)
#endif
         enddo

! provisional precipitation falling through model layer
         do i = 1,ncol
!!$            prprov(i) =  precab(i) + prodprec(i,k)*pdel(i,k)/gravit
! rain produced in this layer not too effective in collection process
            wtthick = max(0._r8,min(0.5_r8,((zi(i,k)-zi(i,k+1))/1000._r8)**2))
            prprov(i) =  precab(i) + wtthick*prodprec(i,k)*pdel(i,k)/gravit
#ifdef DIRIND
            prprovx(i) =  precabx(i) + wtthick*prainx(i,k)*pdel(i,k)/gravit
#endif
         end do

#ifdef DIRIND
!initialize fields for findmcnew
   cxsind  (:ncol,k) = 0._r8
   ntotal  (:ncol,k) = 0._r8
   n1      (:ncol,k) = 0._r8
   n2      (:ncol,k) = 0._r8
   n3      (:ncol,k) = 0._r8
   n4      (:ncol,k) = 0._r8
   n5      (:ncol,k) = 0._r8
   n6      (:ncol,k) = 0._r8
   n7      (:ncol,k) = 0._r8
   n8      (:ncol,k) = 0._r8
   n9      (:ncol,k) = 0._r8
   n10     (:ncol,k) = 0._r8
   n11     (:ncol,k) = 0._r8
   n12     (:ncol,k) = 0._r8
   n13     (:ncol,k) = 0._r8
   n14     (:ncol,k) = 0._r8
   relh    (:ncol,k) = 0._r8
   cdnc    (:ncol,k) = 0._r8
#endif

! calculate conversion of condensate to precipitation by cloud microphysics 
#ifdef DIRIND
         call findmcnew (lchnk   ,ncol    , &
                         k       ,prprov  ,snowab,  t       ,p        , &
                         cwm     ,cldm    ,cldmax  ,fice(1,k),coef    , &
                         fwaut(1,k),fsaut(1,k),fracw(1,k),fsacw(1,k),fsaci(1,k), &
                         seaicef, snowh, pracwo(1,k), psacwo(1,k), psacio(1,k), &
                         fwautx(1,k),fsautx(1,k),fracwx(1,k),fsacwx(1,k),fsacix(1,k), &
                         qm      ,precc   ,cwmx    ,prprovx  ,coefx   , &
                         volrad  ,capna   ,ntotal(1,k), cxsind(1,k)   , &
                         n1(1,k) ,n2(1,k) ,n3(1,k) ,n4(1,k)  ,n5(1,k) , &
                         n6(1,k) ,n7(1,k) ,n8(1,k) ,n9(1,k)  ,n10(1,k), &
                         n11(1,k),n12(1,k),n13(1,k),n14(1,k) ,nummode,  &
			 relh(1,k), sigw, cldn, cldo, omega, nuclx, icnlx, &     
			 ncoefx, supersat)                                
#else
         call findmcnew (lchnk   ,ncol    , &
                         k       ,prprov  ,snowab,  t       ,p        , &
                         cwm     ,cldm    ,cldmax  ,fice(1,k),coef    , &
                         fwaut(1,k),fsaut(1,k),fracw(1,k),fsacw(1,k),fsaci(1,k), &
                         seaicef, snowh, pracwo(1,k), psacwo(1,k), psacio(1,k))
#endif
!
! calculate the precip rate
!
         error_found = .false.
         do i = 1,ncol

#ifdef DIRIND
    	if(cldo(i,k) < 1.e-4_r8) then 
          reffl=1.e-6_r8*(2.1_r8+66.25_r8*(MAX(10._r8,2*nuclx(i,k))**(-1._r8/3))) !fit through simulated values by Choularton & Bower, QJRMS 1992
          dropmass = 1.333_r8*1000*3.1416_r8*(reffl)**3 !use radius which depends on CCN
          !convert dlf/dropmass from #/kg into #/cm3
          nucrat(i,k) = cldm(i)*(nuclx(i,k)/deltat) &
               + 1.e-6_r8*moistrho(i,k)*(1._r8-cmec1(i))*dlf(i,k)*(1._r8/dropmass)   !reduce by evaporating fraction
        else
          reffl=1.e-6_r8*(2.1_r8+66.25_r8*(MAX(10._r8,2*nuclx(i,k))**(-1._r8/3))) !fit through simulated values by Choularton & Bower, QJRMS 1992
          dropmass = 1.333_r8*1000*3.1416_r8*(reffl)**3 !use radius which depends on CCN
          !convert dlf/dropmass from #/kg into #/cm3
          nucrat(i,k) = cldm(i)*max((nuclx(i,k)-icnlx(i))/deltat,0._r8) &
               + 1.e-6_r8*moistrho(i,k)*(1._r8-cmec1(i))*dlf(i,k)*(1._r8/dropmass)
	end if
#endif


#ifdef DIRIND
            if (cldm(i) > 0._r8) then
#else
            if (cldm(i) > 0) then  
#endif
!
! first predict the cloud water
!
               cdt = coef(i)*deltat
               if (cdt > 0.01_r8) then
                  pol = cme(i,k)/coef(i) ! production over loss
                  cwn(i) = max(0._r8,(cwat(i,k)-pol)*exp(-cdt)+ pol)
               else
                  cwn(i) = max(0._r8,(cwat(i,k) + cme(i,k)*deltat)/(1+cdt))
               endif
#ifdef DIRIND
	       cdtx = coefx(i)*deltat
               if (cdtx > 0.01_r8) then
		 polx = cmex(i,k)/coefx(i) ! production over loss
                 cwnx(i) = max(0._r8,(cwatx(i,k)-polx)*exp(-cdtx)+ polx)
               else
                 cwnx(i) = max(0._r8,(cwatx(i,k) + cmex(i,k)*deltat)/(1+cdtx))
               endif

               !corinna: if only a very small amount of water is left, precipitate all
               !if (cwn(i) /cldm(i) < 1.E-6_r8) cwn(i)  = 0._r8 
               !if (cwnx(i)/cldm(i) < 1.E-6_r8) cwnx(i) = 0._r8
               !end corinna

!TS 	Predict the cloud droplet number concentration
!
	rhocgs =  p(i,k)/(287._r8*t(i,k))*1.e-3_r8 
	! Calculating the selfcollection rate for cloud droplets following Beheng (1994):
!CH        selfx(i,k) = -cldm(i)*1.29e10_r8*(rhocgs*icwcx(i))**2
!CH     introduce limit to available droplets
        selfx(i,k) = -MIN(cnlx(i,k)/deltat,cldm(i)*1.29e10_r8*(rhocgs*icwcx(i))**2)
!CH
	! Adding up the loss rates from selfcollection, evaporation and freezing:
        lossx(i,k) = selfx(i,k) + evapx(i,k) + freez(i,k)
	! Making sure that we don't deplete more droplets than we actually have:
	if(abs(lossx(i,k)) > cnlx(i,k)/deltat .and. cnlx(i,k) > 0._r8) then
	lossx(i,k) = - cnlx(i,k)/deltat
	end if 

	ncdtx = ncoefx(i)*deltat
	if(ncdtx > 0.01_r8) then
	npolx  =  (nucrat(i,k)+lossx(i,k))/ncoefx(i) 
	ncwnx(i) = max(0._r8,(cnlx(i,k)-npolx)*exp(-ncdtx)+ npolx)
	else
	ncwnx(i) = max(0._r8,(cnlx(i,k) + (nucrat(i,k)+lossx(i,k))*deltat)/(1+ncdtx))
	end if
#endif

!
! now back out the tendency of net rain production
!
               prodprec(i,k) =  max(0._r8,cme(i,k)-(cwn(i)-cwat(i,k))/deltat)
#ifdef DIRIND
               prainx(i,k) = max(0._r8,cmex(i,k)-(cwnx(i)-cwatx(i,k))/deltat)
               nrainx(i,k) = max(0._r8,(nucrat(i,k)+lossx(i,k))-(ncwnx(i)-cnlx(i,k))/deltat)
#endif
            else
               prodprec(i,k) = 0.0_r8
               cwn(i) = 0._r8
#ifdef DIRIND
               prainx(i,k) = 0.0_r8
               cwnx(i) = 0._r8
               nrainx(i,k) = 0.0_r8
               ncwnx(i) = 0._r8
#endif
            endif

            ! provisional calculation of conversion terms
            ice2pr(i,k) = prodprec(i,k)*(fsaut(i,k)+fsaci(i,k))
            liq2pr(i,k) = prodprec(i,k)*(fwaut(i,k)+fsacw(i,k)+fracw(i,k))
!old        liq2snow(i,k) = prodprec(i,k)*fsacw(i,k)
#ifdef DIRIND
	    liq2prx(i,k) =  prainx(i,k)*(fwautx(i,k)+fsacwx(i,k)+fracwx(i,k))
#endif

!           revision suggested by Jim McCaa
!           it controls the amount of snow hitting the sfc 
!           by forcing a lot of conversion of cloud liquid to snow phase
!           it might be better done later by an explicit representation of 
!           rain accreting ice (and freezing), or by an explicit freezing of raindrops
            liq2snow(i,k) = max(prodprec(i,k)*fsacw(i,k), fsnow(i,k)*liq2pr(i,k))

            ! bounds
            nice2pr = min(ice2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*fice(i,k)/deltat)
            nliq2pr = min(liq2pr(i,k),(cwat(i,k)+cme(i,k)*deltat)*(1._r8-fice(i,k))/deltat)
#ifdef DIRIND
            nliq2prx = min(liq2prx(i,k),(cwatx(i,k)+cmex(i,k)*deltat)*(1._r8-fice(i,k))/deltat)
#endif
!            write(iulog,*) ' prodprec ', i, k, prodprec(i,k)
!            write(iulog,*) ' nliq2pr, nice2pr ', nliq2pr, nice2pr
            if (liq2pr(i,k).ne.0._r8) then
               nliq2snow = liq2snow(i,k)*nliq2pr/liq2pr(i,k)   ! correction
            else
               nliq2snow = liq2snow(i,k)
            endif

!           avoid roundoff problems generating negatives
            nliq2snow = nliq2snow*omsm
            nliq2pr = nliq2pr*omsm
            nice2pr = nice2pr*omsm
#ifdef DIRIND
            nliq2prx = nliq2prx*omsm
#endif
            
!           final estimates of conversion to precip and snow
            prodprec(i,k) = (nliq2pr + nice2pr)
            prodsnow(i,k) = (nice2pr + nliq2snow)

            rcwn(i,l,k) =  cwat(i,k) + (cme(i,k)-   prodprec(i,k))*deltat
            rliq(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)*(1._r8-fice(i,k)) - nliq2pr * deltat
            rice(i,l,k) = (cwat(i,k) + cme(i,k)*deltat)* fice(i,k)      -    nice2pr                     *deltat

#ifdef DIRIND
	    prainx(i,k) = (nliq2prx + nice2pr)
            rcwnx(i,l,k) =  cwatx(i,k) + (cmex(i,k)-   prainx(i,k))*deltat
#endif

!           Save for sanity check later...  
!           Putting sanity checks inside loops 100 and 800 screws up the 
!           IBM compiler for reasons as yet unknown.  TBH
            cwnsave(i,l,k)      = cwn(i)
            cmesave(i,l,k)      = cme(i,k)
            prodprecsave(i,l,k) = prodprec(i,k)
!           End of save for sanity check later...  

!           final version of condensate to precip terms
            liq2pr(i,k) = nliq2pr
            liq2snow(i,k) = nliq2snow
            ice2pr(i,k) = nice2pr

            cwn(i) = rcwn(i,l,k)

#ifdef DIRIND
            liq2prx(i,k) = nliq2prx
            cwnx(i) = rcwnx(i,l,k)
#endif
!
! update any remaining  provisional values
!
            cwm(i) = (cwn(i) + cwat(i,k))*0.5_r8
#ifdef DIRIND
            cwmx(i) = (cwnx(i) + cwatx(i,k))*0.5_r8      !cak
            if(cwmx(i) < 0.) cwmx(i) = 0._r8             !cak
#endif
!
! update in cloud water
!
            if(cldm(i) > mincld) then
               icwc(i) = cwm(i)/cldm(i)
#ifdef DIRIND
               icwcx(i) = cwmx(i)/cldm(i)                !cak
#endif
            else
               icwc(i) = 0.0_r8
#ifdef DIRIND
               icwcx(i) = 0.0_r8                         !cak
#endif
            endif
!PJR the above logic give zero icwc with nonzero cwat, dont like it!
!PJR generates problems with csigma
!PJR set the icwc to a very small number, so we can start from zero cloud cover and make some clouds
!         icwc(i) = max(1.e-8_r8,cwm(i)/cldm(i))

         end do              ! end of do i = 1,ncol

#ifdef DIRIND
      do i = 1,ncol
!       still keeping the old rer3fac parameterization, 
!       for initialization (otherwise: occational division by 0)
        if (landm(i,lchnk).eq.1._r8) then
          rer3fac = 0.87503401_r8
        else
          rer3fac = 0.92831777_r8
        endif
       !corinna: new beta from Rotstayn & Liu (GRL 2009), eq. 2 (OLDBETA)
       if (cwmx(i)*(1.-fice(i,k)).gt.1.e-10_r8.and. &
          cldm(i).gt.mincld.and.ncwnx(i)/cldm(i).ge.1._r8) then
            epsbeta=1._r8-0.7_r8*exp(-0.003_r8*ncwnx(i)/cldm(i))
            rer3fac=((1._r8+2._r8*epsbeta**2)**(2._r8/3))/((1._r8+epsbeta**2)**(1._r8/3))
            rer3fac=1._r8/rer3fac
       endif
       !end corinna
          !corinna: use update mass & number to calculate relca
          volrad(i,k) = 1.e6_r8*1._r8/(1.333_r8*pi)**0.333_r8 * 0.01_r8* &
               ((rhocgs/rhow)*cwmx(i)*(1._r8-fice(i,k))/(ncwnx(i)+1._r8))**0.333_r8
	  relca(i,k) = volrad(i,k)/rer3fac
          relca(i,k) = min(max(relca(i,k),4.0_r8),20._r8)
!         CDNC (only used) for output in history file:  
 	  cdnc(i,k)  = ncwnx(i)/MAX(mincld,cldm(i))  !corinna: diagnose cdnc after microphysics
      end do
#endif

!
! calculate provisional value of cloud water for
! evaporation of precipitate (evapprec) calculation
!
      do i = 1,ncol
         qtmp = qn(i,k) - cme(i,k)*deltat
         ttmp = tn(i,k) + deltat/cp * ( meltheat(i,k)       &
              + (hlatv + hlatf*fice(i,k)) * cme(i,k) )
         esn = estblf(ttmp)
         qsn = min(epsqs*esn/(p(i,k) - omeps*esn),1._r8)
         qtl(i) = max((qsn - qtmp)/deltat,0._r8)
         relhum1(i) = qtmp/qsn
      end do
!
      do i = 1,ncol
#ifdef PERGRO
         evapprec(i,k) = conke*(1._r8 - max(cldm(i),mincld))* &
                         sqrt(precab(i))*(1._r8 - min(relhum1(i),1._r8))
#ifdef DIRIND
 	 evaprx(i,k) = conke*(1._r8 - max(cldm(i),mincld))* &
                      sqrt(precabx(i))*(1._r8 - min(relhum1(i),1._r8))
#endif
#else
         evapprec(i,k) = conke*(1._r8 - cldm(i))*sqrt(precab(i)) &
                         *(1._r8 - min(relhum1(i),1._r8))
#ifdef DIRIND
         evaprx(i,k) = conke*(1._r8 - cldm(i))*sqrt(precabx(i)) &
                      *(1._r8 - min(relhum1(i),1._r8))
#endif
#endif
!
! limit the evaporation to the amount which is entering the box
! or saturates the box
!
         prtmp = precab(i)*gravit/pdel(i,k)
         evapprec(i,k) = min(evapprec(i,k), prtmp, qtl(i))*omsm
#ifdef DIRIND
         prtmpx = precabx(i)*gravit/pdel(i,k)
 	 evaprx(i,k) = min(evaprx(i,k), prtmpx, qtl(i))*omsm
#endif
#ifdef PERGRO
!           zeroing needed for pert growth
         evapprec(i,k) = 0._r8
#ifdef DIRIND
 	 evaprx(i,k) = 0._r8
#endif
#endif
!
! Partition evaporation of precipitate between rain and snow using
! the fraction of snow falling into the box. Determine the heating
! due to evaporation. Note that evaporation is positive (loss of precip,
! gain of vapor) and that heating is negative.
         if (evapprec(i,k) > 0._r8) then
            evapsnow(i,k) = evapprec(i,k) * snowab(i) / precab(i)
            evapheat(i,k) = -hlatv * evapprec(i,k) - hlatf * evapsnow(i,k)
         else 
            evapsnow(i,k) = 0._r8
            evapheat(i,k) = 0._r8
         end if
! Account for the latent heat of fusion for liquid drops collected by falling snow
         prfzheat(i,k) = hlatf * liq2snow(i,k)
      end do

! now remove the residual of any over-saturation. Normally,
! the oversaturated water vapor should have been removed by 
! cme formulation plus constraints by wet bulb tsp/qsp
! as computed above. However, because of non-linearity,
! addition of (cme-evapprec) to update t and q may still cause
! a very small amount of over saturation. It is called a
! residual of over-saturation because theoretically, cme
! should have taken care of all of large scale condensation.
! 

       do i = 1,ncol
          qtmp = qn(i,k)-(cme(i,k)-evapprec(i,k))*deltat
          ttmp = tn(i,k) + deltat/cp * ( meltheat(i,k) + evapheat(i,k) + prfzheat(i,k)      &
              + (hlatv + hlatf*fice(i,k)) * cme(i,k) )
          esn = estblf(ttmp)
          qsn = min(epsqs*esn/(p(i,k) - omeps*esn),1._r8)
          !
          !Upper stratosphere and mesosphere, qsn calculated
          !above may be negative. Here just to skip it instead
          !of resetting it to 1 as in aqsat
          !
          if(qtmp > qsn .and. qsn > 0) then
             !calculate dqsdt, a more precise calculation
             !which taking into account different range of T 
             !can be found in aqsatd.F. Here follows
             !cond.F to calculate it.
             !
             denom = (p(i,k)-omeps*esn)*ttmp*ttmp
             dqsdt = clrh2o*qsn*p(i,k)/denom
             !
             !now extra condensation to bring air to just saturation
             !
             ctmp = (qtmp-qsn)/(1._r8+hlocp*dqsdt)/deltat
             cme(i,k) = cme(i,k)+ctmp
!
! save residual on cmeres to addtion to cme on entering next iteration
! cme exit here contain the residual but overrided if back to iteration
!
             cmeres(i) = ctmp
          else
             cmeres(i) = 0.0_r8
          endif
       end do
              
 100 continue              ! end of do l = 1,iter

!
! precipitation
!
      do i = 1,ncol
         precip(i) = precip(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
         precab(i) = precab(i) + pdel(i,k)/gravit * (prodprec(i,k) - evapprec(i,k))
         if(precab(i).lt.0._r8) precab(i)=0._r8
!         snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodprec(i,k)*fice(i,k) - evapsnow(i,k))
         snowab(i) = snowab(i) + pdel(i,k)/gravit * (prodsnow(i,k) - evapsnow(i,k))
#ifdef DIRIND
         prtmpx = pdel(i,k)/gravit * (prainx(i,k) - evaprx(i,k))
         precabx(i) = precabx(i) + prtmpx
         if(precabx(i).lt.0._r8) precabx(i)=0._r8
#endif
         ! If temperature above freezing, all precip is rain flux.  if temperature below freezing, all precip is snow flux.
         rkflxprc(i,k+1) = precab(i)   !! making this consistent with other precip fluxes.  prc = rain + snow
#ifdef DIRIND
#ifndef AEROFFL
         rkflxprc(i,k+1) = precabx(i)   !! making this consistent with other precip fluxes.  prc = rain + snow
#endif
#endif
         !!rkflxprc(i,k+1) = precab(i) - snowab(i)
         rkflxsnw(i,k+1) = snowab(i)

!!$         if ((precab(i)) < 1.e-10_r8) then      
!!$            precab(i) = 0._r8
!!$            snowab(i) = 0._r8
!!$         endif
      end do
 800 continue                ! level loop (k=1,pver)

! begin sanity checks
   error_found = .false.
   do k = k1mb,pver
      do l = 1,iter
         do i = 1,ncol
	    if (abs(rcwn(i,l,k)).lt.1.e-300_r8) rcwn(i,l,k) = 0._r8
	    if (abs(rliq(i,l,k)).lt.1.e-300_r8) rliq(i,l,k) = 0._r8
	    if (abs(rice(i,l,k)).lt.1.e-300_r8) rice(i,l,k) = 0._r8
            if (rcwn(i,l,k).lt.0._r8) error_found = .true.
            if (rliq(i,l,k).lt.0._r8) error_found = .true.
            if (rice(i,l,k).lt.0._r8) error_found = .true.
         enddo
      enddo
   enddo
   if (error_found) then
      do k = k1mb,pver
         do l = 1,iter
            do i = 1,ncol
               if (rcwn(i,l,k).lt.0._r8) then
                  write(iulog,*) ' prob with neg rcwn1 ', rcwn(i,l,k),  &
                     cwnsave(i,l,k)
                  write(iulog,*) ' cwat, cme*deltat, prodprec*deltat ', &
                     cwat(i,k), cmesave(i,l,k)*deltat,               &
                     prodprecsave(i,l,k)*deltat,                     &
                     (cmesave(i,l,k)-prodprecsave(i,l,k))*deltat
                  call endrun('PCOND')
               endif
               if (rliq(i,l,k).lt.0._r8) then
                  write(iulog,*) ' prob with neg rliq1 ', rliq(i,l,k)
                  call endrun('PCOND')
               endif
               if (rice(i,l,k).lt.0._r8) then
                  write(iulog,*) ' prob with neg rice ', rice(i,l,k)
                  call endrun('PCOND')
               endif
            enddo
         enddo
      enddo
   end if
! end sanity checks

   return
end subroutine pcond

!##############################################################################

#ifdef DIRIND
subroutine findmcnew (lchnk   ,ncol    , &
                      k       ,precab  ,snowab,  t       ,p       , &
                      cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
                      fwaut   ,fsaut   ,fracw   ,fsacw   ,fsaci   , &
                      seaicef ,snowh,  pracwo, psacwo, psacio     , &
                      fwautx  ,fsautx  ,fracwx  ,fsacwx  ,fsacix  , &
                      qm      ,precc   , cwmx   ,precabx ,coefx   , &
                      volrad  ,capnx   , Ntot   ,cxstot  ,n1,n2,n3, &
                      n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,nummode, &
		      relh, sigw1, cldn, cldo, omega, nuclx, icnlx,  & 
		      ncoefx, supersat)
#else
subroutine findmcnew (lchnk   ,ncol    , &
                      k       ,precab  ,snowab,  t       ,p       , &
                      cwm     ,cldm    ,cldmax  ,fice    ,coef    , &
                      fwaut   ,fsaut   ,fracw   ,fsacw   ,fsaci   , &
                      seaicef ,snowh,  pracwo, psacwo, psacio )
#endif
 
!----------------------------------------------------------------------- 
! 
! Purpose: 
! calculate the conversion of condensate to precipitate
! 
! Method: 
! See: Rasch, P. J, and J. E. Kristjansson, A Comparison of the CCM3
!  model climate using diagnosed and 
!  predicted condensate parameterizations, 1998, J. Clim., 11,
!  pp1587---1614.
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
   use phys_grid, only: get_rlat_all_p
   use comsrf,        only: landm
#ifdef DIRIND
   use opttab, only: nmodes, nbmodes
   use constituents, only: pcnst, cnst_get_ind
   use const, only: Mso4, Msv, Ms
   use aerosoldef
#endif
!
! input args
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: k                     ! level index

   real(r8), intent(in) :: precab(pcols)        ! rate of precipitation from above (kg / (m**2 * s))
   real(r8), intent(in) :: t(pcols,pver)        ! temperature       (K)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure          (Pa)
   real(r8), intent(in) :: cldm(pcols)          ! cloud fraction
   real(r8), intent(in) :: cldmax(pcols)        ! max cloud fraction above this level
   real(r8), intent(in) :: cwm(pcols)           ! condensate mixing ratio (kg/kg)
   real(r8), intent(in) :: fice(pcols)          ! fraction of cwat that is ice
   real(r8), intent(in) :: seaicef(pcols)       ! sea ice fraction 
   real(r8), intent(in) :: snowab(pcols)        ! rate of snow from above (kg / (m**2 * s))
   real(r8), intent(in) :: snowh(pcols)         ! Snow depth over land, water equivalent (m)
#ifdef DIRIND
   real(r8), intent(in) :: cldo(pcols,pver)     ! TS Old cloud fraction
   real(r8), intent(in) :: cldn(pcols,pver)     ! TS New cloud fraction
   real(r8), intent(in) :: omega(pcols,pver)    ! TS
   real(r8), intent(in) :: qm(pcols,pver,pcnst) ! Common aerosol (only!) tracers for indirect and direct calc. 
   real(r8), intent(in) :: precc(pcols)         ! convective-scale precipitation rate
   real(r8), intent(in) :: cwmx(pcols)          ! condensate mixing ratio (kg/kg)
   real(r8), intent(in) :: precabx(pcols)       ! rate of precipitation from above (kg / (m**2 * s))
   real(r8), intent(in) :: icnlx(pcols)         ! In-cloud droplet number concentration (cm-3)
   real(r8), intent(in) :: sigw1(pcols)	        ! standard deviation, vertical velocity distribution
#endif

! output arguments
   real(r8), intent(out) :: coef(pcols)          ! conversion rate (1/s)
   real(r8), intent(out) :: fwaut(pcols)         ! relative importance of liquid autoconversion (a diagnostic)
   real(r8), intent(out) :: fsaut(pcols)         ! relative importance of ice autoconversion    (a diagnostic)
   real(r8), intent(out) :: fracw(pcols)         ! relative importance of rain accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsacw(pcols)         ! relative importance of snow accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsaci(pcols)         ! relative importance of snow accreting ice    (a diagnostic)
   real(r8), intent(out) :: pracwo(pcols)        ! accretion of cloud water by rain (1/s)
   real(r8), intent(out) :: psacwo(pcols)        ! accretion of cloud water by snow (1/s)
   real(r8), intent(out) :: psacio(pcols)        ! accretion of cloud ice by snow (1/s)
#ifdef DIRIND
   real(r8), intent(out) :: Ntot(pcols)          ! total aerosol number concentration 
   real(r8), intent(out) :: n1(pcols)
   real(r8), intent(out) :: n2(pcols)
   real(r8), intent(out) :: n3(pcols)
   real(r8), intent(out) :: n4(pcols)
   real(r8), intent(out) :: n5(pcols)
   real(r8), intent(out) :: n6(pcols)
   real(r8), intent(out) :: n7(pcols)
   real(r8), intent(out) :: n8(pcols)
   real(r8), intent(out) :: n9(pcols)
   real(r8), intent(out) :: n10(pcols)
   real(r8), intent(out) :: n11(pcols)
   real(r8), intent(out) :: n12(pcols)
   real(r8), intent(out) :: n13(pcols)
   real(r8), intent(out) :: n14(pcols)
   real(r8), intent(out) :: cxstot(pcols)        ! total excess mass concentration wrt. internal mixing
   real(r8), intent(out) :: nummode(pcols,pcnst) ! cos
   real(r8), intent(out) :: coefx(pcols)         ! conversion rate (1/s)
   real(r8), intent(out) :: ncoefx(pcols)        ! conversion rate (1/s) (Droplets)
   real(r8), intent(out) :: capnx(pcols)         ! local cloud particle number concentration (cm-3)
   real(r8), intent(out) :: volrad(pcols,pver)   ! Mean volume radius (um)
   real(r8), intent(out) :: fwautx(pcols)        ! relative importance of liquid autoconversion (a diagnostic)
   real(r8), intent(out) :: fsautx(pcols)        ! relative importance of ice autoconversion    (a diagnostic)
   real(r8), intent(out) :: fracwx(pcols)        ! relative importance of rain accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsacwx(pcols)        ! relative importance of snow accreting liquid (a diagnostic)
   real(r8), intent(out) :: fsacix(pcols)        ! relative importance of snow accreting ice    (a diagnostic)
#endif

! work variables

   integer i
   integer ii
   integer ind(pcols)
   integer ncols
#ifdef DIRIND
   integer m
   integer n
   integer kcomp               ! aerosol mode index (1 to nmodes=14)
#endif

   real(r8), parameter :: degrad = 57.296_r8 ! divide by this to convert degrees to radians
   real(r8) capn                 ! local cloud particles / cm3
   real(r8) capnoice             ! local cloud particles when not over sea ice / cm3
   real(r8) ciaut                ! coefficient of autoconversion of ice (1/s)
   real(r8) cldloc(pcols)        ! non-zero amount of cloud
   real(r8) cldpr(pcols)         ! assumed cloudy volume occupied by rain and cloud
   real(r8) con1                 ! work constant
   real(r8) con2                 ! work constant
   real(r8) csacx                ! constant used for snow accreting liquid or ice
!!$   real(r8) dtice                ! interval for transition from liquid to ice
   real(r8) icemr(pcols)         ! in-cloud ice mixing ratio
   real(r8) icrit                ! threshold for autoconversion of ice
   real(r8) liqmr(pcols)         ! in-cloud liquid water mixing ratio
   real(r8) pracw                ! rate of rain accreting water
   real(r8) prlloc(pcols)        ! local rain flux in mm/day
   real(r8) prscgs(pcols)        ! local snow amount in cgs units
   real(r8) psaci                ! rate of collection of ice by snow (lin et al 1983)
   real(r8) psacw                ! rate of collection of liquid by snow (lin et al 1983)
   real(r8) psaut                ! rate of autoconversion of ice condensate
   real(r8) ptot                 ! total rate of conversion
   real(r8) pwaut                ! rate of autoconversion of liquid condensate
   real(r8) r3l                  ! volume radius
   real(r8) rainmr(pcols)        ! in-cloud rain mixing ratio
   real(r8) rat1                 ! work constant
   real(r8) rat2                 ! work constant
!!$   real(r8) rdtice               ! recipricol of dtice
   real(r8) rho(pcols)           ! density (mks units)
   real(r8) rhocgs               ! density (cgs units)
   real(r8) rlat(pcols)          ! latitude (radians)
   real(r8) snowfr               ! fraction of precipate existing as snow
   real(r8) totmr(pcols)         ! in-cloud total condensate mixing ratio
   real(r8) vfallw               ! fall speed of precipitate as liquid
   real(r8) wp                   ! weight factor used in calculating pressure dep of autoconversion
   real(r8) wsi                  ! weight factor for sea ice
   real(r8) wt                   ! fraction of ice
   real(r8) wland                ! fraction of land

#ifdef DIRIND
!  local variables
   real(r8) :: liqmrx(pcols)       ! in-cloud liquid water mixing ratio
   real(r8) :: totmrx(pcols)       ! in-cloud total condensate mixing ratio
   real(r8) :: icemrx(pcols)       ! in-cloud ice mixing ratio
   real(r8) :: rat2x               ! work constant
   real(r8) :: con2x               ! work constant
   real(r8) :: r3lx                ! Mean volume radius (m)
   real(r8) :: pwautx              ! rate of autoconversion of liquid condensate
   real(r8) :: psautx              ! rate of autoconversion of ice condensate
   real(r8) :: psacwx              ! rate of collection of liquid by snow (lin et al 1983)
   real(r8) :: prllocx(pcols)      ! local rain flux in mm/day
   real(r8) :: pracwx              ! rate of rain accreting water
   real(r8) :: rainmrx(pcols)      ! in-cloud rain mixing ratio
   real(r8) :: psacix              ! rate of collection of ice by snow (lin et al 1983)
   real(r8) :: ptotx               ! total rate of conversion
   real(r8) :: rtotx               ! total rate of droplet conversion
   real(r8) :: Cnso4(pcols)	   ! so4(n) (ug/m3)
   real(r8) :: Caitso4(pcols)	   ! so4(ait) (ug/m3)
   real(r8) :: Caitbc(pcols)	   ! BC(ait) (ug/m3)
   real(r8) :: Caitocbc(pcols)	   ! OC+BC(ait) (ug/m3)
   real(r8) :: rhoocbc(pcols)	   ! OC+BC rho
   real(r8) :: Cabce(pcols)	   ! BC(ax) (ug/m3)	    
   real(r8) :: Cnbc(pcols)	   ! BC(n) (ug/m3)	
   real(r8) :: Cnoc(pcols)         ! OC(n) (ug/m3)
   real(r8) :: Cnocibc(pcols)      ! OC(ni) (ug/m3)
   real(r8) :: Cnbcioc(pcols)      ! BC(ni) (ug/m3)
   real(r8) :: rhobcocn(pcols)     ! Density BC/OC(n) (ug/m3)
   real(r8) :: Catot(pcols)        ! Total added so4, BC and OC (ug/m3)
   real(r8) :: Cas75(pcols)        ! so4(Ait75) (ug/m3)
   real(r8) :: C_dst2(pcols)       ! Dust mode 2 (ug/m3)
   real(r8) :: C_dst3(pcols)       ! Dust mode 3 (ug/m3)
   real(r8) :: C_ss1(pcols)        ! Sea salt mode 1 (ug/m3)
   real(r8) :: C_ss2(pcols)        ! Sea salt 2 (ug/m3)
   real(r8) :: C_ss3(pcols)        ! Sea salt 3 (ug/m3)
   real(r8) :: fnbc(pcols)         ! = Cbc/(Cbc+Coc) for BC&OC(n)
   real(r8) :: faitbc(pcols)       ! = Cbc/(Cbc+Coc) for BC&OC(Ait)
   real(r8) :: f_c(pcols)          ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
   real(r8) :: f_bc(pcols)         ! = Cbc/(Cbc+Coc)
   real(r8) :: f_aq(pcols)         ! = Cso4a2/(Cso4a1+Cso4a2)
   real(r8) :: Nnatk(pcols,0:nmodes) ! modal number concentration
   real(r8) :: lnsig(pcols,nmodes) ! ln of standard deviation, lognormal fits, method 1.
   real(r8) :: relh(pcols)	   ! relative humidity at supersaturation
   real(r8) :: rhoall(pcols)       ! Densitiy(mks) in all points   ! cos
   real(r8) :: lthick(pcols,nbmodes)   ! coating layer thickness
   real(r8) :: fso4coat(pcols,nbmodes) ! SO4-bound mass fraction in the coating
   real(r8) nuclx(pcols,pver)	   ! Droplet concentration activated over the timestep
   real(r8) supersat(pcols,pver)   ! Supersaturation as calculated by ARG activation scheme
   real(r8) ::	cam(pcols,nbmodes)	! TS: Internally mixed material 
   real(r8) ::	facm(pcols,nbmodes)	! TS: Carbonaceous fraction of cam
   real(r8) ::	fbcm(pcols,nbmodes) 	! TS: BC fraction of cam
   real(r8) ::  faqm(pcols,nbmodes)	! TS: Fraction of so4 from in-cloud oxidation
   real(r8) ::  Massratio(pcols,nbmodes)! TS: Ratio of modified to original mass, Lognormal fit method
   real(r8) ::  logsig(pcols,nbmodes)	! TS: Log (log10) of standard deviation for lognormal modes, method 
!TS++ Variables for ARG activation scheme
   logical top, activ, version2	   ! top= .true. --> cloudtop, activ=.true. --> activate will be called 
				   ! version2= .true. --> version 2 of lognormal fit calculations used
   real(r8) wminf	           ! Lower bound, vertical velocity
   real(r8) wmaxf	           ! Upper bound, vertical velocity
   real(r8) wdiab		   ! Diabatic velocity
   real(r8) wbar		   ! Grid box mean vertical velocity
   real(r8) smax		   ! Max. supersaturation
   real(r8) rhoair, tair	   ! Density air, Temperature air
   real(r8) lnsigman(12)           ! ln of standard deviation for lognormal modes
   real(r8) na(12)                 ! Number concentration pr mode, cm-3
   real(r8) sigw		   ! Standard deviation, vertical velocity
   real(r8) ma(5,12)		   ! Mass of each aerosol type within given mode
   real(r8) rhodry(5,12)	   ! Dry aerosol density
   real(r8) solubl(5,12)	   ! Solubility of aerosol type
   real(r8) nu(5,12)	
   real(r8) phi(5,12)	
   real(r8) mwaer(5,12)	  	   ! Molecular weight
   real(r8) fn(12)		   ! Number fraction activated
   real(r8) fs(12)		   ! Surface fraction activated
   real(r8) fm(12)		   ! Mass fraction activated
   real(r8) fluxn(12)		   ! Number flux (not used)
   real(r8) fluxs(12)		   ! Surface flux (not used)
   real(r8) fluxm(12)		   ! Mass flux (not used)
   real(r8) rhocomp(5)
   real(r8) phicomp(5)		
   real(r8) mwaercomp(5)	
   real(r8) solubcomp(5)	
   integer  nucomp(5)		
   integer ptype,ntype(12),pmode,nmode ! Number of aerosol modes and types
!TS--
#endif

!      real(r8) csaci
!      real(r8) csacw
!      real(r8) cwaut
!      real(r8) efact
!      real(r8) lamdas
!      real(r8) lcrit
!      real(r8) rcwm
!      real(r8) r3lc2
!      real(r8) snowmr(pcols)
!      real(r8) vfalls

   real(8) ftot

!     inline statement functions
   real(r8) heavy, heavym, a1, a2, heavyp, heavymp
   heavy(a1,a2) = max(0._r8,sign(1._r8,a1-a2))  ! heavyside function
   heavym(a1,a2) = max(0.01_r8,sign(1._r8,a1-a2))  ! modified heavyside function
!
! New heavyside functions to perhaps address error growth problems
!
   heavyp(a1,a2) = a1/(a2+a1+1.e-36_r8)
   heavymp(a1,a2) = (a1+0.01_r8*a2)/(a2+a1+1.e-36_r8)

#ifdef DIRIND
!Call ARG activation scheme
! 	activ = .false.
	activ = .true.
!old	version2 = .false. ! old mass-conserving tables for aerosol size distributions   !cak_test
	version2 = .true.  ! new tables with approximate mass and NUMBER conservation
!
! Note: rhocomp and solubcomp for SS originally 1900 and 1, but the product of the two are unchanged
! Dry densities for SO4, BC, POM, Sea salt & dust:
   data rhocomp / 1769._r8, 2000._r8, 1500._r8, 2200._r8, 2600._r8 / 
! Number of ions the salt disassociates into. Fictive values for BC,OC and dust:
   data nucomp / 3, 1, 1, 2, 1 / 
! Osmotic coefficients. Fictive for BC, OC and dust:
   data phicomp / 0.7_r8, 0.093_r8, 0.3_r8, 1.0_r8, 0.5_r8/ 
! Molecular weights of aerosol. Fictive values for BC, OC and dust:
   data mwaercomp / 132._r8, 3750._r8, 11.57_r8, 59._r8, 19.93_r8 /
!corinna (back to Storelvmo et al 1996):  
   data solubcomp /1._r8, 1.e-3_r8, 0.2_r8, 0.864_r8, 0.013_r8 / 
#endif

!
! find all the points where we need to do the microphysics
! and set the output variables to zero
!
   ncols = 0
   do i = 1,ncol
      coef(i) = 0._r8
      fwaut(i) = 0._r8
      fsaut(i) = 0._r8
      fracw(i) = 0._r8
      fsacw(i) = 0._r8
      fsaci(i) = 0._r8
      liqmr(i) = 0._r8
      rainmr(i) = 0._r8
      if (cwm(i) > 1.e-20_r8) then
         ncols = ncols + 1
         ind(ncols) = i
      endif
#ifdef DIRIND
      coefx(i) = 0._r8
      ncoefx(i) = 0._r8
      fwautx(i)= 0._r8
      fsautx(i)= 0._r8
      fracwx(i)= 0._r8
      fsacwx(i)= 0._r8
      fsacix(i)= 0._r8
      liqmrx(i)= 0._r8      
      volrad(i,k) = 0._r8  
!       implemented 30 June 1999 by JEK. Later adapted to new aerosol
!       schemes by Alf Kirkevag, latest in August 2005. Used only when
!       DIAGNCDNC=.true. 
        if(precc(i) > 1.1574e-8_r8) then   ! convective clouds
          relh(i) = landm(i,lchnk)*1.008_r8+(1.-landm(i,lchnk))*1.0025_r8
	else                            ! stratiform clouds
          relh(i) = 1.001_r8
	end if
#endif
   end do

!cdir nodep
!DIR$ CONCURRENT
   do ii = 1,ncols
      i = ind(ii)
!
! the local cloudiness at this level
!
      cldloc(i) = max(cldmin,cldm(i))
!
! a weighted mean between max cloudiness above, and this layer
!
      cldpr(i) = max(cldmin,(cldmax(i)+cldm(i))*0.5_r8)
!
! decompose the suspended condensate into
! an incloud liquid and ice phase component
!
      totmr(i) = cwm(i)/cldloc(i)
      icemr(i) = totmr(i)*fice(i)
      liqmr(i) = totmr(i)*(1._r8-fice(i))
#ifdef DIRIND
      totmrx(i)= cwmx(i)/cldloc(i)
      liqmrx(i)= totmrx(i)*(1._r8-fice(i))
      icemrx(i)= totmrx(i)*fice(i)
#endif
!
! density
!
      rho(i) = p(i,k)/(287._r8*t(i,k))
      rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
!
! decompose the precipitate into a liquid and ice phase
!
      if (t(i,k) > t0) then
         vfallw = convfw/sqrt(rho(i))
         rainmr(i) = precab(i)/(rho(i)*vfallw*cldpr(i))
#ifdef DIRIND
         rainmrx(i)= precabx(i)/(rho(i)*vfallw*cldpr(i))
#endif
         snowfr = 0
!        snowmr(i)
      else
         snowfr = 1
         rainmr(i) = 0._r8
#ifdef DIRIND
         rainmrx(i)= 0._r8
#endif
      endif
!     rainmr(i) = (precab(i)-snowab(i))/(rho(i)*vfallw*cldpr(i))
!
! local snow amount in cgs units
!
      prscgs(i) = precab(i)/cldpr(i)*0.1_r8*snowfr
#ifdef DIRIND
!ikke med i r128!?   prscgs(i) = precabx(i)/cldpr(i)*0.1_r8*snowfr
#endif
!     prscgs(i) = snowab(i)/cldpr(i)*0.1_r8
!
! local rain amount in mm/day
!
      prlloc(i) = precab(i)*86400._r8/cldpr(i)
#ifdef DIRIND
      prllocx(i)= precabx(i)*86400._r8/cldpr(i)
#endif
   end do  ! ii->i

   con1 = 1._r8/(1.333_r8*pi)**0.333_r8 * 0.01_r8 ! meters

#ifdef DIRIND  ! long ifdef here... --> separate subroutine instead? 
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   do n=1,12
      lnsigman(n)=0._r8
      ntype(n)=0
      fluxn(n)=0._r8
      fluxs(n)=0._r8
      fluxm(n)=0._r8
      fn(n)=0._r8
      fs(n)=0._r8
      fm(n)=0._r8
      na(n)=0._r8
      do m=1,5
         mwaer(m,n)=0._r8
         nu(m,n)=0._r8
         phi(m,n)=0._r8	
         solubl(m,n)=0._r8
         rhodry(m,n)=0._r8
         ma(m,n)=0._r8
      end do
   end do

   do i=1,ncol
      rhoall(i)= p(i,k)/(287._r8*t(i,k))  ! cos
      Ntot(i)  = 0.0_r8
      n1(i)    = 0.0_r8
      n2(i)    = 0.0_r8
      n3(i)    = 0.0_r8
      n4(i)    = 0.0_r8
      n5(i)    = 0.0_r8
      n6(i)    = 0.0_r8
      n7(i)    = 0.0_r8
      n8(i)    = 0.0_r8
      n9(i)    = 0.0_r8
      n10(i)   = 0.0_r8
      n11(i)   = 0.0_r8
      n12(i)   = 0.0_r8
      n13(i)   = 0.0_r8
      n14(i)   = 0.0_r8
      capnx(i) = 0.0_r8 
      cxstot(i)= 0.0_r8
      nuclx(i,k)= 0.0_r8
      do kcomp = 1 , nmodes
	 lnsig(i,kcomp) = 0.0_r8
      end do
   enddo ! i
   nummode(:,:)= 0._r8 ! cos

!     Convert aerosol mass mixing ratios to number concentrations
!cos  Replaced rho with rhoall 

      call convaer(lchnk,ncol,k,rhoall,qm,Cnso4,Cas75,Cnbc,Cnoc, &
                   Cabce,Catot,f_c,f_bc,f_aq,Nnatk,fnbc,faitbc, &
	           Caitso4, Caitbc, Caitocbc, rhoocbc, C_dst2, &
		   C_dst3, C_ss1, C_ss2, C_ss3, rhobcocn, Cnbcioc, Cnocibc)

!     Find modified Nnatk due to lumping of excess internally 
!     mixed mass into external Aitken modes, and then find CCN 
!     concentrations by use of look-up tables  

      call parmix_progncdnc(lchnk,ncol,ncols,ind,k,Nnatk,Cnso4,Cnbc,Cnoc,Cabce,&
         Catot,f_c,f_bc,f_aq,fnbc,faitbc,relh,capnx,cxstot, &
         n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14, &
         Caitso4, Caitbc, Caitocbc, rhoocbc, Cas75, C_dst2, C_dst3, &
         C_ss1, C_ss2, C_ss3, cam, facm, fbcm, faqm, lnsig, Massratio, &
         logsig, version2, lthick, fso4coat)

!     Total aerosol number concentration. Cos: added nummode
      do i=1,ncol
        do kcomp = 0 , nmodes
          Ntot(i)=Ntot(i)+Nnatk(i,kcomp)
        enddo 
          nummode(i,l_dst_a2)=Nnatk(i,6)
          nummode(i,l_dst_a3)=Nnatk(i,7)
          nummode(i,l_ss_a1)=Nnatk(i,8) 
          nummode(i,l_ss_a2)=Nnatk(i,9)
          nummode(i,l_ss_a3)=Nnatk(i,10)
      end do

!TS++ Calling aerosol activation code by Abdul-Razzak and Ghan
   do ii = 1,ncols
      i = ind(ii)
	if(t(i,k) > 233.16_r8) then
	! Test wether the cloud actually contains liquid water here!!
	nmode = 0
        !set to 0 again here
	do n=1,12
	 lnsigman(n)=0._r8
	 ntype(n)=0
	 fluxn(n)=0._r8
	 fluxs(n)=0._r8
	 fluxm(n)=0._r8
	 fn(n)=0._r8
	 fs(n)=0._r8
 	 fm(n)=0._r8
  	 na(n)=0._r8
	 do m=1,5
	  mwaer(m,n)=0._r8
	  nu(m,n)=0._r8
	  phi(m,n)=0._r8	
	  solubl(m,n)=0._r8
	  rhodry(m,n)=0._r8
	  ma(m,n)=0._r8
 	 end do
	end do
! if-testen avgjör bare om skytoppen befinner seg i dette modell-nivaaet eller ikke 
! (logisk variabel top = .true. eller .false.). Ghan-koden kalles i begge tilfeller.
! Forklaring: Hvis skyen dannes i det aktuelle tidsskrittet skal top uansett settes til 
! .false., altsaa tester jeg om det var en sky der i tidsskrittet för. Hvis det var en 
! sky der i forrige tidsskritt, tester jeg om det er sky i nivaaet over. Hvis ikke er 
! top =.true. Jeg sjekker ogsaa om det er sky i nivaaet under, for hvis ikke er det 
! aktuelle nivaaet baade skybunn & skytopp, og da skal top=.false. 
        if(k==1.or.k==pver) then
         top = .false.
        else
!TS      test om det er sky over, under eller i forrige tisskritt:
!AK      test om det er sky under men ikke over, samtidig som at det i dette nivået var sky i 
!AK      forrige tidsskritt. Hvorfor ikke i dette tidskrittet?  
!TS      at det faktisk er sky i dette tiddkrittet er testet allerede, lengre oppe....
!TS      NB: top=.false. alltid ved dannelse av ny sky
!+ 	 if(cldn(i,k-1) < 0.01_r8 .and. cldo(i,k) > 0.01_r8 .and. cldn(i,k+1) > 0.01_r8) then
 	 if(cldn(i,k-1) < 1.e-4_r8 .and. cldo(i,k) > 1.e-4_r8 .and. cldn(i,k+1) > 1.e-4_r8) then
	  top = .true.
	 else 
 	  top = .false.  
	 end if
        end if
        wminf = -1000._r8
        wmaxf = 1000._r8
        wdiab = 0._r8
	wbar = (omega(i,k)/(-gravit*rhoall(i)))*100._r8 ! Grid vertical velocity (small)
	!sigw= max(sigw1(i),30._r8)
	sigw= max(sigw1(i),10._r8) !corinna: Morrison and Gettelman (JClim 2008)
        wbar=sigw !corinna: Morrison and Gettelman (JClim 2008)
        sigw=0._r8 !corinna: uniform updraft Morrison and Gettelman (JClim 2008)
!	 Loop over 14 aerosol modes
       do kcomp= 1,14  ! -> nmodes
!       Note: In parmix 1.e-6 is used as lower limit (should be ok)
	if(Nnatk(i,kcomp) > 1.e-3_r8 .and. kcomp /= 12 .and. kcomp /= 3) then
	nmode=nmode+1
	if(kcomp < 5) ntype(nmode)= 2 ! Internally mixed modes
	if(kcomp > 4 .and. kcomp < 11) ntype(nmode)= 4
	if(kcomp > 10) ntype(nmode) = 1  ! Externally mixed sulfate
 	na(nmode)=Nnatk(i,kcomp)
	lnsigman(nmode) = lnsig(i,kcomp)

	if(kcomp==1) then ! Fine sulfate mode
	  ma(1,nmode)=Caitso4(i)*1.e-12_r8*Msv/Mso4
	  rhodry(1,nmode)=rhocomp(1)*1.e-3_r8
	  nu(1,nmode) = nucomp(1)
	  phi(1,nmode) = phicomp(1)
	  mwaer(1,nmode) = mwaercomp(1)
	  solubl(1,nmode) = solubcomp(1)
	elseif(kcomp==2) then ! Fine BC mode
	  ma(1,nmode)=Caitbc(i)*1.e-12_r8
	  rhodry(1,nmode)=rhocomp(2)*1.e-3_r8
	  nu(1,nmode) = nucomp(2)
	  phi(1,nmode) = phicomp(2)
	  mwaer(1,nmode) = mwaercomp(2)
	  solubl(1,nmode) = solubcomp(2)
	elseif(kcomp==4) then ! Fine OC/BC mode
	  ma(1,nmode)=Caitocbc(i)*1.e-12_r8
	  rhodry(1,nmode)=rhoocbc(i)*1.e-3_r8 !corinna: can become 0?
	  !rhodry(1,nmode)=MAX(rhocomp(3),rhoocbc(i))*1.e-3_r8
 	  nu(1,nmode) = faitbc(i)*nucomp(2)+(1.0_r8-faitbc(i))*nucomp(3)
	  phi(1,nmode) =  faitbc(i)*phicomp(2)+(1.0_r8-faitbc(i))*phicomp(3)
	  mwaer(1,nmode) =  faitbc(i)*mwaercomp(2)+(1.0_r8-faitbc(i))*mwaercomp(3)
	  solubl(1,nmode) =  faitbc(i)*solubcomp(2)+(1.0_r8-faitbc(i))*solubcomp(3)
	elseif(kcomp==5) then !SO4 aitken mode
	  ma(1,nmode)=Cas75(i)*1.e-12_r8*Ms/Mso4
	  rhodry(1,nmode)=rhocomp(1)*1.e-3_r8
	  nu(1,nmode) = nucomp(1)
	  phi(1,nmode) = phicomp(1)
	  mwaer(1,nmode) = mwaercomp(1)
	  solubl(1,nmode) = solubcomp(1)
	elseif(kcomp==6 .or. kcomp==7) then
	  rhodry(1,nmode)=rhocomp(5)*1.e-3_r8
	  nu(1,nmode) = nucomp(5)
	  phi(1,nmode) = phicomp(5)
	  mwaer(1,nmode) = mwaercomp(5)
	  solubl(1,nmode) = solubcomp(5)
	  if(kcomp==6) ma(1,nmode)=C_dst2(i)*1.e-12_r8
	  if(kcomp==7) ma(1,nmode)=C_dst3(i)*1.e-12_r8
	elseif(kcomp > 7 .and. kcomp < 11) then
	  rhodry(1,nmode)=rhocomp(4)*1.e-3_r8
	  nu(1,nmode) = nucomp(4)
	  phi(1,nmode) = phicomp(4)
	  mwaer(1,nmode) = mwaercomp(4)
	  solubl(1,nmode) = solubcomp(4)
	  if(kcomp==8) ma(1,nmode)=C_ss1(i)*1.e-12_r8
	  if(kcomp==9) ma(1,nmode)=C_ss2(i)*1.e-12_r8
	  if(kcomp==10) ma(1,nmode)=C_ss3(i)*1.e-12_r8
	elseif(kcomp==11) then
	  rhodry(1,nmode) = rhocomp(1)*1.e-3_r8
	  nu(1,nmode) = nucomp(1)
          phi(1,nmode) = phicomp(1)
          mwaer(1,nmode) = mwaercomp(1)
          solubl(1,nmode) = solubcomp(1)
	  ma(1,nmode)=Cnso4(i)*1.e-12_r8*Msv/Mso4
	elseif(kcomp==13) then
	  rhodry(1,nmode) = rhocomp(3)*1.e-3_r8
	  nu(1,nmode) = nucomp(3)
          phi(1,nmode) = phicomp(3)
          mwaer(1,nmode) = mwaercomp(3)
          solubl(1,nmode) = solubcomp(3)
	  ma(1,nmode)=Cnoc(i)*1.e-12_r8
	elseif(kcomp==14) then
	  rhodry(1,nmode) = rhobcocn(i)*1.e-3_r8 
 	  nu(1,nmode) = fnbc(i)*nucomp(2)+(1.0_r8-fnbc(i))*nucomp(3)
	  phi(1,nmode) =  fnbc(i)*phicomp(2)+(1.0_r8-fnbc(i))*phicomp(3)
	  mwaer(1,nmode) =  fnbc(i)*mwaercomp(2)+(1.0_r8-fnbc(i))*mwaercomp(3)
	  solubl(1,nmode) =  fnbc(i)*solubcomp(2)+(1.0_r8-fnbc(i))*solubcomp(3)
	  ma(1,nmode)=(Cnocibc(i)+Cnbcioc(i))*1.e-12_r8	
	end if

 	if(ntype(nmode) > 1) then           ! i.e., only kcomp = 1 to 10
!	 Sulfate hygroscopic properties
	  rhodry(2,nmode)=rhocomp(1)*1.e-3_r8
	  nu(2,nmode) = nucomp(1)
	  phi(2,nmode) = phicomp(1)
	  mwaer(2,nmode) = mwaercomp(1)
	  solubl(2,nmode) = solubcomp(1)
!	 Black carbon hygroscopic properties
	  rhodry(3,nmode)=rhocomp(2)*1.e-3_r8
	  nu(3,nmode) = nucomp(2)
	  phi(3,nmode) = phicomp(2)
	  mwaer(3,nmode) = mwaercomp(2)
	  solubl(3,nmode) = solubcomp(2)
!	 Organic carbon hygroscopic properties
	  rhodry(4,nmode)=rhocomp(3)*1.e-3_r8
	  nu(4,nmode) = nucomp(3)
	  phi(4,nmode) = phicomp(3)
	  mwaer(4,nmode) = mwaercomp(3)
	  solubl(4,nmode) = solubcomp(3)
          if(kcomp<=3) then ! all added mass is H2SO4
  	    ma(2,nmode)=cam(i,kcomp)*1.e-12_r8*(Msv/Mso4)
            ma(3,nmode)=0.0_r8
  	    ma(4,nmode)=0.0_r8
          elseif(kcomp==4) then ! all added mass is H2SO4 or (NH4)2SO4 (BC&OC is background)
 	    ma(2,nmode)=cam(i,kcomp)*1.e-12_r8*((1.0_r8-faqm(i,kcomp))*Msv+faqm(i,kcomp)*Ms)/Mso4
            ma(3,nmode)=0.0_r8
  	    ma(4,nmode)=0.0_r8
          elseif(kcomp>=5) then ! all added mass is H2SO4, (NH4)2SO4, BC or OC
 	    ma(2,nmode)=cam(i,kcomp)*1.e-12_r8*(1.0_r8-facm(i,kcomp))*((1.0_r8-faqm(i,kcomp))*Msv+faqm(i,kcomp)*Ms)/Mso4
            ma(3,nmode)=cam(i,kcomp)*1.e-12_r8*facm(i,kcomp)*fbcm(i,kcomp)
  	    ma(4,nmode)=cam(i,kcomp)*1.e-12_r8*facm(i,kcomp)*(1.-fbcm(i,kcomp))
          end if          
 !t+  more physically based coating assumptions (coating if layer thickness > 2 nm) 
          if(lthick(i,kcomp)>0.002_r8) then 
           if(kcomp==2) then                             ! for pure SO4 coating on BC
  	    nu(1,nmode) = nucomp(1)
	    phi(1,nmode) = phicomp(1)
	    mwaer(1,nmode) = mwaercomp(1)
	    solubl(1,nmode) = solubcomp(1)
  	    nu(2,nmode)     = nu(1,nmode)
  	    phi(2,nmode)    = phi(1,nmode)
  	    mwaer(2,nmode)  = mwaer(1,nmode)
  	    solubl(2,nmode) = solubl(1,nmode)
  	    nu(3,nmode)     = nu(1,nmode)
  	    phi(3,nmode)    = phi(1,nmode)
  	    mwaer(3,nmode)  = mwaer(1,nmode)
  	    solubl(3,nmode) = solubl(1,nmode)
  	    nu(4,nmode)     = nu(1,nmode)
  	    phi(4,nmode)    = phi(1,nmode)
  	    mwaer(4,nmode)  = mwaer(1,nmode)
  	    solubl(4,nmode) = solubl(1,nmode)
           elseif(kcomp>=4.and.kcomp<=7) then             ! SO4&OC (volume mixed) coating on BC or DU&BC
  	    nu(1,nmode)     = fso4coat(i,kcomp)*nucomp(1)   +(1.0_r8-fso4coat(i,kcomp))*nucomp(3)
  	    phi(1,nmode)    = fso4coat(i,kcomp)*phicomp(1)  +(1.0_r8-fso4coat(i,kcomp))*phicomp(3)
  	    mwaer(1,nmode)  = fso4coat(i,kcomp)*mwaercomp(1)+(1.0_r8-fso4coat(i,kcomp))*mwaercomp(3)
  	    solubl(1,nmode) = fso4coat(i,kcomp)*solubcomp(1)+(1.0_r8-fso4coat(i,kcomp))*solubcomp(3)
  	    nu(2,nmode)     = nu(1,nmode)
  	    phi(2,nmode)    = phi(1,nmode)
  	    mwaer(2,nmode)  = mwaer(1,nmode)
  	    solubl(2,nmode) = solubl(1,nmode)
  	    nu(3,nmode)     = nu(1,nmode)
  	    phi(3,nmode)    = phi(1,nmode)
  	    mwaer(3,nmode)  = mwaer(1,nmode)
  	    solubl(3,nmode) = solubl(1,nmode)
  	    nu(4,nmode)     = nu(1,nmode)
  	    phi(4,nmode)    = phi(1,nmode)
  	    mwaer(4,nmode)  = mwaer(1,nmode)
  	    solubl(4,nmode) = solubl(1,nmode)
           endif
          endif         
	end if  ! ntype(nmode) > 1

	if(version2 .and. kcomp <= 10) then
 	  do m=1,ntype(nmode)
	    ma(m,nmode) = Massratio(i,kcomp)*ma(m,nmode)
	    lnsigman(nmode) = log(10**logsig(i,kcomp))
	  end do
	end if

	end if ! kcomp-test
	end do ! kcomp

	tair = t(i,k)
        smax=0._r8    ! Corinna!

	if(activ .and. nmode > 0) then

       call activate_ny(top,wbar,sigw,wdiab,wminf,wmaxf,tair,rhocgs, &
                          na,ntype,nmode,ma,lnsigman,nu,phi,         &                                   
                          mwaer,rhodry,solubl,fn,fluxn,smax)
	end if

	do n = 1,nmode
 	 nuclx(i,k) = nuclx(i,k) + fn(n)*na(n)		
	end do
	supersat(i,k) = smax 
	end if
       end do ! ii->i
!TS-- End of activation code

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#endif  ! DIRIND

!
! calculate the conversion terms
!
   call get_rlat_all_p(lchnk, ncol, rlat)

!cdir nodep
!DIR$ CONCURRENT
   do ii = 1,ncols
      i = ind(ii)
      rhocgs = rho(i)*1.e-3_r8     ! density in cgs units
!
! exponential temperature factor
!
!        efact = exp(0.025_r8*(t(i,k)-t0))
!
! some temperature dependent constants
!
!!$      wt = min(1._r8,max(0._r8,(t0-t(i,k))*rdtice))
      wt = fice(i)
      icrit = icritc*wt + icritw*(1-wt)
!
! jrm Reworked droplet number concentration algorithm
      ! Start with pressure-dependent value appropriate for continental air
      ! Note: reltab has a temperature dependence here
      capn = capnw + (capnc-capnw) * min(1._r8,max(0._r8,1.0_r8-(p(i,k)-0.8_r8*p(i,pver))/(0.2_r8*p(i,pver))))
      ! Modify for snow depth over land
      capn = capn + (capnc-capn) * min(1.0_r8,max(0.0_r8,snowh(i)*10._r8))
      ! Ramp between polluted value over land to clean value over ocean.
      capn = capn + (capnc-capn) * min(1.0_r8,max(0.0_r8,1.0_r8-landm(i,lchnk)))
      ! Ramp between the resultant value and a sea ice value in the presence of ice.
      capn = capn + (capnsi-capn) * min(1.0_r8,max(0.0_r8,seaicef(i)))
! end jrm
!      
#ifdef DEBUG2
      if ( (lat(i) == latlook(1)) .or. (lat(i) == latlook(2)) ) then
         if (i == ilook(1)) then
            write(iulog,*) ' findmcnew: lat, k, seaicef, landm, wp, capnoice, capn ', &
                 lat(i), k, seaicef(i), landm(i,lat(i)), wp, capnoice, capn
         endif
      endif
#endif

!
! useful terms in following calculations
!
      rat1 = rhocgs/rhow
      rat2 = liqmr(i)/capn
      con2 = (rat1*rat2)**0.333_r8
#ifdef DIRIND
#ifdef AEROFFL
      rat2 = liqmr(i)/capn
#else
      rat2 = liqmr(i)/max(icnlx(i), 1._r8)  ! Note: CDNCmin = 1 cm-3
#endif 
      con2 = (rat1*rat2)**0.333_r8
      rat2x = liqmrx(i)/max(icnlx(i), 1._r8)  ! Note: CDNCmin = 1 cm-3  
      con2x = (rat1*rat2x)**0.333_r8
#endif
!
! volume radius
!
!        r3l = (rhocgs*liqmr(i)/(1.333_r8*pi*capn*rhow))**0.333_r8 * 0.01_r8 ! meters
      r3l = con1*con2
#ifdef DIRIND
       r3lx = con1*con2x
#ifdef AEROFFL
       volrad(i,k) = 1.e6_r8*r3lx
#else
       volrad(i,k) = 1.e6_r8*r3l
#endif
#endif
!
! critical threshold for autoconversion if modified for mixed phase
! clouds to mimic a bergeron findeisen process
! r3lc2 = r3lcrit*(1._r8-0.5_r8*fice(i)*(1-fice(i)))
!
! autoconversion of liquid
!
!        cwaut = 2.e-4_r8
!        cwaut = 1.e-3_r8
!        lcrit = 2.e-4_r8
!        lcrit = 5.e-4_r8
!        pwaut = max(0._r8,liqmr(i)-lcrit)*cwaut
!
! pwaut is following tripoli and cotton (and many others)
! we reduce the autoconversion below critpr, because these are regions where
! the drop size distribution is likely to imply much smaller collector drops than
! those relevant for a cloud distribution corresponding to the value of effc = 0.55
! suggested by cotton (see austin 1995 JAS, baker 1993)

! easy to follow form
!        pwaut = capc*liqmr(i)**2*rhocgs/rhow
!    $           *(liqmr(i)*rhocgs/(rhow*capn))**(.333)
!    $           *heavy(r3l,r3lcrit)
!    $           *max(0.10_r8,min(1._r8,prlloc(i)/critpr))
! somewhat faster form
#define HEAVYNEW
#ifdef HEAVYNEW
!#ifdef PERGRO
      pwaut = capc*liqmr(i)**2*rat1*con2*heavymp(r3l,r3lcrit) * &
              max(0.10_r8,min(1._r8,prlloc(i)/critpr))
#ifdef DIRIND
      pwautx= capc*liqmrx(i)**2*rat1*con2x*heavymp(r3lx,r3lcrit) * &
              max(0.10_r8,min(1._r8,prllocx(i)/critpr))
#endif
#else
      pwaut = capc*liqmr(i)**2*rat1*con2*heavym(r3l,r3lcrit)* &
              max(0.10_r8,min(1._r8,prlloc(i)/critpr))
#ifdef DIRIND
      pwautx = capc*liqmrx(i)**2*rat1*con2x*heavym(r3lx,r3lcrit)* &
              max(0.10_r8,min(1._r8,prllocx(i)/critpr))
#endif
#endif
!
! autoconversion of ice
!
!        ciaut = ciautb*efact
      ciaut = ciautb
!        psaut = capc*totmr(i)**2*rhocgs/rhoi
!     $           *(totmr(i)*rhocgs/(rhoi*capn))**(.333)
!
! autoconversion of ice condensate
!
#ifdef PERGRO
      psaut = heavyp(icemr(i),icrit)*icemr(i)*ciaut
#else
      psaut = max(0._r8,icemr(i)-icrit)*ciaut
#endif
#ifdef DIRIND  ! for off-line use
      psautx = max(0._r8,icemrx(i)-icrit)*ciaut ! autoconversion of ice condensate
#endif
!
! collection of liquid by rain
!
!        pracw = cracw*rho(i)*liqmr(i)*rainmr(i) !(beheng 1994)
      pracw = cracw*rho(i)*sqrt(rho(i))*liqmr(i)*rainmr(i) !(tripoli and cotton)
#ifdef DIRIND  ! for off-line use
      pracwx= cracw*rho(i)*sqrt(rho(i))*liqmrx(i)*rainmrx(i) !(tripoli and cotton)
#endif

      pracwo(i)=pracw

!!      pracw = 0._r8
!
! the following lines calculate the slope parameter and snow mixing ratio
! from the precip rate using the equations found in lin et al 83
! in the most natural form, but it is expensive, so after some tedious
! algebraic manipulation you can use the cheaper form found below
!            vfalls = c*gam4pd/(6._r8*lamdas**d)*sqrt(rhonot/rhocgs)
!     $               *0.01_r8   ! convert from cm/s to m/s
!            snowmr(i) = snowfr*precab(i)/(rho(i)*vfalls*cldpr(i))
!            snowmr(i) = ( prscgs(i)*mcon02 * (rhocgs**mcon03) )**mcon04
!            lamdas = (prhonos/max(rhocgs*snowmr(i),small))**0.25_r8
!            csacw = mcon01*sqrt(rhonot/rhocgs)/(lamdas**thrpd)
!
! coefficient for collection by snow independent of phase
!
      csacx = mcon07*rhocgs**mcon08*prscgs(i)**mcon05

!
! collection of liquid by snow (lin et al 1983)
!
      psacw = csacx*liqmr(i)*esw
#ifdef DIRIND
      psacwx= csacx*liqmrx(i)*esw
#endif
#ifdef PERGRO
! this is necessary for pergro
      psacw = 0._r8
#endif

      psacwo(i)=psacw

!
! collection of ice by snow (lin et al 1983)
!
      psaci = csacx*icemr(i)*esi
#ifdef DIRIND
      psacix= csacx*icemrx(i)*esi
#endif
!
      psacio(i)=psaci

! total conversion of condensate to precipitate
!
      ptot = pwaut + psaut + pracw + psacw + psaci
#ifdef DIRIND
      ptotx= pwautx + psautx + pracwx + psacwx + psacix
#ifdef AEROFFL
      if(liqmrx(i) > 1.e-20_r8) then
        rtotx = (icnlx(i)/liqmrx(i))*(pwautx + pracwx + psacwx)
       else
	rtotx = 0._r8 
       end if
#else
       if(liqmr(i) > 1.e-20_r8) then
        rtotx = (icnlx(i)/liqmr(i))*(pwaut + pracw + psacw)  ! note: pwautx etc originally
       else
	rtotx = 0._r8 
       end if
#endif
#endif
!
! the recipricol of cloud water amnt (or zero if no cloud water)
!
!         rcwm =  totmr(i)/(max(totmr(i),small)**2)
!
! turn the tendency back into a loss rate (1/seconds)
!
      if (totmr(i) > 0._r8) then
         coef(i) = ptot/totmr(i)
      else
         coef(i) = 0._r8
      endif
#ifdef DIRIND
      if (totmrx(i) > 1.e-20_r8) then 
         coefx(i) = ptotx/totmrx(i)
      else
         coefx(i) = 0._r8
      endif
      if (icnlx(i) > 1.e-10_r8) then 
         ncoefx(i) = rtotx/icnlx(i)
      else
         ncoefx(i) = 0._r8
      endif
#endif
      if (ptot.gt.0._r8) then
         fwaut(i) = pwaut/ptot
         fsaut(i) = psaut/ptot
         fracw(i) = pracw/ptot
         fsacw(i) = psacw/ptot
         fsaci(i) = psaci/ptot
      else
         fwaut(i) = 0._r8
         fsaut(i) = 0._r8
         fracw(i) = 0._r8
         fsacw(i) = 0._r8
         fsaci(i) = 0._r8
      endif
#ifdef DIRIND
      if (ptotx.gt.0._r8) then
         fwautx(i) = pwautx/ptotx
         fsautx(i) = psautx/ptotx
         fracwx(i) = pracwx/ptotx
         fsacwx(i) = psacwx/ptotx
         fsacix(i) = psacix/ptotx
      else
         fwautx(i) = 0._r8
         fsautx(i) = 0._r8
         fracwx(i) = 0._r8
         fsacwx(i) = 0._r8
         fsacix(i) = 0._r8
      endif
#endif

      ftot = fwaut(i)+fsaut(i)+fracw(i)+fsacw(i)+fsaci(i)
!      if (abs(ftot-1._r8).gt.1.e-14_r8.and.ftot.ne.0._r8) then
!         write(iulog,*) ' something is wrong in findmcnew ', ftot, &
!              fwaut(i),fsaut(i),fracw(i),fsacw(i),fsaci(i)
!         write(iulog,*) ' unscaled ', ptot, &
!              pwaut,psaut,pracw,psacw,psaci
!         write(iulog,*) ' totmr, liqmr, icemr ', totmr(i), liqmr(i), icemr(i)
!         call endrun()
!      endif
   end do
#ifdef DEBUG
   i = icollook(nlook)
   if (lchnk == lchnklook(nlook) ) then
      write(iulog,*)
      write(iulog,*) '------', k, i, lchnk
      write(iulog,*) ' liqmr, rainmr,precab ', liqmr(i), rainmr(i), precab(i)*8.64e4_r8
      write(iulog,*) ' frac: waut,saut,racw,sacw,saci ', &
           fwaut(i), fsaut(i), fracw(i), fsacw(i), fsaci(i)
   endif
#endif

   return
end subroutine findmcnew

!##############################################################################

subroutine findsp (lchnk, ncol, q, t, p, tsp, qsp)
!----------------------------------------------------------------------- 
! 
! Purpose: 
!     find the wet bulb temperature for a given t and q
!     in a longitude height section
!     wet bulb temp is the temperature and spec humidity that is 
!     just saturated and has the same enthalpy
!     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method: 
! a Newton method is used
! first guess uses an algorithm provided by John Petch from the UKMO
! we exclude points where the physical situation is unrealistic
! e.g. where the temperature is outside the range of validity for the
!      saturation vapor pressure, or where the water vapor pressure
!      exceeds the ambient pressure, or the saturation specific humidity is 
!      unrealistic
! 
! Author: P. Rasch
! 
!-----------------------------------------------------------------------
!
!     input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(pcols,pver)        ! water vapor (kg/kg)
   real(r8), intent(in) :: t(pcols,pver)        ! temperature (K)
   real(r8), intent(in) :: p(pcols,pver)        ! pressure    (Pa)
!
! output arguments
!
   real(r8), intent(out) :: tsp(pcols,pver)      ! saturation temp (K)
   real(r8), intent(out) :: qsp(pcols,pver)      ! saturation mixing ratio (kg/kg)
!
! local variables
!
   integer i                 ! work variable
   integer k                 ! work variable
   logical lflg              ! work variable
   integer iter              ! work variable
   integer l                 ! work variable
   logical :: error_found

   real(r8) omeps                ! 1 minus epsilon
   real(r8) trinv                ! work variable
   real(r8) es                   ! sat. vapor pressure
   real(r8) desdt                ! change in sat vap pressure wrt temperature
!     real(r8) desdp                ! change in sat vap pressure wrt pressure
   real(r8) dqsdt                ! change in sat spec. hum. wrt temperature
   real(r8) dgdt                 ! work variable
   real(r8) g                    ! work variable
   real(r8) weight(pcols)        ! work variable
   real(r8) hlatsb               ! (sublimation)
   real(r8) hlatvp               ! (vaporization)
   real(r8) hltalt(pcols,pver)   ! lat. heat. of vap.
   real(r8) tterm                ! work var.
   real(r8) qs                   ! spec. hum. of water vapor
   real(r8) tc                   ! crit temp of transition to ice

! work variables
   real(r8) t1, q1, dt, dq
   real(r8) dtm, dqm
   real(r8) qvd, a1, tmp
   real(r8) rair
   real(r8) r1b, c1, c2, c3
   real(r8) denom
   real(r8) dttol
   real(r8) dqtol
   integer doit(pcols) 
   real(r8) enin(pcols), enout(pcols)
   real(r8) tlim(pcols)

   omeps = 1.0_r8 - epsqs
   trinv = 1.0_r8/ttrice
   a1 = 7.5_r8*log(10._r8)
   rair =  287.04_r8
   c3 = rair*a1/cp
   dtm = 0._r8    ! needed for iter=0 blowup with f90 -ei
   dqm = 0._r8    ! needed for iter=0 blowup with f90 -ei
   dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
   dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
!  tmin = 173.16_r8 ! the coldest temperature we can deal with
!
! max number of times to iterate the calculation
   iter = 8
!
   do k = k1mb,pver

!
! first guess on the wet bulb temperature
!
      do i = 1,ncol

#ifdef DEBUG
         if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
            write(iulog,*) ' '
            write(iulog,*) ' level, t, q, p', k, t(i,k), q(i,k), p(i,k)
         endif
#endif
! limit the temperature range to that relevant to the sat vap pres tables
#if ( ! defined WACCM_PHYS )
         tlim(i) = min(max(t(i,k),173._r8),373._r8)
#else
         tlim(i) = min(max(t(i,k),128._r8),373._r8)
#endif
         es = estblf(tlim(i))
         denom = p(i,k) - omeps*es
         qs = epsqs*es/denom
         doit(i) = 0
         enout(i) = 1._r8
! make sure a meaningful calculation is possible
         if (p(i,k) > 5._r8*es .and. qs > 0._r8 .and. qs < 0.5_r8) then
!
! Saturation specific humidity
!
             qs = min(epsqs*es/denom,1._r8)
!
! "generalized" analytic expression for t derivative of es
!  accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
             tc     = tlim(i) - t0
             lflg   = (tc >= -ttrice .and. tc < 0.0_r8)
             weight(i) = min(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0_r8*tc
             if (tlim(i) < t0) then
                hltalt(i,k) = hlatsb
             else
                hltalt(i,k) = hlatvp
             end if
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36._r8)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
#ifdef DEBUG
             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                write(iulog,*) ' relative humidity ', q(i,k)/qs
                write(iulog,*) ' first guess ', tsp(i,k)
             endif
#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
          else
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1._r8
          endif
       end do   ! end do i
!
! now iterate on first guess
!
      do l = 1, iter
         dtm = 0
         dqm = 0
         do i = 1,ncol
            if (doit(i) == 0) then
               es = estblf(tsp(i,k))
!
! Saturation specific humidity
!
               qs = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
               tc     = tsp(i,k) - t0
               lflg   = (tc >= -ttrice .and. tc < 0.0_r8)
               weight(i) = min(-tc*trinv,1.0_r8)
               hlatsb = hlatv + weight(i)*hlatf
               hlatvp = hlatv - 2369.0_r8*tc
               if (tsp(i,k) < t0) then
                  hltalt(i,k) = hlatsb
               else
                  hltalt(i,k) = hlatvp
               end if
               if (lflg) then
                  tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
               else
                  tterm = 0.0_r8
               end if
               desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
               dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
!              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
               g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
               dgdt = -(cp + hltalt(i,k)*dqsdt)
               t1 = tsp(i,k) - g/dgdt
               dt = abs(t1 - tsp(i,k))/t1
               tsp(i,k) = max(t1,tmin)
               es = estblf(tsp(i,k))
               q1 = min(epsqs*es/(p(i,k) - omeps*es),1._r8)
               dq = abs(q1 - qsp(i,k))/max(q1,1.e-12_r8)
               qsp(i,k) = q1
#ifdef DEBUG
               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                  write(iulog,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
               endif
#endif
               dtm = max(dtm,dt)
               dqm = max(dqm,dq)
! if converged at this point, exclude it from more iterations
               if (dt < dttol .and. dq < dqtol) then
                  doit(i) = 2
               endif
               enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
! bail out if we are too near the end of temp range
#if ( ! defined WACCM_PHYS )
               if (tsp(i,k) < 174.16_r8) then
#else
               if (tsp(i,k) < 130.16_r8) then
#endif
                  doit(i) = 4
               endif
            else
            endif
         end do              ! do i = 1,ncol

         if (dtm < dttol .and. dqm < dqtol) then
            go to 10
         endif

      end do                 ! do l = 1,iter
10    continue

      error_found = .false.
      if (dtm > dttol .or. dqm > dqtol) then
         do i = 1,ncol
            if (doit(i) == 0) error_found = .true.
         end do
         if (error_found) then
            do i = 1,ncol
               if (doit(i) == 0) then
                  write(iulog,*) ' findsp not converging at point i, k ', i, k
                  write(iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                  write(iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                  call endrun ('FINDSP')
               endif
            end do
         endif
      endif
      do i = 1,ncol
         if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) then
            error_found = .true.
         endif
      end do
      if (error_found) then
         do i = 1,ncol
            if (doit(i) == 2 .and. abs((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) then
               write(iulog,*) ' the enthalpy is not conserved for point ', &
                  i, k, enin(i), enout(i)
               write(iulog,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
               write(iulog,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
               call endrun ('FINDSP')
            endif
         end do
      endif
      
   end do                    ! level loop (k=1,pver)

   return
end subroutine findsp

end module cldwat
