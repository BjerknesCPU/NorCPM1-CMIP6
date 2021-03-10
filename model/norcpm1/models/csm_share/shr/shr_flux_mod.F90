!===============================================================================
! SVN $Id: shr_flux_mod.F90 28565 2011-05-20 19:24:42Z fischer $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/branch_tags/cesm1_0_4_rel_tags/cesm1_0_4_n02_share3_110527/shr/shr_flux_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: flux_mod -- CCSM shared flux calculations.
!
! !DESCRIPTION:
!
!     CCSM shared flux calculations.
!     
! !REVISION HISTORY:
!     2006-Nov-07 - B. Kauffman - first version, code taken/migrated from cpl6
!
! !INTERFACE: ------------------------------------------------------------------

module shr_flux_mod

! !USES:

   use shr_kind_mod    ! shared kinds
   use shr_const_mod   ! shared constants
   use shr_sys_mod     ! shared system routines
   use shr_log_mod, only: s_loglev  => shr_log_Level
   use shr_log_mod, only: s_logunit => shr_log_Unit

   implicit none

   private ! default private

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_flux_atmOcn      ! computes atm/ocn fluxes
   public :: shr_flux_atmIce      ! computes atm/ice fluxes
   public :: shr_flux_MOstability ! boundary layer stability scales/functions
   public :: psit_30
   public :: psiuo

! !PUBLIC DATA MEMBERS:

  integer(SHR_KIND_IN),parameter,public :: shr_flux_MOwScales   = 1 ! w scales  option
  integer(SHR_KIND_IN),parameter,public :: shr_flux_MOfunctions = 2 ! functions option
  real   (SHR_KIND_R8),parameter,public :: shr_flux_MOgammaM = 3.59_SHR_KIND_R8
  real   (SHR_KIND_R8),parameter,public :: shr_flux_MOgammaS = 7.86_SHR_KIND_R8

 ! added by tht 29/01/2014 
!  integer(SHR_KIND_IN),parameter,public :: shr_flux_DOCOARE = 1 ! COARE flux computation option

!EOP

   !--- rename kinds for local readability only ---
   integer,parameter :: R8 = SHR_KIND_R8  ! 8 byte real
   integer,parameter :: IN = SHR_KIND_IN  ! native/default integer

   integer,parameter :: debug = 0 ! internal debug level
   integer, save :: flux_scheme 
   public :: shr_map_docoare 
!===============================================================================
contains
!===============================================================================
!===============================================================================
! !BOP =========================================================================
!
! IROUTINE: shr_flux_atmOcn -- internal atm/ocn flux calculation
!
! !DESCRIPTION:
!
!     Internal atm/ocn flux calculation
!     
! !REVISION HISTORY:
!     2002-Jun-10 - B. Kauffman - code migrated from cpl5 to cpl6
!     2003-Apr-02 - B. Kauffman - taux & tauy now utilize ocn velocity
!     2003-Apr-02 - B. Kauffman - tref,qref,duu10n mods as per Bill Large
!     2006-Nov-07 - B. Kauffman - code migrated from cpl6 to share
!
!     2014-Jan-01 - T. Toniazzo - introduced convergence loop & cleaned up the
!                                 general mess a little bit & added comments;
!                                 altern. option for COAREv3.0 flux computation
!                                 controlled by integer variable shr_flux_DOCOARE
!                                 (=0 -> use old algorithm; any other value -> use 
!                                 COARE algorithm)
!
! !INTERFACE: ------------------------------------------------------------------

SUBROUTINE shr_flux_atmOcn(nMax  ,zbot  ,ubot  ,vbot  ,thbot ,   & 
           &               qbot  ,rbot  ,tbot  ,us    ,vs    ,   &
           &               ts    ,mask  ,sen   ,lat   ,lwup  ,   &
           &               evap  ,taux  ,tauy  ,tref  ,qref  ,   &
           &               duu10n,  ustar_sv   ,re_sv ,ssq_sv,   &
           &               missval    )

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in) ::       nMax  ! data vector length
   integer(IN),intent(in) :: mask (nMax) ! ocn domain mask       0 <=> out of domain
   real(R8)   ,intent(in) :: zbot (nMax) ! atm level height      (m)
   real(R8)   ,intent(in) :: ubot (nMax) ! atm u wind            (m/s)
   real(R8)   ,intent(in) :: vbot (nMax) ! atm v wind            (m/s)
   real(R8)   ,intent(in) :: thbot(nMax) ! atm potential T       (K)
   real(R8)   ,intent(in) :: qbot (nMax) ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in) :: rbot (nMax) ! atm air density       (kg/m^3)
   real(R8)   ,intent(in) :: tbot (nMax) ! atm T                 (K) 
   real(R8)   ,intent(in) :: us   (nMax) ! ocn u-velocity        (m/s)
   real(R8)   ,intent(in) :: vs   (nMax) ! ocn v-velocity        (m/s)
   real(R8)   ,intent(in) :: ts   (nMax) ! ocn temperature       (K)

   !--- output arguments -------------------------------
   real(R8),intent(out)  ::  sen  (nMax) ! heat flux: sensible    (W/m^2)
   real(R8),intent(out)  ::  lat  (nMax) ! heat flux: latent      (W/m^2)
   real(R8),intent(out)  ::  lwup (nMax) ! heat flux: lw upward   (W/m^2)
   real(R8),intent(out)  ::  evap (nMax) ! water flux: evap  ((kg/s)/m^2)
   real(R8),intent(out)  ::  taux (nMax) ! surface stress, zonal      (N)
   real(R8),intent(out)  ::  tauy (nMax) ! surface stress, maridional (N)
   real(R8),intent(out)  ::  tref (nMax) ! diag:  2m ref height T     (K)
   real(R8),intent(out)  ::  qref (nMax) ! diag:  2m ref humidity (kg/kg)
   real(R8),intent(out)  :: duu10n(nMax) ! diag: 10m wind speed squared (m/s)^2

   real(R8),intent(out),optional :: ustar_sv(nMax) ! diag: ustar
   real(R8),intent(out),optional :: re_sv   (nMax) ! diag: sqrt of exchange coefficient (water)
   real(R8),intent(out),optional :: ssq_sv  (nMax) ! diag: sea surface humidity  (kg/kg)
 
   real(R8),intent(in) ,optional :: missval        ! masked value

! !EOP

   !--- local variables --------------------------------
   integer(IN) :: n      ! vector loop index
   integer(IN) :: k      ! interation loop index
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
   real(R8)    :: thvbot ! virtual temperature      (K)
   real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
   real(R8)    :: delt   ! potential T difference   (K)
   real(R8)    :: delq   ! humidity difference      (kg/kg)
   real(R8)    :: stable ! stability factor
   real(R8)    :: rdn    ! sqrt of neutral exchange coeff (momentum) 
   real(R8)    :: rhn    ! sqrt of neutral exchange coeff (heat)     
   real(R8)    :: ren    ! sqrt of neutral exchange coeff (water)    
   real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)         
   real(R8)    :: rh     ! sqrt of exchange coefficient (heat)             
   real(R8)    :: re     ! sqrt of exchange coefficient (water)            
   real(R8)    :: ustar  ! ustar             
   real(R8)    :: qstar  ! qstar             
   real(R8)    :: tstar  ! tstar             
   real(R8)    :: hol    ! H (at zbot) over L
   real(R8)    :: xsq    ! ?
   real(R8)    :: xqq    ! ?
   real(R8)    :: u10n   ! 10m neutral wind 
   real(R8)    :: tau    ! stress at zbot
   real(R8)    :: cp     ! specific heat of moist air
   real(R8)    :: bn     ! exchange coef funct for interpolation
   real(R8)    :: bh     ! exchange coef funct for interpolation
   real(R8)    :: fac    ! vertical interpolation factor
   real(R8)    :: spval  ! local missing value

! non-COARE
   real(R8)    :: psimh  ! stability function at zbot (momentum)
   real(R8)    :: psixh  ! stability function at zbot (heat and water)
   real(R8)    :: psix2  ! stability function at ztref reference height
   real(R8)    :: alz    ! ln(zbot/zref)
   real(R8)    :: al2    ! ln(zref/ztref)

! COARE
   real(R8)    :: zo,zot,zoq      ! roughness lengths
   real(R8)    :: hsb,hlb         ! sens & lat heat flxs at zbot
   real(R8)    :: trf,qrf,urf,vrf ! reference-height quantities


   !--- local constants --------------------------------
   real(R8),parameter :: zref  = 10.0_R8 ! reference height           (m)
   real(R8),parameter :: ztref =  2.0_R8 ! reference height for air T (m)

! non-COARE
   real(R8),parameter :: umin  =  0.5_R8 ! minimum wind speed       (m/s)
   real(R8),parameter :: cexcd = 0.0346_R8  ! ratio Ch(water)/CD 
   real(R8),parameter :: chxcds = 0.018_R8  ! ratio Ch(heat)/CD for stable case
   real(R8),parameter :: chxcdu = 0.0327_R8 ! ratio Ch(heat)/CD for unstable case

! tht 12/11/2013
   integer(IN),parameter :: kmax = 1  ! No. of interations

! COARE
   real(R8),parameter :: zpbl  =700.0_R8 ! PBL depth [m] for gustiness parametriz.
 

   !--- local functions --------------------------------
   real(R8)    :: qsat   ! function: the saturation humididty of air (kg/m^3)
   real(R8)    :: Tk     ! dummy arg ~ temperature (K)
    qsat(Tk)   = 640380.0_R8 / exp(5107.4_R8/Tk)

! non-COARE: Large et al. (1994) formula v*=[c4/U10+c5+c6*U10]*U10
   real(R8)    :: cdn    ! function: neutral drag coeff at 10m
   real(R8)    :: Umps   ! dummy arg ~ wind velocity (m/s)
    cdn(Umps)  =   0.0027_R8 / Umps + 0.000142_R8 + 0.0000764_R8 * Umps

! non-COARE: stability functions 
   real(R8)    :: psimhu ! function: unstable part of psimh
   real(R8)    :: psixhu ! function: unstable part of psimx
   real(R8)    :: xd     ! dummy arg ~ Chi := (1 - 16*zeta)^{1/4} with zeta=(z-d)/L
    psimhu(xd) = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) - 2.0_R8*atan(xd) + 1.571_R8
    psixhu(xd) = 2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)


   !--- formats ----------------------------------------
   character(*),parameter :: subName = '(shr_flux_atmOcn) '
   character(*),parameter ::   F00 = "('(shr_flux_atmOcn) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   computes atm/ocn surface fluxes
!
! NOTES: 
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + sen + lat
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
! 
! ASSUMPTIONS:
!  non-COARE:
!   o Neutral 10m drag coeff: cdn = .0027/U10 + .000142 + .0000764 U10
!   o Neutral 10m stanton number: ctn = .0327 sqrt(cdn), unstable
!                                 ctn = .0180 sqrt(cdn), stable
!   o Neutral 10m dalton number:  cen = .0346 sqrt(cdn)
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!  COARE:
!   o use COAREv3.0 function (tht 22/11/2013)
!-------------------------------------------------------------------------------

   if (debug > 0 .and. s_loglev > 0) write(s_logunit,F00) "enter"

   if (present(missval)) then
      spval = missval
   else
      spval = shr_const_spval
   endif
 
!=============================================================
   IF (flux_scheme .eq.0) THEN 
!----------------------------
!    non-COARE computation
!----------------------------
    al2 = log(zref/ztref)

    DO n=1,nMax
     if (mask(n) /= 0) then
    
        !--- compute some needed quantities ---
        vmag   = max(umin, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) )
        ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)

        thvbot = thbot(n) * (1.0_R8 + shr_const_zvir * qbot(n)) ! virtual temp (K)
        delt   = thbot(n) - ts(n)                  ! pot temp diff (K)
        delq   = qbot(n) - ssq                     ! spec hum dif (kg/kg)
        alz    = log(zbot(n)/zref) 
        cp     = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*ssq) 
   
        !------------------------------------------------------------
        ! first estimate of Z/L and ustar, tstar and qstar
        !------------------------------------------------------------
   
        !--- neutral coefficients, z/L = 0.0 ---
        rdn    = sqrt(cdn(vmag))
!tht: this uses fixed St/Da numbers i.e. ratios of Ch/CD
        ren    = cexcd
!tht: iterate once towards value of Rhn
        stable = 0.5_R8 + sign(0.5_R8 , delt)
        rhn    = (1.0_R8-stable) * chxcdu + stable * chxcds 
        stable = 0.5_R8 + sign(0.5_R8 ,(rhn*delt/thbot(n)+ren*delq/(1.0_R8/shr_const_zvir+qbot(n))))
        rhn    = (1.0_R8-stable) * chxcdu + stable * chxcds 
   
        !--- ustar, tstar, qstar ---
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  

!--- iteration loop tht 12/11/2013
        DO k=1,kmax
!--- iteration loop tht 12/11/2013
   
        !--- compute stability & evaluate all stability functions ---
        hol  = shr_const_karman*shr_const_g*zbot(n)*  &
               (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)

        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?
!       psimh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
!       psixh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?
   
        !--- shift wind speed using old coefficient ---
        rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
        u10n = vmag * rd / rdn 
   
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
!tht: fixed St/Da numbers i.e. ratios of Ch/CD
        rhn    = (1.0_R8-stable) * chxcdu + stable * chxcds 
!       ren    = cexcd ! tht: unchanged
    
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh)) 
   
        !--- update ustar, tstar, qstar using updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 
    
        !------------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !------------------------------------------------------------
    
        !--- compute stability & evaluate all stability functions ---
        hol  = shr_const_karman*shr_const_g*zbot(n)* &
               (tstar/thbot(n)+qstar/(1.0_R8/shr_const_zvir+qbot(n)))/ustar**2
        hol  = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)

        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?
!       psimh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
!       psixh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?
    
        !--- shift wind speed using old coeffs ---
        rd   = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh))
        u10n = vmag * rd/rdn 
    
        !--- update transfer coeffs at 10m and neutral stability ---
        rdn = sqrt(cdn(u10n))
!tht: fixed St/Da numbers i.e. ratios of Ch/CD
        rhn    = (1.0_R8-stable) * chxcdu + stable * chxcds 
!       ren    = cexcd ! tht: unchanged
   
        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8 + rdn/shr_const_karman*(alz-psimh)) 
        rh = rhn / (1.0_R8 + rhn/shr_const_karman*(alz-psixh)) 
        re = ren / (1.0_R8 + ren/shr_const_karman*(alz-psixh)) 
    
        !--- update ustar, tstar, qstar using updated, shifted coeffs ---
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

!--- iteration loop tht12/11/2013
        ENDDO
!--- iteration loop tht12/11/2013

    
        !------------------------------------------------------------
        ! compute the fluxes
        !------------------------------------------------------------
    
        tau = rbot(n) * ustar * ustar 
       
        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag 
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag 
        
        !--- heat flux ---
        sen (n) =                cp * tau * tstar / ustar 
        lat (n) =  shr_const_latvap * tau * qstar / ustar
        lwup(n) = -shr_const_stebol * ts(n)**4 
      
        !--- water flux ---
        evap(n) = lat(n)/shr_const_latvap 
    
        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        hol = hol*ztref/zbot(n)
        xsq = max( 1.0_R8, sqrt(abs(1.0_R8-16.0_R8*hol)) )
        xqq = sqrt(xsq)
!tht 15/11/2013 5->4 as in UKMO UM?
        psix2   = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!       psix2   = -4.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?
        fac     = (rh/shr_const_karman) * (alz + al2 - psixh + psix2 )
        tref(n) = thbot(n) - delt*fac 
        tref(n) = tref(n) - 0.01_R8*ztref   ! pot temp to temp correction
        fac     = (re/shr_const_karman) * (alz + al2 - psixh + psix2 )
        qref(n) =  qbot(n) - delq*fac
    
        duu10n(n) = u10n*u10n ! 10m wind speed squared

        !------------------------------------------------------------
        ! optional diagnostics, needed for water tracer fluxes (dcn)
        !------------------------------------------------------------
        if (present(ustar_sv)) ustar_sv(n) = ustar
        if (present(re_sv   )) re_sv(n)    = re
        if (present(ssq_sv  )) ssq_sv(n)   = ssq

     else
        !------------------------------------------------------------
        ! no valid data here -- out of domain
        !------------------------------------------------------------
        sen   (n) = spval  ! sensible         heat flux  (W/m^2)
        lat   (n) = spval  ! latent           heat flux  (W/m^2)
        lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
        evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
        taux  (n) = spval  ! x surface stress (N)
        tauy  (n) = spval  ! y surface stress (N)
        tref  (n) = spval  !  2m reference height temperature (K)
        qref  (n) = spval  !  2m reference height humidity (kg/kg)
        duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

        if (present(ustar_sv)) ustar_sv(n) = spval
        if (present(re_sv   )) re_sv   (n) = spval
        if (present(ssq_sv  )) ssq_sv  (n) = spval
     endif
    ENDDO 

!=============================================================
   ELSE ! use COARE algorithm
!=============================================================

    DO n=1,nMax
     if (mask(n) /= 0) then
    
        !--- compute some needed quantities ---
        vmag   = max(1e-6, sqrt( (ubot(n)-us(n))**2 + (vbot(n)-vs(n))**2) ) 
        ssq    = 0.98_R8 * qsat(ts(n)) / rbot(n)   ! sea surf hum (kg/kg)

        call cor30a(ubot(n),vbot(n),tbot(n),qbot(n),rbot(n) &  ! in atm params		
                 & ,us(n),vs(n),ts(n),ssq                   &  ! in surf params		
                 & ,zpbl,zbot(n),zbot(n),zref,ztref,ztref   &  ! in heights			
                 & ,tau,hsb,hlb                             &  ! out: fluxes			
                 & ,zo,zot,zoq,hol,ustar,tstar,qstar        &  ! out: ss scales		
                 & ,rd,rh,re                                &  ! out: exch. coeffs		
                 & ,trf,qrf,urf,vrf)			       ! out: reference-height params

! for the sake of maintaining same defs
        hol=zbot(n)/hol
        rd=sqrt(rd)
        rh=sqrt(rh)
        re=sqrt(re)
       
        !--- momentum flux ---
        taux(n) = tau * (ubot(n)-us(n)) / vmag 
        tauy(n) = tau * (vbot(n)-vs(n)) / vmag 
        
        !--- heat flux ---
        sen (n) =  hsb
        lat (n) =  hlb
        lwup(n) = -shr_const_stebol * ts(n)**4 
      
        !--- water flux ---
        evap(n) = lat(n)/shr_const_latvap 
    
        !------------------------------------------------------------
        ! compute diagnositcs: 2m ref T & Q, 10m wind speed squared
        !------------------------------------------------------------
        tref(n) = trf
        qref(n) = qrf
        duu10n(n) = urf**2+vrf**2 

        !------------------------------------------------------------
        ! optional diagnostics, needed for water tracer fluxes (dcn)
        !------------------------------------------------------------
        if (present(ustar_sv)) ustar_sv(n) = ustar
        if (present(re_sv   )) re_sv(n)    = re
        if (present(ssq_sv  )) ssq_sv(n)   = ssq

     else
        !------------------------------------------------------------
        ! no valid data here -- out of domain
        !------------------------------------------------------------
        sen   (n) = spval  ! sensible         heat flux  (W/m^2)
        lat   (n) = spval  ! latent           heat flux  (W/m^2)
        lwup  (n) = spval  ! long-wave upward heat flux  (W/m^2)
        evap  (n) = spval  ! evaporative water flux ((kg/s)/m^2)
        taux  (n) = spval  ! x surface stress (N)
        tauy  (n) = spval  ! y surface stress (N)
        tref  (n) = spval  !  2m reference height temperature (K)
        qref  (n) = spval  !  2m reference height humidity (kg/kg)
        duu10n(n) = spval  ! 10m wind speed squared (m/s)^2

        if (present(ustar_sv)) ustar_sv(n) = spval
        if (present(re_sv   )) re_sv   (n) = spval
        if (present(ssq_sv  )) ssq_sv  (n) = spval
     endif
    ENDDO 

   ENDIF
!=============================================================

END subroutine shr_flux_atmOcn

!BOP ===========================================================================
!
! !IROUTINE: shr_flux_atmIce -- computes atm/ice fluxes
!
! !DESCRIPTION:
!    Computes atm/ice fluxes
!
! !REVISION HISTORY:
!    2006-Jun-12 - B. Kauffman, first version, adapted from dice6 code
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_flux_atmIce(mask  ,zbot  ,ubot  ,vbot  ,thbot  &
               &          ,qbot  ,rbot  ,tbot  ,ts    ,sen    &
               &          ,lat   ,lwup  ,evap  ,taux  ,tauy   &
               &          ,tref  ,qref                        )

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   !--- input arguments --------------------------------
   integer(IN),intent(in)  :: mask (:)    ! 0 <=> cell NOT in model domain
   real(R8)   ,intent(in)  :: zbot (:)    ! atm level height  (m)
   real(R8)   ,intent(in)  :: ubot (:)    ! atm u wind     (m/s)
   real(R8)   ,intent(in)  :: vbot (:)    ! atm v wind     (m/s)
   real(R8)   ,intent(in)  :: thbot(:)    ! atm potential T   (K)
   real(R8)   ,intent(in)  :: qbot (:)    ! atm specific humidity (kg/kg)
   real(R8)   ,intent(in)  :: rbot (:)    ! atm air density   (kg/m^3)
   real(R8)   ,intent(in)  :: tbot (:)    ! atm T       (K) 
   real(R8)   ,intent(in)  :: ts   (:)    ! surface temperature

   !--- output arguments -------------------------------
   real(R8)   ,intent(out) :: sen  (:)    ! sensible      heat flux  (W/m^2)
   real(R8)   ,intent(out) :: lat  (:)    ! latent        heat flux  (W/m^2)
   real(R8)   ,intent(out) :: lwup (:)    ! long-wave upward heat flux  (W/m^2)
   real(R8)   ,intent(out) :: evap (:)    ! evaporative water flux ((kg/s)/m^2)
   real(R8)   ,intent(out) :: taux (:)    ! x surface stress (N)
   real(R8)   ,intent(out) :: tauy (:)    ! y surface stress (N)
   real(R8)   ,intent(out) :: tref (:)    ! 2m reference height temperature
   real(R8)   ,intent(out) :: qref (:)    ! 2m reference height humidity
 
!EOP

   !--- local constants --------------------------------
   real(R8),parameter :: umin   =  1.0_R8            ! minimum wind speed (m/s)
   real(R8),parameter :: zref   = 10.0_R8            ! ref height           ~ m
   real(R8),parameter :: ztref  =  2.0_R8            ! ref height for air T ~ m
   real(R8),parameter :: spval  = shr_const_spval    ! special value
   real(R8),parameter :: g      = shr_const_g        ! gravity
   real(R8),parameter :: cpdair = shr_const_cpdair   ! spec heat of dry air
   real(R8),parameter :: cpvir  = shr_const_cpvir    ! cpwv/cpdair - 1.0
   real(R8),parameter :: zvir   = shr_const_zvir     ! rh2o/rair   - 1.0
   real(R8),parameter :: latvap = shr_const_latvap   ! latent heat of evap
   real(R8),parameter :: latice = shr_const_latice   ! latent heat of fusion
   real(R8),parameter :: stebol = shr_const_stebol   ! Stefan-Boltzmann
   real(R8),parameter :: karman = shr_const_karman   ! Von Karman constant
   real(R8),parameter :: zzsice = 0.0005_R8          ! ice surface roughness

   !--- local variables --------------------------------
   integer(IN) :: lsize  ! array dimensions
   integer(IN) :: n      ! array indicies
   real(R8)    :: vmag   ! surface wind magnitude   (m/s)
   real(R8)    :: thvbot ! virtual temperature      (K)
   real(R8)    :: ssq    ! sea surface humidity     (kg/kg)
   real(R8)    :: dssqdt ! derivative of ssq wrt Ts (kg/kg/K)
   real(R8)    :: delt   ! potential T difference   (K)
   real(R8)    :: delq   ! humidity difference      (kg/kg)
   real(R8)    :: stable ! stability factor
   real(R8)    :: rdn    ! sqrt of neutral exchange coefficient (momentum)
   real(R8)    :: rhn    ! sqrt of neutral exchange coefficient (heat)
   real(R8)    :: ren    ! sqrt of neutral exchange coefficient (water)
   real(R8)    :: rd     ! sqrt of exchange coefficient (momentum)
   real(R8)    :: rh     ! sqrt of exchange coefficient (heat)
   real(R8)    :: re     ! sqrt of exchange coefficient (water)      
   real(R8)    :: ustar  ! ustar
   real(R8)    :: qstar  ! qstar
   real(R8)    :: tstar  ! tstar
   real(R8)    :: hol    ! H (at zbot) over L
   real(R8)    :: xsq    ! temporary variable
   real(R8)    :: xqq    ! temporary variable
   real(R8)    :: psimh  ! stability function at zbot (momentum)
   real(R8)    :: psixh  ! stability function at zbot (heat and water)
   real(R8)    :: alz    ! ln(zbot/z10)
   real(R8)    :: ltheat ! latent heat for surface
   real(R8)    :: tau    ! stress at zbot
   real(R8)    :: cp     ! specific heat of moist air

   real(R8)    :: bn     ! exchange coef funct for interpolation
   real(R8)    :: bh     ! exchange coef funct for interpolation
   real(R8)    :: fac    ! interpolation factor
   real(R8)    :: ln0    ! log factor for interpolation
   real(R8)    :: ln3    ! log factor for interpolation

   !--- local functions --------------------------------
   real(R8)   :: Tk      ! temperature (K)
   real(R8)   :: qsat    ! the saturation humididty of air (kg/m^3)
   real(R8)   :: dqsatdt ! derivivative of qsat wrt surface temperature
   real(R8)   :: xd      ! dummy argument  
   real(R8)   :: psimhu  ! unstable part of psimh
   real(R8)   :: psixhu  ! unstable part of psimx

   qsat(Tk)    = 627572.4_R8 / exp(5107.4_R8/Tk)
   dqsatdt(Tk) = (5107.4_R8 / Tk**2) * 627572.4_R8 / exp(5107.4_R8/Tk)
   psimhu(xd)  = log((1.0_R8+xd*(2.0_R8+xd))*(1.0_R8+xd*xd)/8.0_R8) &
              & -2.0_R8*atan(xd) +1.571_R8
   psixhu(xd)  =  2.0_R8 * log((1.0_R8 + xd*xd)/2.0_R8)

   !--- formats ----------------------------------------
   character(*),parameter :: subName =  "(shr_flux_atmIce) "

!-------------------------------------------------------------------------------
! PURPOSE:
!   using atm & ice state variables, compute atm/ice fluxes 
!   and diagnostic 10m air temperature and humidity
!
! NOTE: 
!   o all fluxes are positive downward
!   o net heat flux = net sw + lw up + lw down + sen + lat
!   o here, tstar = <WT>/U*, and qstar = <WQ>/U*.
!   o wind speeds should all be above a minimum speed (eg. 1.0 m/s)
! 
! ASSUME:
!   o The saturation humidity of air at T(K): qsat(T)  (kg/m^3)
!-------------------------------------------------------------------------------

   lsize = size(tbot)

   do n = 1,lsize

     if (mask(n) == 0) then
        sen  (n) = spval
        lat  (n) = spval
        lwup (n) = spval
        evap (n) = spval
        taux (n) = spval
        tauy (n) = spval
        tref (n) = spval
        qref (n) = spval
     else
        !--- define some needed variables ---
        vmag   = max(umin, sqrt(ubot(n)**2+vbot(n)**2))
        thvbot = thbot(n)*(1.0_R8 + zvir * qbot(n)) ! virtual pot temp (K)
         ssq   =  qsat  (ts(n)) / rbot(n)           ! sea surf hum (kg/kg)
        dssqdt = dqsatdt(ts(n)) / rbot(n)           ! deriv of ssq wrt Ts 
        delt   = thbot(n) - ts(n)                   ! pot temp diff (K)
        delq   = qbot(n) - ssq                        ! spec hum dif (kg/kg)
        alz    = log(zbot(n)/zref) 
        cp     = cpdair*(1.0_R8 + cpvir*ssq) 
        ltheat = latvap + latice

        !----------------------------------------------------------
        ! first estimate of Z/L and ustar, tstar and qstar
        !----------------------------------------------------------

        !--- neutral coefficients, z/L = 0.0 ---
        rdn = karman/log(zref/zzsice)
        rhn = rdn
        ren = rdn

        !--- ustar,tstar,qstar ----
        ustar = rdn * vmag
        tstar = rhn * delt  
        qstar = ren * delq  

        !--- compute stability & evaluate all stability functions ---
        hol    = karman * g * zbot(n) &
        &     * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
        hol    = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
!tht 15/11/2013 5->4 as in UKMO UM?
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!       psimh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
!       psixh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8+rdn/karman*(alz-psimh))
        rh = rhn / (1.0_R8+rhn/karman*(alz-psixh))
        re = ren / (1.0_R8+ren/karman*(alz-psixh))

        !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

        !----------------------------------------------------------
        ! iterate to converge on Z/L, ustar, tstar and qstar
        !----------------------------------------------------------

        !--- compute stability & evaluate all stability functions ---
        hol    = karman * g * zbot(n) &
        &      * (tstar/thvbot+qstar/(1.0_R8/zvir+qbot(n))) / ustar**2
        hol    = sign( min(abs(hol),10.0_R8), hol )
        stable = 0.5_R8 + sign(0.5_R8 , hol)
        xsq    = max(sqrt(abs(1.0_R8 - 16.0_R8*hol)) , 1.0_R8)
        xqq    = sqrt(xsq)
!tht 15/11/2013 5->4 as in UKMO UM?
        psimh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
        psixh  = -5.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!       psimh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psimhu(xqq)
!       psixh  = -4.0_R8*hol*stable + (1.0_R8-stable)*psixhu(xqq)
!tht 15/11/2013 5->4 as in UKMO UM?

        !--- shift all coeffs to measurement height and stability ---
        rd = rdn / (1.0_R8+rdn/karman*(alz-psimh)) 
        rh = rhn / (1.0_R8+rhn/karman*(alz-psixh)) 
        re = ren / (1.0_R8+ren/karman*(alz-psixh)) 

        !--- update ustar, tstar, qstar w/ updated, shifted coeffs --
        ustar = rd * vmag 
        tstar = rh * delt 
        qstar = re * delq 

        !----------------------------------------------------------
        ! compute the fluxes
        !----------------------------------------------------------

        tau = rbot(n) * ustar * ustar 
    
        !--- momentum flux ---
        taux(n) = tau * ubot(n) / vmag 
        tauy(n) = tau * vbot(n) / vmag 
     
        !--- heat flux ---
        sen (n) =   cp * tau * tstar / ustar 
        lat (n) =  ltheat * tau * qstar / ustar
        lwup(n) = -stebol * ts(n)**4 
     
        !--- water flux ---
        evap(n) = lat(n)/ltheat 

        !----------------------------------------------------------
        ! compute diagnostic: 2m reference height temperature
        !----------------------------------------------------------

        !--- Compute function of exchange coefficients. Assume that 
        !--- cn = rdn*rdn, cm=rd*rd and ch=rh*rd, and therefore 
        !--- 1/sqrt(cn(n))=1/rdn and sqrt(cm(n))/ch(n)=1/rh 
        bn = karman/rdn
        bh = karman/rh

        !--- Interpolation factor for stable and unstable cases
        ln0 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn) - 1.0_R8))
        ln3 = log(1.0_R8 + (ztref/zbot(n))*(exp(bn - bh) - 1.0_R8))
        fac = (ln0 - ztref/zbot(n)*(bn - bh))/bh * stable &
        &   + (ln0 - ln3)/bh * (1.0_R8-stable)
        fac = min(max(fac,0.0_R8),1.0_R8)

        !--- actual interpolation
        tref(n) = ts(n) + (tbot(n) - ts(n))*fac
        qref(n) = qbot(n) - delq*fac

     endif
   enddo 

end subroutine shr_flux_atmIce
 
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: shr_flux_MOstability -- Monin-Obukhov BL stability functions
!
! !DESCRIPTION:
!
!    Monin-Obukhov boundary layer stability functions, two options:
!    turbulent velocity scales or gradient and integral functions
!    via option = shr_flux_MOwScales or shr_flux_MOfunctions
!     
! !REVISION HISTORY:
!    2007-Sep-19 - B. Kauffman, Bill Large - first version
!
!===============================================================================
!   THIS ROUTINE NOT USED, AND PROBABLY FOR THE BEST (tht 22/11/2013)
!===============================================================================
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_flux_MOstability(option,arg1,arg2,arg3,arg4,arg5)

! !USES:

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(IN),intent(in)           :: option ! shr_flux_MOwScales or MOfunctions 
   real(R8)   ,intent(in)           :: arg1   ! scales: uStar (in)  funct: zeta (in)
   real(R8)   ,intent(inout)        :: arg2   ! scales: zkB   (in)  funct: phim (out)
   real(R8)   ,intent(out)          :: arg3   ! scales: phim  (out) funct: phis (out)
   real(R8)   ,intent(out)          :: arg4   ! scales: phis  (out) funct: psim (out)
   real(R8)   ,intent(out),optional :: arg5   ! scales:    (unused) funct: psis (out)
 
! !EOP

   !----- local variables -----
   real(R8)           :: zeta  ! z/L
   real(R8)           :: uStar ! friction velocity
   real(R8)           :: zkB   ! (height)*(von Karman)*(surface bouyancy flux)
   real(R8)           :: phim  ! momentum    gradient function or scale
   real(R8)           :: phis  ! temperature gradient function or scale
   real(R8)           :: psim  ! momentum    integral function or scale
   real(R8)           :: psis  ! temperature integral function or scale
   real(R8)           :: temp  ! temporary-variable/partial calculation

   !----- local variables, stable case -----
   real(R8),parameter :: uStarMin = 0.001_R8 ! lower bound on uStar
   real(R8),parameter :: a = 1.000_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: b = 0.667_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: c = 5.000_R8  ! constant from Holtslag & de Bruin, equation 12
   real(R8),parameter :: d = 0.350_R8  ! constant from Holtslag & de Bruin, equation 12

   !----- local variables, unstable case -----
   real(R8),parameter :: a2 = 3.0_R8   ! constant from Wilson, equation 10 

   !----- formats -----
   character(*),parameter :: subName = '(shr_flux_MOstability) '
   character(*),parameter ::   F00 = "('(shr_flux_MOstability) ',4a)"
   character(*),parameter ::   F01 = "('(shr_flux_MOstability) ',a,i5)"

!-------------------------------------------------------------------------------
! Notes::
!   o this could be two routines, but are one to help keep them aligned
!   o the stable calculation is taken from...
!     A.A.M. HoltSlag and H.A.R. de Bruin, 1988:
!     "Applied Modeling of the Nighttime Surface Energy Balance over Land",
!     Journal of Applied Meteorology, Vol. 27, No. 6, June 1988, 659-704
!   o the unstable calculation is taken from...
!     D. Keith Wilson, 2001: "An Alternative Function for the Wind and 
!     Temperature Gradients in Unstable Surface Layers", 
!     Boundary-Layer Meteorology, 99 (2001), 151-158
!-------------------------------------------------------------------------------

   !----- check for consistancy between option and arguments ------------------
   if (debug > 1 .and. s_loglev > 0) then
      if (debug > 2) write(s_logunit,F01) "enter, option = ",option
      if ( option == shr_flux_MOwScales .and. present(arg5) ) then
         write(s_logunit,F01) "ERROR: option1 must have four arguments"
         call shr_sys_abort(subName//"option inconsistant with arguments")
      else if ( option == shr_flux_MOfunctions .and. .not. present(arg5) ) then
         write(s_logunit,F01) "ERROR: option2 must have five arguments"
         call shr_sys_abort(subName//"option inconsistant with arguments")
      else 
         write(s_logunit,F01) "invalid option = ",option
         call shr_sys_abort(subName//"invalid option")
      end if
   end if

   !------ velocity scales option ----------------------------------------------
   if (option == shr_flux_MOwScales) then

      !--- input ---
      uStar = arg1
      zkB   = arg2

      if (zkB >= 0.0_R8) then ! ----- stable -----
         zeta = zkB/(max(uStar,uStarMin)**3)
         temp = exp(-d*zeta)
         phim = uStar/(1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp))
         phis = phim
      else                    ! ----- unstable -----
         temp = (zkB*zkB)**(1.0_R8/a2)   ! note: zkB < 0, zkB*zkB > 0
         phim = sqrt(uStar**2 + shr_flux_MOgammaM*temp)
         phis = sqrt(uStar**2 + shr_flux_MOgammaS*temp)
      end if

      !--- output ---
      arg3 = phim
      arg4 = phis
   !  arg5 = <unused>

   !------ stability function option -------------------------------------------
   else if (option == shr_flux_MOfunctions) then

      !--- input ---
      zeta  = arg1

      if (zeta >= 0.0_R8) then ! ----- stable -----
         temp = exp(-d*zeta)
         phim =        1.0_R8 + zeta*(a + b*(1.0_R8 + c - d*zeta)*temp)
         phis = phim
         psim = -a*zeta - b*(zeta - c/d)*temp - b*c/d
         psis = psim 
      else                    ! ----- unstable ----
         temp = (zeta*zeta)**(1.0_R8/a2)   ! note: zeta < 0, zeta*zeta > 0
         phim = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaM*temp)
         phis = 1.0_R8/sqrt(1.0_R8 + shr_flux_MOgammaS*temp)
         psim = a2*log(0.5_R8 + 0.5_R8/phim)
         psis = a2*log(0.5_R8 + 0.5_R8/phis)
      end if

      !--- output ---
      arg2 = phim
      arg3 = phis
      arg4 = psim
      arg5 = psis
   !----------------------------------------------------------------------------
   else 
      write(s_logunit,F01) "invalid option = ",option
      call shr_sys_abort(subName//"invalid option")
   endif

end subroutine shr_flux_MOstability


!===============================================================================
!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: cor30a 
!
! !DESCRIPTION:
!
!   COARE v3.0 parametrisation
!     
! !REVISION HISTORY:
!   2013-Nov-22: Thomas Toniazzo's adaptation of Chris Fairall's code,
!    downloaded from
!    ftp://ftp1.esrl.noaa.gov/users/cfairall/wcrp_wgsf/computer_programs/cor3_0/
!     * no wave, standard coare 2.6 charnock
!     * skin parametrisation also off (would require radiative fluxes and
!      rainrate in input) 
!     * added diagnostics, comments and references
!===============================================================================
!
! !INTERFACE: ------------------------------------------------------------------

subroutine cor30a(ubt,vbt,tbt,qbt,rbt        &    ! in atm params		
               & ,uss,vss,tss,qss            &    ! in surf params		
               & ,zbl,zbu,zbt,zrfu,zrfq,zrft &    ! in heights			
               & ,tau,hsb,hlb                &    ! out: fluxes			
               & ,zo,zot,zoq,L,usr,tsr,qsr   &    ! out: ss scales		
               & ,Cd,Ch,Ce                   &    ! out: exch. coeffs		
               & ,trf,qrf,urf,vrf)                ! out: reference-height params

! !USES:

IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

real(R8),intent(in) :: ubt,vbt,tbt,qbt,rbt,uss,vss,tss,qss
real(R8),intent(in) :: zbl,zbu,zbt,zrfu,zrfq,zrft
real(R8),intent(out):: tau,hsb,hlb,zo,zot,zoq,L,usr,tsr,qsr,Cd,Ch,Ce &
                    & ,trf,qrf,urf,vrf
! !EOP

real ua,va,ta,q,rb,us,vs,ts,qs,zi,zu,zt,zq,zru,zrq,zrt ! internal vars

real(R8):: cpa,rgas,grav,pi,von,beta ! phys. params
real(R8):: le,rhoa,cpv               ! derived phys. params
real(R8):: t,visa,du,dq,dt           ! params of problem

real(R8):: u10,zo10,zot10,cd10,ch10,ct10,ct,cc,ribu,zetu,l10,charn ! init vars
real(R8):: zet,rr,bf,ug,ut     ! loop iter vars
real(R8):: cdn_10,chn_10,cen_10  ! aux. output vars

integer(IN):: i,nits ! iter loop counters

integer(IN):: jcool                  ! aux. cool-skin vars 
real(R8):: dter,wetc,dqer 


ua=ubt  !wind components (m/s) at height zu (m)
va=vbt
ta=tbt  !bulk air temperature (K), height zt
Q =qbt  !bulk air spec hum (kg/kg), height zq
rb=rbt  ! air density
us=uss  !surface current components (m/s)
vs=vss
ts=tss  !bulk water temperature (K) if jcool=1, interface water T if jcool=0  
qs=qss  !bulk water spec hum (kg/kg) if jcool=1 etc
zi=zbl  !PBL depth (m)
zu=zbu  !wind speed measurement height (m)
zt=zbt  !air T measurement height (m)
zq=zbt  !air q measurement height (m)
zru=zrfu ! reference height for st.diagn.U
zrq=zrfq ! reference height for st.diagn.T,q
zrt=zrft ! reference height for st.diagn.T,q

!**** constants
    Beta= 1.2 
    von = 0.4 
    pi  = 3.141593
    grav= SHR_CONST_G 
    Rgas= SHR_CONST_RGAS
    cpa = SHR_CONST_CPDAIR  

!*** physical parameters
    Le  = SHR_CONST_LATVAP -.00237e6*(ts-273.16)
!   cpv = shr_const_cpdair*(1.0_R8 + shr_const_cpvir*Qs) ! form in NCAR code
    cpv = cpa*(1+0.84*Q) 
!   rhoa= P/(Rgas*ta*(1+0.61*Q)) ! if input were pressure
    rhoa= rb

! parametrisation for air kinematic viscosity (Andreas 1989,p.31)
    t   = ta-273.16
    visa= 1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t) 

    du  = sqrt((ua-us)**2+(va-vs)**2)
    dt  = ts-ta -.0098*zt 
    dq  = Qs-Q 

!*** don't use cool-skin params for now, but assign values to Ter and Qer
jcool=0
    dter=0.3  
    wetc=0.622*Le*Qs/(Rgas*ts**2) 
    dqer=wetc*dter 

!***************** Begin bulk-model calculations ***************
     
!*************** first guess 
    ug=.5 

    ut   = sqrt(du*du+ug*ug) 
    u10  = ut*log(10/1e-4)/log(zu/1e-4) 
    usr  = .035*u10 
    zo10 = 0.011*usr*usr/grav+0.11*visa/usr 
    Cd10 = (von/log(10/zo10))**2 
    Ch10 = 0.00115 
    Ct10 = Ch10/sqrt(Cd10) 
    zot10= 10/exp(von/Ct10) 
    Cd   =(von/log(zu/zo10))**2 
    Ct   = von/log(zt/zot10) 
    CC   = von*Ct/Cd 

! Bulk Richardson number
    Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/ut**2 
! initial guess for stability parameter...
    if (Ribu .LT. 0) then 
    ! pbl-height dependent
        zetu=CC*Ribu/( 1- (.004*Beta**3*zi/zu) * Ribu ) 
    else 
        zetu=CC*Ribu*(1+27/9*Ribu/CC)
    endif 
! ...and MO length
    L10=zu/zetu 

    if (zetu .GT. 50) then 
        nits=1 
    else 
        nits=3 
    endif 

    usr =  ut*von/(log(zu/zo10)-psiuo(zu/L10))
    tsr = (dt-dter*jcool)*von/(log(zt/zot10)-psit_30(zt/L10)) 
    qsr = (dq-dqer*jcool)*von/(log(zq/zot10)-psit_30(zq/L10)) 

! parametrisation for Charney parameter (section 3c of Fairall et al. 2003)
    charn=0.011 
    if (ut .GT. 10) then
      charn=0.011+(ut-10)/(18-10)*(0.018-0.011) 
    endif 
    if (ut .GT. 18) then
      charn=0.018 
    endif 
        
!***************  iteration loop ************
    do i=1, nits 
     
     ! stability parameter
     zet=-von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q) 

     ! momentum roughness length...
     zo = charn*usr*usr/grav+0.11*visa/usr  
     ! ...& MO length
     L  = zu/zet 

     ! tracer roughness length
     rr = zo*usr/visa 
     zoq= min(1.15e-4,5.5e-5/rr**.6) 
     zot= zoq ! N.B. same for vapour and heat

     ! new surface-layer scales
     usr =  ut            *von/(log(zu/zo )-psiuo(zu/L)) 
     tsr = (dt-dter*jcool)*von/(log(zt/zot)-psit_30(zt/L)) 
     qsr = (dq-dqer*jcool)*von/(log(zq/zoq)-psit_30(zq/L)) 

     ! gustiness parametrisation
     Bf=-grav/ta*usr*(tsr+.61*ta*qsr) 
     if (Bf .GT. 0) then
       ug=Beta*(Bf*zi)**.333 
     else
       ug=.2 
     endif
     ut=sqrt(du*du+ug*ug) 

    enddo 
!***************     end loop    ************

   !******** fluxes @ measurement heights zu,zt,zq ********
   tau= rhoa*usr*usr*du/ut                !stress magnitude
   hsb=-rhoa*cpa*usr*tsr                  !heat downwards
   hlb=-rhoa*Le*usr*qsr                   !wv downwards

   !****** transfer coeffs relative to ut @meas. hts ******
   Cd= tau/rhoa/ut/max(.1,du) 
   Ch= usr/ut*tsr/(dt-dter*jcool) 
   Ce= usr/ut*qsr/(dq-dqer*jcool)

   !**********  10-m neutral coeff relative to ut *********
   Cdn_10=von*von/log(10/zo)/log(10/zo) 
   Chn_10=von*von/log(10/zo)/log(10/zot) 
   Cen_10=von*von/log(10/zo)/log(10/zoq) 

   !**********  reference-height values for u,q,T *********
   urf=us+(ua-us)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
   vrf=vs+(va-vs)*(log(zru/zo)-psiuo(zru/L))/(log(zu/zo)-psiuo(zu/L))
   qrf=qs-dq*(log(zrq/zoq)-psit_30(zrq/L))/(log(zq/zoq)-psit_30(zq/L))
   trf=ts-dt*(log(zrt/zot)-psit_30(zrt/L))/(log(zt/zot)-psit_30(zt/L))
   trf=trf+.0098*zrt

end subroutine cor30a


!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: PSIUo
!
! !DESCRIPTION:
!
!   momentum stability functions adopted in COARE v3.0 parametrisation.
!   Chris Fairall's code (see cor30a)
!
! !REVISION HISTORY:
!   22/11/2013: Thomas Toniazzo: comments added
!
! !INTERFACE: ------------------------------------------------------------------
real (SHR_KIND_R8) function psiuo(zet)
! !INPUT/OUTPUT PARAMETERS:
real(R8),intent(in)  :: zet
! !EOP
real(R8) ::c,x,psik,psic,f
!-----------------------------------------------------------------
! N.B.: z0/L always neglected compared to z/L and to 1
!-----------------------------------------------------------------
    if(zet>0)then 
! Beljaars & Holtslag (1991)
     c=min(50.,.35*zet) 
     psiuo=-((1+1.0*zet)**1.0+.667*(zet-14.28)/exp(c)+8.525)
    else 
! Dyer & Hicks (1974) for weak instability
     x=(1.-15.*zet)**.25                   ! 15 instead of 16
     psik=2.*log((1.+x)/2.)+log((1.+x*x)/2.)-2.*atan(x)+2.*atan(1.) 
! Fairall et al. (1996) for strong instability (Eq.(13))
     x=(1.-10.15*zet)**.3333 
     psic= 1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.)) &
         & +4.*atan(1.)/sqrt(3.) 
     f=zet*zet/(1+zet*zet) 
     psiuo=(1-f)*psik+f*psic                                                
    endif 
END FUNCTION psiuo 


!===============================================================================
! !BOP =========================================================================
!
! !IROUTINE: PSIT_30
!
! !DESCRIPTION:
!
!   momentum stability functions adopted in COARE v3.0 parametrisation.
!   Chris Fairall's code (see cor30a)
!
! !REVISION HISTORY:
!   22/11/2013: Thomas Toniazzo: comments added
!
! !INTERFACE: ------------------------------------------------------------------
real (SHR_KIND_R8) function psit_30(zet)
! !INPUT/OUTPUT PARAMETERS:
real(R8),intent(in)  :: zet
! !EOP
real(R8) ::c,x,psik,psic,f
!-----------------------------------------------------------------
! N.B.: z0/L always neglected compared to z/L and to 1
!-----------------------------------------------------------------
    if(zet>0)then 
! Beljaars & Holtslag (1991)
     c=min(50.,.35*zet) 
     psit_30=-((1.+2./3.*zet)**1.5+.667*(zet-14.28)/exp(c)+8.525)
    else 
! Dyer & Hicks (1974) for weak instability
     x=(1.-15.*zet)**.5                    ! 15 instead of 16
     psik=2*log((1+x)/2) 
! Fairall et al. (1996) for strong instability
     x=(1.-(34.15*zet))**.3333 
     psic= 1.5*log((1.+x+x*x)/3.)-sqrt(3.)*atan((1.+2.*x)/sqrt(3.)) &
         & +4.*atan(1.)/sqrt(3.) 
     f=zet*zet/(1+zet*zet) 
     psit_30=(1-f)*psik+f*psic
   endif
end FUNCTION psit_30


!===============================================================================
!===============================================================================
!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_map_docoare-- set docoare flag
!
! !DESCRIPTION:
!     set docoare flag
!     \newline
!     call shr\_map\_setDopole(flag)
!
! !REVISION HISTORY:
!     2009-Jun-22 - T. Craig - first version
!
! !INTERFACE: ------------------------------------------------------------------

subroutine shr_map_docoare(iflag)

  implicit none

! !INPUT/OUTPUT PARAMETERS:

  integer, intent(in) :: iflag


  flux_scheme = iflag

end subroutine shr_map_docoare

end module shr_flux_mod
