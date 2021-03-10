  module stratiform

  !-------------------------------------------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the CAM interface to the prognostic cloud macro and microphysics
  !
  ! Author: Byron Boville  Sept 04, 2002
  !         modified by D.B. Coleman May 2004
  !         modified by Sungsu Park. Dec.2009
  !         modified by J. Kay Jan. 2010 to add COSP simulator info from RK microphysics to physics buffer, Nov 2010 for more COSP vars
  !         modified by A. Gettelman and C. Craig Nov 2010 to remove MG microphysics
  !         Modified by Øyvind Seland and Alf Kirkevåg, September 2011. 
  ! Note for DIAGNCDNC (and DIRIND version) of CAM4-Oslo:
  ! This is old code for diagnostic CCN-calculations. It is kept here as a possibility to address
  ! at a later stage (e.g. for testing ourposes). In that case, the code needs to be re-structured 
  ! and tested, since this option has not been used or checked for some time.
  !-------------------------------------------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use ppgrid,        only: pcols, pver, pverp
  use physconst,     only: gravit, latvap, latice
  use abortutils,    only: endrun
  use chemistry,     only: chem_is
  use phys_control,  only: phys_getopts

  use perf_mod
  use cam_logfile,   only: iulog

  implicit none
  private
  save

  public :: stratiform_register, stratiform_init_cnst, stratiform_implements_cnst
  public :: stratiform_init
  public :: stratiform_tend

  ! ------------------------- !
  ! Private Module Parameters !
  ! ------------------------- !

  ! Choose either 'intermediate' ('inter') or complete ('compl') cloud microphysics 
  ! inter : Microphysics assumes 'liquid stratus frac = ice stratus frac = max( liquid stratus frac, ice stratus frac )'.
  ! compl : Microphysics explicitly treats 'liquid stratus frac .ne. ice stratus frac'  
  ! for CAM5, only 'inter' is functional

    character(len=5), private, parameter :: micro_treatment = 'inter' 

  ! 'cu_det_st' : If .true. (.false.), detrain cumulus liquid condensate into the pre-existing liquid stratus 
  !               (environment) without (with) macrophysical evaporation. If there is no pre-esisting stratus, 
  !               evaporate cumulus liquid condensate. This option only influences the treatment of cumulus
  !               liquid condensate, not cumulus ice condensate.

    logical,          private, parameter :: cu_det_st  = .false.  

  ! -------------------------------- !
  ! End of Private Module Parameters !
  ! -------------------------------- !

! Physics buffer indices 
integer  ::  qcwat_idx          = 0 
integer  ::  lcwat_idx          = 0 
#ifdef DIRIND
integer  ::  lcwatx_idx         = 0 
#endif
integer  ::  iccwat_idx         = 0 
integer  ::  nlwat_idx          = 0 
integer  ::  niwat_idx          = 0 
integer  ::  cc_t_idx           = 0
integer  ::  cc_qv_idx          = 0
integer  ::  cc_ql_idx          = 0
integer  ::  cc_qi_idx          = 0
integer  ::  cc_nl_idx          = 0
integer  ::  cc_ni_idx          = 0
integer  ::  cc_qlst_idx        = 0
integer  ::  tcwat_idx          = 0 
integer  ::  cld_idx            = 0 
integer  ::  cldo_idx           = 0  
integer  ::  ast_idx            = 0 
integer  ::  aist_idx           = 0 
integer  ::  alst_idx           = 0 
integer  ::  qist_idx           = 0 
integer  ::  qlst_idx           = 0 
integer  ::  concld_idx         = 0 
integer  ::  rhdfda_idx         = 0 
integer  ::  rhu00_idx          = 0 
integer  ::  rel2_idx           = 0  
integer  ::  rei2_idx           = 0
integer  ::  concldql_idx       = 0 
integer  ::  fice_idx           = 0 
integer  ::  sh_frac_idx        = 0 
integer  ::  dp_frac_idx        = 0 
integer  ::  qini_idx           = 0 
integer  ::  cldliqini_idx      = 0 
integer  ::  cldiceini_idx      = 0 
integer  ::  tini_idx           = 0 
integer  ::  qme_idx            = 0 
integer  ::  prain_idx          = 0 
integer  ::  nevapr_idx         = 0 
integer  ::  wsedl_idx          = 0 
integer  ::  rei_idx            = 0 
integer  ::  rel_idx            = 0 
integer  ::  rel_fn_idx         = 0 
integer  ::  dei_idx            = 0 
integer  ::  mu_idx             = 0 
integer  ::  lambdac_idx        = 0 
integer  ::  iciwp_idx          = 0 
integer  ::  iclwp_idx          = 0 
integer  ::  deiconv_idx        = 0 
integer  ::  muconv_idx         = 0 
integer  ::  lambdaconv_idx     = 0
integer  ::  iciwpst_idx        = 0
integer  ::  iclwpst_idx        = 0
integer  ::  iciwpconv_idx      = 0
integer  ::  iclwpconv_idx      = 0
integer  ::  des_idx            = 0
integer  ::  icswp_idx          = 0 
integer  ::  cldfsnow_idx       = 0
integer  ::  dgnumwet_idx       = 0
integer  ::  dgnum_idx          = 0 
integer  ::  tke_idx            = 0
integer  ::  turbtype_idx       = 0
integer  ::  smaw_idx           = 0
integer  ::  rate1_cw2pr_st_idx = 0 
integer  ::  shfrc_idx          = 0

#ifdef DIRIND
  integer, parameter :: ncnstmax = 4                    ! Number of constituent
#else
  integer, parameter :: ncnstmax = 2                    ! Number of constituent
#endif
  integer            :: ncnst 	  		        ! Number of constituents (can vary)
#ifdef DIRIND
  character(len=8), dimension(ncnstmax), parameter &    ! Constituent names
    :: cnst_names = (/'CLDLIQ', 'CLDICE','CLDLIX','CDNCIX'/)
#else
  character(len=8), dimension(ncnstmax), parameter &    ! Constituent names
                     :: cnst_names = (/'CLDLIQ', 'CLDICE'/)
#endif
  logical            :: use_shfrc                       ! Local copy of flag from convect_shallow_use_shfrc
  character(len=16)  :: microp_scheme                   ! Microphysics scheme

  integer :: &
     kvh_idx,      &! kvh index in physics buffer
     ixcldliq,     &! cloud liquid amount index
     ixcldice,     &! cloud ice amount index
#ifdef DIRIND
     icldlx,       &! cloud liquid water index offline
     icldnx,       &! Cloud droplet number conc. index
#endif
     ixnumliq,     &! cloud liquid number index
     ixnumice,     &! cloud ice water index
     ls_flxprc_idx,&
     ls_flxsnw_idx


  contains

  ! ===============================================================================

  subroutine stratiform_register

  !---------------------------------------------------------------------- !
  !                                                                       !
  ! Register the constituents (cloud liquid and cloud ice) and the fields !
  ! in the physics buffer.                                                !
  !                                                                       !
  !---------------------------------------------------------------------- !

    use constituents, only: cnst_add, pcnst
    use physconst,    only: mwdry, cpair
    use phys_buffer,  only: pbuf_times, pbuf_add


  !-----------------------------------------------------------------------

    call phys_getopts(microp_scheme_out=microp_scheme)

    if ( microp_scheme .eq. 'RK' ) then
#ifdef DIRIND
	ncnst = 4
#else
       ncnst  = 2
#endif
    end if

  ! Register cloud water and determine index.

    call cnst_add(cnst_names(1), mwdry, cpair, 0._r8, ixcldliq, &
       longname='Grid box averaged cloud liquid amount')
    call cnst_add(cnst_names(2), mwdry, cpair, 0._r8, ixcldice, &
       longname='Grid box averaged cloud ice amount')

#ifdef DIRIND
    call cnst_add(cnst_names(3), mwdry, cpair, 0._r8, icldlx, &
         longname='Aerosol calc averaged liquid condensate amount')
    call cnst_add(cnst_names(4), mwdry, cpair, 0._r8, icldnx, &
         longname='Aerosol calc cloud droplet number concentration')
#endif


  ! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add('QCWAT',   'global',  1, pver, pbuf_times,   qcwat_idx)
    call pbuf_add('LCWAT',   'global',  1, pver, pbuf_times,   lcwat_idx)
#ifdef DIRIND
    call pbuf_add('LCWATX',  'global',  1, pver, pbuf_times,  lcwatx_idx)
#endif
    call pbuf_add('ICCWAT',  'global',  1, pver, pbuf_times,  iccwat_idx)
    call pbuf_add('NLWAT',   'global',  1, pver, pbuf_times,   nlwat_idx)
    call pbuf_add('NIWAT',   'global',  1, pver, pbuf_times,   niwat_idx)
    call pbuf_add('CC_T',    'global',  1, pver, pbuf_times,    cc_t_idx)
    call pbuf_add('CC_qv',   'global',  1, pver, pbuf_times,   cc_qv_idx)
    call pbuf_add('CC_ql',   'global',  1, pver, pbuf_times,   cc_ql_idx)
    call pbuf_add('CC_qi',   'global',  1, pver, pbuf_times,   cc_qi_idx)
    call pbuf_add('CC_nl',   'global',  1, pver, pbuf_times,   cc_nl_idx)
    call pbuf_add('CC_ni',   'global',  1, pver, pbuf_times,   cc_ni_idx)
    call pbuf_add('CC_qlst', 'global',  1, pver, pbuf_times, cc_qlst_idx)
    call pbuf_add('TCWAT',   'global',  1, pver, pbuf_times,   tcwat_idx)
    call pbuf_add('CLD',     'global',  1, pver, pbuf_times,     cld_idx)
    call pbuf_add('CLDO',    'global',  1, pver, pbuf_times,    cldo_idx)
    call pbuf_add('AST',     'global',  1, pver, pbuf_times,     ast_idx)
    call pbuf_add('AIST',    'global',  1, pver, pbuf_times,    aist_idx)
    call pbuf_add('ALST',    'global',  1, pver, pbuf_times,    alst_idx)
    call pbuf_add('QIST',    'global',  1, pver, pbuf_times,    qist_idx)
    call pbuf_add('QLST',    'global',  1, pver, pbuf_times,    qlst_idx)
    call pbuf_add('CONCLD',  'global',  1, pver, pbuf_times,  concld_idx)
    call pbuf_add('RHDFDA',  'global',  1, pver, pbuf_times,  rhdfda_idx)
    call pbuf_add('RHU00',   'global',  1, pver, pbuf_times,   rhu00_idx)
    call pbuf_add('REL2',    'global',  1, pver, pbuf_times,    rel2_idx)
    call pbuf_add('REI2',    'global',  1, pver, pbuf_times,    rei2_idx)

  ! Physics buffer variables for convective cloud properties.

    call pbuf_add('CONCLDQL',   'physpkg', 1, pver, 1, concldql_idx) 
    call pbuf_add('FICE',       'physpkg', 1, pver, 1, fice_idx) 
    call pbuf_add('SH_FRAC',    'physpkg', 1, pver, 1, sh_frac_idx) 
    call pbuf_add('DP_FRAC',    'physpkg', 1, pver, 1, dp_frac_idx) 

    call pbuf_add('QINI'      , 'physpkg', 1,pver, 1, qini_idx)
    call pbuf_add('CLDLIQINI' , 'physpkg', 1,pver, 1, cldliqini_idx)
    call pbuf_add('CLDICEINI' , 'physpkg', 1,pver, 1, cldiceini_idx)
    call pbuf_add('TINI'      , 'physpkg', 1,pver, 1, tini_idx)

    call pbuf_add('QME',        'physpkg', 1, pver, 1, qme_idx)
    call pbuf_add('PRAIN' ,     'physpkg', 1, pver, 1, prain_idx)
    call pbuf_add('NEVAPR' ,    'physpkg', 1, pver, 1, nevapr_idx)

    call pbuf_add('WSEDL',      'physpkg', 1, pver, 1, wsedl_idx)

    call pbuf_add('REI',        'physpkg', 1, pver, 1, rei_idx)
    call pbuf_add('REL',        'physpkg', 1, pver, 1, rel_idx)
    call pbuf_add('REL_FN',     'physpkg', 1, pver, 1, rel_fn_idx)          ! REL at fixed number for indirect rad forcing

    call pbuf_add('DEI',        'physpkg', 1, pver, 1, dei_idx)          ! Mitchell ice effective diameter for radiation
    call pbuf_add('MU',         'physpkg', 1, pver, 1, mu_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('LAMBDAC',    'physpkg', 1, pver, 1, lambdac_idx)          ! Size distribution shape parameter for radiation
    call pbuf_add('ICIWP',      'physpkg', 1, pver, 1, iciwp_idx)          ! In cloud ice water path for radiation
    call pbuf_add('ICLWP',      'physpkg', 1, pver, 1, iclwp_idx)          ! In cloud liquid water path for radiation

    call pbuf_add('DEICONV',    'physpkg', 1, pver, 1, deiconv_idx)          ! Convective ice effective diameter for radiation
    call pbuf_add('MUCONV',     'physpkg', 1, pver, 1, muconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('LAMBDACONV', 'physpkg', 1, pver, 1, lambdaconv_idx)          ! Convective size distribution shape parameter for radiation
    call pbuf_add('ICIWPST',    'physpkg', 1, pver, 1, iciwpst_idx)          ! Stratiform only in cloud ice water path for radiation
    call pbuf_add('ICLWPST',    'physpkg', 1, pver, 1, iclwpst_idx)          ! Stratiform in cloud liquid water path for radiation
    call pbuf_add('ICIWPCONV',  'physpkg', 1, pver, 1, iciwpconv_idx)          ! Convective only in cloud ice water path for radiation
    call pbuf_add('ICLWPCONV',  'physpkg', 1, pver, 1, iclwpconv_idx)          ! Convective in cloud liquid water path for radiation

    call pbuf_add('DES',        'physpkg', 1, pver, 1, des_idx)          ! Snow effective diameter for radiation
    call pbuf_add('ICSWP',      'physpkg', 1, pver, 1, icswp_idx)          ! In cloud snow water path for radiation
    call pbuf_add('CLDFSNOW',   'physpkg', 1, pver ,pbuf_times, cldfsnow_idx) ! Cloud fraction for liquid drops + snow

    call pbuf_add('LS_FLXPRC',  'physpkg', 1, pverp, 1, ls_flxprc_idx)
    call pbuf_add('LS_FLXSNW',  'physpkg', 1, pverp, 1, ls_flxsnw_idx)

  end subroutine stratiform_register

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  function stratiform_implements_cnst(name)

  !----------------------------------------------------------------------------- ! 
  !                                                                              !    
  ! Purpose: return true if specified constituent is implemented by this package !
  !                                                                              !
  ! Author: B. Eaton                                                             !
  !                                                                              ! 
  !----------------------------------------------------------------------------- !
     implicit none
  !-----------------------------Arguments---------------------------------
     character(len=*), intent(in) :: name      ! constituent name
     logical :: stratiform_implements_cnst     ! return value
  !---------------------------Local workspace-----------------------------
     integer :: m
  !-----------------------------------------------------------------------

     stratiform_implements_cnst = .false.

     do m = 1, ncnst
        if (name == cnst_names(m)) then
           stratiform_implements_cnst = .true.
           return
        end if
     end do
  end function stratiform_implements_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_init_cnst(name, q, gcid)

  !----------------------------------------------------------------------- !
  !                                                                        !
  ! Initialize the cloud water mixing ratios (liquid and ice), if they are !
  ! not read from the initial file                                         ! 
  !                                                                        !
  !----------------------------------------------------------------------- !
    implicit none
  !---------------------------- Arguments ---------------------------------
    character(len=*), intent(in)  :: name     ! constituent name
    real(r8),         intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)
    integer,          intent(in)  :: gcid(:)  ! global column id
  !-----------------------------------------------------------------------

    if ( name == 'CLDLIQ' ) then
       q = 0.0_r8
       return
    else if ( name == 'CLDICE' ) then
       q = 0.0_r8
       return
#ifdef DIRIND
    else if ( name == 'CLDLIX' ) then
       q = 0.0_r8
       return
    else if ( name == 'CDNCIX' ) then
       q = 0.0_r8
       return
    end if
#else
    end if
#endif DIRIND



  end subroutine stratiform_init_cnst

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_init

  !-------------------------------------------- !
  !                                             !
  ! Initialize the cloud water parameterization !
  !                                             ! 
  !-------------------------------------------- !

    use microp_aero,     only: ini_microp_aero
    use cldwat2m_micro,  only: ini_micro
    use cldwat2m_macro,  only: ini_macro
    use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, sflxnam, apcnst, bpcnst
    use cam_history,     only: addfld, add_default, phys_decomp
    use physconst,       only: tmelt, rh2o, rhodair
    use phys_buffer,     only: pbuf_get_fld_idx
    use convect_shallow, only: convect_shallow_use_shfrc
    use dycore,          only: dycore_is
    use phys_control,    only: cam_physpkg_is

    integer              :: m, mm
    logical              :: history_aerosol      ! Output the MAM aerosol tendencies
    logical              :: history_microphysics ! Output variables for microphysics diagnostics package
    logical              :: history_budget       ! Output tendencies and state variables for CAM4
                                                 ! temperature, water vapor, cloud ice and cloud
                                                 ! liquid budgets.
    integer              :: history_budget_histfile_num ! output history file number for budget fields
  !-----------------------------------------------------------------------

    call phys_getopts( history_aerosol_out        = history_aerosol      , &
                       history_microphysics_out   = history_microphysics , & 
                       history_budget_out         = history_budget       , &
                       history_budget_histfile_num_out = history_budget_histfile_num)

  ! Find out whether shfrc from convect_shallow will be used in cldfrc

    if( convect_shallow_use_shfrc() ) then
        use_shfrc = .true.
        shfrc_idx =     pbuf_get_fld_idx('shfrc')
    else 
        use_shfrc = .false.
    endif

  ! Register history variables

    do m = 1, ncnst
       call cnst_get_ind( cnst_names(m), mm )
       call addfld( cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm)                   , phys_decomp )
       call addfld( sflxnam  (mm), 'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp )
       call add_default( cnst_name(mm), 1, ' ' )
       call add_default( sflxnam  (mm), 1, ' ' )
    enddo

    call addfld (apcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' after physics'  , phys_decomp)
    call addfld (apcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' after physics'  , phys_decomp)
    call addfld (bpcnst(ixcldliq), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldliq))//' before physics' , phys_decomp)
    call addfld (bpcnst(ixcldice), 'kg/kg   ', pver, 'A', trim(cnst_name(ixcldice))//' before physics' , phys_decomp)

    if( history_budget) then
       call add_default (cnst_name(ixcldliq), history_budget_histfile_num, ' ')
       call add_default (cnst_name(ixcldice), history_budget_histfile_num, ' ')
       call add_default (apcnst   (ixcldliq), history_budget_histfile_num, ' ')
       call add_default (apcnst   (ixcldice), history_budget_histfile_num, ' ')
       call add_default (bpcnst   (ixcldliq), history_budget_histfile_num, ' ')
       call add_default (bpcnst   (ixcldice), history_budget_histfile_num, ' ')
    end if

    call addfld ('FWAUT    ', 'fraction', pver, 'A', 'Relative importance of liquid autoconversion'            ,phys_decomp)
    call addfld ('FSAUT    ', 'fraction', pver, 'A', 'Relative importance of ice autoconversion'               ,phys_decomp)
    call addfld ('FRACW    ', 'fraction', pver, 'A', 'Relative importance of rain accreting liquid'            ,phys_decomp)
    call addfld ('FSACW    ', 'fraction', pver, 'A', 'Relative importance of snow accreting liquid'            ,phys_decomp)
    call addfld ('FSACI    ', 'fraction', pver, 'A', 'Relative importance of snow accreting ice'               ,phys_decomp)
    call addfld ('CME      ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap within the cloud'                      ,phys_decomp)
    call addfld ('CMEICE   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of ice within the cloud'               ,phys_decomp)
    call addfld ('CMELIQ   ', 'kg/kg/s ', pver, 'A', 'Rate of cond-evap of liq within the cloud'               ,phys_decomp)
    call addfld ('ICE2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of ice to precip'                     ,phys_decomp)
    call addfld ('LIQ2PR   ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of liq to precip'                     ,phys_decomp)
    call addfld ('ZMDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from ZM convection'               ,phys_decomp)
    call addfld ('SHDLF    ', 'kg/kg/s ', pver, 'A', 'Detrained liquid water from shallow convection'          ,phys_decomp)

    call addfld ('PRODPREC ', 'kg/kg/s ', pver, 'A', 'Rate of conversion of condensate to precip'              ,phys_decomp)
    call addfld ('EVAPPREC ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling precip'                   ,phys_decomp)
    call addfld ('EVAPSNOW ', 'kg/kg/s ', pver, 'A', 'Rate of evaporation of falling snow'                     ,phys_decomp)
    call addfld ('HPROGCLD ', 'W/kg'    , pver, 'A', 'Heating from prognostic clouds'                          ,phys_decomp)
    call addfld ('HCME     ', 'W/kg'    , pver, 'A', 'Heating from cond-evap within the cloud'                 ,phys_decomp)
    call addfld ('HEVAP    ', 'W/kg'    , pver, 'A', 'Heating from evaporation of falling precip'              ,phys_decomp)
    call addfld ('HFREEZ   ', 'W/kg'    , pver, 'A', 'Heating rate due to freezing of precip'                  ,phys_decomp)
    call addfld ('HMELT    ', 'W/kg'    , pver, 'A', 'Heating from snow melt'                                  ,phys_decomp)
    call addfld ('HREPART  ', 'W/kg'    , pver, 'A', 'Heating from cloud ice/liquid repartitioning'            ,phys_decomp)
    call addfld ('REPARTICE', 'kg/kg/s' , pver, 'A', 'Cloud ice tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('REPARTLIQ', 'kg/kg/s' , pver, 'A', 'Cloud liq tendency from cloud ice/liquid repartitioning' ,phys_decomp)
    call addfld ('FICE     ', 'fraction', pver, 'A', 'Fractional ice content within cloud'                     ,phys_decomp)
    call addfld ('ICWMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud water mixing ratio'                  ,phys_decomp)
    call addfld ('ICIMR    ', 'kg/kg   ', pver, 'A', 'Prognostic in-cloud ice mixing ratio'                    ,phys_decomp)
    call addfld ('ICWMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus water mixing ratio'                ,phys_decomp)
    call addfld ('ICIMRST  ', 'kg/kg   ', pver, 'A', 'Prognostic in-stratus ice mixing ratio'                  ,phys_decomp)
    call addfld ('PCSNOW   ', 'm/s     ', 1   , 'A', 'Snow fall from prognostic clouds'                        ,phys_decomp)
 
    call addfld ('DQSED    ', 'kg/kg/s ', pver, 'A', 'Water vapor tendency from cloud sedimentation'           ,phys_decomp)
    call addfld ('DLSED    ', 'kg/kg/s ', pver, 'A', 'Cloud liquid tendency from sedimentation'                ,phys_decomp)
    call addfld ('DISED    ', 'kg/kg/s ', pver, 'A', 'Cloud ice tendency from sedimentation'                   ,phys_decomp)
    call addfld ('HSED     ', 'W/kg    ', pver, 'A', 'Heating from cloud sediment evaporation'                 ,phys_decomp)
    call addfld ('SNOWSED  ', 'm/s     ', 1   , 'A', 'Snow from cloud ice sedimentation'                       ,phys_decomp)
    call addfld ('RAINSED  ', 'm/s     ', 1   , 'A', 'Rain from cloud liquid sedimentation'                    ,phys_decomp)
    call addfld ('PRECSED  ', 'm/s     ', 1   , 'A', 'Precipitation from cloud sedimentation'                  ,phys_decomp)

    call add_default ('FICE    ', 1, ' ')

    call addfld ('CNVCLD   ', 'fraction', 1,    'A', 'Vertically integrated convective cloud amount'           ,phys_decomp)
    call addfld ('CLDST    ', 'fraction', pver, 'A', 'Stratus cloud fraction'                                  ,phys_decomp)
    call addfld ('CONCLD   ', 'fraction', pver, 'A', 'Convective cloud cover'                                  ,phys_decomp)
	
    call add_default ('CONCLD  ', 1, ' ')

    call addfld ('AST','fraction',pver, 'A','Stratus cloud fraction',phys_decomp)

    call addfld ('LIQCLDF  ', 'fraction', pver, 'A', 'Stratus Liquid cloud fraction'                           ,phys_decomp)
    call addfld ('ICECLDF  ', 'fraction', pver, 'A', 'Stratus ICE cloud fraction'                              ,phys_decomp)
    call addfld ('IWC      ', 'kg/m3   ', pver, 'A', 'Grid box average ice water content'                      ,phys_decomp)
    call addfld ('LWC      ', 'kg/m3   ', pver, 'A', 'Grid box average liquid water content'                   ,phys_decomp)
    call addfld ('ICWNC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud water number conc'                   ,phys_decomp)
    call addfld ('ICINC    ', 'm-3     ', pver, 'A', 'Prognostic in-cloud ice number conc'                     ,phys_decomp)
    call addfld ('EFFLIQ   ', 'Micron  ', pver, 'A', 'Prognostic droplet effective radius'                     ,phys_decomp)
    call addfld ('EFFLIQ_IND','Micron  ', pver, 'A', 'Prognostic droplet effective radius (indirect effect)'   ,phys_decomp)
    call addfld ('EFFICE   ', 'Micron  ', pver, 'A', 'Prognostic ice effective radius'                         ,phys_decomp)

#ifdef DIRIND
    call addfld ('ICWMRX  ','kg/kg   ',pver, 'A','Prognostic in-cloud water mixing ratio'       ,phys_decomp)
    call add_default ('ICWMRX  ', 1, ' ')
#endif
    if ( history_budget ) then

       call add_default ('EVAPSNOW ', history_budget_histfile_num, ' ')
       call add_default ('EVAPPREC ', history_budget_histfile_num, ' ')
       call add_default ('CMELIQ   ', history_budget_histfile_num, ' ')

       if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then

          call add_default ('ZMDLF    ', history_budget_histfile_num, ' ')
          call add_default ('CME      ', history_budget_histfile_num, ' ')
          call add_default ('DQSED    ', history_budget_histfile_num, ' ')
          call add_default ('DISED    ', history_budget_histfile_num, ' ')
          call add_default ('DLSED    ', history_budget_histfile_num, ' ')
          call add_default ('HSED     ', history_budget_histfile_num, ' ')
          call add_default ('CMEICE   ', history_budget_histfile_num, ' ')
          call add_default ('LIQ2PR   ', history_budget_histfile_num, ' ')
          call add_default ('ICE2PR   ', history_budget_histfile_num, ' ')
          call add_default ('HCME     ', history_budget_histfile_num, ' ')
          call add_default ('HEVAP    ', history_budget_histfile_num, ' ')
          call add_default ('HFREEZ   ', history_budget_histfile_num, ' ')
          call add_default ('HMELT    ', history_budget_histfile_num, ' ')
          call add_default ('HREPART  ', history_budget_histfile_num, ' ')
          call add_default ('HPROGCLD ', history_budget_histfile_num, ' ')
          call add_default ('REPARTLIQ', history_budget_histfile_num, ' ')
          call add_default ('REPARTICE', history_budget_histfile_num, ' ')

       end if

    end if

    call add_default ('ICWMR', 1, ' ')
    call add_default ('ICIMR', 1, ' ')

    ! History Variables for COSP/CFMIP
    call addfld ('LS_FLXPRC', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface rain+snow flux', phys_decomp)
    call addfld ('LS_FLXSNW', 'kg/m2/s', pverp, 'A', 'ls stratiform gbm interface snow flux', phys_decomp)
    call addfld ('PRACWO', '1/s', pver, 'A', 'Accretion of cloud water by rain', phys_decomp)
    call addfld ('PSACWO', '1/s', pver, 'A', 'Accretion of cloud water by snow', phys_decomp)
    call addfld ('PSACIO', '1/s', pver, 'A', 'Accretion of cloud ice by snow', phys_decomp)

    if( microp_scheme .eq. 'RK' ) then
       call addfld ('CLDLIQSTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDLIQ'                                  ,phys_decomp)
       call addfld ('CLDICESTR   ', 'kg/kg', pver, 'A', 'Stratiform CLDICE'                                  ,phys_decomp)
       call addfld ('CLDLIQCON   ', 'kg/kg', pver, 'A', 'Convective CLDLIQ'                                  ,phys_decomp)
       call addfld ('CLDICECON   ', 'kg/kg', pver, 'A', 'Convective CLDICE'                                  ,phys_decomp)
    end if

    kvh_idx =       pbuf_get_fld_idx('kvh')
    tke_idx =       pbuf_get_fld_idx('tke')
    turbtype_idx =  pbuf_get_fld_idx('turbtype')
    smaw_idx =      pbuf_get_fld_idx('smaw')

    return
  end subroutine stratiform_init

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

  subroutine stratiform_tend(                             &
             state, ptend_all, dtime, icefrac, landfrac,  &
             ocnfrac, landm, snowh,                       &
             dlf, dlf2, rliq, cmfmc, cmfmc2, ts,          &
             sst, zdu, prec_str, snow_str, prec_sed,      &
             snow_sed, prec_pcw, snow_pcw,                & 
#ifdef DIRIND
#ifdef DIAGNCDNC
             pbuf ,state_eq,precc, relca, totcwx)
#else
             pbuf ,state_eq,precc, relca, totcwx, cdncint, fochana)
#endif
#else 
             pbuf, state_eq )
#endif


  !-------------------------------------------------------- !  
  !                                                         ! 
  ! Purpose:                                                !
  !                                                         !
  ! Interface to sedimentation, detrain, cloud fraction and !
  !        cloud macro - microphysics subroutines           !
  !                                                         ! 
  ! Author: D.B. Coleman                                    !
  ! Date: Apr 2004                                          !
  !                                                         !
  !-------------------------------------------------------- !

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use ppgrid
  use cloud_fraction,   only: cldfrc
  use physics_types,    only: physics_state, physics_ptend, physics_tend
  use physics_types,    only: physics_ptend_init, physics_update, physics_tend_init
  use physics_types,    only: physics_ptend_sum,  physics_state_copy
  use cam_history,      only: outfld
  use phys_buffer,      only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
  use constituents,     only: cnst_get_ind, pcnst
  use pkg_cld_sediment, only: cld_sediment_vel, cld_sediment_tend
  use cldwat,           only: pcond, cldwat_fice
  use cldwat2m_micro,   only: mmicro_pcond
  use microp_aero,      only: microp_aero_ts 
  use cldwat2m_macro,   only: mmacro_pcond
  use physconst,        only: cpair
  use time_manager,     only: is_first_step, get_nstep
  use pkg_cldoptics,    only: cldefr
  use phys_control,     only: cam_physpkg_is
  ! Debug
    use phys_debug_util,  only: phys_debug_col
  ! Debug

  implicit none

  ! Debug
    integer icol
  ! Debug

!
! Parameters
!
  real(r8) pnot                  ! Reference pressure
  parameter (pnot = 1.e5_r8)

  !
  ! Input arguments
  !

  type(physics_state), intent(in)    :: state       ! State variables
  type(physics_ptend), intent(out)   :: ptend_all   ! Package tendencies
  type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

  real(r8), intent(in)  :: dtime                    ! Timestep
  real(r8), intent(in)  :: icefrac (pcols)          ! Sea ice fraction (fraction)
  real(r8), intent(in)  :: landfrac(pcols)          ! Land fraction (fraction)
  real(r8), intent(in)  :: ocnfrac (pcols)          ! Ocean fraction (fraction)
  real(r8), intent(in)  :: landm(pcols)             ! Land fraction ramped over water
  real(r8), intent(in)  :: snowh(pcols)             ! Snow depth over land, water equivalent (m)

  real(r8), intent(in)  :: dlf(pcols,pver)          ! Detrained water from convection schemes
  real(r8), intent(in)  :: dlf2(pcols,pver)         ! Detrained water from shallow convection scheme
  real(r8), intent(in)  :: rliq(pcols)              ! Vertical integral of liquid not yet in q(ixcldliq)
  real(r8), intent(in)  :: cmfmc(pcols,pverp)       ! Deep + Shallow Convective mass flux [ kg /s/m^2 ]
  real(r8), intent(in)  :: cmfmc2(pcols,pverp)      ! Shallow convective mass flux [ kg/s/m^2 ]

  real(r8), intent(in)  :: ts(pcols)                ! Surface temperature
  real(r8), intent(in)  :: sst(pcols)               ! Sea surface temperature
  real(r8), intent(in)  :: zdu(pcols,pver)          ! Detrainment rate from deep convection
#ifdef DIRIND
  real(r8), intent(in) :: precc(pcols)                ! convective-scale precipitation rate
#endif
  real(r8), intent(out) :: prec_str(pcols)          ! [Total] Sfc flux of precip from stratiform [ m/s ] 
  real(r8), intent(out) :: snow_str(pcols)          ! [Total] Sfc flux of snow from stratiform   [ m/s ]
  real(r8), intent(out) :: prec_sed(pcols)          ! Surface flux of total cloud water from sedimentation
  real(r8), intent(out) :: snow_sed(pcols)          ! Surface flux of cloud ice from sedimentation
  real(r8), intent(out) :: prec_pcw(pcols)          ! Sfc flux of precip from microphysics [ m/s ]
  real(r8), intent(out) :: snow_pcw(pcols)          ! Sfc flux of snow from microphysics [ m/s ]

  ! Equilibrium state variables at the end of macrophysics
  ! Below 'state_eq' is for future use as the input of radiation'PBL scheme

  type(physics_state), intent(out) :: state_eq   ! Equilibrium state variables at the end of macrophysics

  !
  ! Local variables
  !

  type(physics_state)   :: state1                   ! Local copy of the state variable
  type(physics_tend )   :: tend                     ! Physics tendencies (empty, needed for physics_update call)
  type(physics_ptend)   :: ptend_loc                ! Package tendencies

  integer i,k,m
  integer :: lchnk                                  ! Chunk identifier
  integer :: ncol                                   ! Number of atmospheric columns

  ! Physics buffer fields

  integer itim, ifld
  real(r8), pointer, dimension(:,:) :: rhdfda       !
  real(r8), pointer, dimension(:,:) :: rhu00        ! 
  real(r8), pointer, dimension(:,:) :: qcwat        ! Cloud water old q
  real(r8), pointer, dimension(:,:) :: tcwat        ! Cloud water old temperature
  real(r8), pointer, dimension(:,:) :: lcwat        ! Cloud liquid water old q
#ifdef DIRIND
  real(r8), pointer, dimension(:,:) :: lcwatx       ! Offline cloud liquid water old q
#endif
  real(r8), pointer, dimension(:,:) :: iccwat       ! Cloud ice water old q
  real(r8), pointer, dimension(:,:) :: nlwat        ! Cloud liquid droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: niwat        ! Cloud ice    droplet number condentration. old.
  real(r8), pointer, dimension(:,:) :: CC_T         ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qv        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ql        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qi        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_nl        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_ni        ! Grid-mean microphysical tendency
  real(r8), pointer, dimension(:,:) :: CC_qlst      ! In-liquid stratus microphysical tendency
  real(r8), pointer, dimension(:,:) :: cld          ! Total cloud fraction
  real(r8), pointer, dimension(:,:) :: ast          ! Relative humidity cloud fraction
  real(r8), pointer, dimension(:,:) :: aist         ! Physical ice stratus fraction
  real(r8), pointer, dimension(:,:) :: alst         ! Physical liquid stratus fraction
  real(r8), pointer, dimension(:,:) :: qist         ! Physical in-cloud IWC
  real(r8), pointer, dimension(:,:) :: qlst         ! Physical in-cloud LWC
  real(r8), pointer, dimension(:,:) :: concld       ! Convective cloud fraction
  real(r8), pointer, dimension(:,:) :: qme
  real(r8), pointer, dimension(:,:) :: prain        ! Total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: nevapr       ! Evaporation of total precipitation (rain + snow)
  real(r8), pointer, dimension(:,:) :: rel          ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei          ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: rel2         ! Liquid effective drop radius (microns)
  real(r8), pointer, dimension(:,:) :: rei2         ! Ice effective drop size (microns)
  real(r8), pointer, dimension(:,:) :: cldo         ! Old cloud fraction
  real(r8), pointer, dimension(:,:) :: kkvh         ! Vertical eddy diffusivity
  real(r8), pointer, dimension(:,:) :: wsedl        ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
  real(r8) :: shfrc(pcols,pver)                     ! Cloud fraction from shallow convection scheme
  real(r8), pointer, dimension(:,:) :: rel_fn       ! Ice effective drop size at fixed number (indirect effect) (microns)

  real(r8)  rate1cld(pcols,pver)                    ! array to hold rate1ord_cw2pr_st from microphysics

  ! physics buffer fields for COSP simulator (RK only)

  real(r8), pointer, dimension(:,:) :: rkflxprc     ! RK grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
  real(r8), pointer, dimension(:,:) :: rkflxsnw     ! RK grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)

  ! physics buffer fields for radiation

  real(r8), pointer, dimension(:,:) :: dei          ! Ice effective diameter (meters) (AG: microns?)
  real(r8), pointer, dimension(:,:) :: mu           ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdac      ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwp        ! In-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwp        ! In-cloud liquid water path for radiation
  
  ! For rrtm optics. specificed distribution.

  real(r8) :: mucon                                 ! Convective size distribution shape parameter
  real(r8) :: dcon                                  ! Convective size distribution effective radius (meters)
  real(r8) :: lamcon                                ! Convective size distribution slope parameter (meters-1)
  real(r8) :: deicon                                ! Convective ice effective diameter (meters)

  ! Physics buffer fields

  real(r8), pointer, dimension(:,:) :: deiconv      ! Ice effective diameter (microns)
  real(r8), pointer, dimension(:,:) :: muconv       ! Size distribution shape parameter for radiation
  real(r8), pointer, dimension(:,:) :: lambdaconv   ! Size distribution slope parameter for radiation
  real(r8), pointer, dimension(:,:) :: iciwpst      ! Stratiform in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpst      ! Stratiform in-cloud liquid water path for radiation
  real(r8), pointer, dimension(:,:) :: iciwpconv    ! Convective in-cloud ice water path for radiation
  real(r8), pointer, dimension(:,:) :: iclwpconv    ! Convective in-cloud liquid water path for radiation

  real(r8), pointer, dimension(:,:) :: tke          ! TKE from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: turbtype     ! Turbulence type from the moist PBL scheme
  real(r8), pointer, dimension(:,:) :: smaw         ! Instability function of momentum from the moist PBL scheme

  ! Convective cloud to the physics buffer for purposes of ql contrib. to radn.

  real(r8), pointer, dimension(:,:) :: concld_ql    ! Convective cloud
  real(r8), pointer, dimension(:,:) :: fice_ql      ! Cloud ice/water partitioning ratio.

  ! Local variables for in-cloud water quantities adjusted for convective water

  real(r8)  allcld_ice (pcols,pver)                 ! All-cloud cloud ice
  real(r8)  allcld_liq (pcols,pver)                 ! All-cloud liquid

  ! Snow

  real(r8), pointer, dimension(:,:) :: cldfsnow     ! Cloud fraction for liquid+snow
  real(r8), pointer, dimension(:,:) :: icswp        ! In-cloud snow water path
  real(r8), pointer, dimension(:,:) :: des          ! Snow effective diameter (m)
  real(r8)  qsout(pcols,pver)                       ! Snow mixing ratio
 
  ! Local variables for stratiform_sediment

  real(r8)  rain(pcols)                             ! Surface flux of cloud liquid
  real(r8)  pvliq(pcols,pver+1)                     ! Vertical velocity of cloud liquid drops (Pa/s)
  real(r8)  pvice(pcols,pver+1)                     ! Vertical velocity of cloud ice particles (Pa/s)

  ! Local variables for cldfrc

  real(r8)  cldst(pcols,pver)                       ! Stratus cloud fraction
  real(r8)  rhcloud(pcols,pver)                     ! Relative humidity cloud (last timestep)
  real(r8)  rhcloud2(pcols,pver)                    ! Relative humidity cloud (perturbation)
  real(r8)  clc(pcols)                              ! Column convective cloud amount
  real(r8)  relhum(pcols,pver)                      ! RH, output to determine drh/da
  real(r8)  rhu002(pcols,pver)                      ! Same as rhu00 but for perturbed rh 
  real(r8)  cld2(pcols,pver)                        ! Same as cld but for perturbed rh
  real(r8)  concld2(pcols,pver)                     ! Same as concld but for perturbed rh 
  real(r8)  cldst2(pcols,pver)                      ! Same as cldst but for perturbed rh 
  real(r8)  relhum2(pcols,pver)                     ! RH after  perturbation            
  real(r8)  icecldf(pcols,pver)                     ! Ice cloud fraction
  real(r8)  liqcldf(pcols,pver)                     ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf_out(pcols,pver)                 ! Ice cloud fraction
  real(r8)  liqcldf_out(pcols,pver)                 ! Liquid cloud fraction (combined into cloud)
  real(r8)  icecldf2(pcols,pver)                    ! Ice cloud fraction
  real(r8)  liqcldf2(pcols,pver)                    ! Liquid cloud fraction (combined into cloud)

  ! Local variables for microphysics

  real(r8)  rdtime                                  ! 1./dtime
  real(r8)  qtend(pcols,pver)                       ! Moisture tendencies
  real(r8)  ttend(pcols,pver)                       ! Temperature tendencies
  real(r8)  ltend(pcols,pver)                       ! Cloud liquid water tendencies
#ifdef DIRIND
  real(r8) :: ltendx(pcols,pver)                    ! Offline cloud liquid water tendencies
#endif
  real(r8)  evapheat(pcols,pver)                    ! Heating rate due to evaporation of precip
  real(r8)  evapsnow(pcols,pver)                    ! Local evaporation of snow
  real(r8)  prfzheat(pcols,pver)                    ! Heating rate due to freezing of precip (W/kg)
  real(r8)  meltheat(pcols,pver)                    ! Heating rate due to phase change of precip
  real(r8)  cmeheat (pcols,pver)                    ! Heating rate due to phase change of precip
  real(r8)  prodsnow(pcols,pver)                    ! Local production of snow
  real(r8)  totcw(pcols,pver)                       ! Total cloud water mixing ratio
  real(r8)  fice(pcols,pver)                        ! Fractional ice content within cloud
  real(r8)  fsnow(pcols,pver)                       ! Fractional snow production
  real(r8)  repartht(pcols,pver)                    ! Heating rate due to phase repartition of input precip
  real(r8)  icimr(pcols,pver)                       ! In cloud ice mixing ratio
  real(r8)  icwmr(pcols,pver)                       ! In cloud water mixing ratio
  real(r8)  icimrst(pcols,pver)                     ! In stratus ice mixing ratio
  real(r8)  icwmrst(pcols,pver)                     ! In stratus water mixing ratio
  real(r8)  icimrst_out(pcols,pver)                 ! In stratus ice mixing ratio
  real(r8)  icwmrst_out(pcols,pver)                 ! In stratus water mixing ratio
  real(r8)  fwaut(pcols,pver)              
  real(r8)  fsaut(pcols,pver)              
  real(r8)  fracw(pcols,pver)              
  real(r8)  fsacw(pcols,pver)              
  real(r8)  fsaci(pcols,pver)              
  real(r8)  cmeice(pcols,pver)                      ! Rate of cond-evap of ice within the cloud
  real(r8)  cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud
  real(r8)  ice2pr(pcols,pver)                      ! Rate of conversion of ice to precip
  real(r8)  liq2pr(pcols,pver)                      ! Rate of conversion of liquid to precip
  real(r8)  liq2snow(pcols,pver)                    ! Rate of conversion of liquid to snow
  real(r8)  temp(pcols)
  real(r8)  res(pcols,pver)
  real(r8)  droprad                                 ! Radius of droplets detrained from cumulus (m)
  real(r8)  invdropmass                             ! Inverse of mean droplet mass (#/kg)

  ! Local variables for CFMIP calculations
  real(r8) :: mr_lsliq(pcols,pver)                     ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
  real(r8) :: mr_lsice(pcols,pver)                     ! mixing_ratio_large_scale_cloud_ice (kg/kg)
  real(r8) :: mr_ccliq(pcols,pver)                     ! mixing_ratio_convective_cloud_liquid (kg/kg)
  real(r8) :: mr_ccice(pcols,pver)                     ! mixing_ratio_convective_cloud_ice (kg/kg)

  real(r8) :: pracwo(pcols,pver)                       ! RK accretion of cloud water by rain (1/s)
  real(r8) :: psacwo(pcols,pver)                       ! RK accretion of cloud water by snow (1/s)
  real(r8) :: psacio(pcols,pver)                       ! RK accretion of cloud ice by snow (1/s)

  ! Variables for MG microphysics

  real(r8)  dum1,dum2
  real(r8)  qc(pcols,pver)
  real(r8)  qi(pcols,pver)
  real(r8)  nc(pcols,pver)
  real(r8)  ni(pcols,pver)
  real(r8)  icinc(pcols,pver)                       ! In cloud ice number conc
  real(r8)  cdnumc(pcols)                           ! Vertically-integrated droplet concentration
  real(r8)  icwnc(pcols,pver)                       ! In cloud water number conc
  real(r8)  iwc(pcols,pver)                         ! Grid box average ice water content
  real(r8)  lwc(pcols,pver)                         ! Grid box average liquid water content  
  real(r8)  effliq(pcols,pver)                      ! In cloud liq eff rad
  real(r8)  effice(pcols,pver)                      ! In cloud ice eff rad
  real(r8)  effliq_fn(pcols,pver)                   ! In cloud liq eff rad at fixed number concentration	
  real(r8)  wsub(pcols,pver)                        ! Sub-grid vertical velocity (m/s)
  real(r8)  wsubi(pcols,pver)                       ! Sub-grid vertical velocity for ice (m/s)

  ! Output from mmicro_pcond

  real(r8)  tlat(pcols,pver)
  real(r8)  qvlat(pcols,pver)
  real(r8)  qcten(pcols,pver)
  real(r8)  qiten(pcols,pver)
  real(r8)  ncten(pcols,pver)
  real(r8)  niten(pcols,pver)
  real(r8)  effc(pcols,pver)
  real(r8)  effc_fn(pcols,pver)                     ! Liquid effective radius at fixed number (for indirect calc)
  real(r8)  effi(pcols,pver)
  real(r8)  prect(pcols)
  real(r8)  preci(pcols)

  ! Output from microp_aero_ts for aerosol actication
  real(r8)  naai(pcols,pver)      !ice nucleation number
  real(r8)  npccn(pcols,pver)     !liquid activation number
  real(r8)  rndst(pcols,pver,4)
  real(r8)  nacon(pcols,pver,4)

  ! Output from mmacro_pcond

  real(r8)  qvadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
  real(r8)  qladj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
  real(r8)  qiadj(pcols,pver)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
  real(r8)  qllim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
  real(r8)  qilim(pcols,pver)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

  ! Averaging arrays for effective radius and number....

  real(r8)  efiout(pcols,pver)
  real(r8)  efcout(pcols,pver)
  real(r8)  ncout(pcols,pver)
  real(r8)  niout(pcols,pver)

  real(r8)  freqi(pcols,pver)
  real(r8)  freql(pcols,pver)

  ! Average cloud top radius & number

  real(r8)  ctrel(pcols)
  real(r8)  ctrei(pcols)
  real(r8)  ctnl(pcols)
  real(r8)  ctni(pcols)
  real(r8)  fcti(pcols)
  real(r8)  fctl(pcols)

  ! Gather mass mixing ratio for all aerosols affecting the climate

  integer               :: naer_all
  real(r8), pointer     :: aermmr1(:,:)
  real(r8), allocatable :: aer_mmr(:,:,:)           ! Aerosol mass mixing ratio


#ifdef DIRIND  ! Local variables used for calculating aerosol optics and direct and indirect forcings
  real(r8), intent(out) :: totcwx(pcols,pver)! liquid cloud water mixing ratio
  real(r8) :: totcwl(pcols,pver)         ! liquid cloud water mixing ratio
  real(r8) :: totcwlx(pcols,pver)        ! liquid cloud water mixing ratio
  real(r8) :: cmex(pcols,pver)           ! local condensation - evaporation of cloud water
#ifndef DIAGNCDNC
  real(r8) :: relh(pcols,pver)
  real(r8) :: totnlx(pcols,pver)         ! total cloud droplet concentration
#endif
  real(r8), intent(out) :: relca(pcols,pver)          ! computed effective radius of liquid clouds
  real(r8) :: relcao(pcols,pver)         ! computed effective radius of liquid clouds
  real(r8) :: cdnc(pcols,pver)       	 ! cloud droplet number concentration
  real(r8) :: cdnco(pcols,pver)       	 ! cloud droplet number concentration

  real(r8) :: cldfoc(pcols,pver)         ! cloud droplet freq. of occurence
  real(r8) :: rehanx(pcols)              ! effective radius of "Han cloud droplets"
  real(r8) :: rehana(pcols)              ! rehanx*fochana
#ifdef DIAGNCDNC
  real(r8) :: fochana(pcols)             ! frequency of occurence of Han droplets
#else
  real(r8),intent(out) :: fochana(pcols)             ! frequency of occurence of Han droplets
#endif
  real(r8) :: wpatha(pcols,pver)              ! liquid water path for "Han cloud droplets"
  real(r8) :: pdelm1(pcols,pver)              ! delta_p "Han cloud droplets"
#ifndef DIAGNCDNC
  real(r8) :: ccdnc(pcols,pver)          ! column cdnc per level
  real(r8),intent(out) :: cdncint(pcols)             ! column cdnc
#endif
  real(r8) :: ntotal(pcols,pver)         ! total aerosol number concentration (3D)
  real(r8) :: cxsind(pcols,pver)         ! total excess mass concentration wrt. internal mixing
  real(r8) :: n1(pcols,pver)
  real(r8) :: n2(pcols,pver)
  real(r8) :: n3(pcols,pver)
  real(r8) :: n4(pcols,pver)
  real(r8) :: n5(pcols,pver)
  real(r8) :: n6(pcols,pver)
  real(r8) :: n7(pcols,pver)
  real(r8) :: n8(pcols,pver)
  real(r8) :: n9(pcols,pver)
  real(r8) :: n10(pcols,pver)
  real(r8) :: n11(pcols,pver)
  real(r8) :: n12(pcols,pver)
  real(r8) :: n13(pcols,pver)
  real(r8) :: n14(pcols,pver)
  real(r8) :: nrmodes(pcols,pver,pcnst) ! number concentration in each mode 
  real(r8) :: icwmrx(pcols,pver)        ! in cloud water mixing ratio
#ifndef DIAGNCDNC
  real(r8) :: sigw1(pcols,pver)         ! standard deviation of the vertical velocity
  real(r8) :: deltaz(pcols,pver)        ! level thickness (in km)
  real(r8) :: selfx(pcols,pver)         ! selfcollection rate
  real(r8) :: selfxo(pcols,pver)        ! selfcollection rate
  real(r8) :: nrainx(pcols,pver)        ! loss of cloud droplets due to precip
  real(r8) :: nrainxo(pcols,pver)       ! loss of cloud droplets due to precip
  real(r8) :: supersat(pcols,pver)      ! supersaturation
  real(r8) :: moistrho(pcols,pver)      ! 
  real(r8) :: liq2prx(pcols,pver)       ! rate of conversion of liquid to precip
  real(r8) :: nucrat(pcols,pver)
  real(r8) :: nucrato(pcols,pver)
  real(r8) :: lossx(pcols,pver)
  real(r8) :: evapx(pcols,pver)         ! Evaporation of cloud droplets (cm-3/s)
  real(r8) :: evapxo(pcols,pver)        ! Evaporation of cloud droplets (cm-3/s)
  real(r8) :: freez(pcols,pver)         ! rate of freezing of cloud droplets (simple) (cm-3/s)
  real(r8) :: freezo(pcols,pver)        ! rate of freezing of cloud droplets (simple) (cm-3/s)
#endif
#endif

  ! For revised macophysics, mmacro_pcond

  real(r8)  itend(pcols,pver)
  real(r8)  lmitend(pcols,pver)
  real(r8)  zeros(pcols,pver)
  real(r8)  t_inout(pcols,pver)
  real(r8)  qv_inout(pcols,pver)
  real(r8)  ql_inout(pcols,pver)
  real(r8)  qi_inout(pcols,pver)
  real(r8)  prsed(pcols,pver)
  real(r8)  pssed(pcols,pver)
  real(r8)  ersed(pcols,pver)
  real(r8)  essed(pcols,pver)
  real(r8)  alst_mic(pcols,pver)
  real(r8)  aist_mic(pcols,pver)
  real(r8)  concld_old(pcols,pver)

  real(r8)  nl_inout(pcols,pver)
  real(r8)  ni_inout(pcols,pver)
  real(r8)  dum1D(pcols)
  real(r8)  nltend(pcols,pver)
  real(r8)  nitend(pcols,pver)

  real(r8)  zero1D(pcols)
  real(r8)  t_out(pcols,pver)
  real(r8)  qv_out(pcols,pver)
  real(r8)  ql_out(pcols,pver)
  real(r8)  qi_out(pcols,pver)
  real(r8)  nl_out(pcols,pver)
  real(r8)  ni_out(pcols,pver)
  real(r8)  QQw(pcols,pver)
  real(r8)  QQi(pcols,pver)
  real(r8)  QQnl(pcols,pver)
  real(r8)  QQni(pcols,pver)

  ! For detraining cumulus condensate into the 'stratus' without evaporation
  ! This is for use in mmacro_pcond

  real(r8)  dlf_T(pcols,pver)
  real(r8)  dlf_qv(pcols,pver)
  real(r8)  dlf_ql(pcols,pver)
  real(r8)  dlf_qi(pcols,pver)
  real(r8)  dlf_nl(pcols,pver)
  real(r8)  dlf_ni(pcols,pver)

  real(r8)  rel_detcu
  real(r8)  rei_detcu
  ! ======================================================================

  lchnk = state%lchnk
  ncol  = state%ncol

  call physics_state_copy(state,state1)             ! Copy state to local state1.
  call physics_ptend_init(ptend_loc)                ! Initialize local ptend type
  call physics_ptend_init(ptend_all)                ! Initialize output ptend type
  call physics_tend_init(tend)                      ! tend here is just a null place holder

  ! Associate pointers with physics buffer fields

  itim = pbuf_old_tim_idx()

  qcwat => pbuf(qcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rhdfda => pbuf(rhdfda_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rhu00 => pbuf(rhu00_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  tcwat => pbuf(tcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  lcwat => pbuf(lcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

#ifdef DIRIND
  lcwatx => pbuf(lcwatx_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
#endif

  iccwat => pbuf(iccwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  nlwat => pbuf(nlwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  niwat => pbuf(niwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_T => pbuf(cc_t_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qv => pbuf(cc_qv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_ql => pbuf(cc_ql_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qi => pbuf(cc_qi_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_nl => pbuf(cc_nl_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_ni => pbuf(cc_ni_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  CC_qlst => pbuf(cc_qlst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  cld => pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  ast => pbuf(ast_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  aist => pbuf(aist_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  alst => pbuf(alst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  qist => pbuf(qist_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  qlst => pbuf(qlst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  concld => pbuf(concld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  cldo => pbuf(cldo_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rel2 => pbuf(rel2_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

  rei2 => pbuf(rei2_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

! For purposes of convective ql.

  concld_ql => pbuf(concldql_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)                                      

  fice_ql => pbuf(fice_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

  qme  => pbuf(qme_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  prain  => pbuf(prain_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  nevapr  => pbuf(nevapr_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rel  => pbuf(rel_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rei  => pbuf(rei_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  rel_fn  => pbuf(rel_fn_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  kkvh => pbuf(kvh_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  dei  => pbuf(dei_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  mu  => pbuf(mu_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  lambdac  => pbuf(lambdac_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwp  => pbuf(iciwp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwp  => pbuf(iclwp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  deiconv  => pbuf(deiconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  muconv  => pbuf(muconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  lambdaconv  => pbuf(lambdaconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwpst  => pbuf(iciwpst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwpst  => pbuf(iclwpst_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iciwpconv  => pbuf(iciwpconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  iclwpconv  => pbuf(iclwpconv_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  des  => pbuf(des_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  icswp  => pbuf(icswp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  cldfsnow  => pbuf(cldfsnow_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, itim)

  tke => pbuf(tke_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  turbtype => pbuf(turbtype_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  smaw => pbuf(smaw_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk, 1)

  wsedl  => pbuf(wsedl_idx)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

  ! If first timestep, initialize heatflux....in pbuf at all time levels.

  if( is_first_step() ) then
      kkvh(:,:)     = 0._r8
      tke(:,:)      = 0._r8
      turbtype(:,:) = 0._r8
      smaw(:,:)     = 0._r8
  endif

  ! Assign default size distribution parameters for no-stratiform clouds (convection only)
  ! Also put into physics buffer for possible separate use by radiation

  dcon   = 25.e-6_r8
  mucon  = 5.3_r8
  deicon = 50._r8

  muconv(:,:)     = mucon
  lambdaconv(:,:) = (mucon + 1._r8)/dcon
  deiconv(:,:)    = deicon

  ! Initialize convective detrainment tendency

  dlf_T(:,:)  = 0._r8
  dlf_qv(:,:) = 0._r8
  dlf_ql(:,:) = 0._r8
  dlf_qi(:,:) = 0._r8
  dlf_nl(:,:) = 0._r8
  dlf_ni(:,:) = 0._r8

   ! ------------------------------------- !
   ! From here, process computation begins ! 
   ! ------------------------------------- !

   ! ------------- !
   ! Sedimentation !
   ! ------------- !

   if( microp_scheme .eq. 'RK' ) then

     ! Allow the cloud liquid drops and ice particles to sediment.
     ! This is done before adding convectively detrained cloud water, 
     ! because the phase of the detrained water is unknown.

       call t_startf('stratiform_sediment')

       ptend_loc%name         = 'pcwsediment'
       ptend_loc%ls           = .TRUE.
       ptend_loc%lq(1)        = .TRUE.
       ptend_loc%lq(ixcldice) = .TRUE.
       ptend_loc%lq(ixcldliq) = .TRUE.

#ifdef DIRIND
       ptend_loc%lq(icldlx) = .TRUE. 
       ptend_loc%lq(icldnx) = .FALSE. !not implemented (yet)
#endif

       call cld_sediment_vel( ncol,                                                           &
                              icefrac, landfrac, ocnfrac, state1%pmid, state1%pdel, state1%t, &
                              cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice),            & 
                              pvliq, pvice, landm, snowh )

       wsedl(:ncol,:pver) = pvliq(:ncol,:pver)/gravit/(state1%pmid(:ncol,:pver)/(287.15_r8*state1%t(:ncol,:pver)))
#ifdef DIRIND
       call cld_sediment_tend( ncol, dtime ,                                                             &
                               state1%pint, state1%pmid, state1%pdel, state1%t,                          &
                               cld, state1%q(:,:,ixcldliq), state1%q(:,:,icldlx), state1%q(:,:,ixcldice), pvliq, pvice,        &
                               ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,icldlx),ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,1), &
                               ptend_loc%s, rain, snow_sed )
#else
       call cld_sediment_tend( ncol, dtime ,                                                             &
                               state1%pint, state1%pmid, state1%pdel, state1%t,                          &
                               cld, state1%q(:,:,ixcldliq), state1%q(:,:,ixcldice), pvliq, pvice,        &
                               ptend_loc%q(:,:,ixcldliq), ptend_loc%q(:,:,ixcldice), ptend_loc%q(:,:,1), &
                               ptend_loc%s, rain, snow_sed )
#endif

     ! Convert rain and snow fluxes at the surface from [kg/m2/s] to [m/s]
     ! Compute total precipitation flux at the surface in [m/s]

       snow_sed(:ncol) = snow_sed(:ncol)/1000._r8
       rain(:ncol)     = rain(:ncol)/1000._r8
       prec_sed(:ncol) = rain(:ncol) + snow_sed(:ncol)

     ! Record history variables
       lchnk = state1%lchnk
       call outfld( 'DQSED'   ,ptend_loc%q(:,:,1)       , pcols,lchnk )
       call outfld( 'DISED'   ,ptend_loc%q(:,:,ixcldice), pcols,lchnk )
       call outfld( 'DLSED'   ,ptend_loc%q(:,:,ixcldliq), pcols,lchnk )
       call outfld( 'HSED'    ,ptend_loc%s              , pcols,lchnk )
       call outfld( 'PRECSED' ,prec_sed                 , pcols,lchnk )
       call outfld( 'SNOWSED' ,snow_sed                 , pcols,lchnk )
       call outfld( 'RAINSED' ,rain                     , pcols,lchnk )

     ! Add tendency from this process to tend from other processes here
       call physics_ptend_sum( ptend_loc, ptend_all, state )

     ! Update physics state type state1 with ptend_loc 
       call physics_update( state1, tend, ptend_loc, dtime )
       call physics_ptend_init( ptend_loc )

       call t_stopf('stratiform_sediment')

     ! Accumulate prec and snow flux at the surface [ m/s ]
       prec_str(:ncol) = prec_sed(:ncol)
       snow_str(:ncol) = snow_sed(:ncol)

   endif  ! End of 'Sediment'

   ! ----------------------------------------------------------------------------- !
   ! Detrainment of convective condensate into the environment or stratiform cloud !
   ! ----------------------------------------------------------------------------- !

   ptend_loc%name = 'pcwdetrain'

   if ( microp_scheme .eq. 'RK' ) then

     ! Put all of the detraining cloud water from convection into the large scale cloud.
     ! It all goes in liquid for the moment.
     ! Strictly speaking, this approach is detraining all the cconvective water into 
     ! the environment, not the large-scale cloud.

       ptend_loc%lq(ixcldliq) = .TRUE.
#ifdef DIRIND
       ptend_loc%lq(icldlx) = .TRUE.
#endif
       do k = 1, pver
       do i = 1, state1%ncol
          ptend_loc%q(i,k,ixcldliq) = dlf(i,k)
#ifdef DIRIND
         ptend_loc%q(i,k,icldlx) = dlf(i,k)
#endif
       end do
       end do

   end if

   call outfld( 'ZMDLF', dlf, pcols, state1%lchnk )
   call outfld( 'SHDLF', dlf2, pcols, state1%lchnk )

 ! Add hie detrainment tendency to tend from the other prior processes

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc )

 ! Accumulate prec and snow, reserved liquid has now been used.

   if( microp_scheme .eq. 'RK' ) then
       prec_str(:ncol) = prec_str(:ncol) - rliq(:ncol)  ! ( snow contribution is zero )
   endif

   ! -------------------------------------- !
   ! Computation of Various Cloud Fractions !
   ! -------------------------------------- !

   ! ----------------------------------------------------------------------------- !
   ! Treatment of cloud fraction in CAM4 and CAM5 differs                          !  
   ! (1) CAM4                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( empirical fcn of mass flux )        !
   !     . Stratus AMT = max( RH stratus AMT, Stability Stratus AMT )              !
   !     . Cumulus and Stratus are 'minimally' overlapped without hierarchy.       !
   !     . Cumulus LWC,IWC is assumed to be the same as Stratus LWC,IWC            !
   ! (2) CAM5                                                                      !
   !     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      !
   !                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) !
   !     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     !
   !     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  !
   !     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              !
   ! ----------------------------------------------------------------------------- ! 

   concld_old(:ncol,:pver) = concld(:ncol,:pver)

   if( use_shfrc ) then
       shfrc(:pcols,:pver) = pbuf(shfrc_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   else 
       shfrc(:,:) = 0._r8
   endif

#ifdef DIRIND
#ifndef DIAGNCDNC
   !save old cloud fraction
   cldo(:,:)=cld(:,:)
#endif
#endif

   ! CAM5 only uses 'concld' output from the below subroutine. 
   ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
   ! will be computed using this updated 'concld' in the stratiform macrophysics 
   ! scheme (mmacro_pcond) later below. 
   ! Note 'shfrc' and' deep convective cloud fraction' will be saved into the 
   ! physical buffer (SH_FRAC,DP_FRAC) within cldfrc.

   call t_startf("cldfrc")
   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld, rhcloud, clc, state1%pdel,                                    &
                cmfmc, cmfmc2, landfrac,snowh, concld, cldst,                      &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                relhum, 0 )    

   ! Re-calculate cloud with perturbed rh add call cldfrc  

   call cldfrc( lchnk, ncol, pbuf,                                                 &
                state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                shfrc, use_shfrc,                                                  &
                cld2, rhcloud2, clc, state1%pdel,                                  &
                cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                relhum2, 1 )              

   call t_stopf("cldfrc")

   ! Add following to estimate rhdfda. Below block is only for CAM4                       

   rhu00(:ncol,1) = 2.0_r8

   do k = 1, pver
   do i = 1, ncol
      if( relhum(i,k) < rhu00(i,k) ) then
          rhdfda(i,k) = 0.0_r8
      elseif( relhum(i,k) >= 1.0_r8 ) then
          rhdfda(i,k) = 0.0_r8
      else
         ! Under certain circumstances, rh+ cause cld not to changed
         ! when at an upper limit, or w/ strong subsidence
         if( ( cld2(i,k) - cld(i,k) ) < 1.e-4_r8 ) then
               rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8  
         else
               rhdfda(i,k) = 0.01_r8*relhum(i,k)/(cld2(i,k)-cld(i,k))
         endif
      endif
   enddo
   enddo

   ! ---------------------------------------------- !
   ! Stratiform Cloud Macrophysics and Microphysics !
   ! ---------------------------------------------- !

   call t_startf('stratiform_microphys')

   lchnk  = state1%lchnk
   ncol   = state1%ncol
   rdtime = 1._r8/dtime

 ! Define fractional amount of stratus condensate and precipitation in ice phase.
 ! This uses a ramp ( -30 ~ -10 for fice, -5 ~ 0 for fsnow ). 
 ! The ramp within convective cloud may be different

   call cldwat_fice( ncol, state1%t, fice, fsnow )

   if( microp_scheme .eq. 'RK' ) then

     ! Perform repartitioning of stratiform condensate.    
     ! Corresponding heating tendency will be added later. 

       totcw(:ncol,:pver)     = state1%q(:ncol,:pver,ixcldice) + state1%q(:ncol,:pver,ixcldliq)

#ifdef DIRIND
#ifdef DIAGNCDNC
!temporary (online) solution:  (without the additional 'prognostic' cloud water) 
       totcwx(:ncol,:pver) = totcw(:ncol,:pver)
#else
       totcwx(:ncol,:pver) = state1%q(:ncol,:pver,ixcldice) + state1%q(:ncol,:pver,icldlx)
#endif
       totcwl(:ncol,:pver) = state1%q(:ncol,:pver,ixcldliq)
       totcwlx(:ncol,:pver) = state1%q(:ncol,:pver,icldlx)
#ifndef DIAGNCDNC
       totnlx(:ncol,:pver) = state1%q(:ncol,:pver,icldnx)
#endif
#endif

       repartht(:ncol,:pver)  = state1%q(:ncol,:pver,ixcldice)
       ptend_loc%q(:ncol,:pver,ixcldice) = rdtime * ( totcw(:ncol,:pver)*fice(:ncol,:pver)          - state1%q(:ncol,:pver,ixcldice) )
       ptend_loc%q(:ncol,:pver,ixcldliq) = rdtime * ( totcw(:ncol,:pver)*(1.0_r8-fice(:ncol,:pver)) - state1%q(:ncol,:pver,ixcldliq) )
#ifdef DIRIND
   ptend_loc%q(:ncol,:pver,icldlx) = &
         ( totcwx(:ncol,:pver) * (1.0_r8 - fice(:ncol,:pver)) - state1%q(:ncol,:pver,icldlx) ) * rdtime
#endif

       call outfld( 'REPARTICE', ptend_loc%q(:,:,ixcldice), pcols, lchnk )
       call outfld( 'REPARTLIQ', ptend_loc%q(:,:,ixcldliq), pcols, lchnk )

       ptend_loc%name         = 'cldwat-repartition'
       ptend_loc%lq(ixcldice) = .true.
       ptend_loc%lq(ixcldliq) = .true.
#ifdef DIRIND
       ptend_loc%lq(icldlx) = .true.
#endif

       call physics_ptend_sum( ptend_loc, ptend_all, state )
       call physics_update( state1, tend, ptend_loc, dtime )
       call physics_ptend_init( ptend_loc )

   endif

   ptend_loc%name         = 'cldwat'
   ptend_loc%ls           = .true.
   ptend_loc%lq(1)        = .true.
   ptend_loc%lq(ixcldice) = .true.
   ptend_loc%lq(ixcldliq) = .true.
#ifdef DIRIND
   ptend_loc%lq(icldlx) = .true.
   ptend_loc%lq(icldnx) = .true.
#endif

   if( microp_scheme .eq. 'RK' ) then

     ! Determine repartition heating from change in cloud ice.

       repartht(:ncol,:pver) = (latice/dtime) * ( state1%q(:ncol,:pver,ixcldice) - repartht(:ncol,:pver) )

     ! Non-micro and non-macrophysical external advective forcings to compute net condensation rate. 
     ! Note that advective forcing of condensate is aggregated into liquid phase.

       qtend(:ncol,:pver) = ( state1%q(:ncol,:pver,1) - qcwat(:ncol,:pver) ) * rdtime
       ttend(:ncol,:pver) = ( state1%t(:ncol,:pver)   - tcwat(:ncol,:pver) ) * rdtime
       ltend(:ncol,:pver) = ( totcw   (:ncol,:pver)   - lcwat(:ncol,:pver) ) * rdtime
#ifdef DIRIND
   ltendx(:ncol,:pver) = (totcwx (:ncol,:pver)   - lcwatx(:ncol,:pver))*rdtime
#endif

#ifdef DIRIND
#ifndef DIAGNCDNC
   ! calculate the standard deviation of the vertical velocity
   do i=1,ncol		
      do k=1,pver
         deltaz(i,k) = 16._r8*log10(state%pint(i,k+1)/state%pint(i,k))	!TS
         !sigw1(i,k) = ((sqrt(2*pi)*kvq(i,k))/(1000._r8*deltaz(i,k)))*100._r8	!TS
         !corinna: Morrison and Gettelman (JClim 2008) subgrid vertical velocity
         sigw1(i,k) = 100._r8*0.5_r8*(kkvh(i,k+1)+kkvh(i,k))/30._r8
       end do
   end do
#endif
#endif

     ! Compute Stratiform Macro-Microphysical Tendencies

       ! Add rain and snow fluxes as output variables from pcond, and into physics buffer
       rkflxprc  => pbuf(ls_flxprc_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)
       rkflxsnw  => pbuf(ls_flxsnw_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,1)

       call t_startf('pcond')

#ifdef DIRIND
       call pcond( lchnk, ncol,                                                &
                   state1%t, ttend, state1%q(1,1,1), qtend, state1%omega,      &
                   totcw, state1%pmid , state1%pdel, cld, fice, fsnow,         &
                   qme, prain, prodsnow, nevapr, evapsnow, evapheat, prfzheat, &
                   meltheat, prec_pcw, snow_pcw, dtime, fwaut,                 &
                   fsaut, fracw, fsacw, fsaci, ltend,                          &
                   rhdfda, rhu00, icefrac, state1%zi, ice2pr, liq2pr,          &
                   liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio,&
                   ltendx, state1%q, precc, totcwx, sigw1, cldo, totnlx, dlf,  &       ! in
                   liq2prx, relca, ntotal, cdnc, cxsind, n1, n2, n3, n4, n5,   &       ! out
                   n6, n7, n8, n9, n10, n11, n12, n13, n14, nrmodes, relh,     &       ! out
                   cmex, selfx, lossx, nucrat, nrainx, supersat, evapx, freez )        ! out
#else
       call pcond( lchnk, ncol,                                                &
                   state1%t, ttend, state1%q(1,1,1), qtend, state1%omega,      &
                   totcw, state1%pmid , state1%pdel, cld, fice, fsnow,         &
                   qme, prain, prodsnow, nevapr, evapsnow, evapheat, prfzheat, &
                   meltheat, prec_pcw, snow_pcw, dtime, fwaut,                 &
                   fsaut, fracw, fsacw, fsaci, ltend,                          &
                   rhdfda, rhu00, icefrac, state1%zi, ice2pr, liq2pr,          &
                   liq2snow, snowh, rkflxprc, rkflxsnw, pracwo, psacwo, psacio )
#endif
       call t_stopf('pcond')


#ifdef DIRIND

!     calculate cdnc concentrations and effective droplet radii in liquid phase clouds.

      do i=1,ncol
       
        do k=1,pver
 	  pdelm1(i,k) = state1%pint(i,k+1) - state1%pint(i,k)
#ifdef DIAGNCDNC
!old  	  wpatha(i,k) = totcwx(i,k)*pdelm1(i,k)/gravit*1000._r8 &
  	  wpatha(i,k) = totcwl(i,k)*pdelm1(i,k)/gravit*1000._r8 &
                    /max(0.01_r8,cld(i,k))
#else   
  	  wpatha(i,k) = totcwlx(i,k)*pdelm1(i,k)/gravit*1000._r8 &
                    /max(0.01_r8,cld(i,k))
#endif
        end do
      end do

  do i=1,ncol
       rehanx(i) = 0._r8
       fochana(i) = 0._r8
#ifndef DIAGNCDNC
       cdncint(i) = 0._r8
#endif
  end do
  do k=1,pver
       do i=1,ncol
        if(cld(i,k) .gt. 0.01_r8 .and. state%t(i,k) .gt. 273.16_r8 .and. &
          wpatha(i,k).gt.30._r8 .and. fochana(i).eq.0._r8) then
          fochana(i) = 1._r8
          rehanx(i) = relca(i,k)
        endif
       end do
  end do

  do i=1,ncol
       rehana(i) = rehanx(i) * fochana(i)
  end do

  do k=1,pver
        do i=1,ncol
! Flag for liquid cloud presence. On a time averaged h-tape this
! will contain the frequency of occurence of cloudiness > 1%
          !if(cld(i,k) .gt. 0.01_r8 .and. state1%t(i,k) .gt. 273.16_r8) then
#ifdef DIAGNCDNC
          if(cld(i,k) .gt. 0.01_r8 .and. state1%t(i,k) .gt. 273.16_r8 .and. state1%q(i,k,ixcldliq)/cld(i,k) .gt. 1.e-6_r8) then
#else
          if(cld(i,k) .gt. 0.01_r8 .and. state1%t(i,k) .gt. 273.16_r8 .and. state1%q(i,k,icldlx)/cld(i,k) .gt. 1.e-6_r8) then
#endif
            cldfoc(i,k) = 1._r8
          else
            cldfoc(i,k) = 0._r8
          endif
          relcao(i,k) = relca(i,k) * cldfoc(i,k)
          cdnco(i,k) = cdnc(i,k) * cldfoc(i,k)
#ifndef DIAGNCDNC
          supersat(i,k) = supersat(i,k) * cldfoc(i,k)
          ccdnc(i,k) = cld(i,k)*cdnc(i,k)*deltaz(i,k)*1.e5_r8
          if(state1%t(i,k) > 273.16_r8) cdncint(i) = cdncint(i) + ccdnc(i,k)
          !nucrato(i,k) = nucrat(i,k) * cldfoc(i,k)
          nucrato(i,k) = nucrat(i,k)/MAX(1.e-4_r8,cld(i,k)) * cldfoc(i,k)
          evapxo(i,k) = evapx(i,k) * cldfoc(i,k)
          selfxo(i,k) = selfx(i,k) * cldfoc(i,k)
          nrainxo(i,k) = nrainx(i,k) * cldfoc(i,k)
          freezo(i,k) = freez(i,k) * cldfoc(i,k)
#endif
        end do
  end do

  !do i=1,ncol !corinna: needs cltot -> radiation
  !     cdncint(i) = (cdncint(i)/max(cltot(i),0.01_r8)) * fochana(i)
  !end do               

#endif


   end if

   if( microp_scheme .eq. 'RK' ) then

       do k = 1, pver
       do i = 1, ncol
          ptend_loc%s(i,k)          =   qme(i,k)*( latvap + latice*fice(i,k) ) + &
                                        evapheat(i,k) + prfzheat(i,k) + meltheat(i,k) + repartht(i,k)
          ptend_loc%q(i,k,1)        = - qme(i,k) + nevapr(i,k)
          ptend_loc%q(i,k,ixcldice) =   qme(i,k)*fice(i,k)         - ice2pr(i,k)
          ptend_loc%q(i,k,ixcldliq) =   qme(i,k)*(1._r8-fice(i,k)) - liq2pr(i,k)
#ifdef DIRIND
#ifdef DIAGNCDNC
          ptend_loc%q(i,k,icldlx) = 0.0_r8
          ptend_loc%q(i,k,icldnx) = 0.0_r8
#else
          ptend_loc%q(i,k,icldlx) = cmex(i,k)*(1._r8-fice(i,k)) - liq2prx(i,k)
!TS++	Tendency for cloud droplet number concentration
          moistrho(i,k) = state1%pmid(i,k)/(287._r8*(1+0.6_r8*state1%q(i,k,1))*state1%t(i,k))
          ptend_loc%q(i,k,icldnx) = 1.e6_r8*(nucrat(i,k) + lossx(i,k) - nrainx(i,k))*(moistrho(i,k)**(-1))
!TS--
#endif
#endif
       end do
       end do
 
       do k = 1, pver
       do i = 1, ncol
          aist(i,k)  = cld(i,k)
          alst(i,k)  = cld(i,k)
          ast(i,k)   = cld(i,k)
          icimr(i,k) = (state1%q(i,k,ixcldice) + dtime*ptend_loc%q(i,k,ixcldice)) / max(0.01_r8,aist(i,k))
          icwmr(i,k) = (state1%q(i,k,ixcldliq) + dtime*ptend_loc%q(i,k,ixcldliq)) / max(0.01_r8,alst(i,k))
       end do
       end do

    ! Convert precipitation from [ kg/m2 ] to [ m/s ]
      snow_pcw(:ncol) = snow_pcw(:ncol)/1000._r8
      prec_pcw(:ncol) = prec_pcw(:ncol)/1000._r8

      do k = 1, pver
      do i = 1, ncol
         cmeheat(i,k) = qme(i,k) * ( latvap + latice*fice(i,k) )
         cmeice (i,k) = qme(i,k) *   fice(i,k)
         cmeliq (i,k) = qme(i,k) * ( 1._r8 - fice(i,k) )
      end do
      end do

    ! Record history variables

      call outfld( 'FWAUT'   , fwaut,       pcols, lchnk )
      call outfld( 'FSAUT'   , fsaut,       pcols, lchnk )
      call outfld( 'FRACW'   , fracw,       pcols, lchnk )
      call outfld( 'FSACW'   , fsacw,       pcols, lchnk )
      call outfld( 'FSACI'   , fsaci,       pcols, lchnk )

      call outfld( 'PCSNOW'  , snow_pcw,    pcols, lchnk )
      call outfld( 'FICE'    , fice,        pcols, lchnk )
      call outfld( 'CMEICE'  , cmeice,      pcols, lchnk )
      call outfld( 'CMELIQ'  , cmeliq,      pcols, lchnk )
      call outfld( 'ICE2PR'  , ice2pr,      pcols, lchnk )
      call outfld( 'LIQ2PR'  , liq2pr,      pcols, lchnk )
      call outfld( 'HPROGCLD', ptend_loc%s, pcols, lchnk )
      call outfld( 'HEVAP   ', evapheat,    pcols, lchnk )
      call outfld( 'HMELT'   , meltheat,    pcols, lchnk )
      call outfld( 'HCME'    , cmeheat ,    pcols, lchnk )
      call outfld( 'HFREEZ'  , prfzheat,    pcols, lchnk )
      call outfld( 'HREPART' , repartht,    pcols, lchnk )

#ifdef DIRIND
   call outfld('CDNC    ',cdnco   ,pcols,lchnk)
   call outfld('CLDFOC  ',cldfoc  ,pcols,lchnk)
   call outfld('REHANA  ',rehana  ,pcols,lchnk)
   call outfld('FOCHANA ',fochana ,pcols,lchnk)
   call outfld('REFFL   ',relcao  ,pcols,lchnk)
#ifndef DIAGNCDNC
   call outfld('RELH    ',relh    ,pcols,lchnk)
   call outfld('SUPERSAT',supersat,pcols,lchnk)
   call outfld('NUCRAT  ',nucrato ,pcols,lchnk)
   call outfld('EVAPX   ',evapxo  ,pcols,lchnk)
   call outfld('SELFX   ',selfxo  ,pcols,lchnk)
   call outfld('NRAINX  ',nrainxo ,pcols,lchnk)
   call outfld('FREEZ   ',freezo  ,pcols,lchnk)
#endif
#endif

      call outfld('LS_FLXPRC', rkflxprc,    pcols, lchnk )
      call outfld('LS_FLXSNW', rkflxsnw,    pcols, lchnk )
      call outfld('PRACWO'   , pracwo,      pcols, lchnk )
      call outfld('PSACWO'   , psacwo,      pcols, lchnk )
      call outfld('PSACIO'   , psacio,      pcols, lchnk )

      ! initialize local variables
      mr_ccliq(1:ncol,1:pver)=0._r8
      mr_ccice(1:ncol,1:pver)=0._r8
      mr_lsliq(1:ncol,1:pver)=0._r8
      mr_lsice(1:ncol,1:pver)=0._r8

      do k=1,pver
      do i=1,ncol
       if (cld(i,k) .gt. 0._r8) then
         mr_ccliq(i,k) = (state%q(i,k,ixcldliq)/cld(i,k))*concld(i,k)
         mr_ccice(i,k) = (state%q(i,k,ixcldice)/cld(i,k))*concld(i,k)
         mr_lsliq(i,k) = (state%q(i,k,ixcldliq)/cld(i,k))*(cld(i,k)-concld(i,k))
         mr_lsice(i,k) = (state%q(i,k,ixcldice)/cld(i,k))*(cld(i,k)-concld(i,k))
       else
         mr_ccliq(i,k) = 0._r8
         mr_ccice(i,k) = 0._r8
         mr_lsliq(i,k) = 0._r8
         mr_lsice(i,k) = 0._r8
       end if
      end do
      end do

      call outfld( 'CLDLIQSTR  ', mr_lsliq,    pcols, lchnk )
      call outfld( 'CLDICESTR  ', mr_lsice,    pcols, lchnk )
      call outfld( 'CLDLIQCON  ', mr_ccliq,    pcols, lchnk )
      call outfld( 'CLDICECON  ', mr_ccice,    pcols, lchnk )

   endif

   ! ------------------------------- !
   ! Update microphysical tendencies !
   ! ------------------------------- !

   call physics_ptend_sum( ptend_loc, ptend_all, state )
   ptend_all%name = 'stratiform'
   call physics_update( state1, tend, ptend_loc, dtime )
   call physics_ptend_init( ptend_loc ) 

   if ( microp_scheme .eq. 'RK' .and. .not. cam_physpkg_is('cam3')) then

      call t_startf("cldfrc")
      call cldfrc( lchnk, ncol, pbuf,                                                 &
                   state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                   shfrc, use_shfrc,                                                  &
                   cld, rhcloud, clc, state1%pdel,                                    &
                   cmfmc, cmfmc2, landfrac, snowh, concld, cldst,                     &
                   ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu00,                &
                   state1%q(:,:,ixcldice), icecldf, liqcldf,                          &
                   relhum, 0 )    
      call cldfrc( lchnk, ncol, pbuf,                                                 &
                   state1%pmid, state1%t, state1%q(:,:,1), state1%omega, state1%phis, &
                   shfrc, use_shfrc,                                                  &
                   cld2, rhcloud2, clc, state1%pdel,                                  &
                   cmfmc, cmfmc2, landfrac, snowh, concld2, cldst2,                   &
                   ts, sst, state1%pint(:,pverp), zdu, ocnfrac, rhu002,               &
                   state1%q(:,:,ixcldice), icecldf2, liqcldf2,                        &
                   relhum2, 1  )              

      call t_stopf("cldfrc")

      do k = 1, pver
         do i = 1, ncol
            if( relhum(i,k) < rhu00(i,k) ) then
               rhdfda(i,k)=0.0_r8
            elseif( relhum(i,k) >= 1.0_r8 ) then
               rhdfda(i,k)=0.0_r8
            else
           ! Under certain circumstances, rh+ causes cld not to be changed
           ! when at an upper limit, or w/ strong subsidence
               if( ( rhcloud2(i,k) - rhcloud(i,k) ) < 1.e-4_r8 ) then
                  rhdfda(i,k) = 0.01_r8*relhum(i,k)*1.e+4_r8 
               else
                  rhdfda(i,k) = 0.01_r8*relhum(i,k)/(rhcloud2(i,k)-rhcloud(i,k))
               endif
            endif
         enddo
      enddo

   endif

 ! Copy of concld/fice to put in physics buffer
 ! Below are used only for convective cloud.

   concld_ql(:ncol,:pver) = concld(:ncol,:pver)
   fice_ql(:ncol,:pver)   = fice(:ncol,:pver)

   if( micro_treatment .eq. 'inter' ) then
       icecldf(:ncol,:pver) = ast(:ncol,:pver)
       liqcldf(:ncol,:pver) = ast(:ncol,:pver)
   elseif( micro_treatment .eq. 'compl' ) then
       icecldf(:ncol,:pver) = aist(:ncol,:pver)
       liqcldf(:ncol,:pver) = alst(:ncol,:pver)
   endif

   call outfld( 'CONCLD  ', concld, pcols, lchnk )
   call outfld( 'CLDST   ', cldst,  pcols, lchnk )
   call outfld( 'CNVCLD  ', clc,    pcols, lchnk )

   call outfld( 'ICECLDF ', aist,   pcols, lchnk )
   call outfld( 'LIQCLDF', alst,   pcols, lchnk )
   call outfld( 'AST',      ast,    pcols, lchnk )   

   if( microp_scheme .eq. 'RK' ) then

     do k = 1, pver
     do i = 1, ncol
	iwc(i,k)   = state1%q(i,k,ixcldice)*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
	lwc(i,k)   = state1%q(i,k,ixcldliq)*state1%pmid(i,k)/(287.15_r8*state1%t(i,k))
        icimr(i,k) = state1%q(i,k,ixcldice) / max(0.01_r8,rhcloud(i,k))
        icwmr(i,k) = state1%q(i,k,ixcldliq) / max(0.01_r8,rhcloud(i,k))
     end do
     end do

   endif ! RK,MG microphysics

   ! --------------------------------------------- !
   ! Common outfield calls for either microphysics !
   ! --------------------------------------------- !

   call outfld( 'IWC'      , iwc,         pcols, lchnk )
   call outfld( 'LWC'      , lwc,         pcols, lchnk )
   call outfld( 'ICIMR'    , icimr,       pcols, lchnk )
   call outfld( 'ICWMR'    , icwmr,       pcols, lchnk )
#ifdef DIRIND
   call outfld('ICWMRX',icwmrx,  pcols,lchnk)
#endif
   call outfld( 'ICIMRST'  , icimrst_out, pcols, lchnk )
   call outfld( 'ICWMRST'  , icwmrst_out, pcols, lchnk )
   call outfld( 'CME'      , qme,         pcols, lchnk )
   call outfld( 'PRODPREC' , prain,       pcols, lchnk )
   call outfld( 'EVAPPREC' , nevapr,      pcols, lchnk )
   call outfld( 'EVAPSNOW' , evapsnow,    pcols, lchnk )

   call t_stopf('stratiform_microphys')

   if( microp_scheme .eq. 'RK' ) then

      prec_str(:ncol) = prec_str(:ncol) + prec_pcw(:ncol)
      snow_str(:ncol) = snow_str(:ncol) + snow_pcw(:ncol)

      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
#ifdef DIRIND
      call cnst_get_ind('CLDLIX', icldlx)
#endif

    ! Save variables for use in the macrophysics at the next time step

      do k = 1, pver
         qcwat(:ncol,k) = state1%q(:ncol,k,1)
         tcwat(:ncol,k) = state1%t(:ncol,k)
         lcwat(:ncol,k) = state1%q(:ncol,k,ixcldice) + state1%q(:ncol,k,ixcldliq)
#ifdef DIRIND
      lcwatx(:ncol,k) = state1%q(:ncol,k,ixcldice) + state1%q(:ncol,k,icldlx)
#endif
      end do
  
    ! Cloud water and ice particle sizes, saved in physics buffer for radiation

      call cldefr( lchnk, ncol, landfrac, state1%t, rel, rei, state1%ps, state1%pmid, landm, icefrac, snowh )
      rel2(:ncol,:pver) = rel(:ncol,:pver)      
      rei2(:ncol,:pver) = rei(:ncol,:pver)      

   end if

   end subroutine stratiform_tend

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_1(state1,ptend,i,k, &
        dtime,qme,fice,snow_pcw,prec_pcw, &
        prain,nevapr,prodsnow, evapsnow, &
        ice2pr,liq2pr,liq2snow)

     use physics_types, only: physics_state, physics_ptend
     use physconst,     only: tmelt

     implicit none
     
     integer, intent(in) :: i,k
     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     type(physics_ptend), intent(in) :: ptend  ! local copy of the ptend variable
     real(r8), intent(in)  :: dtime                ! timestep
     real(r8), intent(in) :: qme(pcols,pver)          ! local condensation - evaporation of cloud water

     real(r8), intent(in) :: prain(pcols,pver)          ! local production of precipitation
     real(r8), intent(in) :: nevapr(pcols,pver)          ! local evaporation of precipitation
     real(r8), intent(in) :: prodsnow(pcols,pver)          ! local production of snow
     real(r8), intent(in) :: evapsnow(pcols,pver)          ! local evaporation of snow
     real(r8), intent(in) :: ice2pr(pcols,pver)   ! rate of conversion of ice to precip
     real(r8), intent(in) :: liq2pr(pcols,pver)   ! rate of conversion of liquid to precip
     real(r8), intent(in) :: liq2snow(pcols,pver)   ! rate of conversion of liquid to snow
     real(r8), intent(in) :: fice    (pcols,pver)          ! Fractional ice content within cloud
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: prec_pcw(pcols)

     real(r8) hs1, qv1, ql1, qi1, qs1, qr1, fice2, pr1, w1, w2, w3, fliq, res
     real(r8) w4, wl, wv, wi, wlf, wvf, wif, qif, qlf, qvf

     pr1 = 0
     hs1 = 0
     qv1 = 0
     ql1 = 0
     qi1 = 0
     qs1 = 0
     qr1 = 0
     w1 = 0
     wl = 0
     wv = 0
     wi = 0
     wlf = 0
     wvf = 0 
     wif = 0


     write(iulog,*) 
     write(iulog,*) ' input state, t, q, l, i ', k, state1%t(i,k), state1%q(i,k,1), state1%q(i,k,ixcldliq),  state1%q(i,k,ixcldice)
     write(iulog,*) ' rain, snow, total from components before accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip before accumulation                      ', k, pr1

     wv = wv + state1%q(i,k,1       )*state1%pdel(i,k)/gravit
     wl = wl + state1%q(i,k,ixcldliq)*state1%pdel(i,k)/gravit
     wi = wi + state1%q(i,k,ixcldice)*state1%pdel(i,k)/gravit

     qvf = state1%q(i,k,1) + ptend%q(i,k,1)*dtime
     qlf = state1%q(i,k,ixcldliq) + ptend%q(i,k,ixcldliq)*dtime
     qif = state1%q(i,k,ixcldice) + ptend%q(i,k,ixcldice)*dtime

     if (qvf.lt.0._r8) then
        write(iulog,*) ' qvf is negative *******', qvf
     endif
     if (qlf.lt.0._r8) then
        write(iulog,*) ' qlf is negative *******', qlf
     endif
     if (qif.lt.0._r8) then
        write(iulog,*) ' qif is negative *******', qif
     endif
     write(iulog,*) ' qvf, qlf, qif ', qvf, qlf, qif

     wvf = wvf + qvf*state1%pdel(i,k)/gravit
     wlf = wlf + qlf*state1%pdel(i,k)/gravit
     wif = wif + qif*state1%pdel(i,k)/gravit

     hs1 = hs1 + ptend%s(i,k)*state1%pdel(i,k)/gravit
     pr1 = pr1 + state1%pdel(i,k)/gravit*(prain(i,k)-nevapr(i,k))
     qv1 = qv1 - (qme(i,k)-nevapr(i,k))*state1%pdel(i,k)/gravit    ! vdot
     w1  = w1  + (qme(i,k)-prain(i,k))*state1%pdel(i,k)/gravit    ! cdot
     qi1 = qi1 + ((qme(i,k))*fice(i,k)        -ice2pr(i,k) )*state1%pdel(i,k)/gravit   ! idot
     ql1 = ql1 + ((qme(i,k))*(1._r8-fice(i,k))-liq2pr(i,k) )*state1%pdel(i,k)/gravit   ! ldot

     qr1 = qr1 &
          + ( liq2pr(i,k)-liq2snow(i,k)   &     ! production of rain
          -(nevapr(i,k)-evapsnow(i,k)) &     ! rain evaporation
          )*state1%pdel(i,k)/gravit
     qs1 = qs1 &
          + ( ice2pr(i,k) + liq2snow(i,k) &     ! production of snow.Note last term has phase change
          -evapsnow(i,k)               &     ! snow evaporation
          )*state1%pdel(i,k)/gravit

     if (state1%t(i,k).gt.tmelt) then
        qr1 = qr1 + qs1
        qs1 = 0._r8
     endif
     write(iulog,*) ' rain, snow, total after accumulation ', qr1, qs1, qr1+qs1
     write(iulog,*) ' total precip after accumulation      ', k, pr1
     write(iulog,*)
     write(iulog,*) ' layer prain, nevapr, pdel ', prain(i,k), nevapr(i,k), state1%pdel(i,k)
     write(iulog,*) ' layer prodsnow, ice2pr+liq2snow ', prodsnow(i,k), ice2pr(i,k)+liq2snow(i,k)
     write(iulog,*) ' layer prain-prodsnow, liq2pr-liq2snow ', prain(i,k)-prodsnow(i,k), liq2pr(i,k)-liq2snow(i,k)
     write(iulog,*) ' layer evapsnow, evaprain ', k, evapsnow(i,k), nevapr(i,k)-evapsnow(i,k)
     write(iulog,*) ' layer ice2pr, liq2pr, liq2snow ', ice2pr(i,k), liq2pr(i,k), liq2snow(i,k)
     write(iulog,*) ' layer ice2pr+liq2pr, prain ', ice2pr(i,k)+liq2pr(i,k), prain(i,k)
     write(iulog,*)
     write(iulog,*) ' qv1 vapor removed from col after accum  (vdot)   ', k, qv1
     write(iulog,*) ' - (precip produced - vapor removed) after accum  ', k, -pr1-qv1
     write(iulog,*) ' condensate produce after accum                   ', k, w1
     write(iulog,*) ' liq+ice tends accum                              ', k, ql1+qi1
     write(iulog,*) ' change in total water after accum                ', k, qv1+ql1+qi1
     write(iulog,*) ' imbalance in colum after accum                   ', k, qs1+qr1+qv1+ql1+qi1
     write(iulog,*) ' fice at this lev ', fice(i,k)
     write(iulog,*)

     res = abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1),abs(ql1),abs(qi1),abs(qs1),abs(qr1),1.e-36_r8))
     write(iulog,*) ' relative residual in column method 1             ', k, res

     write(iulog,*) ' relative residual in column method 2             ',&
	 k, abs((qs1+qr1+qv1+ql1+qi1)/max(abs(qv1+ql1+qi1),1.e-36_r8))
     !            if (abs((qs1+qr1+qv1+ql1+qi1)/(qs1+qr1+1.e-36_r8)).gt.1.e-14_r8) then
     if (res.gt.1.e-14_r8) then
        call endrun ('STRATIFORM_TEND')
     endif

     !             w3  = qme(i,k) * (latvap + latice*fice(i,k)) &
     !               + evapheat(i,k) + prfzheat(i,k) + meltheat(i,k)

     res = qs1+qr1-pr1
     w4 = max(abs(qs1),abs(qr1),abs(pr1)) 
     if (w4.gt.0._r8)  then
        if (res/w4.gt.1.e-14_r8) then
           write(iulog,*) ' imbalance in precips calculated two ways '
           write(iulog,*) ' res/w4, pr1, qr1, qs1, qr1+qs1 ', &
                res/w4, pr1, qr1, qs1, qr1+qs1
           !                   call endrun()
        endif
     endif
     if (k.eq.pver) then
        write(iulog,*) ' pcond returned precip, rain and snow rates ', prec_pcw(i), prec_pcw(i)-snow_pcw(i), snow_pcw(i)
        write(iulog,*) ' I calculate ', pr1, qr1, qs1
        !               call endrun
        write(iulog,*) ' byrons water check ', wv+wl+wi-pr1*dtime, wvf+wlf+wif
     endif
     write(iulog,*)


   end subroutine debug_microphys_1

  !============================================================================ !
  !                                                                             !
  !============================================================================ !

   subroutine debug_microphys_2(state1,&
        snow_pcw,fsaut,fsacw ,fsaci, meltheat)

     use ppgrid,        only: pver
     use physconst,     only: tmelt
     use physics_types, only: physics_state
     
     implicit none

     type(physics_state), intent(in) :: state1   ! local copy of the state variable
     real(r8), intent(in) :: snow_pcw(pcols)
     real(r8), intent(in) :: fsaut(pcols,pver)              
     real(r8), intent(in) :: fsacw(pcols,pver)              
     real(r8), intent(in) :: fsaci(pcols,pver)              
     real(r8), intent(in) :: meltheat(pcols,pver)          ! heating rate due to phase change of precip


     integer  i,ncol,lchnk


     ncol = state1%ncol
     lchnk = state1%lchnk
     
     do i = 1,ncol
        if (snow_pcw(i) .gt. 0.01_r8/8.64e4_r8  .and.  state1%t(i,pver) .gt. tmelt) then
           write(iulog,*) ' stratiform: snow, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
        
        if (snow_pcw(i)*8.64e4_r8 .lt. -1.e-5_r8) then
           write(iulog,*) ' neg snow ', snow_pcw(i)*8.64e4_r8
           write(iulog,*) ' stratiform: snow_pcw, temp, ', i, lchnk, &
                snow_pcw(i), state1%t(i,pver)
           write(iulog,*) ' t ', state1%t(i,:)
           write(iulog,*) ' fsaut ', fsaut(i,:)
           write(iulog,*) ' fsacw ', fsacw(i,:)
           write(iulog,*) ' fsaci ', fsaci(i,:)
           write(iulog,*) ' meltheat ', meltheat(i,:)
           call endrun ('STRATIFORM_TEND')
        endif
     end do
     
   end subroutine debug_microphys_2

  end module stratiform
