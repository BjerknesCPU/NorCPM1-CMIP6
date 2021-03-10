module cam_diagnostics

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

use shr_kind_mod,  only: r8 => shr_kind_r8
use infnan,        only: nan
use camsrfexch_types, only: cam_in_t, cam_out_t
use physics_types, only: physics_state, physics_tend
use ppgrid,        only: pcols, pver, pverp, begchunk, endchunk
use cam_history,   only: outfld, write_inithist, hist_fld_active
use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, dmetendnam, apcnst, bpcnst, &
                         cnst_get_ind
use chemistry,     only: chem_is
use abortutils,    only: endrun
use scamMod,       only: single_column,wfld
use dycore,        only: dycore_is
use phys_control,  only: phys_getopts
use wv_saturation, only: aqsat, aqsat_water, polysvp

implicit none
private
save

! Public interfaces

public :: &
   diag_register,      &! register pbuf space
   diag_defaultopts,   &! set default values of namelist variables
   diag_setopts,       &! get namelist input
   diag_init,          &! initialization
   diag_allocate,      &! allocate memory for module variables
   diag_deallocate,    &! deallocate memory for module variables
   diag_conv_tend_ini, &! initialize convective tendency calcs
   diag_phys_writeout, &! output diagnostics of the dynamics
   diag_phys_tend_writeout, & ! output physics tendencies
   diag_state_b4_phys_write,& ! output state before physics execution
   diag_conv,          &! output diagnostics of convective processes
   diag_surf,          &! output diagnostics of the surface
   diag_export,        &! output export state
   diag_physvar_ic

logical, public :: inithist_all = .false. ! Flag to indicate set of fields to be 
                                          ! included on IC file
                                          !  .false.  include only required fields
                                          !  .true.   include required *and* optional fields

! Private data

real(r8), allocatable :: &
   dtcond(:,:,:),  &! temperature tendency due to convection
   dqcond(:,:,:,:)  ! constituent tendencies due to convection

character(len=8) :: diag_cnst_conv_tend = 'q_only' ! output constituent tendencies due to convection
                                                   ! 'none', 'q_only' or 'all'
logical          :: history_budget                 ! output tendencies and state variables for CAM4
                                                   ! temperature, water vapor, cloud ice and cloud
                                                   ! liquid budgets.
integer          :: history_budget_histfile_num    ! output history file number for budget fields

character(len=16), public :: dcconnam(pcnst)       ! names of convection tendencies

!Physics buffer indices
integer  ::      qcwat_idx  = 0 
integer  ::      tcwat_idx  = 0 
integer  ::      lcwat_idx  = 0 
integer  ::      cld_idx    = 0 
integer  ::      concld_idx = 0 
integer  ::      tke_idx    = 0 
integer  ::      kvm_idx    = 0 
integer  ::      kvh_idx    = 0 
integer  ::      cush_idx   = 0 
integer  ::      t_ttend_idx 

contains

  ! ===============================================================================

  subroutine diag_register

    use phys_buffer,  only: pbuf_times, pbuf_add

  ! Request physics buffer space for fields that persist across timesteps.

    call pbuf_add('T_TTEND',   'global',  1, pver, pbuf_times,   t_ttend_idx)

  end subroutine diag_register

!===============================================================================

subroutine diag_defaultopts(diag_cnst_conv_tend_out)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(out), optional :: diag_cnst_conv_tend_out
!-----------------------------------------------------------------------
   if ( present(diag_cnst_conv_tend_out) ) then
      if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then
         diag_cnst_conv_tend_out = 'none'
      else
         diag_cnst_conv_tend_out = diag_cnst_conv_tend
      end if
   endif
end subroutine diag_defaultopts

!================================================================================================

subroutine diag_setopts(diag_cnst_conv_tend_in)
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
!-----------------------------------------------------------------------

   character(len=*), intent(in), optional :: diag_cnst_conv_tend_in
!-----------------------------------------------------------------------
   if ( present(diag_cnst_conv_tend_in) ) then
      diag_cnst_conv_tend = diag_cnst_conv_tend_in
   endif
end subroutine diag_setopts

!================================================================================================

subroutine diag_init

  ! Declare the history fields for which this module contains outfld calls.

   use cam_history,        only: addfld, add_default, phys_decomp
   use constituent_burden, only: constituent_burden_init
   use cam_control_mod,    only: moist_physics, ideal_phys
   use tidal_diag,         only: tidal_diag_init 
   use phys_buffer,        only: pbuf_get_fld_idx

   integer :: k, m
   ! Note - this is a duplication of information in ice_constants 
   ! Cannot put in a use statement if want to swap ice models to cice4
   integer, parameter :: plevmx = 4       ! number of subsurface levels
   character(len=8), parameter :: tsnam(plevmx) = (/ 'TS1', 'TS2', 'TS3', 'TS4' /)
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.

   ! outfld calls in diag_phys_writeout

   call addfld ('NSTEP   ','timestep',1,    'A','Model timestep',phys_decomp)
   call addfld ('PHIS    ','m2/s2   ',1,    'I','Surface geopotential',phys_decomp)

   call addfld ('PS      ','Pa      ',1,    'A','Surface pressure',phys_decomp)
   call addfld ('T       ','K       ',pver, 'A','Temperature',phys_decomp)
   call addfld ('U       ','m/s     ',pver, 'A','Zonal wind',phys_decomp)
   call addfld ('V       ','m/s     ',pver, 'A','Meridional wind',phys_decomp)
   call addfld (cnst_name(1),'kg/kg ',pver, 'A',cnst_longname(1),phys_decomp)

   ! State before physics
   call addfld ('TBP     ','K       ',pver, 'A','Temperature (before physics)'       ,phys_decomp)
   call addfld (bpcnst(1) ,'kg/kg   ',pver, 'A',cnst_longname(1)//' (before physics)',phys_decomp)
   ! State after physics
   call addfld ('TAP     ','K       ',pver, 'A','Temperature (after physics)'       ,phys_decomp)
   call addfld ('UAP     ','m/s     ',pver, 'A','Zonal wind (after physics)'        ,phys_decomp)
   call addfld ('VAP     ','m/s     ',pver, 'A','Meridional wind (after physics)'   ,phys_decomp)
   call addfld (apcnst(1) ,'kg/kg   ',pver, 'A',cnst_longname(1)//' (after physics)',phys_decomp)
   if ( dycore_is('LR') ) then
      call addfld ('TFIX    ','K/s     ',1,    'A'     ,'T fixer (T equivalent of Energy correction)',phys_decomp)
      call addfld ('PTTEND_RESID','K/s ',pver, 'A'     ,&
                   'T-tendency due to BAB kluge at end of tphysac (diagnostic not part of T-budget)' ,phys_decomp)
   end if
   call addfld ('TTEND_TOT   ','K/s' ,pver, 'A','Total temperature tendency'   ,phys_decomp)
  
   ! column burdens for all constituents except water vapor
   call constituent_burden_init

   call addfld ('Z3      ','m       ',pver, 'A','Geopotential Height (above sea level)',phys_decomp)
   call addfld ('Z700    ','m       ',1,    'A','Geopotential Z at 700 mbar pressure surface',phys_decomp)
   call addfld ('Z500    ','m       ',1,    'A','Geopotential Z at 500 mbar pressure surface',phys_decomp)
   call addfld ('Z300    ','m       ',1,    'A','Geopotential Z at 300 mbar pressure surface',phys_decomp)
   call addfld ('Z250    ','m       ',1,    'A','Geopotential Z at 250 mbar pressure surface',phys_decomp)
   call addfld ('Z200    ','m       ',1,    'A','Geopotential Z at 200 mbar pressure surface',phys_decomp)
   call addfld ('Z100    ','m       ',1,    'A','Geopotential Z at 100 mbar pressure surface',phys_decomp)
   call addfld ('Z050    ','m       ',1,    'A','Geopotential Z at 50 mbar pressure surface',phys_decomp)
   call addfld ('Z010    ','m       ',1,    'A','Geopotential Z at 10 mbar pressure surface',phys_decomp)

   call addfld ('ZZ      ','m2      ',pver, 'A','Eddy height variance' ,phys_decomp)
   call addfld ('VZ      ','m2/s    ',pver, 'A','Meridional transport of geopotential energy',phys_decomp)
   call addfld ('VT      ','K m/s   ',pver, 'A','Meridional heat transport',phys_decomp)
   call addfld ('VU      ','m2/s2   ',pver, 'A','Meridional flux of zonal momentum' ,phys_decomp)
   call addfld ('VV      ','m2/s2   ',pver, 'A','Meridional velocity squared' ,phys_decomp)
   call addfld ('VQ      ','m/skg/kg',pver, 'A','Meridional water transport',phys_decomp)
   call addfld ('QQ      ','kg2/kg2 ',pver, 'A','Eddy moisture variance',phys_decomp)
   call addfld ('OMEGAV  ','m Pa/s2 ',pver ,'A','Vertical flux of meridional momentum' ,phys_decomp)
   call addfld ('OMGAOMGA','Pa2/s2  ',pver ,'A','Vertical flux of vertical momentum' ,phys_decomp)
   call addfld ('OMEGAQ  ','kgPa/kgs',pver ,'A','Vertical water transport' ,phys_decomp)

   call addfld ('UU      ','m2/s2   ',pver, 'A','Zonal velocity squared' ,phys_decomp)
   call addfld ('WSPEED  ','m/s     ',pver, 'X','Horizontal total wind speed' ,phys_decomp)

   call addfld ('OMEGA   ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)
   call addfld ('OMEGAT  ','K Pa/s  ',pver, 'A','Vertical heat flux' ,phys_decomp)
   call addfld ('OMEGAU  ','m Pa/s2 ',pver, 'A','Vertical flux of zonal momentum' ,phys_decomp)
   call addfld ('OMEGA850','Pa/s    ',1,    'A','Vertical velocity at 850 mbar pressure surface',phys_decomp)
   call addfld ('OMEGA500','Pa/s    ',1,    'A','Vertical velocity at 500 mbar pressure surface',phys_decomp)

   call addfld ('MQ      ','kg/m2   ',pver, 'A','Water vapor mass in layer',phys_decomp)
   call addfld ('TMQ     ','kg/m2   ',1,    'A','Total (vertically integrated) precipitatable water',phys_decomp)
   call addfld ('RELHUM  ','percent ',pver, 'A','Relative humidity',phys_decomp)
   call addfld ('RHW  ','percent '   ,pver, 'A','Relative humidity with respect to liquid',phys_decomp)
   call addfld ('RHI  ','percent '   ,pver, 'A','Relative humidity with respect to ice',phys_decomp)
   call addfld ('RHCFMIP','percent ' ,pver, 'A','Relative humidity with respect to water above 273 K, ice below 273 K',phys_decomp)
   call addfld ('PSL     ','Pa      ',1,    'A','Sea level pressure',phys_decomp)

   call addfld ('T850    ','K       ',1,    'A','Temperature at 850 mbar pressure surface',phys_decomp)
   call addfld ('T500    ','K       ',1,    'A','Temperature at 500 mbar pressure surface',phys_decomp)
   call addfld ('T300    ','K       ',1,    'A','Temperature at 300 mbar pressure surface',phys_decomp)
   call addfld ('T200    ','K       ',1,    'A','Temperature at 200 mbar pressure surface',phys_decomp)
   call addfld ('T010    ','K       ',1,    'A','Temperature at 10 mbar pressure surface',phys_decomp)
   call addfld ('Q850    ','kg/kg   ',1,    'A','Specific Humidity at 850 mbar pressure surface',phys_decomp)
   call addfld ('Q200    ','kg/kg   ',1,    'A','Specific Humidity at 200 mbar pressure surface',phys_decomp)
   call addfld ('U850    ','m/s     ',1,    'A','Zonal wind at 850 mbar pressure surface',phys_decomp)
   call addfld ('U200    ','m/s     ',1,    'A','Zonal wind at 200 mbar pressure surface',phys_decomp)
   call addfld ('U010    ','m/s     ',1,    'A','Zonal wind at 10 mbar pressure surface',phys_decomp)
   call addfld ('U250    ','m/s     ',1,    'A','Zonal wind at 250 mbar pressure surface',phys_decomp)
   call addfld ('U500    ','m/s     ',1,    'A','Zonal wind at 500 mbar pressure surface',phys_decomp)
   call addfld ('V850    ','m/s     ',1,    'A','Meridional wind at 850 mbar pressure surface',phys_decomp)
   call addfld ('V200    ','m/s     ',1,    'A','Meridional wind at 200 mbar pressure surface',phys_decomp)
   call addfld ('V250    ','m/s     ',1,    'A','Meridional wind at 250 mbar pressure surface',phys_decomp)
   call addfld ('V500    ','m/s     ',1,    'A','Meridional wind at 500 mbar pressure surface',phys_decomp)

   call addfld ('TT      ','K2      ',pver, 'A','Eddy temperature variance' ,phys_decomp)

   call addfld ('UBOT    ','m/s     ',1,    'A','Lowest model level zonal wind',phys_decomp)
   call addfld ('VBOT    ','m/s     ',1,    'A','Lowest model level meridional wind',phys_decomp)
   call addfld ('QBOT    ','kg/kg   ',1,    'A','Lowest model level water vapor mixing ratio',phys_decomp)
   call addfld ('ZBOT    ','m       ',1,    'A','Lowest model level height', phys_decomp)

   call addfld ('ATMEINT  ','J/m2    ',1, 'A','Vertically integrated total atmospheric energy ',phys_decomp)

   call addfld ('T1000      ','K     ',1,   'A','Temperature at 1000 mbar pressure surface',phys_decomp)
   call addfld ('T925       ','K     ',1,   'A','Temperature at 925 mbar pressure surface',phys_decomp)   
   call addfld ('T700       ','K     ',1,   'A','Temperature at 700 mbar pressure surface',phys_decomp)
   call addfld ('Q1000      ','kg/kg ',1,   'A','Specific Humidity at 1000 mbar pressure surface',phys_decomp)   
   call addfld ('Q925       ','kg/kg ',1,   'A','Specific Humidity at 925 mbar pressure surface',phys_decomp)

   call addfld ('T7001000   ','K     ',1,   'A','Temperature difference 700 mb - 1000 mb',phys_decomp)
   call addfld ('TH7001000  ','K     ',1,   'A','Theta difference 700 mb - 1000 mb',phys_decomp)
   call addfld ('THE7001000 ','K     ',1,   'A','ThetaE difference 700 mb - 1000 mb',phys_decomp)

   call addfld ('T8501000   ','K     ',1,   'A','Temperature difference 850 mb - 1000 mb',phys_decomp)
   call addfld ('TH8501000  ','K     ',1,   'A','Theta difference 850 mb - 1000 mb',phys_decomp)   
   call addfld ('THE8501000 ','K     ',1,   'A','ThetaE difference 850 mb - 1000 mb',phys_decomp)
   call addfld ('T9251000   ','K     ',1,   'A','Temperature difference 925 mb - 1000 mb',phys_decomp) 
   call addfld ('TH9251000  ','K     ',1,   'A','Theta difference 925 mb - 1000 mb',phys_decomp)   
   call addfld ('THE9251000 ','K     ',1,   'A','ThetaE difference 925 mb - 1000 mb',phys_decomp) 

   ! defaults
   call add_default ('PHIS    '  , 1, ' ')
   call add_default ('PS      '  , 1, ' ')
   call add_default ('T       '  , 1, ' ')
   call add_default ('U       '  , 1, ' ')
   call add_default ('V       '  , 1, ' ')
   call add_default (cnst_name(1), 1, ' ')
   call add_default ('Z3      '  , 1, ' ')
   call add_default ('OMEGA   '  , 1, ' ')
   call add_default ('T700    '  , 1, ' ')
   call add_default ('T850    '  , 1, ' ')
   call add_default ('ATMEINT '  , 1, ' ')

   if (moist_physics) call add_default ('RELHUM  ', 1, ' ')

   call phys_getopts(history_budget_out = history_budget, history_budget_histfile_num_out = history_budget_histfile_num)
   if ( history_budget ) then
      call add_default ('PHIS    '  , history_budget_histfile_num, ' ')
      call add_default ('PS      '  , history_budget_histfile_num, ' ')
      call add_default ('T       '  , history_budget_histfile_num, ' ')
      call add_default ('U       '  , history_budget_histfile_num, ' ')
      call add_default ('V       '  , history_budget_histfile_num, ' ')
      call add_default (cnst_name(1), history_budget_histfile_num, ' ')
      ! State before physics (FV)
      call add_default ('TBP     '  , history_budget_histfile_num, ' ')
      call add_default (bpcnst(1)   , history_budget_histfile_num, ' ')
      ! State after physics (FV)
      call add_default ('TAP     '  , history_budget_histfile_num, ' ')
      call add_default ('UAP     '  , history_budget_histfile_num, ' ')
      call add_default ('VAP     '  , history_budget_histfile_num, ' ')
      call add_default (apcnst(1)   , history_budget_histfile_num, ' ')
      if ( dycore_is('LR') ) then
         call add_default ('TFIX    '    , history_budget_histfile_num, ' ')
         call add_default ('PTTEND_RESID', history_budget_histfile_num, ' ')
      end if
   end if

#ifdef AERLIFE

   call addfld ('C_DMS    ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_SO2    ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_SO4    ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_BC     ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_POM    ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_SS     ','kg m-2   ',1,    'A','Column burden',phys_decomp)
   call addfld ('C_DUST   ','kg m-2   ',1,    'A','Column burden',phys_decomp)

   call addfld ('DMSCO    ','mol/mol  ',pver,    'A','Concentration',phys_decomp)
   call addfld ('SO2CO    ','mol/mol  ',pver,    'A','Concentration',phys_decomp)
   call addfld ('SO4    ','kg m-3   ',pver,    'A','Concentration',phys_decomp)
   call addfld ('BC     ','kg m-3   ',pver,    'A','Concentration',phys_decomp)
   call addfld ('POM    ','kg m-3   ',pver,    'A','Concentration',phys_decomp)
   call addfld ('SS     ','kg m-3   ',pver,    'A','Concentration',phys_decomp)
   call addfld ('DUST   ','kg m-3   ',pver,    'A','Concentration',phys_decomp)

   call addfld ('EMI_DMS    ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_SO2    ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_SO4    ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC     ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_POM    ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_SS     ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_DUST   ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)

#ifdef CMIP6
   call addfld ('EMI_SO2_VOLC','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_SO2_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_SO4_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_AIR ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_N   ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_NI  ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_AX  ','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_N_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_NI_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_BC_AX_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
   call addfld ('EMI_POM_CMIP6','kg m-2 s-1  ',1,    'A','Emissions',phys_decomp)
#endif

!   call addfld ('WET_MSA    ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_SO2    ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_SO4    ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_BC     ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_POM    ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_SS     ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)
   call addfld ('WET_DUST   ','kg m-2 s-1   ',1,    'A','Wet deposition',phys_decomp)

!   call addfld ('DRY_MSA    ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_SO2    ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_SO4    ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_BC     ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_POM    ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_SS     ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   call addfld ('DRY_DUST   ','kg m-2 s-1   ',1,    'A','Dry deposition',phys_decomp)
   
   call addfld ('KHET    ','s-1     ',pver,    'A','Aqueous-phase reaction rate',phys_decomp)
   call addfld ('CH2O2   ','kg/kg   ',pver,    'A','H2O2 MMR from file ',phys_decomp)
   call addfld ('NUSO4N  ','kg/kg s-1',pver,    'A','SO4 nucleation rate' ,phys_decomp)

#ifdef DIRIND
#ifdef CMIP6
   call addfld ('AODVVOLC','unitless',1,    'A','CMIP6 volcanic aerosol optical depth at 0.442-0.625um',phys_decomp) ! CAM4-Oslo: 0.35-0.64um
   call addfld ('ABSVVOLC','unitless',1,    'A','CMIP6 volcanic aerosol absorptive optical depth at 0.442-0.625um',phys_decomp) ! CAM4-Oslo: 0.35-0.64um
#endif 
   call addfld ('TAUE550 ','unitless',1,    'A','Aerosol optical depth at 550nm',phys_decomp)  
   call addfld ('TAUA550 ','unitless',1,    'A','Aerosol absorptive optical depth at 550nm',phys_decomp)
   call addfld ('CDODVIS ','unitless',1,    'A','Clear air aerosol optical depth',phys_decomp)  
   call addfld ('CABSVIS ','unitless',1,    'A','Clear air aerosol absorptive optical depth',phys_decomp)
   call addfld ('CLDFREE ','unitless',1,    'A','Cloud free fraction wrt CDODVIS and CABSVIS',phys_decomp)
   call addfld ('DAYFOC  ','unitless',1,    'A','Daylight fraction',phys_decomp)
   call addfld ('N_AER   ','unitless',pver, 'A','Aerosol number concentration',phys_decomp)
   call addfld ('N_AERORG','unitless',pver, 'A','Aerosol number concentration',phys_decomp)
   call addfld ('WAK     ','unitless',pver, 'A','Aerosol single scattering albedo',phys_decomp)    
   call addfld ('GAK     ','unitless',pver, 'A','Aerosol assymetry factor',phys_decomp)    
   call addfld ('BAK     ','1/km    ',pver, 'A','Aerosol extinction',phys_decomp)     
   call addfld ('AKCXS   ','mg/m2   ',1,    'A','Scheme excess aerosol mass burden',phys_decomp)     
   call addfld ('CDNC    ','1/CM3   ',pver, 'A','Cloud Droplet Number Concentration',phys_decomp)
   call addfld ('CLDFOC  ','FRACTION',pver, 'A','Frequency of Warm Cloud Occurence',phys_decomp)
   call addfld ('REFFL   ','uM      ',pver, 'A','Effective Radius of Cloud Droplets',phys_decomp)
   call addfld ('REHANA  ','uM      ',1,    'A','Effective radius as seen from satellite',phys_decomp)
   call addfld ('FOCHANA ','FRACTION',1,    'A','Frequency of Occurrence of Clouds with REHANA /= 0',phys_decomp)
   call addfld ('RELH    ','unitless',pver, 'A','Fictive relative humidity',phys_decomp)
!#define DIAGNCDNC
#ifndef DIAGNCDNC
   call addfld ('SUPERSAT','percent ',pver, 'A','Supersaturation from AGR schem',phys_decomp)
   call addfld ('NUCRAT  ','cm-3s-1 ',pver, 'A','Nucleation of cloud droplets',phys_decomp)
   call addfld ('SELFX   ','cm-3s-1 ',pver, 'A','Selfcollection of cloud droplets',phys_decomp)
   call addfld ('EVAPX   ','cm-3s-1 ',pver, 'A','Evaporation of cloud droplets ',phys_decomp)
   call addfld ('FREEZ   ','cm-3s-1 ',pver, 'A',' rate of freezing of cloud droplets',phys_decomp)
   call addfld ('NRAINX  ','cm-3s-1 ',pver, 'A','Loss of cloud droplets due to precip '   ,phys_decomp)
   call addfld ('CDNCINT ','1/CM2   ',1,    'A','Cloud Droplet Number Column',phys_decomp)
#endif
!#define AEROCOM
!#define AEROFFL 
#ifdef AEROFFL
   call addfld ('FSNT_DRF','W/m^2   ',1,    'A','Total column absorbed solar flux (DIRind)' ,phys_decomp)
   call addfld ('FSNTCDRF','W/m^2   ',1,    'A','Clear sky total column absorbed solar flux (DIRind)' ,phys_decomp)
   call addfld ('FSNS_DRF','W/m^2   ',1,    'A','Surface absorbed solar flux (DIRind)' ,phys_decomp)
   call addfld ('FSNSCDRF','W/m^2   ',1,    'A','Clear sky surface absorbed solar flux (DIRind)' ,phys_decomp)
   call addfld ('QRS_DRF ','K/s     ',pver, 'A','Solar heating rate (DIRind)' ,phys_decomp)
   call addfld ('QRSC_DRF','K/s     ',pver, 'A','Clearsky solar heating rate (DIRind)' ,phys_decomp)
   call addfld ('FSNT_AIE','W/m^2   ',1,    'A','Total column absorbed solar flux (dirIND)' ,phys_decomp)
   call addfld ('FLNT_AIE','W/m^2   ',1,    'A','Total column longwave flux (dirIND)' ,phys_decomp)
   call addfld ('FSNTCAIE','W/m^2   ',1,    'A','Clear sky total column absorbed solar flux (dirIND)' ,phys_decomp)
   call addfld ('FSNS_AIE','W/m^2   ',1,    'A','Surface absorbed solar flux (dirIND)' ,phys_decomp)
   call addfld ('FSNSCAIE','W/m^2   ',1,    'A','Clear sky surface absorbed solar flux (dirIND)' ,phys_decomp)
   call addfld ('QRS_AIE ','K/s     ',pver, 'A','Solar heating rate (dirIND)' ,phys_decomp)
   call addfld ('QRSC_AIE','K/s     ',pver, 'A','Clearsky solar heating rate (dirIND)' ,phys_decomp)
#ifdef AEROCOM 
   call addfld ('FSUTADRF','W/m^2   ',1,    'A','SW upwelling flux at TOA',phys_decomp)
   call addfld ('FSDS_DRF','W/m^2   ',1,    'A','SW downelling flux at surface',phys_decomp)
   call addfld ('FSUS_DRF','W/m^2   ',1,    'A','SW upwelling flux at surface',phys_decomp)
   call addfld ('FSDSCDRF','W/m^2   ',1,    'A','SW downwelling clear sky flux at surface',phys_decomp)
   call addfld ('FLUS    ','W/m^2   ',1,    'A','LW surface upwelling flux',phys_decomp)
#endif  ! aerocom
#endif  ! aeroffl
#ifdef AEROCOM 
      call addfld ('GRIDAREA','m2      ',1,    'A','Grid area for 1.9x2.5 horizontal resolution',phys_decomp)
      call addfld ('DAERH2O ','mg/m2   ',1,    'A','Aerosol water load',phys_decomp)
      call addfld ('MMR_AH2O','kg/kg   ',pver, 'A','Aerosol water mmr',phys_decomp)
      call addfld ('ECDRYAER','m-1     ',pver, 'A','Dry aerosol extinction at 550nm',phys_decomp)  
      call addfld ('ABSDRYAE','m-1     ',pver, 'A','Dry aerosol absorption at 550nm',phys_decomp)  
      call addfld ('OD550DRY','unitless',1,    'A','Dry aerosol optical depth at 550nm',phys_decomp)  
      call addfld ('AB550DRY','unitless',1,    'A','Dry aerosol absorptive optical depth at 550nm',phys_decomp)   
      call addfld ('DERLT05 ','um      ',1,    'A','Effective aerosol dry radius<0.5um',phys_decomp)
      call addfld ('DERGT05 ','um      ',1,    'A','Effective aerosol dry radius>0.5um',phys_decomp) 
      call addfld ('DER     ','um      ',1,    'A','Effective aerosol dry radius',phys_decomp) 
      call addfld ('DOD440  ','unitless',1,    'A','Aerosol optical depth at 440nm',phys_decomp)  
!      call addfld ('ABS440  ','unitless',1,    'A','Aerosol absorptive optical depth at 440nm',phys_decomp)   
      call addfld ('DOD500  ','unitless',1,    'A','Aerosol optical depth at 500nm',phys_decomp)  
!      call addfld ('ABS500  ','unitless',1,    'A','Aerosol absorptive optical depth at 500nm',phys_decomp)   
      call addfld ('DOD550  ','unitless',1,    'A','Aerosol optical depth at 550nm',phys_decomp)  
      call addfld ('ABS550  ','unitless',1,    'A','Aerosol absorptive optical depth at 550nm',phys_decomp)   
      call addfld ('ABS550AL','unitless',1,    'A','Alt. aerosol absorptive optical depth at 550nm',phys_decomp)   
      call addfld ('DOD670  ','unitless',1,    'A','Aerosol optical depth at 670nm',phys_decomp)  
!      call addfld ('ABS670  ','unitless',1,    'A','Aerosol absorptive optical depth at 870nm',phys_decomp)   
      call addfld ('DOD870  ','unitless',1,    'A','Aerosol optical depth at 670nm',phys_decomp)  
!      call addfld ('ABS870  ','unitless',1,    'A','Aerosol absorptive optical depth at 870nm',phys_decomp)   
      call addfld ('DLOAD_MI','mg/m2   ',1,    'A','mineral aerosol load',phys_decomp)     
      call addfld ('DLOAD_SS','mg/m2   ',1,    'A','sea-salt aerosol load',phys_decomp)     
      call addfld ('EC550AER','m-1     ',pver, 'A','aerosol extinction coefficient',phys_decomp)     
      call addfld ('ABS550_A','m-1     ',pver, 'A','aerosol absorption coefficient',phys_decomp)     
      call addfld ('BS550AER','m-1 sr-1',pver, 'A','aerosol backscatter coefficient',phys_decomp)     
      call addfld ('CDOD550  ','unitless',1,    'A','Clear air Aerosol optical depth at 550nm',phys_decomp)  
      call addfld ('CABS550  ','unitless',1,    'A','Clear air Aerosol abs optical depth at 550nm',phys_decomp)  
      call addfld ('CABS550A ','unitless',1,    'A','Clear air Aerosol abs optical depth at 550nm',phys_decomp)  
      call addfld ('CDOD870 ','unitless',1,    'A','Clear air Aerosol optical depth at 870nm',phys_decomp)  
      call addfld ('A550_DU ','unitless',1,    'A','mineral abs. aerosol optical depth 550nm',phys_decomp)     
      call addfld ('A550_SS ','unitless',1,    'A','sea-salt abs aerosol optical depth 550nm',phys_decomp)     
      call addfld ('A550_SO4','unitless',1,    'A','SO4 aerosol abs. optical depth 550nm',phys_decomp)     
      call addfld ('A550_POM','unitless',1,    'A','OC abs. aerosol optical depth 550nm',phys_decomp)     
      call addfld ('A550_BC ','unitless',1,    'A','BC abs. aerosol optical depth 550nm',phys_decomp)
      call addfld ('ABS5503D','unitless',pver,    'A','aerosol 3d abs. optical depth 550nm',phys_decomp)
!-      call addfld ('D443_DU ','unitless',pver,    'A','mineral aerosol 3d optical depth 440nm',phys_decomp)     
!-      call addfld ('D443_SS ','unitless',pver,    'A','sea-salt aerosol 3d optical depth 440nm',phys_decomp)     
!-      call addfld ('D443_SO4','unitless',pver,    'A','SO4 aerosol 3d optical depth 440nm',phys_decomp)     
!-      call addfld ('D443_POM','unitless',pver,    'A','OC aerosol 3d optical depth 440nm',phys_decomp)     
!-      call addfld ('D443_BC ','unitless',pver,    'A','BC aerosol 3d optical depth 440nm',phys_decomp)
!-      call addfld ('D503_DU ','unitless',pver,    'A','mineral aerosol 3d optical depth 500nm',phys_decomp)     
!-      call addfld ('D503_SS ','unitless',pver,    'A','sea-salt aerosol 3d optical depth 500nm',phys_decomp)     
!-      call addfld ('D503_SO4','unitless',pver,    'A','SO4 aerosol 3d optical depth 500nm',phys_decomp)     
!-      call addfld ('D503_POM','unitless',pver,    'A','OC aerosol 3d optical depth 500nm',phys_decomp)     
!-      call addfld ('D503_BC ','unitless',pver,    'A','BC aerosol 3d optical depth 500nm',phys_decomp)
      call addfld ('D553_DU ','unitless',pver,    'A','mineral aerosol 3d optical depth 550nm',phys_decomp)     
      call addfld ('D553_SS ','unitless',pver,    'A','sea-salt aerosol 3d optical depth 550nm',phys_decomp)     
      call addfld ('D553_SO4','unitless',pver,    'A','SO4 aerosol 3d optical depth 550nm',phys_decomp)     
      call addfld ('D553_POM','unitless',pver,    'A','OC aerosol 3d optical depth 550nm',phys_decomp)     
      call addfld ('D553_BC ','unitless',pver,    'A','BC aerosol 3d optical depth 550nm',phys_decomp)
!-      call addfld ('D673_DU ','unitless',pver,    'A','mineral aerosol 3d optical depth 670nm',phys_decomp)     
!-      call addfld ('D673_SS ','unitless',pver,    'A','sea-salt aerosol 3d optical depth 670nm',phys_decomp)     
!-      call addfld ('D673_SO4','unitless',pver,    'A','SO4 aerosol 3d optical depth 670nm',phys_decomp)     
!-      call addfld ('D673_POM','unitless',pver,    'A','OC aerosol 3d optical depth 670nm',phys_decomp)     
!-      call addfld ('D673_BC ','unitless',pver,    'A','BC aerosol 3d optical depth 670nm',phys_decomp)
!-      call addfld ('D873_DU ','unitless',pver,    'A','mineral aerosol 3d optical depth 870nm',phys_decomp)     
!-      call addfld ('D873_SS ','unitless',pver,    'A','sea-salt aerosol 3d optical depth 870nm',phys_decomp)     
!-      call addfld ('D873_SO4','unitless',pver,    'A','SO4 aerosol 3d optical depth 870nm',phys_decomp)     
!-      call addfld ('D873_POM','unitless',pver,    'A','OC aerosol 3d optical depth 870nm',phys_decomp)     
!-      call addfld ('D873_BC ','unitless',pver,    'A','BC aerosol 3d optical depth 870nm',phys_decomp)
      call addfld ('D440_DU ','unitless',1,    'A','mineral aerosol optical depth 440nm',phys_decomp)     
      call addfld ('D440_SS ','unitless',1,    'A','sea-salt aerosol optical depth 440nm',phys_decomp)     
      call addfld ('D440_SO4','unitless',1,    'A','SO4 aerosol optical depth 440nm',phys_decomp)     
      call addfld ('D440_POM','unitless',1,    'A','OC aerosol optical depth 440nm',phys_decomp)     
      call addfld ('D440_BC ','unitless',1,    'A','BC aerosol optical depth 440nm',phys_decomp)     
      call addfld ('D500_DU ','unitless',1,    'A','mineral aerosol optical depth 500nm',phys_decomp)     
      call addfld ('D500_SS ','unitless',1,    'A','sea-salt aerosol optical depth 500nm',phys_decomp)     
      call addfld ('D500_SO4','unitless',1,    'A','SO4 aerosol optical depth 500nm',phys_decomp)     
      call addfld ('D500_POM','unitless',1,    'A','OC aerosol optical depth 500nm',phys_decomp)     
      call addfld ('D500_BC ','unitless',1,    'A','BC aerosol optical depth 500nm',phys_decomp)     
      call addfld ('D550_DU ','unitless',1,    'A','mineral aerosol optical depth 550nm',phys_decomp)     
      call addfld ('D550_SS ','unitless',1,    'A','sea-salt aerosol optical depth 550nm',phys_decomp)     
      call addfld ('D550_SO4','unitless',1,    'A','SO4 aerosol optical depth 550nm',phys_decomp)     
      call addfld ('D550_POM','unitless',1,    'A','OC aerosol optical depth 550nm',phys_decomp)     
      call addfld ('D550_BC ','unitless',1,    'A','BC aerosol optical depth 550nm',phys_decomp)     
      call addfld ('D670_DU ','unitless',1,    'A','mineral aerosol optical depth 670nm',phys_decomp)     
      call addfld ('D670_SS ','unitless',1,    'A','sea-salt aerosol optical depth 670nm',phys_decomp)     
      call addfld ('D670_SO4','unitless',1,    'A','SO4 aerosol optical depth 670nm',phys_decomp)     
      call addfld ('D670_POM','unitless',1,    'A','OC aerosol optical depth 670nm',phys_decomp)     
      call addfld ('D670_BC ','unitless',1,    'A','BC aerosol optical depth 670nm',phys_decomp)     
      call addfld ('D870_DU ','unitless',1,    'A','mineral aerosol optical depth 870nm',phys_decomp)     
      call addfld ('D870_SS ','unitless',1,    'A','sea-salt aerosol optical depth 870nm',phys_decomp)     
      call addfld ('D870_SO4','unitless',1,    'A','SO4 aerosol optical depth 870nm',phys_decomp)     
      call addfld ('D870_POM','unitless',1,    'A','OC aerosol optical depth 870nm',phys_decomp)     
      call addfld ('D870_BC ','unitless',1,    'A','BC aerosol optical depth 870nm',phys_decomp)     
      call addfld ('DLT_DUST','unitless',1,    'A','mineral aerosol optical depth 550nm lt05',phys_decomp)     
      call addfld ('DLT_SS  ','unitless',1,    'A','sea-salt aerosol optical depth 550nm lt05',phys_decomp)     
      call addfld ('DLT_SO4 ','unitless',1,    'A','SO4 aerosol optical depth 550nm lt05',phys_decomp)     
      call addfld ('DLT_POM ','unitless',1,    'A','OC aerosol optical depth 550nm lt05',phys_decomp)     
      call addfld ('DLT_BC  ','unitless',1,    'A','BC aerosol optical depth 550nm lt05',phys_decomp)     
      call addfld ('DGT_DUST','unitless',1,    'A','mineral aerosol optical depth 550nm gt05',phys_decomp)     
      call addfld ('DGT_SS  ','unitless',1,    'A','sea-salt aerosol optical depth 550nm gt05',phys_decomp)     
      call addfld ('DGT_SO4 ','unitless',1,    'A','SO4 aerosol optical depth 550nm gt05',phys_decomp)     
      call addfld ('DGT_POM ','unitless',1,    'A','OC aerosol optical depth 550nm gt05',phys_decomp)     
      call addfld ('DGT_BC  ','unitless',1,    'A','BC aerosol optical depth 550nm gt05',phys_decomp)     
      call addfld ('NNAT_0  ','1/cm3   ',pver, 'A','Aerosol mode 0 number concentration',phys_decomp)     
      call addfld ('NNAT_1  ','1/cm3   ',pver, 'A','Aerosol mode 1 number concentration',phys_decomp)     
      call addfld ('NNAT_2  ','1/cm3   ',pver, 'A','Aerosol mode 2 number concentration',phys_decomp)     
      call addfld ('NNAT_4  ','1/cm3   ',pver, 'A','Aerosol mode 4 number concentration',phys_decomp)     
      call addfld ('NNAT_5  ','1/cm3   ',pver, 'A','Aerosol mode 5 number concentration',phys_decomp)     
      call addfld ('NNAT_6  ','1/cm3   ',pver, 'A','Aerosol mode 6 number concentration',phys_decomp)     
      call addfld ('NNAT_7  ','1/cm3   ',pver, 'A','Aerosol mode 7 number concentration',phys_decomp)     
      call addfld ('NNAT_8  ','1/cm3   ',pver, 'A','Aerosol mode 8 number concentration',phys_decomp)     
      call addfld ('NNAT_9  ','1/cm3   ',pver, 'A','Aerosol mode 9 number concentration',phys_decomp)     
      call addfld ('NNAT_10 ','1/cm3   ',pver, 'A','Aerosol mode 10 number concentration',phys_decomp)     
      call addfld ('NNAT_11 ','1/cm3   ',pver, 'A','Aerosol mode 11 number concentration',phys_decomp)     
      call addfld ('NNAT_12 ','1/cm3   ',pver, 'A','Aerosol mode 12 number concentration',phys_decomp)     
      call addfld ('NNAT_13 ','1/cm3   ',pver, 'A','Aerosol mode 13 number concentration',phys_decomp)     
      call addfld ('NNAT_14 ','1/cm3   ',pver, 'A','Aerosol mode 14 number concentration',phys_decomp)     
      call addfld ('AIRMASS ','kg/m3   ',pver, 'A','Layer airmass',phys_decomp)     
!-      call addfld ('C_BCPM  ','ug C/m3 ',pver, 'A',' BC concentration' ,phys_decomp)    
!-      call addfld ('C_BC05  ','ug C/m3 ',pver, 'A',' BC concentration < 0.5um' ,phys_decomp)    
!-      call addfld ('C_BC125 ','ug C/m3 ',pver, 'A',' BC concentration > 1.25um' ,phys_decomp)    
!-      call addfld ('C_OCPM  ','ug C/m3 ',pver, 'A',' OC concentration' ,phys_decomp)    
!-      call addfld ('C_OC05  ','ug C/m3 ',pver, 'A',' OC concentration < 0.5um' ,phys_decomp)    
!-      call addfld ('C_OC125 ','ug C/m3 ',pver, 'A',' OC concentration > 1.25um' ,phys_decomp)    
!-      call addfld ('C_S4PM  ','ug C/m3 ',pver, 'A',' SO4 concentration' ,phys_decomp)    
!-      call addfld ('C_S405  ','ug C/m3 ',pver, 'A',' SO4 concentration < 0.5um' ,phys_decomp)    
!-      call addfld ('C_S4125 ','ug C/m3 ',pver, 'A',' SO4 concentration > 1.25um' ,phys_decomp)    
!-      call addfld ('C_MIPM  ','ug C/m3 ',pver, 'A',' dust concentration' ,phys_decomp)    
!-      call addfld ('C_MI05  ','ug C/m3 ',pver, 'A',' dust concentration < 0.5um' ,phys_decomp)    
!-      call addfld ('C_MI125 ','ug C/m3 ',pver, 'A',' dust concentration > 1.25um' ,phys_decomp)    
!-      call addfld ('C_SSPM  ','ug C/m3 ',pver, 'A',' sea-salt concentration' ,phys_decomp)    
!-      call addfld ('C_SS05  ','ug C/m3 ',pver, 'A',' sea-salt concentration < 0.5um' ,phys_decomp)    
!-      call addfld ('C_SS125 ','ug C/m3 ',pver, 'A',' sea-salt concentration > 1.25um' ,phys_decomp)    
      call addfld ('BETOT550','unitless',pver, 'A','Aerosol 3d optical depth at 550nm',phys_decomp)  
      call addfld ('BATOT550','unitless',pver, 'A','Aerosol 3d absorptive optical depth at 550nm',phys_decomp)   
#endif  ! aerocom
#endif  ! dirind

#endif  ! aerlife
#ifdef AERLIFE

   call add_default('DMSCO    ', 1,'  ')
!   call add_default('MSA    ', 1,'  ')
   call add_default('SO2CO    ', 1,'  ')
   call add_default('SO4    ', 1,'  ')
   call add_default('BC     ', 1,'  ')
   call add_default('POM    ', 1,'  ')
   call add_default('SS     ', 1,'  ')
   call add_default('DUST   ', 1,'  ')

   call add_default('C_DMS    ', 1,'  ')
!   call add_default('C_MSA    ', 1,'  ')
   call add_default('C_SO2    ', 1,'  ')
   call add_default('C_SO4    ', 1,'  ')
   call add_default('C_BC     ', 1,'  ')
   call add_default('C_POM    ', 1,'  ')
   call add_default('C_SS     ', 1,'  ')
   call add_default('C_DUST   ', 1,'  ')

   call add_default('EMI_DMS    ', 1,'  ')
   call add_default('EMI_SO2    ', 1,'  ')
   call add_default('EMI_SO4    ', 1,'  ')
   call add_default('EMI_BC     ', 1,'  ')
   call add_default('EMI_POM    ', 1,'  ')
   call add_default('EMI_SS     ', 1,'  ')
   call add_default('EMI_DUST   ', 1,'  ')

!   call add_default('WET_MSA    ', 1,'  ')
   call add_default('WET_SO2    ', 1,'  ')
   call add_default('WET_SO4    ', 1,'  ')
   call add_default('WET_BC     ', 1,'  ')
   call add_default('WET_POM    ', 1,'  ')
   call add_default('WET_SS     ', 1,'  ')
   call add_default('WET_DUST   ', 1,'  ')

!   call add_default('DRY_MSA    ', 1,'  ')
   call add_default('DRY_SO2    ', 1,'  ')
   call add_default('DRY_SO4    ', 1,'  ')
   call add_default('DRY_BC     ', 1,'  ')
   call add_default('DRY_POM    ', 1,'  ')
   call add_default('DRY_SS     ', 1,'  ')
   call add_default('DRY_DUST   ', 1,'  ')

   call add_default ('KHET    ', 1, ' ')
   call add_default ('CH2O2   ', 1, ' ')
   call add_default ('NUSO4N  ', 1, ' ')

#ifdef DIRIND
   call add_default ('TAUE550 ', 1, ' ')
   call add_default ('TAUA550 ', 1, ' ')
#ifdef CMIP6
   call add_default ('AODVVOLC', 1, ' ')
   call add_default ('ABSVVOLC', 1, ' ')
#endif
   call add_default ('DAYFOC  ', 1, ' ')
   call add_default ('CDODVIS ', 1, ' ')
   call add_default ('CABSVIS ', 1, ' ')
   call add_default ('CLDFREE ', 1, ' ')
   call add_default ('N_AER   ', 1, ' ')
   call add_default ('N_AERORG', 1, ' ')
   call add_default ('WAK     ', 1, ' ')
   call add_default ('GAK     ', 1, ' ')
   call add_default ('BAK     ', 1, ' ')
   call add_default ('AKCXS   ', 1, ' ')
   call add_default ('CDNC    ', 1, ' ')
   call add_default ('CLDFOC  ', 1, ' ')
   call add_default ('REFFL   ', 1, ' ')
   call add_default ('REHANA  ', 1, ' ')
   call add_default ('FOCHANA ', 1, ' ')
   call add_default ('RELH    ', 1, ' ')
!#define DIAGNCDNC
#ifndef DIAGNCDNC
   call add_default ('SUPERSAT', 1, ' ')
   call add_default ('NUCRAT  ', 1, ' ')
   call add_default ('SELFX   ', 1, ' ')
   call add_default ('EVAPX   ', 1, ' ')
   call add_default ('NRAINX  ', 1, ' ')
   call add_default ('FREEZ   ', 1, ' ')
   call add_default ('CDNCINT ', 1, ' ')
#endif
#ifdef AEROFFL
   call add_default ('FSNT_DRF', 1, ' ')
   call add_default ('FSNTCDRF', 1, ' ')
   call add_default ('FSNS_DRF', 1, ' ')
   call add_default ('FSNSCDRF', 1, ' ')
   call add_default ('QRS_DRF ', 1, ' ')
   call add_default ('QRSC_DRF', 1, ' ')
   call add_default ('FSNT_AIE', 1, ' ')
   call add_default ('FLNT_AIE', 1, ' ')
   call add_default ('FSNTCAIE', 1, ' ')
   call add_default ('FSNS_AIE', 1, ' ')
   call add_default ('FSNSCAIE', 1, ' ')
   call add_default ('QRS_AIE ', 1, ' ')
   call add_default ('QRSC_AIE', 1, ' ')
#ifdef AEROCOM 
   call add_default ('FSUTADRF', 1, ' ')
   call add_default ('FSDS_DRF', 1, ' ')
   call add_default ('FSUS_DRF', 1, ' ')
   call add_default ('FSDSCDRF', 1, ' ')
   call add_default ('FLUS    ', 1, ' ')
#endif  ! aerocom
#endif  ! aeroffl
#ifdef AEROCOM 
      call add_default ('GRIDAREA', 1, ' ')
      call add_default ('DAERH2O ', 1, ' ')
      call add_default ('MMR_AH2O', 1, ' ')
      call add_default ('ECDRYAER', 1, ' ')
      call add_default ('ABSDRYAE', 1, ' ')
      call add_default ('OD550DRY', 1, ' ')
      call add_default ('AB550DRY', 1, ' ')
      call add_default ('DERLT05 ', 1, ' ')
      call add_default ('DERGT05 ', 1, ' ')
      call add_default ('DER     ', 1, ' ')
      call add_default ('DOD440  ', 1, ' ')
!      call add_default ('ABS440  ', 1, ' ')
      call add_default ('DOD500  ', 1, ' ')
!      call add_default ('ABS500  ', 1, ' ')
      call add_default ('DOD550  ', 1, ' ')
      call add_default ('ABS550  ', 1, ' ')
      call add_default ('ABS550AL', 1, ' ')
      call add_default ('DOD670  ', 1, ' ')
!      call add_default ('ABS670  ', 1, ' ')
      call add_default ('DOD870  ', 1, ' ')
!      call add_default ('ABS870  ', 1, ' ')
      call add_default ('DLOAD_MI', 1, ' ')
      call add_default ('DLOAD_SS', 1, ' ')
      call add_default ('EC550AER', 1, ' ')
      call add_default ('ABS550_A', 1, ' ')
      call add_default ('BS550AER', 1, ' ')
      call add_default ('CDOD550 ', 1, ' ')
      call add_default ('CABS550 ', 1, ' ')
      call add_default ('CABS550A', 1, ' ')
      call add_default ('CDOD870 ', 1, ' ')
      call add_default ('A550_DU ', 1, ' ')
      call add_default ('A550_SS ', 1, ' ')
      call add_default ('A550_SO4', 1, ' ')
      call add_default ('A550_POM', 1, ' ')
      call add_default ('A550_BC ', 1, ' ')
      call add_default ('ABS5503D', 1, ' ')
!-      call add_default ('D443_DU ', 1, ' ')
!-      call add_default ('D443_SS ', 1, ' ')
!-      call add_default ('D443_SO4', 1, ' ')
!-      call add_default ('D443_POM', 1, ' ')
!-      call add_default ('D443_BC ', 1, ' ')
!-      call add_default ('D503_DU ', 1, ' ')
!-      call add_default ('D503_SS ', 1, ' ')
!-      call add_default ('D503_SO4', 1, ' ')
!-      call add_default ('D503_POM', 1, ' ')
!-      call add_default ('D503_BC ', 1, ' ')
      call add_default ('D553_DU ', 1, ' ')
      call add_default ('D553_SS ', 1, ' ')
      call add_default ('D553_SO4', 1, ' ')
      call add_default ('D553_POM', 1, ' ')
      call add_default ('D553_BC ', 1, ' ')
!-      call add_default ('D673_DU ', 1, ' ')
!-      call add_default ('D673_SS ', 1, ' ')
!-      call add_default ('D673_SO4', 1, ' ')
!-      call add_default ('D673_POM', 1, ' ')
!-      call add_default ('D673_BC ', 1, ' ')
!-      call add_default ('D873_DU ', 1, ' ')
!-      call add_default ('D873_SS ', 1, ' ')
!-      call add_default ('D873_SO4', 1, ' ')
!-      call add_default ('D873_POM', 1, ' ')
!-      call add_default ('D873_BC ', 1, ' ')
      call add_default ('D440_DU ', 1, ' ')
      call add_default ('D440_SS ', 1, ' ')
      call add_default ('D440_SO4', 1, ' ')
      call add_default ('D440_POM', 1, ' ')
      call add_default ('D440_BC ', 1, ' ')
      call add_default ('D500_DU ', 1, ' ')
      call add_default ('D500_SS ', 1, ' ')
      call add_default ('D500_SO4', 1, ' ')
      call add_default ('D500_POM', 1, ' ')
      call add_default ('D500_BC ', 1, ' ')
      call add_default ('D550_DU ', 1, ' ')
      call add_default ('D550_SS ', 1, ' ')
      call add_default ('D550_SO4', 1, ' ')
      call add_default ('D550_POM', 1, ' ')
      call add_default ('D550_BC ', 1, ' ')
      call add_default ('D670_DU ', 1, ' ')
      call add_default ('D670_SS ', 1, ' ')
      call add_default ('D670_SO4', 1, ' ')
      call add_default ('D670_POM', 1, ' ')
      call add_default ('D670_BC ', 1, ' ')
      call add_default ('D870_DU ', 1, ' ')
      call add_default ('D870_SS ', 1, ' ')
      call add_default ('D870_SO4', 1, ' ')
      call add_default ('D870_POM', 1, ' ')
      call add_default ('D870_BC ', 1, ' ')
      call add_default ('DLT_DUST', 1, ' ')
      call add_default ('DLT_SS  ', 1, ' ')
      call add_default ('DLT_SO4 ', 1, ' ')
      call add_default ('DLT_POM ', 1, ' ')
      call add_default ('DLT_BC  ', 1, ' ')
      call add_default ('DGT_DUST', 1, ' ')
      call add_default ('DGT_SS  ', 1, ' ')
      call add_default ('DGT_SO4 ', 1, ' ')
      call add_default ('DGT_POM ', 1, ' ')
      call add_default ('DGT_BC  ', 1, ' ')
      call add_default ('NNAT_0  ', 1, ' ')
      call add_default ('NNAT_1  ', 1, ' ')
      call add_default ('NNAT_2  ', 1, ' ')
      call add_default ('NNAT_4  ', 1, ' ')
      call add_default ('NNAT_5  ', 1, ' ')
      call add_default ('NNAT_6  ', 1, ' ')
      call add_default ('NNAT_7  ', 1, ' ')
      call add_default ('NNAT_8  ', 1, ' ')
      call add_default ('NNAT_9  ', 1, ' ')
      call add_default ('NNAT_10 ', 1, ' ')
      call add_default ('NNAT_11 ', 1, ' ')
      call add_default ('NNAT_12 ', 1, ' ')
      call add_default ('NNAT_13 ', 1, ' ')
      call add_default ('NNAT_14 ', 1, ' ')
      call add_default ('AIRMASS ', 1, ' ')
!-      call add_default ('C_BCPM  ', 1, ' ')
!-      call add_default ('C_BC05  ', 1, ' ')
!-      call add_default ('C_BC125 ', 1, ' ')
!-      call add_default ('C_OCPM  ', 1, ' ')
!-      call add_default ('C_OC05  ', 1, ' ')
!-      call add_default ('C_OC125 ', 1, ' ')
!-      call add_default ('C_S4PM  ', 1, ' ')
!-      call add_default ('C_S405  ', 1, ' ')
!-      call add_default ('C_S4125 ', 1, ' ')
!-      call add_default ('C_MIPM  ', 1, ' ')
!-      call add_default ('C_MI05  ', 1, ' ')
!-      call add_default ('C_MI125 ', 1, ' ')
!-      call add_default ('C_SSPM  ', 1, ' ')
!-      call add_default ('C_SS05  ', 1, ' ')
!-      call add_default ('C_SS125 ', 1, ' ')
      call add_default ('BETOT550', 1, ' ')
      call add_default ('BATOT550', 1, ' ')
#endif  ! aerocom
#endif  ! dirind

#endif  ! aerlife

   ! This field is added by radiation when full physics is used
   if ( ideal_phys )then
      call addfld('QRS     ', 'K/s     ', pver, 'A', 'Solar heating rate', phys_decomp)
      call add_default('QRS     ', 1, ' ')
   end if

   call add_default ('VT      ', 1, ' ')
   call add_default ('VU      ', 1, ' ')
   call add_default ('VV      ', 1, ' ')
   call add_default ('VQ      ', 1, ' ')
   call add_default ('UU      ', 1, ' ')
   call add_default ('OMEGAT  ', 1, ' ')
   call add_default ('TMQ     ', 1, ' ')
   call add_default ('PSL     ', 1, ' ')

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Exit here for adiabatic/ideal physics cases !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (.not. moist_physics) return


   call addfld ('PDELDRY ','Pa      ',pver, 'A','Dry pressure difference between levels',phys_decomp)
   call addfld ('PSDRY   ','Pa      ',1,    'A','Surface pressure',phys_decomp)

   if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then
      call add_default ('PS      ', 2, ' ')
      call add_default ('T       ', 2, ' ')
   end if

   ! outfld calls in diag_conv

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call addfld ('DTCOND  ','K/s     ',pver, 'A','T tendency - moist processes',phys_decomp)
   call addfld ('DTCOND_24_COS','K/s',pver, 'A','T tendency - moist processes 24hr. cos coeff.',phys_decomp)
   call addfld ('DTCOND_24_SIN','K/s',pver, 'A','T tendency - moist processes 24hr. sin coeff.',phys_decomp)
   call addfld ('DTCOND_12_COS','K/s',pver, 'A','T tendency - moist processes 12hr. cos coeff.',phys_decomp)
   call addfld ('DTCOND_12_SIN','K/s',pver, 'A','T tendency - moist processes 12hr. sin coeff.',phys_decomp)
   do m = 1, pcnst
      dcconnam(m) = 'DC'//cnst_name(m)
   end do
   if (diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' .or. history_budget) then
      call addfld (dcconnam(1), 'kg/kg/s',pver,'A',trim(cnst_name(1))//' tendency due to moist processes',phys_decomp)
      if ( diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' ) then
         call add_default (dcconnam(1),                           1, ' ')
      end if
      if( history_budget ) then
         call add_default (dcconnam(1), history_budget_histfile_num, ' ')
      end if
      if (diag_cnst_conv_tend == 'all' .or. history_budget) then
         do m = 2, pcnst
            call addfld (dcconnam(m), 'kg/kg/s',pver,'A',trim(cnst_name(m))//' tendency due to moist processes',phys_decomp)
            if( diag_cnst_conv_tend == 'all' ) then
               call add_default (dcconnam(m),                           1, ' ')
            end if
            if( history_budget .and. (m == ixcldliq .or. m == ixcldice) ) then
               call add_default (dcconnam(m), history_budget_histfile_num, ' ')
            end if
         end do
      end if
   end if
   call addfld ('PRECL   ','m/s     ',1,    'A','Large-scale (stable) precipitation rate (liq + ice)'                ,phys_decomp)
   call addfld ('PRECC   ','m/s     ',1,    'A','Convective precipitation rate (liq + ice)'                          ,phys_decomp)
   call addfld ('PRECT   ','m/s     ',1,    'A','Total (convective and large-scale) precipitation rate (liq + ice)'  ,phys_decomp)
   call addfld ('PRECTMX ','m/s     ',1,    'X','Maximum (convective and large-scale) precipitation rate (liq+ice)'  ,phys_decomp)
   call addfld ('PRECSL  ','m/s     ',1,    'A','Large-scale (stable) snow rate (water equivalent)'                  ,phys_decomp)
   call addfld ('PRECSC  ','m/s     ',1,    'A','Convective snow rate (water equivalent)'                            ,phys_decomp)
   call addfld ('PRECCav ','m/s     ',1,    'A','Average large-scale precipitation (liq + ice)'                      ,phys_decomp)
   call addfld ('PRECLav ','m/s     ',1,    'A','Average convective precipitation  (liq + ice)'                      ,phys_decomp)

   ! defaults
   call add_default ('DTCOND  ', 1, ' ')
   call add_default ('PRECL   ', 1, ' ')
   call add_default ('PRECC   ', 1, ' ')
   call add_default ('PRECT   ', 1, ' ')
   call add_default ('PRECSL  ', 1, ' ')
   call add_default ('PRECSC  ', 1, ' ')

   ! outfld calls in diag_surf

   call addfld ('SHFLX   ','W/m2    ',1,    'A','Surface sensible heat flux',phys_decomp)
   call addfld ('LHFLX   ','W/m2    ',1,    'A','Surface latent heat flux',phys_decomp)
   call addfld ('QFLX    ','kg/m2/s ',1,    'A','Surface water flux',phys_decomp)

   call addfld ('TAUX    ','N/m2    ',1,    'A','Zonal surface stress',phys_decomp)
   call addfld ('TAUY    ','N/m2    ',1,    'A','Meridional surface stress',phys_decomp)
   call addfld ('TREFHT  ','K       ',1,    'A','Reference height temperature',phys_decomp)
   call addfld ('TREFHTMN','K       ',1,    'M','Minimum reference height temperature over output period',phys_decomp)
   call addfld ('TREFHTMX','K       ',1,    'X','Maximum reference height temperature over output period',phys_decomp)
   call addfld ('QREFHT  ','kg/kg   ',1,    'A','Reference height humidity',phys_decomp)
   call addfld ('U10     ','m/s     ',1,    'A','10m wind speed',phys_decomp)
   call addfld ('UAS     ','m/s     ',1,    'A','Eastward 10m wind',phys_decomp)
   call addfld ('VAS     ','m/s     ',1,    'A','Northward 10m wind',phys_decomp)
   call addfld ('RHREFHT ','fraction',1,    'A','Reference height relative humidity',phys_decomp)

   call addfld ('LANDFRAC','fraction',1,    'A','Fraction of sfc area covered by land',phys_decomp)
   call addfld ('ICEFRAC ','fraction',1,    'A','Fraction of sfc area covered by sea-ice',phys_decomp)
   call addfld ('OCNFRAC ','fraction',1,    'A','Fraction of sfc area covered by ocean',phys_decomp)

   call addfld ('TREFMNAV','K       ',1,    'A','Average of TREFHT daily minimum',phys_decomp)
   call addfld ('TREFMXAV','K       ',1,    'A','Average of TREFHT daily maximum',phys_decomp)
   call addfld ('U10MXAV ','K       ',1,    'A','Average of U10 daily maximum',phys_decomp)

   call addfld ('TS      ','K       ',1,    'A','Surface temperature (radiative)',phys_decomp)
   call addfld ('TSMN    ','K       ',1,    'M','Minimum surface temperature over output period',phys_decomp)
   call addfld ('TSMX    ','K       ',1,    'X','Maximum surface temperature over output period',phys_decomp)
   call addfld ('SNOWHLND','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
   call addfld ('SNOWHICE','m       ',1,    'A','Water equivalent snow depth',phys_decomp)
   call addfld ('TBOT    ','K       ',1,    'A','Lowest model level temperature', phys_decomp)

   call addfld ('ASDIR',   '1',       1,    'A','albedo: shortwave, direct', phys_decomp)
   call addfld ('ASDIF',   '1',       1,    'A','albedo: shortwave, diffuse', phys_decomp)
   call addfld ('ALDIR',   '1',       1,    'A','albedo: longwave, direct', phys_decomp)
   call addfld ('ALDIF',   '1',       1,    'A','albedo: longwave, diffuse', phys_decomp)
   call addfld ('SST',     'K',       1,    'A','sea surface temperature', phys_decomp)

   ! defaults
   call add_default ('ASDIR   ', 1, ' ')
   call add_default ('ASDIF   ', 1, ' ')
   call add_default ('ALDIR   ', 1, ' ')
   call add_default ('ALDIF   ', 1, ' ')
   call add_default ('SHFLX   ', 1, ' ')
   call add_default ('LHFLX   ', 1, ' ')
   call add_default ('QFLX    ', 1, ' ')
   call add_default ('TAUX    ', 1, ' ')
   call add_default ('TAUY    ', 1, ' ')
   call add_default ('TREFHT  ', 1, ' ')
   call add_default ('TREFMNAV', 1, ' ') 
   call add_default ('TREFMXAV', 1, ' ')  
   call add_default ('U10MXAV', 1, ' ')  
   call add_default ('LANDFRAC', 1, ' ')
   call add_default ('OCNFRAC ', 1, ' ')
   call add_default ('QREFHT  ', 1, ' ')
   call add_default ('U10     ', 1, ' ')
   call add_default ('UAS     ', 1, ' ')
   call add_default ('VAS     ', 1, ' ')
   call add_default ('RHREFHT ', 1, ' ')
   call add_default ('ICEFRAC ', 1, ' ')
   call add_default ('TS      ', 1, ' ')
   call add_default ('TSMN    ', 1, ' ')
   call add_default ('TSMX    ', 1, ' ')
   call add_default ('SNOWHLND', 1, ' ')
   call add_default ('SNOWHICE', 1, ' ')

   ! outfld calls in diag_phys_tend_writeout

   call addfld ('PTTEND  '   ,'K/s     ',pver, 'A','T total physics tendency'                             ,phys_decomp)
   call addfld (ptendnam(       1),  'kg/kg/s ',pver, 'A',trim(cnst_name(       1))//' total physics tendency '      ,phys_decomp)
   call addfld (ptendnam(ixcldliq),  'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldliq))//' total physics tendency '      ,phys_decomp)
   call addfld (ptendnam(ixcldice),  'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldice))//' total physics tendency '      ,phys_decomp)
   if ( dycore_is('LR') )then
      call addfld (dmetendnam(       1),'kg/kg/s ',pver, 'A',trim(cnst_name(       1))//' dme adjustment tendency (FV) ',phys_decomp)
      call addfld (dmetendnam(ixcldliq),'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldliq))//' dme adjustment tendency (FV) ',phys_decomp)
      call addfld (dmetendnam(ixcldice),'kg/kg/s ',pver, 'A',trim(cnst_name(ixcldice))//' dme adjustment tendency (FV) ',phys_decomp)
   end if

   if ( history_budget ) then
      call add_default ('PTTEND'          , history_budget_histfile_num, ' ')
      call add_default (ptendnam(       1), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldliq), history_budget_histfile_num, ' ')
      call add_default (ptendnam(ixcldice), history_budget_histfile_num, ' ')
      if ( dycore_is('LR') )then
         call add_default(dmetendnam(1)       , history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldliq), history_budget_histfile_num, ' ')
         call add_default(dmetendnam(ixcldice), history_budget_histfile_num, ' ')
      end if
      if( history_budget_histfile_num > 1 ) then
         call add_default ('DTCOND  '         , history_budget_histfile_num, ' ')
      end if
   end if

   ! outfld calls in diag_physvar_ic

   call addfld ('QCWAT&IC   ','kg/kg   ',pver, 'I','q associated with cloud water'                   ,phys_decomp)
   call addfld ('TCWAT&IC   ','kg/kg   ',pver, 'I','T associated with cloud water'                   ,phys_decomp)
   call addfld ('LCWAT&IC   ','kg/kg   ',pver, 'I','Cloud water (ice + liq'                          ,phys_decomp)
   call addfld ('CLOUD&IC   ','fraction',pver, 'I','Cloud fraction'                                  ,phys_decomp)
   call addfld ('CONCLD&IC   ','fraction',pver, 'I','Convective cloud fraction'                      ,phys_decomp)
   call addfld ('TKE&IC     ','m2/s2   ',pverp,'I','Turbulent Kinetic Energy'                        ,phys_decomp)
   call addfld ('CUSH&IC    ','m       ',1,    'I','Convective Scale Height'                         ,phys_decomp)
   call addfld ('KVH&IC     ','m2/s    ',pverp,'I','Vertical diffusion diffusivities (heat/moisture)',phys_decomp)
   call addfld ('KVM&IC     ','m2/s    ',pverp,'I','Vertical diffusion diffusivities (momentum)'     ,phys_decomp)
   call addfld ('PBLH&IC    ','m       ',1,    'I','PBL height'                                      ,phys_decomp)
   call addfld ('TPERT&IC   ','K       ',1,    'I','Perturbation temperature (eddies in PBL)'        ,phys_decomp)
   call addfld ('QPERT&IC   ','kg/kg   ',1,    'I','Perturbation specific humidity (eddies in PBL)'  ,phys_decomp)
   call addfld ('TBOT&IC    ','K       ',1,    'I','Lowest model level temperature'                  ,phys_decomp)

   ! Initial file - Optional fields

   if (inithist_all) then
      call add_default ('CONCLD&IC  ',0, 'I')
      call add_default ('QCWAT&IC   ',0, 'I')
      call add_default ('TCWAT&IC   ',0, 'I')
      call add_default ('LCWAT&IC   ',0, 'I')
      call add_default ('PBLH&IC    ',0, 'I')
      call add_default ('TPERT&IC   ',0, 'I')
      call add_default ('QPERT&IC   ',0, 'I')
      call add_default ('CLOUD&IC   ',0, 'I')
      call add_default ('TKE&IC     ',0, 'I')
      call add_default ('CUSH&IC    ',0, 'I')
      call add_default ('KVH&IC     ',0, 'I')
      call add_default ('KVM&IC     ',0, 'I')
      call add_default ('TBOT&IC    ',0, 'I')
   end if

   ! CAM export state 
   call addfld('a2x_BCPHIWET', 'kg/m2/s', 1, 'A', 'wetdep of hydrophilic black carbon',   phys_decomp)
   call addfld('a2x_BCPHIDRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophilic black carbon',   phys_decomp)
   call addfld('a2x_BCPHODRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophobic black carbon',   phys_decomp)
   call addfld('a2x_OCPHIWET', 'kg/m2/s', 1, 'A', 'wetdep of hydrophilic organic carbon', phys_decomp)
   call addfld('a2x_OCPHIDRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophilic organic carbon', phys_decomp)
   call addfld('a2x_OCPHODRY', 'kg/m2/s', 1, 'A', 'drydep of hydrophobic organic carbon', phys_decomp)
   call addfld('a2x_DSTWET1',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin1)',                phys_decomp)
   call addfld('a2x_DSTDRY1',  'kg/m2/s', 1, 'A', 'drydep of dust (bin1)',                phys_decomp)
   call addfld('a2x_DSTWET2',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin2)',                phys_decomp)
   call addfld('a2x_DSTDRY2',  'kg/m2/s', 1, 'A', 'drydep of dust (bin2)',                phys_decomp)
   call addfld('a2x_DSTWET3',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin3)',                phys_decomp)
   call addfld('a2x_DSTDRY3',  'kg/m2/s', 1, 'A', 'drydep of dust (bin3)',                phys_decomp)
   call addfld('a2x_DSTWET4',  'kg/m2/s', 1, 'A', 'wetdep of dust (bin4)',                phys_decomp)
   call addfld('a2x_DSTDRY4',  'kg/m2/s', 1, 'A', 'drydep of dust (bin4)',                phys_decomp)

   !---------------------------------------------------------
   ! CAM history fields for CAM-DOM/CAM-CSIM 
   !---------------------------------------------------------

   ! CAM-DOM history fields
#ifdef COUP_DOM
   call addfld ('TSOCN&IC   ','m       ',1,    'I','Ocean tempertare',phys_decomp)
   call add_default ('TSOCN&IC   ',0, 'I')
#endif

  ! CAM-CSIM history fields

  do k=1,plevmx
     call addfld (tsnam(k),'K       ',1,'A',tsnam(k)//' subsoil temperature',phys_decomp)
  end do
  call addfld ('SICTHK  '   ,'m       ',1,'A','Sea ice thickness',phys_decomp)
  call addfld ('TSICE   '   ,'K       ',1,'A','Ice temperature',phys_decomp)
  do k = 1,plevmx
     call addfld (trim(tsnam(k))//'&IC','K       ',1,'I',tsnam(k)//' subsoil temperature',phys_decomp)
  end do
  call addfld ('SICTHK&IC  ','m       ',1,'I','Sea ice thickness'                      ,phys_decomp)
  call addfld ('TSICE&IC   ','K       ',1,'I','Ice temperature'                        ,phys_decomp)
  call addfld ('SNOWHICE&IC','m       ',1,'I','Water equivalent snow depth'            ,phys_decomp)
  call addfld ('ICEFRAC&IC ','fraction',1,'I','Fraction of sfc area covered by sea-ice',phys_decomp)
  call addfld ('TSICERAD&IC','K       ',1,'I','Radiatively equivalent ice temperature' ,phys_decomp)
  do k = 1,plevmx
     call add_default(trim(tsnam(k))//'&IC',0, 'I')
  end do
  call add_default ('SICTHK&IC  ',0, 'I')
  call add_default ('TSICE&IC   ',0, 'I')
  call add_default ('SNOWHICE&IC',0, 'I')
  call add_default ('ICEFRAC&IC ',0, 'I')
  if (inithist_all) then
     call add_default ('TSICERAD&IC',0, 'I')
  end if

  !---------------------------------------------------------
  ! WACCM diagnostic history fields 
  !---------------------------------------------------------

  if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then

    ! create history variables for fourier coefficients of the diurnal 
    ! and semidiurnal tide in T, U, V, and Z3

    call tidal_diag_init() 

  endif

      qcwat_idx  = pbuf_get_fld_idx('QCWAT')
      tcwat_idx  = pbuf_get_fld_idx('TCWAT')
      lcwat_idx  = pbuf_get_fld_idx('LCWAT')
      cld_idx    = pbuf_get_fld_idx('CLD')
      concld_idx = pbuf_get_fld_idx('CONCLD')
      tke_idx    = pbuf_get_fld_idx('tke')
      kvm_idx    = pbuf_get_fld_idx('kvm')
      kvh_idx    = pbuf_get_fld_idx('kvh')
      cush_idx   = pbuf_get_fld_idx('cush')
  

end subroutine diag_init

!===============================================================================

subroutine diag_allocate()

! Allocate memory for module variables.
! Done at the begining of a physics step at same point as the pbuf allocate for
! variables with "physpkg" scope.

! Local variables
   character(len=*), parameter :: sub = 'diag_allocate'
   integer :: istat

   allocate(dtcond(pcols,pver,begchunk:endchunk),        &
            dqcond(pcols,pver,pcnst,begchunk:endchunk), &
            stat=istat)
   if ( istat /= 0 ) then
      call endrun (sub//': ERROR: allocate failed')
   end if
   dtcond = nan
   dqcond = nan

end subroutine diag_allocate

!===============================================================================

subroutine diag_deallocate()

! Deallocate memory for module variables.
! Done at the end of a physics step at same point as the pbuf deallocate for
! variables with "physpkg" scope.

! Local variables
   character(len=*), parameter :: sub = 'diag_deallocate'
   integer :: istat

   deallocate(dtcond, dqcond, stat=istat)
   if ( istat /= 0 ) then
      call endrun (sub//': ERROR: deallocate failed')
   end if
end subroutine diag_deallocate
!===============================================================================

subroutine diag_conv_tend_ini(state,pbuf)

! Initialize convective tendency calcs.

   use phys_buffer,       only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx

! Argument:

   type(physics_state), intent(in) :: state
   type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf

! Local variables:

   integer :: i, k, m, lchnk, ncol
   real(r8), pointer, dimension(:,:) :: t_ttend

   lchnk = state%lchnk
   ncol  = state%ncol

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = state%s(i,k)
      end do
   end do

   do m = 1, pcnst
      do k = 1, pver
         do i = 1, ncol
            dqcond(i,k,m,lchnk) = state%q(i,k,m)
         end do
      end do
   end do

   !! initialize to pbuf T_TTEND to temperature at first timestep
   t_ttend => pbuf(t_ttend_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1) 
   t_ttend(:ncol,:pver) = state%t(:ncol,:pver)

end subroutine diag_conv_tend_ini
!===============================================================================

  subroutine diag_phys_writeout(state, psl)

!----------------------------------------------------------------------- 
! 
! Purpose: record dynamics variables on physics grid
!
!-----------------------------------------------------------------------
    use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
    use time_manager,       only: get_nstep
    use interpolate_data,   only: vertinterp
    use constituent_burden, only: constituent_burden_comp
    use cam_control_mod,    only: moist_physics
    use co2_cycle,          only: c_i, co2_transport

    use tidal_diag,         only: tidal_diag_write
!-----------------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(inout) :: state
   real(r8), optional , intent(out)   :: psl(pcols) 
!
!---------------------------Local workspace-----------------------------
!
    real(r8) ftem(pcols,pver) ! temporary workspace
    real(r8) ftem1(pcols,pver) ! another temporary workspace
    real(r8) ftem2(pcols,pver) ! another temporary workspace
    real(r8) psl_tmp(pcols)   ! Sea Level Pressure
    real(r8) z3(pcols,pver)   ! geo-potential height
    real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t1(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_t2(pcols)    ! data interpolated to a pressure surface
    real(r8) p_surf_q1(pcols)    ! data interpolated to a pressure surface    
    real(r8) p_surf_q2(pcols)    ! data interpolated to a pressure surface        
    real(r8) tem2(pcols,pver) ! temporary workspace
    real(r8) timestep(pcols)  ! used for outfld call
    real(r8) esl(pcols,pver)   ! saturation vapor pressures 
    real(r8) esi(pcols,pver)   ! 
    real(r8) dlon(pcols)      ! width of grid cell (meters)
    integer  plon             ! number of longitudes

    integer i, k, m, lchnk, ncol, nstep
!
!-----------------------------------------------------------------------
!
    lchnk = state%lchnk
    ncol  = state%ncol

    ! Output NSTEP for debugging
    nstep = get_nstep()
    timestep(:ncol) = nstep
    call outfld ('NSTEP   ',timestep, pcols, lchnk)

    call outfld('T       ',state%t , pcols   ,lchnk   )
    call outfld('PS      ',state%ps, pcols   ,lchnk   )
    call outfld('U       ',state%u , pcols   ,lchnk   )
    call outfld('V       ',state%v , pcols   ,lchnk   )
    do m=1,pcnst
       if ( cnst_cam_outfld(m) ) then
          call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
       end if
    end do

    if (co2_transport()) then
       do m = 1,4
          call outfld(trim(cnst_name(c_i(m)))//'_BOT', state%q(1,pver,c_i(m)), pcols, lchnk)
       end do
    end if

    ! column burdens of all constituents except water vapor
    call constituent_burden_comp(state)

    if ( moist_physics) then
       call outfld('PDELDRY ',state%pdeldry, pcols, lchnk)
       call outfld('PSDRY',   state%psdry,   pcols, lchnk) 
    end if

    call outfld('PHIS    ',state%phis,    pcols,   lchnk     )



#if (defined BFB_CAM_SCAM_IOP )
    call outfld('phis    ',state%phis,    pcols,   lchnk     )
#endif

!
! Add height of surface to midpoint height above surface 
!
    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do
    call outfld('Z3      ',z3,pcols,lchnk)
!           
! Output Z3 on pressure surfaces
!
    if (hist_fld_active('Z700')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
       call outfld('Z700    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
       call outfld('Z500    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, z3, p_surf)
       call outfld('Z250    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
       call outfld('Z300    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, z3, p_surf)
       call outfld('Z200    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z100')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 10000._r8, z3, p_surf)
       call outfld('Z100    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z050')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  5000._r8, z3, p_surf)
       call outfld('Z050    ', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('Z010')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  1000._r8, z3, p_surf)
       call outfld('Z010    ', p_surf, pcols, lchnk)
    end if
!
! Quadratic height fiels Z3*Z3
!
    ftem(:ncol,:) = z3(:ncol,:)*z3(:ncol,:)
    call outfld('ZZ      ',ftem,pcols,lchnk)

    ftem(:ncol,:) = z3(:ncol,:)*state%v(:ncol,:)*gravit
    call outfld('VZ      ',ftem,  pcols,lchnk)
!
! Meridional advection fields
!
    ftem(:ncol,:) = state%v(:ncol,:)*state%t(:ncol,:)
    call outfld ('VT      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)
    call outfld ('VQ      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%q(:ncol,:,1)*state%q(:ncol,:,1)
    call outfld ('QQ      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)**2
    call outfld ('VV      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:) * state%u(:ncol,:)
    call outfld ('VU      ',ftem    ,pcols   ,lchnk     )

! zonal advection

    ftem(:ncol,:) = state%u(:ncol,:)**2
    call outfld ('UU      ',ftem    ,pcols   ,lchnk     )

! Wind speed
    ftem(:ncol,:) = sqrt( state%u(:ncol,:)**2 + state%v(:ncol,:)**2)
    call outfld ('WSPEED  ',ftem    ,pcols   ,lchnk     )

! Vertical velocity and advection

    if (single_column) then
       call outfld('OMEGA   ',wfld,    pcols,   lchnk     )
    else
       call outfld('OMEGA   ',state%omega,    pcols,   lchnk     )
    endif

#if (defined BFB_CAM_SCAM_IOP )
    call outfld('omega   ',state%omega,    pcols,   lchnk     )
#endif

    ftem(:ncol,:) = state%omega(:ncol,:)*state%t(:ncol,:)
    call outfld('OMEGAT  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%u(:ncol,:)
    call outfld('OMEGAU  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%v(:ncol,:)
    call outfld('OMEGAV  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%q(:ncol,:,1)
    call outfld('OMEGAQ  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%omega(:ncol,:)
    call outfld('OMGAOMGA',ftem,    pcols,   lchnk     )
!
! Output omega at 850 and 500 mb pressure levels
!
    if (hist_fld_active('OMEGA850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
       call outfld('OMEGA850', p_surf, pcols, lchnk)
    end if
    if (hist_fld_active('OMEGA500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%omega, p_surf)
       call outfld('OMEGA500', p_surf, pcols, lchnk)
    end if
!     
! Mass of q, by layer and vertically integrated
!
    ftem(:ncol,:) = state%q(:ncol,:,1) * state%pdel(:ncol,:) * rga
    call outfld ('MQ      ',ftem    ,pcols   ,lchnk     )

    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TMQ     ',ftem, pcols   ,lchnk     )

    if (moist_physics) then

       ! Relative humidity
       if (hist_fld_active('RELHUM')) then
          call aqsat (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
             ncol ,pver  ,1       ,pver    )
          ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
          call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
       end if

       if (hist_fld_active('RHW') .or. hist_fld_active('RHI') .or. hist_fld_active('RHCFMIP') ) then
	  
          ! RH w.r.t liquid (water)
          call aqsat_water (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
             ncol ,pver  ,1       ,pver    )
          ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
          call outfld ('RHW  ',ftem    ,pcols   ,lchnk     )

          ! Convert to RHI (ice)
          do i=1,ncol
             do k=1,pver
                esl(i,k)=polysvp(state%t(i,k),0)
                esi(i,k)=polysvp(state%t(i,k),1)
                ftem1(i,k)=ftem(i,k)*esl(i,k)/esi(i,k)
             end do
          end do
          call outfld ('RHI  ',ftem1    ,pcols   ,lchnk     )

	  ! use temperature to decide if you populate with ftem (liquid, above 0 C) or ftem1 (ice, below 0 C)

	  ftem2(:ncol,:)=ftem(:ncol,:)

          do i=1,ncol
             do k=1,pver
		if (state%t(i,k) .gt. 273) then
                   ftem2(i,k)=ftem(i,k)  !!wrt water
 		else
                   ftem2(i,k)=ftem1(i,k) !!wrt ice
		end if
             end do
          end do
          
          call outfld ('RHCFMIP  ',ftem2    ,pcols   ,lchnk     )

       end if

    end if
!
! Sea level pressure
!
    if (present(psl) .or. hist_fld_active('PSL')) then
       call cpslec (ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair) 
       call outfld ('PSL     ',psl_tmp  ,pcols, lchnk     )
       if (present(psl)) then	
          psl(:ncol) = psl_tmp(:ncol)
       end if
    end if
!
! Output T,q,u,v fields on pressure surfaces
!
    if (hist_fld_active('T850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
       call outfld('T850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%t, p_surf)
       call outfld('T500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T300')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
       call outfld('T300    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%t, p_surf)
       call outfld('T200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('T010')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  1000._r8, state%t, p_surf)
       call outfld('T010    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
       call outfld('Q850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('Q200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
       call outfld('Q200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
       call outfld('U850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%u, p_surf)
       call outfld('U500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%u, p_surf)
       call outfld('U250    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
       call outfld('U200    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('U010')) then
       call vertinterp(ncol, pcols, pver, state%pmid,  1000._r8, state%u, p_surf)
       call outfld('U010    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V850')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
       call outfld('V850    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V500')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%v, p_surf)
       call outfld('V500    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V250')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 25000._r8, state%v, p_surf)
       call outfld('V250    ', p_surf, pcols, lchnk )
    end if
    if (hist_fld_active('V200')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
       call outfld('V200    ', p_surf, pcols, lchnk )
    end if

    ftem(:ncol,:) = state%t(:ncol,:)*state%t(:ncol,:)
    call outfld('TT      ',ftem    ,pcols   ,lchnk   )
!
! Output U, V, T, Q, P and Z at bottom level
!
    call outfld ('UBOT    ', state%u(1,pver)  ,  pcols, lchnk)
    call outfld ('VBOT    ', state%v(1,pver)  ,  pcols, lchnk)
    call outfld ('QBOT    ', state%q(1,pver,1),  pcols, lchnk)
    call outfld ('ZBOT    ', state%zm(1,pver) , pcols, lchnk)

! Total energy of the atmospheric column for atmospheric heat storage calculations

    !! temporary variable to get surface geopotential in dimensions of (ncol,pver)
    do k=1,pver
      ftem1(:ncol,k)=state%phis(:ncol)  !! surface geopotential in units (m2/s2)
    end do

    !! calculate sum of sensible, kinetic, latent, and surface geopotential energy
    !! E=CpT+PHIS+Lv*q+(0.5)*(u^2+v^2)
    ftem(:ncol,:) = (cpair*state%t(:ncol,:) +  ftem1(:ncol,:) + latvap*state%q(:ncol,:,1) + &
         0.5_r8*(state%u(:ncol,:)**2+state%v(:ncol,:)**2))*(state%pdel(:ncol,:)/gravit)
    !! vertically integrate
    do k=2,pver       
	ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('ATMEINT   ',ftem(:ncol,1)  ,pcols   ,lchnk     )

!! Boundary layer atmospheric stability, temperature, water vapor diagnostics

    if (hist_fld_active('T1000')      .or. &
        hist_fld_active('T9251000')   .or. & 
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('T8501000')   .or. &
        hist_fld_active('TH8501000')  .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('T7001000')   .or. &
        hist_fld_active('TH7001000')  .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%t, p_surf_t1)
    end if

    if (hist_fld_active('T925')       .or. &
        hist_fld_active('T9251000')   .or. & 
        hist_fld_active('TH9251000')  .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%t, p_surf_t2)
    end if

    if (hist_fld_active('Q1000')      .or. &
        hist_fld_active('THE9251000') .or. &
        hist_fld_active('THE8501000') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 100000._r8, state%q(1,1,1), p_surf_q1)
    end if

    if (hist_fld_active('Q925')       .or. &
        hist_fld_active('THE9251000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 92500._r8, state%q(1,1,1), p_surf_q2)
    end if

    !!! at 1000 mb and 925 mb
    if (hist_fld_active('T1000')) then
       call outfld('T1000    ', p_surf_t1, pcols, lchnk )
    end if

    if (hist_fld_active('T925')) then
       call outfld('T925    ', p_surf_t2, pcols, lchnk )
    end if

    if (hist_fld_active('Q1000')) then
       call outfld('Q1000    ', p_surf_q1, pcols, lchnk ) 
    end if

    if (hist_fld_active('Q925')) then
       call outfld('Q925    ', p_surf_q2, pcols, lchnk )  
    end if

    if (hist_fld_active('T9251000')) then
       p_surf = p_surf_t2-p_surf_t1  
       call outfld('T9251000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('TH9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)   
       call outfld('TH9251000    ', p_surf, pcols, lchnk )    
    end if

    if (hist_fld_active('THE9251000')) then
       p_surf = (p_surf_t2*(1000.0_r8/925.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))-(p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))  
       call outfld('THE9251000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('T8501000')  .or. &
        hist_fld_active('TH8501000') .or. &
        hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf_t2)
    end if

    !!! at 1000 mb and 850 mb
    if (hist_fld_active('T8501000')) then
       p_surf = p_surf_t2-p_surf_t1  
       call outfld('T8501000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('TH8501000')) then
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)   
       call outfld('TH8501000    ', p_surf, pcols, lchnk )   
    end if

    if (hist_fld_active('THE8501000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/850.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))-(p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))  
       call outfld('THE8501000    ', p_surf, pcols, lchnk ) 
    end if

    if (hist_fld_active('T7001000')  .or. &
        hist_fld_active('TH7001000') .or. &
        hist_fld_active('T700') .or. &
        hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%t, p_surf_t2)
    end if

   !!! at 700 mb
    if (hist_fld_active('T700')) then
       call outfld('T700    ', p_surf_t2, pcols, lchnk )
    end if

    !!! at 1000 mb and 700 mb
    if (hist_fld_active('T7001000')) then
       p_surf = p_surf_t2-p_surf_t1
       call outfld('T7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('TH7001000')) then
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)-(p_surf_t1*(1.0_r8)**cappa)
       call outfld('TH7001000    ', p_surf, pcols, lchnk )
    end if

    if (hist_fld_active('THE7001000')) then
       call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, state%q(1,1,1), p_surf_q2)
       p_surf = (p_surf_t2*(1000.0_r8/700.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q2)/(1004.0_r8*p_surf_t2))-(p_surf_t1*(1.0_r8)**cappa)*exp((2500000.0_r8*p_surf_q1)/(1004.0_r8*p_surf_t1))
       call outfld('THE7001000    ', p_surf, pcols, lchnk )
    end if

  !---------------------------------------------------------
  ! WACCM tidal diagnostics 
  !---------------------------------------------------------

  if (chem_is('waccm_ghg') .or. chem_is('waccm_mozart')) then

    call tidal_diag_write(state) 

  endif

    return
  end subroutine diag_phys_writeout
!===============================================================================

subroutine diag_conv(state, ztodt,   &
     prec_zmc, snow_zmc, prec_cmf, snow_cmf, prec_sed, snow_sed, prec_pcw, snow_pcw)

!----------------------------------------------------------------------- 
! 
! Output diagnostics associated with all convective processes.
!
!-----------------------------------------------------------------------
   use physconst,     only: cpair
   use tidal_diag,    only: get_tidal_coeffs

! Arguments:

   real(r8),            intent(in) :: ztodt   ! timestep for computing physics tendencies
   type(physics_state), intent(in) :: state
! convective precipitation variables
   real(r8), intent(in) :: prec_zmc(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: snow_zmc(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: prec_cmf(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: snow_cmf(pcols)                ! snow from   Hack   convection
   real(r8), intent(in) :: prec_sed(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: snow_sed(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: prec_pcw(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: snow_pcw(pcols)                ! snow from Hack   convection

! Local variables:
   
   integer :: i, k, m, lchnk, ncol
   real(r8) :: rtdt

   real(r8):: precc(pcols)                ! convective precip rate
   real(r8):: precl(pcols)                ! stratiform precip rate
   real(r8):: snowc(pcols)                ! convective snow rate
   real(r8):: snowl(pcols)                ! stratiform snow rate
   real(r8):: prect(pcols)                ! total (conv+large scale) precip rate
   real(r8) :: dcoef(4)                   ! for tidal component of T tend

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1._r8/ztodt

! Precipitation rates (multi-process)
   precc(:ncol) = prec_zmc(:ncol) + prec_cmf(:ncol)
   precl(:ncol) = prec_sed(:ncol) + prec_pcw(:ncol)
   snowc(:ncol) = snow_zmc(:ncol) + snow_cmf(:ncol)
   snowl(:ncol) = snow_sed(:ncol) + snow_pcw(:ncol)
   prect(:ncol) = precc(:ncol) + precl(:ncol)

   call outfld('PRECC   ', precc, pcols, lchnk )
   call outfld('PRECL   ', precl, pcols, lchnk )
   call outfld('PRECSC  ', snowc, pcols, lchnk )
   call outfld('PRECSL  ', snowl, pcols, lchnk )
   call outfld('PRECT   ', prect, pcols, lchnk )
   call outfld('PRECTMX ', prect, pcols, lchnk )

   call outfld('PRECLav ', precl, pcols, lchnk )
   call outfld('PRECCav ', precc, pcols, lchnk )

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ' , prect, pcols, lchnk )
#endif

! Total convection tendencies.

   do k = 1, pver
      do i = 1, ncol
         dtcond(i,k,lchnk) = (state%s(i,k) - dtcond(i,k,lchnk))*rtdt / cpair
      end do
   end do
   call outfld('DTCOND  ', dtcond(:,:,lchnk), pcols, lchnk)

!  output tidal coefficients
   call get_tidal_coeffs( dcoef )
   call outfld( 'DTCOND_24_SIN', dtcond(:ncol,:,lchnk)*dcoef(1), ncol, lchnk )
   call outfld( 'DTCOND_24_COS', dtcond(:ncol,:,lchnk)*dcoef(2), ncol, lchnk )
   call outfld( 'DTCOND_12_SIN', dtcond(:ncol,:,lchnk)*dcoef(3), ncol, lchnk )
   call outfld( 'DTCOND_12_COS', dtcond(:ncol,:,lchnk)*dcoef(4), ncol, lchnk )

   do m = 1, pcnst
      do k = 1, pver
         do i = 1, ncol
            dqcond(i,k,m,lchnk) = (state%q(i,k,m) - dqcond(i,k,m,lchnk))*rtdt
         end do
      end do
   end do

   if (diag_cnst_conv_tend == 'q_only' .or. diag_cnst_conv_tend == 'all' .or. history_budget) then
      if ( cnst_cam_outfld(1) ) then
         call outfld(dcconnam(1), dqcond(:,:,1,lchnk), pcols, lchnk)
      end if
      if (diag_cnst_conv_tend == 'all' .or. history_budget) then
         do m = 2, pcnst
            if ( cnst_cam_outfld(m) ) then
               call outfld(dcconnam(m), dqcond(:,:,m,lchnk), pcols, lchnk)
            end if
         end do
      end if
   end if

end subroutine diag_conv

!===============================================================================

subroutine diag_surf (cam_in, cam_out, ps, trefmxav, trefmnav, u10mxav)

!----------------------------------------------------------------------- 
! 
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use time_manager,     only: is_end_curr_day
   use co2_cycle,        only: c_i, co2_transport
   use constituents,     only: sflxnam

!-----------------------------------------------------------------------
!
! Input arguments
!
    type(cam_in_t),  intent(in) :: cam_in
    type(cam_out_t), intent(in) :: cam_out

    real(r8), intent(inout) :: trefmnav(pcols) ! daily minimum tref  
    real(r8), intent(inout) :: trefmxav(pcols) ! daily maximum tref
    real(r8), intent(inout) :: u10mxav(pcols) ! daily maximum 10m wind speed

    real(r8), intent(in)    :: ps(pcols)       ! Surface pressure.
!
!---------------------------Local workspace-----------------------------
!
    integer :: i, k, m      ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
    real(r8) tem2(pcols)    ! temporary workspace
    real(r8) ftem(pcols)    ! temporary workspace
!
!-----------------------------------------------------------------------
!
    lchnk = cam_in%lchnk
    ncol  = cam_in%ncol

    call outfld('SHFLX',    cam_in%shf,       pcols, lchnk)
    call outfld('LHFLX',    cam_in%lhf,       pcols, lchnk)
    call outfld('QFLX',     cam_in%cflx(1,1), pcols, lchnk)

    call outfld('TAUX',     cam_in%wsx,       pcols, lchnk)
    call outfld('TAUY',     cam_in%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMX', cam_in%tref,      pcols, lchnk)
    call outfld('TREFHTMN', cam_in%tref,      pcols, lchnk)
    call outfld('QREFHT',   cam_in%qref,      pcols, lchnk)
    call outfld('U10',      cam_in%u10,       pcols, lchnk)
    call outfld('UAS',      cam_in%uas,       pcols, lchnk)
    call outfld('VAS',      cam_in%vas,       pcols, lchnk)
! 
! Calculate and output reference height RH (RHREFHT)

   call aqsat (cam_in%tref   ,ps  ,tem2    ,ftem    ,pcols   , &
          ncol ,1  ,1       ,1    )
       ftem(:ncol) = cam_in%qref(:ncol)/ftem(:ncol)*100._r8

      
    call outfld('RHREFHT',   ftem,      pcols, lchnk)


#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',cam_in%shf,   pcols,   lchnk)
    call outfld('lhflx   ',cam_in%lhf,   pcols,   lchnk)
    call outfld('trefht  ',cam_in%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', cam_in%landfrac, pcols, lchnk)
    call outfld('ICEFRAC',  cam_in%icefrac,  pcols, lchnk)
    call outfld('OCNFRAC',  cam_in%ocnfrac,  pcols, lchnk)
!
! Compute daily minimum and maximum of TREF
!
    do i = 1,ncol
       trefmxav(i) = max(cam_in%tref(i),trefmxav(i))
       trefmnav(i) = min(cam_in%tref(i),trefmnav(i))
       u10mxav(i) = max(cam_in%u10(i),u10mxav(i))
    end do
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmxav,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmnav,pcols,   lchnk     )
       call outfld('U10MXAV', u10mxav,pcols,   lchnk     )
       trefmxav(:ncol) = -1.0e36_r8
       trefmnav(:ncol) =  1.0e36_r8
       u10mxav(:ncol) = -1.0e36_r8
    endif

    call outfld('TBOT',     cam_out%tbot,     pcols, lchnk)
    call outfld('TS',       cam_in%ts,        pcols, lchnk)
    call outfld('TSMN',     cam_in%ts,        pcols, lchnk)
    call outfld('TSMX',     cam_in%ts,        pcols, lchnk)
    call outfld('SNOWHLND', cam_in%snowhland, pcols, lchnk)
    call outfld('SNOWHICE', cam_in%snowhice,  pcols, lchnk)
    call outfld('ASDIR',    cam_in%asdir,     pcols, lchnk)
    call outfld('ASDIF',    cam_in%asdif,     pcols, lchnk)
    call outfld('ALDIR',    cam_in%aldir,     pcols, lchnk)
    call outfld('ALDIF',    cam_in%aldif,     pcols, lchnk)
    call outfld('SST',      cam_in%sst,       pcols, lchnk)

    if (co2_transport()) then
       do m = 1,4
          call outfld(sflxnam(c_i(m)), cam_in%cflx(:,c_i(m)), pcols, lchnk)
       end do
    end if

end subroutine diag_surf

!===============================================================================

subroutine diag_export(cam_out)

!----------------------------------------------------------------------- 
! 
! Purpose: Write export state to history file
!
!-----------------------------------------------------------------------

   ! arguments
   type(cam_out_t), intent(inout) :: cam_out

   ! Local variables:
   integer :: lchnk        ! chunk identifier
   logical :: atm_dep_flux ! true ==> sending deposition fluxes to coupler.
                           ! Otherwise, set them to zero.
   !-----------------------------------------------------------------------

   lchnk = cam_out%lchnk

   call phys_getopts(atm_dep_flux_out=atm_dep_flux)

   if (.not. atm_dep_flux) then
      ! set the fluxes to zero before outfld and sending them to the
      ! coupler
      cam_out%bcphiwet = 0.0_r8
      cam_out%bcphidry = 0.0_r8
      cam_out%bcphodry = 0.0_r8
      cam_out%ocphiwet = 0.0_r8
      cam_out%ocphidry = 0.0_r8
      cam_out%ocphodry = 0.0_r8
      cam_out%dstwet1  = 0.0_r8
      cam_out%dstdry1  = 0.0_r8
      cam_out%dstwet2  = 0.0_r8
      cam_out%dstdry2  = 0.0_r8
      cam_out%dstwet3  = 0.0_r8
      cam_out%dstdry3  = 0.0_r8
      cam_out%dstwet4  = 0.0_r8
      cam_out%dstdry4  = 0.0_r8
   end if

   call outfld('a2x_BCPHIWET', cam_out%bcphiwet, pcols, lchnk)
   call outfld('a2x_BCPHIDRY', cam_out%bcphidry, pcols, lchnk)
   call outfld('a2x_BCPHODRY', cam_out%bcphodry, pcols, lchnk)
   call outfld('a2x_OCPHIWET', cam_out%ocphiwet, pcols, lchnk)
   call outfld('a2x_OCPHIDRY', cam_out%ocphidry, pcols, lchnk)
   call outfld('a2x_OCPHODRY', cam_out%ocphodry, pcols, lchnk)
   call outfld('a2x_DSTWET1',  cam_out%dstwet1,  pcols, lchnk)
   call outfld('a2x_DSTDRY1',  cam_out%dstdry1,  pcols, lchnk)
   call outfld('a2x_DSTWET2',  cam_out%dstwet2,  pcols, lchnk)
   call outfld('a2x_DSTDRY2',  cam_out%dstdry2,  pcols, lchnk)
   call outfld('a2x_DSTWET3',  cam_out%dstwet3,  pcols, lchnk)
   call outfld('a2x_DSTDRY3',  cam_out%dstdry3,  pcols, lchnk)
   call outfld('a2x_DSTWET4',  cam_out%dstwet4,  pcols, lchnk)
   call outfld('a2x_DSTDRY4',  cam_out%dstdry4,  pcols, lchnk)

end subroutine diag_export

!#######################################################################

   subroutine diag_physvar_ic (lchnk, pbuf, cam_out, cam_in)
!
!---------------------------------------------
!
! Purpose: record physics variables on IC file
!
!---------------------------------------------
!
   use phys_buffer,       only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx 
   use buffer     ,       only: pblht, tpert, qpert
!
! Arguments
!
   integer       , intent(in) :: lchnk  ! chunk identifier
   type(pbuf_fld), intent(in), dimension(pbuf_size_max) :: pbuf
   type(cam_out_t), intent(inout) :: cam_out
   type(cam_in_t),  intent(inout) :: cam_in 
!
!---------------------------Local workspace-----------------------------
!
   integer  :: k                 ! indices
   integer  :: itim              ! indices
   real(r8), pointer, dimension(:,:) :: cwat_var
   real(r8), pointer, dimension(:,:) :: conv_var_3d
   real(r8), pointer, dimension(:  ) :: conv_var_2d
!
!-----------------------------------------------------------------------
!
   if( write_inithist() ) then

      !following line added temporarily in cam3_5_45 as workaround to 
      !bluevista compiler problems
      nullify(cwat_var, conv_var_3d, conv_var_2d)

      !
      ! Associate pointers with physics buffer fields
      !
      itim = pbuf_old_tim_idx()

      cwat_var => pbuf(qcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

      !following if block added temporarily in cam3_5_45 as workaround to 
      !bluevista compiler problems
      if(.not. associated(cwat_var)) then
         call endrun('Could not associate to pbuf fld_ptr')
      end if

      call outfld('QCWAT&IC   ',cwat_var, pcols,lchnk)

      cwat_var => pbuf(tcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('TCWAT&IC   ',cwat_var, pcols,lchnk)

      cwat_var => pbuf(lcwat_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('LCWAT&IC   ',cwat_var, pcols,lchnk)

      cwat_var => pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('CLOUD&IC   ',cwat_var, pcols,lchnk)

      cwat_var => pbuf(concld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
      call outfld('CONCLD&IC   ',cwat_var, pcols,lchnk)

      conv_var_3d => pbuf(tke_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,itim)
      call outfld('TKE&IC    ',conv_var_3d, pcols,lchnk)

      conv_var_3d => pbuf(kvm_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,itim)
      call outfld('KVM&IC    ',conv_var_3d, pcols,lchnk)

      conv_var_3d => pbuf(kvh_idx)%fld_ptr(1,1:pcols,1:pverp,lchnk,itim)
      call outfld('KVH&IC    ',conv_var_3d, pcols,lchnk)

      conv_var_2d => pbuf(cush_idx)%fld_ptr(1,1:pcols,1,lchnk,itim)
      call outfld('CUSH&IC   ',conv_var_2d, pcols,lchnk)

      call outfld('PBLH&IC    ', pblht(1,  lchnk), pcols, lchnk)
      call outfld('TPERT&IC   ', tpert(1,  lchnk), pcols, lchnk)
      call outfld('QPERT&IC   ', qpert(1,1,lchnk), pcols, lchnk)

      ! The following is only needed for cam-csim
      call outfld('TBOT&IC    ', cam_out%tbot, pcols, lchnk)
   end if

   end subroutine diag_physvar_ic


!#######################################################################

   subroutine diag_phys_tend_writeout (state, pbuf, tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
                                       tmp_t, qini, cldliqini, cldiceini)
!
!---------------------------------------------------------------
!
! Purpose:  Dump physics tendencies for moisture and temperature
!
!---------------------------------------------------------------
!
   use check_energy,    only: check_energy_get_integrals
   use physconst,       only: cpair
   use phys_buffer,     only: pbuf_fld,pbuf_size_max,pbuf_get_fld_idx

!
! Arguments
!
   type(physics_state), intent(in   ) :: state 
   type(pbuf_fld),      intent(in), dimension(pbuf_size_max) :: pbuf
   type(physics_tend ), intent(in   ) :: tend
   real(r8)           , intent(in   ) :: ztodt                  ! physics timestep
   real(r8)           , intent(inout) :: tmp_q     (pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_cldliq(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_cldice(pcols,pver) ! As input, holds pre-adjusted tracers (FV)
   real(r8)           , intent(inout) :: tmp_t     (pcols,pver) ! holds last physics_updated T (FV)
   real(r8)           , intent(in   ) :: qini      (pcols,pver) ! tracer fields at beginning of physics
   real(r8)           , intent(in   ) :: cldliqini (pcols,pver) ! tracer fields at beginning of physics
   real(r8)           , intent(in   ) :: cldiceini (pcols,pver) ! tracer fields at beginning of physics

!
!---------------------------Local workspace-----------------------------
!
   integer  :: m      ! constituent index
   integer  :: lchnk  ! chunk index
   integer  :: ncol   ! number of columns in chunk
   real(r8) :: ftem2(pcols     ) ! Temporary workspace for outfld variables
   real(r8) :: ftem3(pcols,pver) ! Temporary workspace for outfld variables
   real(r8) :: rtdt
   real(r8) :: heat_glob         ! global energy integral (FV only)
   integer  :: ixcldice, ixcldliq! constituent indices for cloud liquid and ice water.
   ! CAM pointers to get variables from the physics buffer
   real(r8), pointer, dimension(:,:) :: t_ttend  
   integer  :: ifld

!
!-----------------------------------------------------------------------
!
   lchnk = state%lchnk
   ncol  = state%ncol
   rtdt  = 1._r8/ztodt
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
!
! Dump out post-physics state (FV only)
!
   if (dycore_is('LR')) then
      tmp_t(:ncol,:pver) = (tmp_t(:ncol,:pver) - state%t(:ncol,:pver))/ztodt
      call outfld('PTTEND_RESID', tmp_t, pcols, lchnk   )
   end if
   call outfld('TAP', state%t, pcols, lchnk   )
   call outfld('UAP', state%u, pcols, lchnk   )
   call outfld('VAP', state%v, pcols, lchnk   )

   if ( cnst_cam_outfld(       1) ) call outfld (apcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (apcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (apcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)
!
! T-tendency due to FV Energy fixer (remove from total physics tendency diagnostic)
!
   if (dycore_is('LR')) then
      call check_energy_get_integrals( heat_glob_out=heat_glob )
      ftem2(:ncol)  = heat_glob/cpair
      call outfld('TFIX', ftem2, pcols, lchnk   )
      ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver) - heat_glob/cpair
   else
      ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver)
   end if
!
! Total physics tendency for Temperature
!
   call outfld('PTTEND',ftem3, pcols, lchnk )
!
! Tendency for dry mass adjustment of q (valid for FV only)
!
   if (dycore_is('LR')) then
      tmp_q     (:ncol,:pver) = (state%q(:ncol,:pver,       1) - tmp_q     (:ncol,:pver))*rtdt
      tmp_cldliq(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - tmp_cldliq(:ncol,:pver))*rtdt
      tmp_cldice(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - tmp_cldice(:ncol,:pver))*rtdt
      if ( cnst_cam_outfld(       1) ) call outfld (dmetendnam(       1), tmp_q     , pcols, lchnk)
      if ( cnst_cam_outfld(ixcldliq) ) call outfld (dmetendnam(ixcldliq), tmp_cldliq, pcols, lchnk)
      if ( cnst_cam_outfld(ixcldice) ) call outfld (dmetendnam(ixcldice), tmp_cldice, pcols, lchnk)
   end if
!
! Total physics tendency for moisture and other tracers
!
   if ( cnst_cam_outfld(       1) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,       1) - qini     (:ncol,:pver) )*rtdt
      call outfld (ptendnam(       1), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldliq) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldliq) - cldliqini(:ncol,:pver) )*rtdt
      call outfld (ptendnam(ixcldliq), ftem3, pcols, lchnk)
   end if
   if ( cnst_cam_outfld(ixcldice) ) then
      ftem3(:ncol,:pver) = (state%q(:ncol,:pver,ixcldice) - cldiceini(:ncol,:pver) )*rtdt
      call outfld (ptendnam(ixcldice), ftem3, pcols, lchnk)
   end if

!
! Total (physics+dynamics, everything!) tendency for Temperature
!
   !! get temperature stored in physics buffer
   ifld = pbuf_get_fld_idx('T_TTEND')
   t_ttend  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)

   !! calculate and outfld the total temperature tendency
   ftem3(:ncol,:pver) = (t_ttend(:ncol,:pver) - state%t(:ncol,:pver))/ztodt
   call outfld('TTEND_TOT',ftem3, pcols, lchnk )

   !! update physics buffer with this time-step's temperature
   t_ttend(:ncol,:pver) = state%t(:ncol,:pver)

   end subroutine diag_phys_tend_writeout

!#######################################################################

   subroutine diag_state_b4_phys_write (state)
!
!---------------------------------------------------------------
!
! Purpose:  Dump state just prior to executing physics
!
!---------------------------------------------------------------
!
! Arguments
!
   type(physics_state), intent(in) :: state 
!
!---------------------------Local workspace-----------------------------
!
   integer :: ixcldice, ixcldliq ! constituent indices for cloud liquid and ice water.
   integer :: lchnk              ! chunk index
!
!-----------------------------------------------------------------------
!
   lchnk = state%lchnk

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   call outfld('TBP', state%t, pcols, lchnk   )
   if ( cnst_cam_outfld(       1) ) call outfld (bpcnst(       1), state%q(1,1,       1), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldliq) ) call outfld (bpcnst(ixcldliq), state%q(1,1,ixcldliq), pcols, lchnk)
   if ( cnst_cam_outfld(ixcldice) ) call outfld (bpcnst(ixcldice), state%q(1,1,ixcldice), pcols, lchnk)

   end subroutine diag_state_b4_phys_write

end module cam_diagnostics
