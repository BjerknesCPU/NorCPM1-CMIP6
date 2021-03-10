!================================================================================================
! This is the 'none' chemistry module.
! Most of the routines return without doing anything.
!================================================================================================

module chemistry
  use shr_kind_mod,        only: r8 => shr_kind_r8,shr_kind_cl
  use physics_types,       only: physics_state, physics_ptend
  use ppgrid,              only: begchunk, endchunk, pcols
  use phys_buffer,         only: pbuf_size_max, pbuf_fld
  use spmd_utils,          only: masterproc
#ifdef CMIP6
  use constituents,        only: pcnst
  use mo_extfrc,           only: extfrc_inti
  use mo_srf_emissions,    only: srf_emissions_inti
#endif

  implicit none
  private
  save
  !
  ! Public interfaces
  !
  public chem_is                        ! identify which chemistry is being used
  public chem_register                  ! register consituents
  public chem_is_active                 ! returns true if this package is active (ghg_chem=.true.)
  public chem_implements_cnst           ! returns true if consituent is implemented by this package
  public chem_init_cnst                 ! initialize mixing ratios if not read from initial file
  public chem_init                      ! initialize (history) variables
  public chem_timestep_init             ! time interpolate chemical loss frequencies
  public chem_timestep_tend             ! interface to tendency computation
  public chem_final
  public chem_write_restart
  public chem_read_restart
  public chem_init_restart
  public chem_readnl                    ! read chem namelist 

  interface chem_write_restart
     module procedure chem_write_restart_bin
     module procedure chem_write_restart_pio
  end interface
  interface chem_read_restart
     module procedure chem_read_restart_bin
     module procedure chem_read_restart_pio
  end interface


#ifdef WACCM_GHG  
  integer, save, public :: imozart=-1

  character(len=shr_kind_cl) :: solar_parms_file = & 
   '/work/shared/noresm/inputdata/atm/waccm/phot/wa_smax_c100517.nc' 
#endif 
#ifdef CMIP6
  character(len=shr_kind_cl) :: ext_frc_specifier(pcnst) = ''
  character(len=24)  :: ext_frc_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer            :: ext_frc_cycle_yr  = 0
  integer            :: ext_frc_fixed_ymd = 0
  integer            :: ext_frc_fixed_tod = 0

  character(len=shr_kind_cl) :: srf_emis_specifier(pcnst) = ''
  character(len=24)  :: srf_emis_type = 'CYCLICAL' ! 'CYCLICAL' | 'SERIAL' | 'INTERP_MISSING_MONTHS'
  integer            :: srf_emis_cycle_yr  = 0
  integer            :: srf_emis_fixed_ymd = 0
  integer            :: srf_emis_fixed_tod = 0
#endif
  ! Private data

!================================================================================================
contains
!================================================================================================

  logical function chem_is (name)

    character(len=*), intent(in) :: name

    chem_is = .false.
    if (name == 'none' ) then
       chem_is = .true.
    end if

  end function chem_is

!================================================================================================

  subroutine chem_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: register advected constituents for parameterized greenhouse gas chemistry
    ! 
    !-----------------------------------------------------------------------

  end subroutine chem_register

!================================================================================================

  subroutine chem_readnl(nlfile)

    use abortutils,      only: endrun
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
  
    ! args
 
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    integer :: unitn, ierr
#ifdef CMIP6
    namelist /chem_inparm/ & 
      & ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod, & 
      & srf_emis_specifier, srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod
#elif WACCM_GHG
    namelist /chem_inparm/ solar_parms_file 
#endif
#if defined(CMIP6) || defined(WACCM_GHG)
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'chem_inparm', status=ierr)
       if (ierr == 0) then
          read(unitn, chem_inparm, iostat=ierr)
          if (ierr /= 0) then
             call endrun('chem_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if
#endif 

#ifdef WACCM_GHG
    call mpibcast (solar_parms_file,  len(solar_parms_file),           mpichar, 0, mpicom)
#endif 

#ifdef CMIP6
    call mpibcast (ext_frc_specifier, len(ext_frc_specifier(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (ext_frc_type,      len(ext_frc_type),               mpichar, 0, mpicom)
    call mpibcast (ext_frc_cycle_yr,  1,                               mpiint, 0, mpicom)
    call mpibcast (ext_frc_fixed_ymd, 1,                               mpiint, 0, mpicom)
    call mpibcast (ext_frc_fixed_tod, 1,                               mpiint, 0, mpicom)

    call mpibcast (srf_emis_specifier, len(srf_emis_specifier(1))*pcnst, mpichar, 0, mpicom)
    call mpibcast (srf_emis_type,      len(srf_emis_type),               mpichar, 0, mpicom)
    call mpibcast (srf_emis_cycle_yr,  1,                               mpiint, 0, mpicom)
    call mpibcast (srf_emis_fixed_ymd, 1,                               mpiint, 0, mpicom)
    call mpibcast (srf_emis_fixed_tod, 1,                               mpiint, 0, mpicom)  
#endif
    WRITE(*,*) 'DEBUG: leaving chem_readnl'

  end subroutine chem_readnl

!================================================================================================

  function chem_is_active()
    !-----------------------------------------------------------------------
    logical :: chem_is_active
    !-----------------------------------------------------------------------
    chem_is_active = .false.
  end function chem_is_active

!================================================================================================

  function chem_implements_cnst(name)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: return true if specified constituent is implemented by this package
    ! 
    ! Author: B. Eaton
    ! 
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------Arguments---------------------------------

    character(len=*), intent(in) :: name   ! constituent name
    logical :: chem_implements_cnst        ! return value

    chem_implements_cnst = .false.

  end function chem_implements_cnst

!===============================================================================

  subroutine chem_init(phys_state)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: initialize parameterized greenhouse gas chemistry
    !          (declare history variables)
    ! 
    !-----------------------------------------------------------------------
    use cam_history,    only: addfld, add_default, phys_decomp
    use mo_solar_parms,    only: solar_parms_init
#ifdef CMIP6
    use mo_tracname, only: solsym
#endif 
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)
#ifdef WACCM_GHG 
    call solar_parms_init('/work/shared/noresm/inputdata/atm/waccm/phot/wa_smax_c100517.nc')
#endif 
#ifdef CMIP6
    write(*,*) 'DEBUG chem_init: call extfrc',ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod 
    solsym=(/ 'SO2     ','SO4_PR  ','BC_N    ','BC_AX   ', 'BC_NI   ','OM_NI   ' /)
    call extfrc_inti(ext_frc_specifier, ext_frc_type, ext_frc_cycle_yr, ext_frc_fixed_ymd, ext_frc_fixed_tod ) 
    call srf_emissions_inti ( srf_emis_specifier, srf_emis_type, srf_emis_cycle_yr, srf_emis_fixed_ymd, srf_emis_fixed_tod)
#endif 
  end subroutine chem_init

!===============================================================================

  subroutine chem_timestep_init(phys_state)
#ifdef CMIP6
    use mo_extfrc,    only: extfrc_timestep_init
    use mo_srf_emissions,  only : set_srf_emissions_time
#endif 
    use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
         is_perpetual
    type(physics_state), intent(in):: phys_state(begchunk:endchunk)                 
#ifdef CMIP6
    call set_srf_emissions_time( phys_state )
    call extfrc_timestep_init( phys_state )
#endif 

  end subroutine chem_timestep_init

!===============================================================================

  subroutine chem_timestep_tend( state, ptend, cam_in, cam_out, dt, pbuf, fh2o, fsds, pblh )
    use cam_history,      only: outfld
    use camsrfexch_types, only: cam_in_t, cam_out_t
    !-----------------------------------------------------------------------
    !
    ! Arguments:
    !
    real(r8),            intent(in)    :: dt              ! time step
    type(physics_state), intent(in)    :: state           ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend           ! indivdual parameterization tendencies
    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    real(r8),            intent(out)   :: fh2o(pcols)     ! h2o flux to balance source from chemistry
    type(pbuf_fld),      intent(in)    :: pbuf(pbuf_size_max)
    real(r8),            intent(in)    :: pblh(pcols)     ! pbl height (m)
    real(r8),            intent(in)    :: fsds(pcols)     ! longwave down at sfc

    return
  end subroutine chem_timestep_tend

!===============================================================================

  subroutine chem_init_cnst(name, q, gcid)

    character(len=*), intent(in) :: name         ! constituent name
    real(r8), intent(out) :: q(:,:)   !  mass mixing ratio (gcol, plev)
    integer, intent(in) :: gcid(:)    !  global column id

    return
  end subroutine chem_init_cnst

!===============================================================================
  subroutine chem_final
    return
  end subroutine chem_final
!===============================================================================
  subroutine chem_write_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_write_restart_bin
!===============================================================================
  subroutine chem_read_restart_bin( nrg )
    implicit none
    integer,intent(in) :: nrg     ! Unit number
    return
  end subroutine chem_read_restart_bin
!===============================================================================
  subroutine chem_write_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_write_restart_pio
!===============================================================================
  subroutine chem_read_restart_pio( File )
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_read_restart_pio
!===============================================================================
  subroutine chem_init_restart(File)
    use pio, only : file_desc_t
    type(file_desc_t) :: File
    return
  end subroutine chem_init_restart

end module chemistry
