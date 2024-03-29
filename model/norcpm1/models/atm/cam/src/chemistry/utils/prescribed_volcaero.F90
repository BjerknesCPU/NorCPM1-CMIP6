!-------------------------------------------------------------------
! manages reading and interpolation of prescribed volcanic aerosol
! Created by: Francis Vitt
!-------------------------------------------------------------------
module prescribed_volcaero

  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
  use spmd_utils,   only : masterproc
  use tracer_data,  only : trfld, trfile
  use cam_logfile,  only : iulog

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: prescribed_volcaero_readnl
  public :: prescribed_volcaero_register
  public :: prescribed_volcaero_init
  public :: prescribed_volcaero_adv
  public :: write_prescribed_volcaero_restart
  public :: read_prescribed_volcaero_restart
  public :: has_prescribed_volcaero
#ifdef CMIP6
  public :: has_prescribed_volcaero_cmip6,solar_bands,terrestrial_bands
#endif 
  public :: init_prescribed_volcaero_restart


  logical :: has_prescribed_volcaero = .false.
  character(len=8), parameter :: volcaero_name = 'VOLC_MMR'
  character(len=13), parameter :: volcrad_name = 'VOLC_RAD_GEOM'
  character(len=9), parameter :: volcmass_name = 'VOLC_MASS'
  character(len=11), parameter :: volcmass_column_name = 'VOLC_MASS_C'

  ! These variables are settable via the namelist (with longer names)
  character(len=16)  :: fld_name = 'MMRVOLC'
  character(len=256) :: filename = ''
  character(len=256) :: filelist = ''
  character(len=256) :: datapath = ''
  character(len=32)  :: data_type = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr  = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  integer            :: radius_ndx

#ifdef CMIP6
  logical, save :: has_prescribed_volcaero_cmip6 = .false.
  integer, parameter :: solar_bands=14, terrestrial_bands=16
  character(len=256) :: locfn
#endif 

contains

!-------------------------------------------------------------------
!-------------------------------------------------------------------
subroutine prescribed_volcaero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand
#ifdef CMIP6
   use cam_pio_utils,   only : cam_pio_openfile, init_pio_subsystem
   use pio,             only : pio_inquire, file_desc_t, pio_inq_dimname
   use pio,             only : pio_nowrite, pio_closefile, pio_noerr
#endif 

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'prescribed_volcaero_readnl'

#ifdef CMIP6
   integer :: dimid,ndims
   type(file_desc_t) :: ncid
   character(len=80) :: dimname
#endif 

   character(len=16)  :: prescribed_volcaero_name
   character(len=256) :: prescribed_volcaero_file
   character(len=256) :: prescribed_volcaero_filelist
   character(len=256) :: prescribed_volcaero_datapath
   character(len=32)  :: prescribed_volcaero_type
   logical            :: prescribed_volcaero_rmfile
   integer            :: prescribed_volcaero_cycle_yr
   integer            :: prescribed_volcaero_fixed_ymd
   integer            :: prescribed_volcaero_fixed_tod

   namelist /prescribed_volcaero_nl/ &
      prescribed_volcaero_name,      &
      prescribed_volcaero_file,      &
      prescribed_volcaero_filelist,  &
      prescribed_volcaero_datapath,  &
      prescribed_volcaero_type,      &
      prescribed_volcaero_rmfile,    &
      prescribed_volcaero_cycle_yr,  &
      prescribed_volcaero_fixed_ymd, &
      prescribed_volcaero_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   prescribed_volcaero_name     = fld_name
   prescribed_volcaero_file     = filename
   prescribed_volcaero_filelist = filelist
   prescribed_volcaero_datapath = datapath
   prescribed_volcaero_type     = data_type
   prescribed_volcaero_rmfile   = rmv_file
   prescribed_volcaero_cycle_yr = cycle_yr
   prescribed_volcaero_fixed_ymd= fixed_ymd
   prescribed_volcaero_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'prescribed_volcaero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, prescribed_volcaero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(prescribed_volcaero_name,     len(prescribed_volcaero_name),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_file,     len(prescribed_volcaero_file),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_filelist, len(prescribed_volcaero_filelist), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_datapath, len(prescribed_volcaero_datapath), mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_type,     len(prescribed_volcaero_type),     mpichar, 0, mpicom)
   call mpibcast(prescribed_volcaero_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(prescribed_volcaero_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(prescribed_volcaero_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   fld_name   = prescribed_volcaero_name
   filename   = prescribed_volcaero_file
   filelist   = prescribed_volcaero_filelist
   datapath   = prescribed_volcaero_datapath
   data_type  = prescribed_volcaero_type
   rmv_file   = prescribed_volcaero_rmfile
   cycle_yr   = prescribed_volcaero_cycle_yr
   fixed_ymd  = prescribed_volcaero_fixed_ymd
   fixed_tod  = prescribed_volcaero_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 ) has_prescribed_volcaero = .true.

#ifdef CMIP6
   ! Check if input file contains CMIP6 forcing   
   if (has_prescribed_volcaero) then
      if (len_trim(datapath) > 0 ) then
        locfn=trim(datapath)//'/'//trim(filename)
      else
        locfn=trim(filename)
      endif
      call init_pio_subsystem('atm_in')
      call cam_pio_openfile(ncid,locfn,PIO_NOWRITE)
      ierr = pio_inquire(ncid,ndimensions=ndims)
      do dimid=1,ndims
         ierr = pio_inq_dimname(ncid,dimid,dimname)
         if ( trim(dimname) == 'altitude' ) then
           has_prescribed_volcaero = .false.
           has_prescribed_volcaero_cmip6 = .true.
         endif
      enddo
      call pio_closefile(ncid)
   endif
#endif

end subroutine prescribed_volcaero_readnl

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_register()
    use ppgrid,         only: pver
    use phys_buffer,    only: pbuf_add
    integer :: idx

#ifdef CMIP6
    integer :: band
    character(len=3) :: c3
#endif 

    if (has_prescribed_volcaero) then
       call pbuf_add(volcaero_name,'physpkg',1,pver,1,idx)
       call pbuf_add(volcrad_name, 'physpkg',1,pver,1,idx)

    endif

#ifdef CMIP6
    if (has_prescribed_volcaero_cmip6) then
       do band=1,solar_bands
         write(c3,'(i3)') band
         call pbuf_add('ext_sun'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
         call pbuf_add('omega_sun'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
         call pbuf_add('g_sun'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
       enddo 
       do band=1,terrestrial_bands
         write(c3,'(i3)') band
         call pbuf_add('ext_earth'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
         call pbuf_add('omega_earth'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
         call pbuf_add('g_earth'//trim(adjustl(c3)),'physpkg',1,pver,1,idx)
       enddo 
    endif
#endif 

  endsubroutine prescribed_volcaero_register

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, phys_decomp
    use ppgrid,      only : pver
    use error_messages, only: handle_err
    use ppgrid,         only: pcols, pver, begchunk, endchunk
    use phys_buffer,  only : pbuf, pbuf_get_fld_idx

    implicit none

    integer :: ndx, istat
    character(len=32) :: specifier(1)

#ifdef CMIP6
    integer :: band, n
    character(len=3) :: c3
    character(len=32) :: specifier_cmip6(3*(solar_bands+terrestrial_bands))
#endif 
    
#ifdef CMIP6
    if ( has_prescribed_volcaero .or. has_prescribed_volcaero_cmip6 ) then
#else
    if ( has_prescribed_volcaero ) then
#endif 
       if ( masterproc ) then
          write(iulog,*) 'volcanic aerosol is prescribed in :'//trim(filename)
       endif
    else
       return
    endif

#ifdef CMIP6
    ! non-CMIP6
    if ( has_prescribed_volcaero ) then
#endif 

    specifier(1) = trim(volcaero_name)//':'//trim(fld_name)

    file%in_pbuf = .true.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type )


    call addfld(volcaero_name,'kg/kg', pver, 'I', 'prescribed volcanic aerosol dry mass mixing ratio', phys_decomp )
    call addfld(volcrad_name,'m', pver, 'I', 'volcanic aerosol geometric-mean radius', phys_decomp )
    call addfld(volcmass_name,'kg/m^2', pver, 'I', 'volcanic aerosol vertical mass path in layer', phys_decomp )
    call addfld(volcmass_column_name,'kg/m^2', 1, 'I', 'volcanic aerosol column mass', phys_decomp )

    radius_ndx = pbuf_get_fld_idx(volcrad_name, failcode=-1)

#ifdef CMIP6
    else

    do band=1,solar_bands
       write(c3,'(i3)') band
       specifier_cmip6(band*3-2) = 'ext_sun'//trim(adjustl(c3))//':'//'ext_sun'//trim(adjustl(c3))
       specifier_cmip6(band*3-1) = 'omega_sun'//trim(adjustl(c3))//':'//'omega_sun'//trim(adjustl(c3))
       specifier_cmip6(band*3-0) = 'g_sun'//trim(adjustl(c3))//':'//'g_sun'//trim(adjustl(c3))
       call addfld('ext_sun'//trim(adjustl(c3)),'1/km', pver, 'I', 'Extinction coefficient of solar bands', phys_decomp )
       call addfld('omega_sun'//trim(adjustl(c3)),'1', pver, 'I', 'Single scattering albedo of solar bands', phys_decomp )
       call addfld('g_sun'//trim(adjustl(c3)),'1', pver, 'I', 'Asymmetry factor of solar bands', phys_decomp )
    enddo
    do band=1,terrestrial_bands
       write(c3,'(i3)') band
       specifier_cmip6((solar_bands+band)*3-2) = 'ext_earth'//trim(adjustl(c3))//':'//'ext_earth'//trim(adjustl(c3))
       specifier_cmip6((solar_bands+band)*3-1) = 'omega_earth'//trim(adjustl(c3))//':'//'omega_earth'//trim(adjustl(c3))
       specifier_cmip6((solar_bands+band)*3-0) = 'g_earth'//trim(adjustl(c3))//':'//'g_earth'//trim(adjustl(c3))
       call addfld('ext_earth'//trim(adjustl(c3)),'1/km', pver, 'I', 'Extinction coefficient of terrestrial bands', phys_decomp )
       call addfld('omega_earth'//trim(adjustl(c3)),'1', pver, 'I', 'Single scattering albedo of terrestrial bands', phys_decomp )
       call addfld('g_earth'//trim(adjustl(c3)),'1', pver, 'I', 'Asymmetry factor of terrestrial bands', phys_decomp )
    enddo

    file%in_pbuf = .true.
    call trcdata_init( specifier_cmip6, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)

    endif 
#endif 

  end subroutine prescribed_volcaero_init

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine prescribed_volcaero_adv( state, pbuf )

    use tracer_data,  only : advance_trcdata
    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use ppgrid,       only : pcols, pver
    use string_utils, only : to_lower, GLC
    use cam_history,  only : outfld
    use physconst,    only : mwdry                ! molecular weight dry air ~ kg/kmole
    use physconst,    only : boltz, gravit        ! J/K/molecule
    use tropopause,   only : tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE
    use phys_buffer,  only : pbuf_size_max, pbuf_fld
    
    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(pbuf_fld),      intent(inout) :: pbuf(pbuf_size_max)

    integer :: c,ncol,i,k
    real(r8) :: to_mmr(pcols,pver)
    real(r8), parameter :: molmass = 47.9981995_r8
    real(r8) :: ptrop
    real(r8) :: concvolc ! micrograms of wetted aerosol per cubic centimeter
    real(r8) :: volcmass(pcols,pver)
    real(r8) :: columnmass(pcols)
    real(r8) :: mmrvolc
    integer  :: tropLev(pcols)

    real(r8) :: outdata(pcols,pver)
    real(r8), pointer :: data(:,:)
    real(r8), pointer :: radius(:,:)

    !WACCM-derived relation between mass concentration and wet aerosol radius in meters
    real(r8),parameter :: radius_conversion = 1.9e-4_r8

#ifdef CMIP6
    ! CMIP6
    integer :: band 
    character(len=3) :: c3

    if ( .not. (has_prescribed_volcaero .or. has_prescribed_volcaero_cmip6) ) return
#else 
    if( .not. has_prescribed_volcaero ) return
#endif

    call advance_trcdata( fields, file, state, pbuf=pbuf )

#ifdef CMIP6
    ! non-CMIP6
    if ( has_prescribed_volcaero ) then
#endif

    ! copy prescribed tracer fields into state svariable with the correct units
    do c = begchunk,endchunk
       radius => pbuf(radius_ndx)%fld_ptr(1,:,:,c,1)
       radius(:,:) = 0._r8
       ncol = state(c)%ncol
       select case ( to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))) )
       case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
          to_mmr(:ncol,:) = (molmass*1.e6_r8*boltz*state(c)%t(:ncol,:))/(mwdry*state(c)%pmiddry(:ncol,:))
       case ('kg/kg','mmr','kg kg-1')
          to_mmr(:ncol,:) = 1._r8
       case ('mol/mol','mole/mole','vmr','fraction')
          to_mmr(:ncol,:) = molmass/mwdry
       case default
          write(iulog,*) 'prescribed_volcaero_adv: units = ',trim(fields(1)%units) ,' are not recognized'
          call endrun('prescribed_volcaero_adv: units are not recognized')
       end select

       data => pbuf(fields(1)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
       data(:ncol,:) = to_mmr(:ncol,:) * data(:ncol,:) ! mmr

       call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
       do i = 1,ncol
          do k = 1,pver
             ! set to zero below tropopause
             if ( k >= tropLev(i) ) then
                data(i,k) = 0._r8
             endif
             mmrvolc = data(i,k)
             if (mmrvolc > 0._r8) then
                concvolc = (mmrvolc * state(c)%pdel(i,k))/(gravit * state(c)%zm(i,k))
                radius(i,k) = radius_conversion*(concvolc**(1._r8/3._r8))
             endif
          enddo
       enddo

       volcmass(:ncol,:) = data(:ncol,:)*state(c)%pdel(:ncol,:)/gravit
       columnmass(:ncol) = sum(volcmass(:ncol,:), 2)

       call outfld( volcaero_name,        data(:,:),     pcols, state(c)%lchnk)
       call outfld( volcrad_name,         radius(:,:),   pcols, state(c)%lchnk)
       call outfld( volcmass_name,        volcmass(:,:), pcols, state(c)%lchnk)
       call outfld( volcmass_column_name, columnmass(:), pcols, state(c)%lchnk)

    enddo

#ifdef CMIP6
    ! CMIP6 
    else 

    do c = begchunk,endchunk
       call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
       ncol = state(c)%ncol
       do band=1,solar_bands
          write(c3,'(i3)') band
          data => pbuf(fields(band*3-2)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0._r8
             enddo
          enddo
          call outfld('ext_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
          data => pbuf(fields(band*3-1)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0.999_r8
             enddo
          enddo
          call outfld('omega_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
          data => pbuf(fields(band*3-0)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0.5_r8
             enddo
          enddo
          call outfld('g_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
       enddo
       do band=1,terrestrial_bands
          write(c3,'(i3)') band
          data => pbuf(fields((solar_bands+band)*3-2)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0._r8
             enddo
          enddo
          call outfld('ext_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
          data => pbuf(fields((solar_bands+band)*3-1)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0.999_r8
             enddo
          enddo
          call outfld('omega_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
          data => pbuf(fields((solar_bands+band)*3-0)%pbuf_ndx)%fld_ptr(1,:,:,c,1)
          do i = 1,ncol
             data(i,:)=data(i,pver:1:-1) ! flip data as workaround for bug in tracer_data.F90  
             do k = 1,pver
                if ( k >= tropLev(i) ) data(i,k) = 0.5_r8
             enddo
          enddo
          call outfld('g_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
       enddo 
    enddo 

    endif 
#endif

  end subroutine prescribed_volcaero_adv

!-------------------------------------------------------------------
  subroutine init_prescribed_volcaero_restart( piofile )
    use pio, only : file_desc_t
    use tracer_data, only : init_trc_restart
    implicit none
    type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer

    call init_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine init_prescribed_volcaero_restart
!-------------------------------------------------------------------
  subroutine write_prescribed_volcaero_restart( piofile )
    use tracer_data, only : write_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call write_trc_restart( piofile, file )

  end subroutine write_prescribed_volcaero_restart
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine read_prescribed_volcaero_restart( pioFile )
    use tracer_data, only : read_trc_restart
    use pio, only : file_desc_t
    implicit none

    type(file_desc_t) :: piofile

    call read_trc_restart( 'prescribed_volcaero', piofile, file )

  end subroutine read_prescribed_volcaero_restart

end module prescribed_volcaero
