module drydep_so2

  !---------------------------------------------------------------------
  !       ... Dry deposition velocity input data and code for netcdf input
  ! Due to troubles with implementing mozart in cam-oslo the module is 
  ! hardcoded for SO2 at present time
  !---------------------------------------------------------------------
  !
  !---------------------------------------------------------------------

  use shr_kind_mod, only : r8 => shr_kind_r8
!  use chem_mods,    only : gas_pcnst
  use spmd_utils,   only : masterproc, iam
  use commap,       only : clat, clon
  use ppgrid,       only : pcols, begchunk, endchunk
  use abortutils,   only : endrun
#ifdef SPMD
  use mpishorthand, only : mpicom, mpir8, mpiint, mpilog
#endif
  use cam_logfile,  only : iulog
  use dyn_grid,     only : get_dyn_grid_parm
  use wrap_nf
  use ioFilemod,    only : getfil
  use pmgrid,       only : plev, plevp
!  use aerosoldef,   only : ndvel_gas
  implicit none

  save
  private
  public dvel_so2_inti
  public dvel_so2
  public interpso2_map
  public set_soilwso2 
  real(r8)              :: dels
  real(r8), allocatable :: days(:)              ! day of year for soilw
  real(r8), allocatable :: dvel(:,:,:,:), &   ! depvel array interpolated to model grid
       dvel_interp(:,:,:) ! depvel array interpolated to grid and time
  integer :: last, next                     ! day indicies
  integer :: ndays                          ! # of days in soilw file

  integer :: so2_ndx



!!$  integer, parameter :: n_land_type = 11
  integer, parameter :: ndvel_gas = 1     ! Total number of dry deposited gases
  integer :: mapping(ndvel_gas) = -99
  real(r8), parameter    :: small_value = 1.e-36_r8
  real(r8), parameter    :: large_value = 1.e36_r8
  real(r8), parameter    :: diffm       = 1.789e-5_r8
  real(r8), parameter    :: diffk       = 1.461e-5_r8
  real(r8), parameter    :: difft       = 2.060e-5_r8
  real(r8), parameter    :: vonkar      = 0.378_r8
  real(r8), parameter    :: ric         = 0.2_r8
  real(r8), parameter    :: r           = 287.04_r8
  real(r8), parameter    :: cp          = 1004._r8
  real(r8), parameter    :: grav        = 9.81_r8
  real(r8), parameter    :: p00         = 100000._r8
  real(r8), parameter    :: wh2o        = 18.0153_r8
  real(r8), parameter    :: ph          = 1.e-5_r8
  real(r8), parameter    :: ph_inv      = 1._r8/ph
  real(r8), parameter    :: rovcp = r/cp

  integer, allocatable :: index_season_lai(:,:,:)

  real(r8)    :: crb = 0.0_r8
  real(r8)    :: foxd(ndvel_gas)     = small_value
  real(r8)    :: drat(ndvel_gas)     = small_value

!  real(r8)                         :: soilw_3d(plon,plat,12)
  real(r8) , allocatable            :: soilw_3d(:,:,:)
  logical :: has_dvel(ndvel_gas) = .true.

  real(r8), public,allocatable  :: fraction_landuse(:,:,:)
  real(r8), allocatable, dimension(:,:,:) :: dep_ra ! [s/m] aerodynamic resistance
  real(r8), allocatable, dimension(:,:,:) :: dep_rb ! [s/m] resistance across sublayer
  integer, parameter :: n_land_type = 11



  integer :: pan_ndx, mpan_ndx, no2_ndx, hno3_ndx, o3_ndx, &
             h2o2_ndx, onit_ndx, onitr_ndx, ch4_ndx, ch2o_ndx, &
             ch3ooh_ndx, pooh_ndx, ch3coooh_ndx, c2h5ooh_ndx, &
             c3h7ooh_ndx, rooh_ndx, ch3cocho_ndx, co_ndx, ch3coch3_ndx, &
             no_ndx, ho2no2_ndx, glyald_ndx, hyac_ndx, ch3oh_ndx, c2h5oh_ndx, &
             hydrald_ndx, h2_ndx, Pb_ndx, o3s_ndx, o3inert_ndx, macrooh_ndx, &
             xooh_ndx, ch3cho_ndx, isopooh_ndx
  integer :: alkooh_ndx, mekooh_ndx, tolooh_ndx, terpooh_ndx, ch3cooh_ndx
  integer :: soa_ndx, so4_ndx, cb1_ndx, cb2_ndx, oc1_ndx, oc2_ndx, nh3_ndx, nh4no3_ndx, &
             sa1_ndx, sa2_ndx, sa3_ndx, sa4_ndx, nh4_ndx



contains

  subroutine dvel_so2_inti( depvel_lnd_file, clim_soilw_file, season_wes_file)
    !-------------------------------------------------------------------------------------
    ! 	... intialize interactive drydep
    !-------------------------------------------------------------------------------------

    use mo_constants,  only : r2d
    use commap,        only : clat, clon
!    use chem_mods,     only : adv_mass
    use mo_drydep_tables
!    use mo_chem_utls,  only : get_spc_ndx

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    character*80, intent(in) :: depvel_lnd_file, clim_soilw_file, season_wes_file 
!    character(len=*), intent(in) :: drydep_list(:)

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer :: i, j, ii, jj, jl, ju
    integer :: nlon_veg, nlat_veg, npft_veg
    integer :: nlat_lai, npft_lai, pos_min, imin
    integer :: ncid, vid, dimid
    integer :: m, n, l, id
    integer :: length1, astat
    integer, allocatable :: wk_lai(:,:,:)
    integer :: k, num_max, k_max
    integer :: num_seas(5)
    integer :: plon, plat
    integer :: ierr

    real(r8)              :: spc_mass
    real(r8)              :: diff_min, target_lat
    real(r8), allocatable :: vegetation_map(:,:,:)
    real(r8), pointer     :: soilw_map(:,:,:)
    real(r8), allocatable :: work(:,:)
    real(r8), allocatable :: landmask(:,:)
    real(r8), allocatable :: urban(:,:)
    real(r8), allocatable :: lake(:,:)
    real(r8), allocatable :: wetland(:,:)
    real(r8), allocatable :: lon_veg(:)
    real(r8), allocatable :: lon_veg_edge(:)
    real(r8), allocatable :: lat_veg(:)
    real(r8), allocatable :: lat_veg_edge(:)
    real(r8), allocatable :: lat_lai(:)

    character(len=32) :: test_name
    character(len=256) locfn      ! local filename
    logical :: do_soilw
! Include ndx even though I do not use it in order to maintain the 
!dry deposition  routine with as small changes as possible

    ch4_ndx      = -1
    h2_ndx       = -1 
    co_ndx       = -1 
    Pb_ndx       = -1 
    pan_ndx      = -1
    mpan_ndx     = -1
    o3_ndx       = -1
    if( o3_ndx < 0 ) then
       o3_ndx  = -1
    end if
!
    so2_ndx     = 1
!
    alkooh_ndx  = -1
    mekooh_ndx  = -1
    tolooh_ndx  = -1
    terpooh_ndx = -1
    ch3cooh_ndx = -1
    soa_ndx     = -1
    so4_ndx     = -1
    cb1_ndx     = -1
    cb2_ndx     = -1
    oc1_ndx     = -1
    oc2_ndx     = -1
    nh3_ndx     = -1
    nh4no3_ndx  = -1
    sa1_ndx     = -1
    sa2_ndx     = -1
    sa3_ndx     = -1
    sa4_ndx     = -1
    nh4_ndx     = -1
!

    plon = get_dyn_grid_parm('plon')
    plat = get_dyn_grid_parm('plat')

          spc_mass = 32.07_r8
          drat(1) = sqrt( spc_mass/wh2o )
    !---------------------------------------------------------------------------
    ! 	... allocate module variables
    !---------------------------------------------------------------------------
    allocate( dep_ra(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_ra; error = ',astat
       call endrun
    end if
    allocate( dep_rb(pcols,n_land_type,begchunk:endchunk),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate dep_rb; error = ',astat
       call endrun
    end if
    allocate( fraction_landuse(plon,plat,n_land_type),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate fraction_landuse; error = ',astat
       call endrun
    end if
    allocate( index_season_lai(plat,n_land_type,12),stat=astat )
    if( astat /= 0 ) then
       write(iulog,*) 'dvel_inti: failed to allocate index_season_lai; error = ',astat
       call endrun
    end if

    do_soilw=.false.
    if(do_soilw) then
       allocate(soilw_3d(plon,plat,12))
    end if

    Masterproc_only :if( masterproc ) then

       !---------------------------------------------------------------------------
       ! 	... read landuse map
       !---------------------------------------------------------------------------

       call getfil (depvel_lnd_file, locfn, 0)


       call wrap_open (trim(locfn), NF_NOWRITE, ncid)

       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       call wrap_inq_dimid( ncid, 'lon', dimid )
       call wrap_inq_dimlen( ncid, dimid, nlon_veg )
       call wrap_inq_dimid( ncid, 'lat', dimid )
       call wrap_inq_dimlen( ncid, dimid, nlat_veg )
       call wrap_inq_dimid( ncid, 'pft', dimid )
       call wrap_inq_dimlen( ncid, dimid, npft_veg )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( vegetation_map(nlon_veg,nlat_veg,npft_veg), work(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( urban(nlon_veg,nlat_veg), lake(nlon_veg,nlat_veg), &
            landmask(nlon_veg,nlat_veg), wetland(nlon_veg,nlat_veg), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( lon_veg(nlon_veg), lat_veg(nlat_veg), &
            lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read the vegetation map and landmask
       !---------------------------------------------------------------------------
       call wrap_inq_varid( ncid, 'PCT_PFT', vid )
       call wrap_get_var_realx( ncid, vid, vegetation_map )

       call wrap_inq_varid( ncid, 'LANDMASK', vid )
       call wrap_get_var_realx( ncid, vid, landmask )

       call wrap_inq_varid( ncid, 'PCT_URBAN', vid )
       call wrap_get_var_realx( ncid, vid, urban )

       call wrap_inq_varid( ncid, 'PCT_LAKE', vid )
       call wrap_get_var_realx( ncid, vid, lake )

       call wrap_inq_varid( ncid, 'PCT_WETLAND', vid )
       call wrap_get_var_realx( ncid, vid, wetland )

       call wrap_close( ncid )

       !---------------------------------------------------------------------------
       ! scale vegetation, urban, lake, and wetland to fraction
       !---------------------------------------------------------------------------
       vegetation_map(:,:,:) = .01_r8 * vegetation_map(:,:,:)
       wetland(:,:)          = .01_r8 * wetland(:,:)
       lake(:,:)             = .01_r8 * lake(:,:)
!       write(iulog,*) 'vegetation_map ', so2_ndx,vegetation_map(:,:,:)
       urban(:,:)            = .01_r8 * urban(:,:)
#ifdef DEBUG
       write(iulog,*) 'minmax vegetation_map ',minval(vegetation_map),maxval(vegetation_map)
       write(iulog,*) 'minmax wetland        ',minval(wetland),maxval(wetland)
       write(iulog,*) 'minmax landmask       ',minval(landmask),maxval(landmask)
#endif
       !---------------------------------------------------------------------------
       ! 	... define lat-lon of vegetation map (1x1)
       !---------------------------------------------------------------------------
       lat_veg(:)      = (/ (-89.5_r8 + (i-1),i=1,nlat_veg  ) /)
       lon_veg(:)      = (/ (  0.5_r8 + (i-1),i=1,nlon_veg  ) /)
       lat_veg_edge(:) = (/ (-90.0_r8 + (i-1),i=1,nlat_veg+1) /)
       lon_veg_edge(:) = (/ (  0.0_r8 + (i-1),i=1,nlon_veg+1) /)
       !---------------------------------------------------------------------------
       ! 	... read soilw table if necessary
       !---------------------------------------------------------------------------
       if( do_soilw ) then
          call soilwso2_inti( clim_soilw_file, nlon_veg, nlat_veg, soilw_map )
       end if

       !---------------------------------------------------------------------------
       ! 	... regrid to model grid
       !---------------------------------------------------------------------------

!#ifdef feil
       call interpso2_map( plon, plat, nlon_veg, nlat_veg, npft_veg, & 
                        lat_veg, lat_veg_edge, &
                        lon_veg, lon_veg_edge, landmask, urban, lake, &
                        wetland, vegetation_map, soilw_map, do_soilw )
!#endif
       deallocate( vegetation_map, work, stat=astat )
       deallocate( lon_veg, lat_veg, lon_veg_edge, lat_veg_edge, stat=astat )
       deallocate( landmask, urban, lake, wetland, stat=astat )
       if( do_soilw ) then
          deallocate( soilw_map, stat=astat )
       end if

       !---------------------------------------------------------------------------
       ! 	... read LAI based season indeces
       !---------------------------------------------------------------------------
       call getfil (season_wes_file, locfn, 0)
       call wrap_open (trim(locfn), NF_NOWRITE, ncid)

       !---------------------------------------------------------------------------
       ! 	... get the dimensions
       !---------------------------------------------------------------------------
       call wrap_inq_dimid( ncid, 'lat', dimid )
       call wrap_inq_dimlen( ncid, dimid, nlat_lai )
       call wrap_inq_dimid( ncid, 'pft', dimid )
       call wrap_inq_dimlen( ncid, dimid, npft_lai )
       !---------------------------------------------------------------------------
       ! 	... allocate arrays
       !---------------------------------------------------------------------------
       allocate( lat_lai(nlat_lai), wk_lai(nlat_lai,npft_lai,12), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation_map; error = ',astat
          call endrun
       end if
       allocate( lon_veg(nlon_veg), lat_veg(nlat_veg), lon_veg_edge(nlon_veg+1), lat_veg_edge(nlat_veg+1), stat=astat )
       if( astat /= 0 ) then
          write(iulog,*) 'dvel_inti: failed to allocate vegation lon, lat arrays; error = ',astat
          call endrun
       end if
       !---------------------------------------------------------------------------
       ! 	... read the latitude and the season indicies
       !---------------------------------------------------------------------------
       call wrap_inq_varid( ncid, 'lat', vid )
       call wrap_get_var_realx( ncid, vid, lat_lai )

       call wrap_inq_varid( ncid, 'season_wes', vid )
       call wrap_get_var_int( ncid, vid, wk_lai )

       call wrap_close( ncid )

       jl = 1
       ju = plat
       imin = 1
       do j = 1,plat
          diff_min = 10._r8
          pos_min  = -99
          target_lat = clat(j)*r2d
          do i = imin,nlat_lai
             if( abs(lat_lai(i) - target_lat) < diff_min ) then
                diff_min = abs(lat_lai(i) - target_lat)
                pos_min  = i
             end if
          end do
          if( pos_min < 0 ) then
             write(iulog,*) 'dvel_inti: cannot find ',target_lat,' at j,pos_min,diff_min = ',j,pos_min,diff_min
             write(iulog,*) 'dvel_inti: imin,nlat_lai = ',imin,nlat_lai
             write(iulog,*) 'dvel_inti: lat_lai'
             write(iulog,'(1p,10g12.5)') lat_lai(:)
             call endrun
          end if
          imin = pos_min
          index_season_lai(j,:,:) = wk_lai(pos_min,:,:)

          !---------------------------------------------------------------------------
          ! specify the season as the most frequent in the 11 vegetation classes
          ! this was done to remove a banding problem in dvel (JFL Oct 04)
          !---------------------------------------------------------------------------
          do m = 1,12
             num_seas = 0
             do l = 1,11
                do k = 1,5
                   if( index_season_lai(j,l,m) == k ) then
                      num_seas(k) = num_seas(k) + 1
                      exit
                   end if
                end do
             end do

             num_max = -1
             do k = 1,5
                if( num_seas(k) > num_max ) then
                   num_max = num_seas(k)
                   k_max = k
                endif
             end do

             index_season_lai(j,:,m) = k_max
          end do
       end do

       deallocate( lat_lai, wk_lai )

    end if Masterproc_only

#ifdef SPMD

    call mpibarrier( mpicom )
    call mpibcast( index_season_lai, plat*n_land_type*12, mpiint, 0, mpicom )
    call mpibcast( fraction_landuse, plon*plat*n_land_type, mpir8, 0, mpicom )
    call mpibcast( do_soilw, 1, mpilog, 0, mpicom )

    if ( do_soilw ) then
       call mpibcast( soilw_3d, plon*plat*12, mpir8, 0,  mpicom )
       call mpibcast( ndays, 1, mpiint, 0, mpicom )
       if ( .not. masterproc ) then
          allocate( days(ndays),stat=ierr )
          if( ierr /= 0 ) then
             write(iulog,*) 'soilwso2_inti: days allocation error = ',ierr
             call endrun
          end if
       endif
       call mpibcast( days, ndays, mpir8, 0, mpicom )
    endif

#endif

    rgss(:,:) = max( 1._r8,rgss(:,:) )
    rac(:,:)  = max( small_value,rac(:,:) )

  end subroutine dvel_so2_inti




    subroutine interpso2_map( plon, plat, nlon_veg, nlat_veg, npft_veg, &
                         lat_veg, lat_veg_edge, &
                         lon_veg, lon_veg_edge, landmask, urban, lake, &
                         wetland, vegetation_map, soilw_map, do_soilw )

    use mo_constants, only : r2d
    use commap,       only : clat, clon
    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)      ::  plon, plat, nlon_veg, nlat_veg, npft_veg
    real(r8), pointer            :: soilw_map(:,:,:)
    real(r8), intent(in)         :: landmask(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: urban(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: lake(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: wetland(nlon_veg,nlat_veg)
    real(r8), intent(in)         :: vegetation_map(nlon_veg,nlat_veg,npft_veg)
    real(r8), intent(in)         :: lon_veg(nlon_veg)
    real(r8), intent(in)         :: lon_veg_edge(nlon_veg+1)
    real(r8), intent(in)         :: lat_veg(nlat_veg)
    real(r8), intent(in)         :: lat_veg_edge(nlat_veg+1)
    logical, intent(in)      :: do_soilw

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    integer, parameter           :: veg_ext = 10

    integer                      :: i, j, ii, jj, jl, ju, i_ndx, n
    integer, dimension(plon+1)   :: ind_lon
    integer, dimension(plat+1)  :: ind_lat
    real(r8)                         :: total_land
    real(r8), dimension(plon+1)      :: lon_edge
    real(r8), dimension(plat+1)     :: lat_edge
    real(r8)                         :: lat1, lat2, lon1, lon2
    real(r8)                         :: x1, x2, y1, y2, dx, dy
    real(r8)                         :: area, total_area
    real(r8), dimension(npft_veg+3)  :: fraction
    real(r8)                         :: total_soilw_area
    real(r8)                         :: fraction_soilw
    real(r8)                         :: total_soilw(12)

    real(r8),    dimension(-veg_ext:nlon_veg+veg_ext) :: lon_veg_edge_ext
    integer, dimension(-veg_ext:nlon_veg+veg_ext) :: mapping_ext

    real(r8), dimension(plon) :: lam
    real(r8), dimension(plat) :: phi

    logical, parameter :: has_npole = .true.

    lam(:) = clon(:,1)
    phi(:) = clat(:)

    jl = 1
    ju = plon

    do i = 1,plon
       lon_edge(i) = lam(i) * r2d - .5_r8*(lam(2) - lam(1)) * r2d
    end do
    lon_edge(plon+1) = lon_edge(plon) + (lam(2) - lam(1)) * r2d
    if( .not. has_npole ) then
       do j = 1,plat+1
          lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
       end do
    else
       do j = 1,plat
          lat_edge(j) = phi(j) * r2d - .5_r8*(phi(2) - phi(1)) * r2d
       end do
       lat_edge(plat+1) = lat_edge(plat) + (phi(2) - phi(1)) * r2d
    end if
    do j = 1,plat+1
       lat_edge(j) = min( lat_edge(j), 90._r8 )
       lat_edge(j) = max( lat_edge(j),-90._r8 )
    end do

    !-------------------------------------------------------------------------------------
    ! wrap around the longitudes
    !-------------------------------------------------------------------------------------
    do i = -veg_ext,0
       lon_veg_edge_ext(i) = lon_veg_edge(nlon_veg+i) - 360._r8
       mapping_ext     (i) =              nlon_veg+i
    end do
    do i = 1,nlon_veg
       lon_veg_edge_ext(i) = lon_veg_edge(i)
       mapping_ext     (i) =              i
    end do
    do i = nlon_veg+1,nlon_veg+veg_ext
       lon_veg_edge_ext(i) = lon_veg_edge(i-nlon_veg) + 360._r8
       mapping_ext     (i) =              i-nlon_veg
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : lon_edge ',lon_edge
    write(iulog,*) 'interp_map : lat_edge ',lat_edge
    write(iulog,*) 'interp_map : mapping_ext ',mapping_ext
#endif
    do j = 1,plon+1
       lon1 = lon_edge(j) 
       do i = -veg_ext,nlon_veg+veg_ext
          dx = lon_veg_edge_ext(i  ) - lon1
          dy = lon_veg_edge_ext(i+1) - lon1
          if( dx*dy <= 0._r8 ) then
             ind_lon(j) = i
             exit
          end if
       end do
    end do

    do j = 1,plat+1
       lat1 = lat_edge(j)
       do i = 1,nlat_veg
          dx = lat_veg_edge(i  ) - lat1
          dy = lat_veg_edge(i+1) - lat1
          if( dx*dy <= 0._r8 ) then
             ind_lat(j) = i
             exit
          end if
       end do
    end do
#ifdef DEBUG
    write(iulog,*) 'interp_map : ind_lon ',ind_lon
    write(iulog,*) 'interp_map : ind_lat ',ind_lat
#endif
    lat_loop : do j = 1,plat
       lon_loop : do i = 1,plon
          total_area       = 0._r8
          fraction         = 0._r8
          total_soilw(:)   = 0._r8
          total_soilw_area = 0._r8
          do jj = ind_lat(j),ind_lat(j+1)
             y1 = max( lat_edge(j),lat_veg_edge(jj) )
             y2 = min( lat_edge(j+1),lat_veg_edge(jj+1) ) 
             dy = (y2 - y1)/(lat_veg_edge(jj+1) - lat_veg_edge(jj))
             do ii =ind_lon(i),ind_lon(i+1)
                i_ndx = mapping_ext(ii)
                x1 = max( lon_edge(i),lon_veg_edge_ext(ii) )
                x2 = min( lon_edge(i+1),lon_veg_edge_ext(ii+1) ) 
                dx = (x2 - x1)/(lon_veg_edge_ext(ii+1) - lon_veg_edge_ext(ii))
                area = dx * dy
                total_area = total_area + area
                !-----------------------------------------------------------------
                ! 	... special case for ocean grid point 
                !-----------------------------------------------------------------
                if( nint(landmask(i_ndx,jj)) == 0 ) then
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area
                else
                   do n = 1,npft_veg
                      fraction(n) = fraction(n) + vegetation_map(i_ndx,jj,n) * area
                   end do
                   fraction(npft_veg+1) = fraction(npft_veg+1) + area * lake   (i_ndx,jj)
                   fraction(npft_veg+2) = fraction(npft_veg+2) + area * wetland(i_ndx,jj)
                   fraction(npft_veg+3) = fraction(npft_veg+3) + area * urban  (i_ndx,jj)
                   !-----------------------------------------------------------------
                   ! 	... check if land accounts for the whole area.
                   !           If not, the remaining area is in the ocean
                   !-----------------------------------------------------------------
                   total_land = sum(vegetation_map(i_ndx,jj,:)) &
                              + urban  (i_ndx,jj) &
                              + lake   (i_ndx,jj) &
                              + wetland(i_ndx,jj)
                   if( total_land < 1._r8 ) then
                      fraction(npft_veg+1) = fraction(npft_veg+1) + (1._r8 - total_land) * area
                   end if
                   !-------------------------------------------------------------------------------------
                   ! 	... compute weighted average of soilw over grid (non-water only)
                   !-------------------------------------------------------------------------------------
                   if( do_soilw ) then
                      fraction_soilw = total_land  - (lake(i_ndx,jj) + wetland(i_ndx,jj))
                      total_soilw_area = total_soilw_area + fraction_soilw * area
                      total_soilw(:)   = total_soilw(:) + fraction_soilw * area * soilw_map(i_ndx,jj,:)
                   end if
                end if
             end do
          end do
          !-------------------------------------------------------------------------------------
          ! 	... divide by total area of grid box
          !-------------------------------------------------------------------------------------
          fraction(:) = fraction(:)/total_area
          !-------------------------------------------------------------------------------------
          ! 	... make sure we don't have too much or too little
          !-------------------------------------------------------------------------------------
          if( abs( sum(fraction) - 1._r8) > .001_r8 ) then
             fraction(:) = fraction(:)/sum(fraction)
          end if
          !-------------------------------------------------------------------------------------
          ! 	... map to Wesely land classification
          !-------------------------------------------------------------------------------------
          fraction_landuse(i,j, 1) =     fraction(20)
          fraction_landuse(i,j, 2) = sum(fraction(16:17))
          fraction_landuse(i,j, 3) = sum(fraction(13:15))
          fraction_landuse(i,j, 4) = sum(fraction( 5: 9))
          fraction_landuse(i,j, 5) = sum(fraction( 2: 4))
          fraction_landuse(i,j, 6) =     fraction(19)
          fraction_landuse(i,j, 7) =     fraction(18)
          fraction_landuse(i,j, 8) =     fraction( 1)
          fraction_landuse(i,j, 9) = 0._r8
          fraction_landuse(i,j,10) = 0._r8
          fraction_landuse(i,j,11) = sum(fraction(10:12))
          if( do_soilw ) then
             if( total_soilw_area > 0._r8 ) then
                soilw_3d(i,j,:) = total_soilw(:)/total_soilw_area
             else
                soilw_3d(i,j,:) = -99._r8
             end if
          end if
       end do lon_loop
    end do lat_loop
    !-------------------------------------------------------------------------------------
    ! 	... make sure there are no negative values
    !-------------------------------------------------------------------------------------
    if( any( fraction_landuse(:,:,:) < 0._r8 ) ) then
       write(iulog,*) 'something is wrong in interpolation of landuse map'
       write(iulog,*) 'minval(fraction_landuse) = ',minval(fraction_landuse)
       call endrun
    end if

    !-------------------------------------------------------------------------------------
    ! 	... reshape according to lat-lon blocks
    !-------------------------------------------------------------------------------------

  end subroutine interpso2_map
  

  subroutine soilwso2_inti( ncfile, nlon_veg, nlat_veg, soilw_map )
    !------------------------------------------------------------------
    !	... read primary soil moisture table
    !------------------------------------------------------------------

    use time_manager,  only : get_calday

    implicit none

    !------------------------------------------------------------------
    !	... dummy args
    !------------------------------------------------------------------
    integer, intent(in) :: &
         nlon_veg, &
         nlat_veg
    real(r8), pointer :: soilw_map(:,:,:)
    character(len=*), intent(in) :: ncfile ! file name of netcdf file containing data

    !------------------------------------------------------------------
    !	... local variables
    !------------------------------------------------------------------
    integer :: gndx = 0
    integer :: nlat, &             ! # of lats in soilw file
               nlon                ! # of lons in soilw file
    integer :: i, ip, k, m
    integer :: j, jl, ju
    integer :: lev, day, ierr
    integer :: ncid, vid
    integer :: dimid_lat, dimid_lon, dimid_time
    integer :: dates(12) = (/ 116, 214, 316, 415,  516,  615, &
                              716, 816, 915, 1016, 1115, 1216 /)

    character(len=256) :: locfn

    !-----------------------------------------------------------------------
    !       ... open netcdf file
    !-----------------------------------------------------------------------
    call getfil (ncfile, locfn, 0)
    call wrap_open (trim(locfn), NF_NOWRITE, ncid)

    !-----------------------------------------------------------------------
    !       ... get longitudes
    !-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'lon', dimid_lon )
    call wrap_inq_dimlen( ncid, dimid_lon, nlon )
    if( nlon /= nlon_veg ) then
       write(iulog,*) 'soilwso2_inti: soil and vegetation lons differ; ',nlon, nlon_veg
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... get latitudes
    !-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'lat', dimid_lat )
    call wrap_inq_dimlen( ncid, dimid_lat, nlat )
    if( nlat /= nlat_veg ) then
       write(iulog,*) 'soilwso2_inti: soil and vegetation lats differ; ',nlat, nlat_veg
       call endrun
    end if
    !-----------------------------------------------------------------------
    !       ... set times (days of year)
    !-----------------------------------------------------------------------
    call wrap_inq_dimid( ncid, 'time', dimid_time )
    call wrap_inq_dimlen( ncid, dimid_time, ndays )
    if( ndays /= 12 ) then
       write(iulog,*) 'soilwso2_inti: dataset not a cyclical year'
       call endrun
    end if
    allocate( days(ndays),stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soilwso2_inti: days allocation error = ',ierr
       call endrun
    end if
    do m = 1,min(12,ndays)
       days(m) = get_calday( dates(m), 0 )
    end do

    !------------------------------------------------------------------
    !	... allocate arrays
    !------------------------------------------------------------------
    allocate( soilw_map(nlon,nlat,ndays), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soilwso2_inti: soilw_map allocation error = ',ierr
       call endrun
    end if

    !------------------------------------------------------------------
    !	... read in the soil moisture
    !------------------------------------------------------------------
    call wrap_inq_varid( ncid, 'SOILW', vid )
    call wrap_get_var_realx( ncid, vid, soilw_map )
    !------------------------------------------------------------------
    !	... close file
    !------------------------------------------------------------------
    call wrap_close( ncid )

  end subroutine soilwso2_inti
  
  subroutine chk_soilwso2( calday )
    !--------------------------------------------------------------------
    !	... check timing for ub values
    !--------------------------------------------------------------------

    use mo_constants, only : dayspy

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    real(r8), intent(in)    :: calday

    !--------------------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------------------
    integer  ::  m, upper
    real(r8)     ::  numer, denom

    !--------------------------------------------------------
    !	... setup the time interpolation
    !--------------------------------------------------------
    if( calday < days(1) ) then
       next = 1
       last = ndays
    else
       if( days(ndays) < dayspy ) then
          upper = ndays
       else
          upper = ndays - 1
       end if
       do m = upper,1,-1
          if( calday >= days(m) ) then
             exit
          end if
       end do
       last = m
       next = mod( m,ndays ) + 1
    end if
    numer = calday - days(last)
    denom = days(next) - days(last)
    if( numer < 0._r8 ) then
       numer = dayspy + numer
    end if
    if( denom < 0._r8 ) then
       denom = dayspy + denom
    end if
    dels = max( min( 1._r8,numer/denom ),0._r8 )

  end subroutine chk_soilwso2

  subroutine set_soilwso2( soilw, ncol , lonndx, latndx, calday )
    !--------------------------------------------------------------------
    !	... set the soil moisture
    !--------------------------------------------------------------------

    implicit none

    !--------------------------------------------------------------------
    !	... dummy args
    !--------------------------------------------------------------------
    integer,  intent(in)    :: ncol
    real(r8), intent(inout) :: soilw(pcols)
    integer,  intent(in)    :: latndx(pcols)           ! chunk latitude indicies
    integer,  intent(in)    :: lonndx(pcols)           ! chunk longitude indicies
    real(r8), intent(in)    :: calday


    integer :: i, ilon,ilat

    call chk_soilwso2( calday )

    do i =1,ncol
       ilon = lonndx(i)
       ilat = latndx(i)
       soilw(i) = soilw_3d(ilon,ilat,last) + dels * ( soilw_3d(ilon,ilat,next) - soilw_3d(ilon,ilat,last) )
    enddo

  end subroutine set_soilwso2



  subroutine dvel_so2( ncdate, sfc_temp, pressure_sfc,  &
                             wind_speed, spec_hum, air_temp, pressure_10m, rain, &
                             snow, solar_flux, dvel, &
                             tv, soilw, rh, ncol, lonndx, latndx, lchnk )
    !-------------------------------------------------------------------------------------
    !   code based on wesely (atmospheric environment, 1989, vol 23, p. 1293-1304) for
    !   calculation of r_c, and on walcek et. al. (atmospheric enviroment, 1986,
    !   vol. 20, p. 949-964) for calculation of r_a and r_b
    !
    !   as suggested in walcek (u_i)(u*_i) = (u_a)(u*_a)
    !   is kept constant where i represents a subgrid environment and a the
    !   grid average environment. thus the calculation proceeds as follows:
    !   va the grid averaged wind is calculated on dots
    !   z0(i) the grid averaged roughness coefficient is calculated
    !   ri(i) the grid averaged richardson number is calculated
    !   --> the grid averaged (u_a)(u*_a) is calculated
    !   --> subgrid scale u*_i is calculated assuming (u_i) given as above
    !   --> final deposotion velocity is weighted average of subgrid scale velocities
    !
    ! code written by P. Hess, rewritten in fortran 90 by JFL (August 2000)
    ! modified by JFL to be used in MOZART-2 (October 2002)
    !-------------------------------------------------------------------------------------

!!$    use mo_seasalt,   only : seasalt_sett_vel, has_seasalt, nbin
    use mo_drydep_tables
    use physconst,    only : tmelt

    implicit none

    !-------------------------------------------------------------------------------------
    ! 	... dummy arguments
    !-------------------------------------------------------------------------------------
    integer, intent(in)   :: ncol
    integer, intent(in)   :: ncdate                   ! present date (yyyymmdd)
    real(r8), intent(in)      :: sfc_temp(pcols)          ! surface temperature (K)
    real(r8), intent(in)      :: pressure_sfc(pcols)      ! surface pressure (Pa)
    real(r8), intent(in)      :: wind_speed(pcols)        ! 10 meter wind speed (m/s)
    real(r8), intent(in)      :: spec_hum(pcols)          ! specific humidity (kg/kg)
    real(r8), intent(in)      :: rh(pcols,1)              ! relative humidity
    real(r8), intent(in)      :: air_temp(pcols)          ! surface air temperature (K)
    real(r8), intent(in)      :: pressure_10m(pcols)      ! 10 meter pressure (Pa)
    real(r8), intent(in)      :: rain(pcols)              
    real(r8), intent(in)      :: snow(pcols)              ! snow height (m)
    real(r8), intent(in)      :: soilw(pcols)             ! soil moisture fraction
    real(r8), intent(in)      :: solar_flux(pcols)        ! direct shortwave radiation at surface (W/m^2)
    real(r8), intent(in)      :: tv(pcols)                ! potential temperature
!    real(r8), intent(in)      :: mmr(pcols,plev,ndvel_gas)    ! constituent concentration (kg/kg)
    real(r8), intent(out)     :: dvel(pcols,ndvel_gas)        ! deposition velocity (cm/s)
!    real(r8), intent(inout)   :: dflx(pcols,ndvel_gas)        ! deposition flux (/cm^2/s)

    integer, intent(in)     ::   latndx(pcols)           ! chunk latitude indicies
    integer, intent(in)     ::   lonndx(pcols)           ! chunk longitude indicies
    integer, intent(in)     ::   lchnk                   ! chunk number

    !-------------------------------------------------------------------------------------
    ! 	... local variables
    !-------------------------------------------------------------------------------------
    real(r8), parameter :: scaling_to_cm_per_s = 100._r8
    real(r8), parameter :: rain_threshold      = 1.e-7_r8  ! of the order of 1cm/day expressed in m/s

    integer :: i, ispec, lt, m
    integer :: sndx
    integer :: month

    real(r8) :: slope = 0._r8
    real(r8) :: z0water ! revised z0 over water
    real(r8) :: p       ! pressure at midpoint first layer
    real(r8) :: pg      ! surface pressure
    real(r8) :: tc      ! temperature in celsius
    real(r8) :: es      ! saturation vapor pressure
    real(r8) :: ws      ! saturation mixing ratio
    real(r8) :: hvar    ! constant to compute xmol
    real(r8) :: h       ! constant to compute xmol
    real(r8) :: psih    ! stability correction factor
    real(r8) :: cts     ! correction to rlu rcl and rgs for frost
    real(r8) :: rs      ! constant for calculating rsmx
    real(r8) :: rmx     ! resistance by vegetation
    real(r8) :: zovl    ! ratio of z to  m-o length
    real(r8) :: cvarb   ! cvar averaged over landtypes
    real(r8) :: bb      ! b averaged over landtypes
    real(r8) :: ustarb  ! ustar averaged over landtypes

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(pcols,ndvel_gas) :: heff

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location only
    !-------------------------------------------------------------------------------------
    integer                :: index_season(pcols,n_land_type)
    real(r8), dimension(pcols) :: tha     ! atmospheric virtual potential temperature
    real(r8), dimension(pcols) :: thg     ! ground virtual potential temperature
    real(r8), dimension(pcols) :: z       ! height of lowest level
    real(r8), dimension(pcols) :: va      ! magnitude of v on cross points
    real(r8), dimension(pcols) :: ribn    ! richardson number
    real(r8), dimension(pcols) :: qs      ! saturation specific humidity
    real(r8), dimension(pcols) :: dewm    ! multiplier for rs when dew occurs
    real(r8), dimension(pcols) :: crs     ! multiplier to calculate crs
    real(r8), dimension(pcols) :: rdc     ! part of lower canopy resistance
    real(r8), dimension(pcols) :: uustar  ! u*ustar (assumed constant over grid)
    real(r8), dimension(pcols) :: z0b     ! average roughness length over grid
    real(r8), dimension(pcols) :: wrk     ! work array
    real(r8), dimension(pcols) :: term    ! work array
    real(r8), dimension(pcols) :: resc    ! work array
    real(r8), dimension(pcols) :: lnd_frc ! work array
    logical, dimension(pcols) :: unstable
    logical, dimension(pcols) :: has_rain

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location and landtype
    !-------------------------------------------------------------------------------------
    real(r8), dimension(pcols,n_land_type) :: rds   ! resistance for deposition of sulfate
    real(r8), dimension(pcols,n_land_type) :: b     ! buoyancy parameter for unstable conditions
    real(r8), dimension(pcols,n_land_type) :: cvar  ! height parameter
    real(r8), dimension(pcols,n_land_type) :: ustar ! friction velocity
    real(r8), dimension(pcols,n_land_type) :: xmol  ! monin-obukhov length

    !-------------------------------------------------------------------------------------
    ! local arrays: dependent on location, landtype and species
    !-------------------------------------------------------------------------------------
    real(r8), dimension(pcols,n_land_type,ndvel_gas) :: rsmx  ! vegetative resistance (plant mesophyll)
    real(r8), dimension(pcols,n_land_type,ndvel_gas) :: rclx  ! lower canopy resistance
    real(r8), dimension(pcols,n_land_type,ndvel_gas) :: rlux  ! vegetative resistance (upper canopy)
    real(r8), dimension(pcols,n_land_type) :: rlux_o3  ! vegetative resistance (upper canopy)
    real(r8), dimension(pcols,n_land_type,ndvel_gas) :: rgsx  ! ground resistance
    real(r8) :: pmid(pcols,1)                             ! for seasalt aerosols
    real(r8) :: tfld(pcols,1)                             ! for seasalt aerosols
!!$    real(r8) :: settling_velocity(pcols,1,nbin)           ! for seasalt aerosols
    real(r8) :: fact, vds
    real(r8) :: rc                                        ! combined surface resistance
    real(r8) :: var_soilw, dv_soil_h2, fact_h2            ! h2 dvel wrking variables
    logical :: fr_lnduse(pcols,n_land_type)           ! wrking array

    !-------------------------------------------------------------------------------------
    ! jfl : mods for PAN
    !-------------------------------------------------------------------------------------
    real(r8)                  :: dv_pan
    real(r8) :: c0_pan(11) = (/ 0.000_r8, 0.006_r8, 0.002_r8, 0.009_r8, 0.015_r8, &
                                0.006_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.002_r8, 0.002_r8 /)
    real(r8) :: k_pan (11) = (/ 0.000_r8, 0.010_r8, 0.005_r8, 0.004_r8, 0.003_r8, &
                                0.005_r8, 0.000_r8, 0.000_r8, 0.000_r8, 0.075_r8, 0.002_r8 /)

    !-------------------------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------------------------
    do m = 1,ndvel_gas
       dvel(:,m) = 0._r8
    end do

    if( all( .not. has_dvel(:) ) ) then
       return
    end if

    !-------------------------------------------------------------------------------------
    ! define species-dependent parameters (temperature dependent)
    !-------------------------------------------------------------------------------------
!    call set_hcoeff( sfc_temp, heff, ncol )
!   Since dry deposition is only calculated for SO2, it is assumed that the amount of water in a lake/oceans is enough to dissolve all of it.
! Heff is just given a large value. 
    heff(:,:)=1.e5_r8
    do lt = 1,n_land_type
       dep_ra (:,lt,lchnk)   = 0._r8
       dep_rb (:,lt,lchnk)   = 0._r8
       rds(:,lt)   = 0._r8
    end do

    !-------------------------------------------------------------------------------------
    ! 	... set month
    !-------------------------------------------------------------------------------------
    month = mod( ncdate,10000 )/100

    !-------------------------------------------------------------------------------------
    ! define which season (relative to Northern hemisphere climate)
    !-------------------------------------------------------------------------------------

    !-------------------------------------------------------------------------------------
    ! define season index based on fixed LAI
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       index_season(i,:) = index_season_lai(latndx(i),:,month)
    end do

    !-------------------------------------------------------------------------------------
    ! special case for snow covered terrain
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       if( snow(i) > .01_r8 ) then
          index_season(i,:) = 4
       end if
    end do
    !-------------------------------------------------------------------------------------
    ! scale rain and define logical arrays
    !-------------------------------------------------------------------------------------
    has_rain(:ncol) = rain(:ncol) > rain_threshold

    !-------------------------------------------------------------------------------------
    ! loop over longitude points
    !-------------------------------------------------------------------------------------
    col_loop :  do i = 1,ncol
       p   = pressure_10m(i)
       pg  = pressure_sfc(i)
       !-------------------------------------------------------------------------------------
       ! potential temperature
       !-------------------------------------------------------------------------------------
       tha(i) = air_temp(i) * (p00/p )**rovcp * (1._r8 + .61_r8*spec_hum(i))
       thg(i) = sfc_temp(i) * (p00/pg)**rovcp * (1._r8 + .61_r8*spec_hum(i))
       !-------------------------------------------------------------------------------------
       ! height of 1st level
       !-------------------------------------------------------------------------------------
       z(i) = - r/grav * air_temp(i) * (1._r8 + .61_r8*spec_hum(i)) * log(p/pg)
       !-------------------------------------------------------------------------------------
       ! wind speed
       !-------------------------------------------------------------------------------------
       va(i) = max( .01_r8,wind_speed(i) )
       !-------------------------------------------------------------------------------------
       ! Richardson number
       !-------------------------------------------------------------------------------------
       ribn(i) = z(i) * grav * (tha(i) - thg(i))/thg(i) / (va(i)*va(i))
       ribn(i) = min( ribn(i),ric )
       unstable(i) = ribn(i) < 0._r8
       !-------------------------------------------------------------------------------------
       ! saturation vapor pressure (Pascals)
       ! saturation mixing ratio
       ! saturation specific humidity
       !-------------------------------------------------------------------------------------
       es    = 611._r8*exp( 5414.77_r8*(sfc_temp(i) - tmelt)/(tmelt*sfc_temp(i)) )
       ws    = .622_r8*es/(pg - es)
       qs(i) = ws/(1._r8 + ws)
       !-------------------------------------------------------------------------------------
       ! multiplier for rs if rain or dew
       !-------------------------------------------------------------------------------------
       dewm(i) = 1._r8
       if( qs(i) <= spec_hum(i) .or. has_rain(i) ) then
          dewm(i) = 3._r8
       end if
       !-------------------------------------------------------------------------------------
       ! no dew if < 0C, effect of frost later
       !-------------------------------------------------------------------------------------
       if( sfc_temp(i) < tmelt ) then
          dewm(i) = 1._r8
       end if
       !-------------------------------------------------------------------------------------
       ! constant in determining rs
       !-------------------------------------------------------------------------------------
       if( sfc_temp(i) > tmelt .and. sfc_temp(i) < 313.15_r8 ) then
          tc = sfc_temp(i) - tmelt
          crs(i) = (1._r8 + (200._r8/(solar_flux(i) + .1_r8))**2) * (400._r8/(tc*(40._r8 - tc)))
       else
          crs(i) = large_value
       end if
       !-------------------------------------------------------------------------------------
       ! rdc (lower canopy res)
       !-------------------------------------------------------------------------------------
       rdc(i) = 100._r8*(1. + 1000._r8/(solar_flux(i) + 10.))/(1. + 1000._r8*slope)
    end do col_loop

    !-------------------------------------------------------------------------------------
    ! 	... form working array
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i=1,ncol
          fr_lnduse(i,lt) = fraction_landuse(lonndx(i),latndx(i),lt) /= 0._r8
       enddo
    end do

    !-------------------------------------------------------------------------------------
    ! find grid averaged z0: z0bar (the roughness length) z_o=exp[S(f_i*ln(z_oi))]
    ! this is calculated so as to find u_i, assuming u*u=u_i*u_i
    !-------------------------------------------------------------------------------------
    z0b(:) = 0._r8
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             z0b(i) = z0b(i) + fraction_landuse(lonndx(i),latndx(i),lt) * log( z0(index_season(i,lt),lt) )
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! find the constant velocity uu*=(u_i)(u*_i)
    !-------------------------------------------------------------------------------------
    do i = 1,ncol
       z0b(i) = exp( z0b(i) )
       cvarb  = vonkar/log( z(i)/z0b(i) )
       !-------------------------------------------------------------------------------------
       ! unstable and stable cases
       !-------------------------------------------------------------------------------------
       if( unstable(i) ) then
          bb = 9.4_r8*(cvarb**2)*sqrt( abs(ribn(i))*z(i)/z0b(i) )
          ustarb = cvarb * va(i) * sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*bb)) )
       else
          ustarb = cvarb * va(i)/(1._r8 + 4.7_r8*ribn(i))
       end if
       uustar(i) = va(i)*ustarb
    end do

    !-------------------------------------------------------------------------------------
    ! calculate the friction velocity for each land type u_i=uustar/u*_i
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( unstable(i) ) then
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)* sqrt( abs(ribn(i))*z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)*sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8 + 7.4_r8*b(i,lt))) ) )
             else
                cvar(i,lt)  = vonkar/log( z(i)/z0(index_season(i,lt),lt) )
                ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
             end if
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! revise calculation of friction velocity and z0 over water
    !-------------------------------------------------------------------------------------
    lt = 7    
    do i = 1,ncol
       if( fr_lnduse(i,lt) ) then
          if( unstable(i) ) then
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log( z(i)/z0water ))
             b(i,lt)     = 9.4_r8*(cvar(i,lt)**2)*sqrt( abs(ribn(i))*z(i)/z0water )
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)* sqrt( 1._r8 - (9.4_r8*ribn(i)/(1._r8+ 7.4_r8*b(i,lt))) ) )
          else
             z0water     = (.016_r8*(ustar(i,lt)**2)/grav) + diffk/(9.1_r8*ustar(i,lt))
             cvar(i,lt)  = vonkar/(log(z(i)/z0water))
             ustar(i,lt) = sqrt( cvar(i,lt)*uustar(i)/(1._r8 + 4.7_r8*ribn(i)) )
          end if
       end if
    end do

    !-------------------------------------------------------------------------------------
    ! compute monin-obukhov length for unstable and stable conditions/ sublayer resistance
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             hvar = (va(i)/0.74_r8) * (tha(i) - thg(i)) * (cvar(i,lt)**2)
             if( unstable(i) ) then                      ! unstable
                h = hvar*(1._r8 - (9.4_r8*ribn(i)/(1._r8 + 5.3_r8*b(i,lt))))
             else
                h = hvar/((1._r8+4.7_r8*ribn(i))**2)
             end if
             xmol(i,lt) = thg(i) * ustar(i,lt) * ustar(i,lt) / (vonkar * grav * h)
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! psih
    !-------------------------------------------------------------------------------------
    do lt = 1,n_land_type
       do i = 1,ncol
          if( fr_lnduse(i,lt) ) then
             if( xmol(i,lt) < 0._r8 ) then
                zovl = z(i)/xmol(i,lt)
                zovl = max( -1._r8,zovl )
                psih = exp( .598_r8 + .39_r8*log( -zovl ) - .09_r8*(log( -zovl ))**2 )
                vds  = 2.e-3_r8*ustar(i,lt) * (1._r8 + (300/(-xmol(i,lt)))**0.666_r8)
             else
                zovl = z(i)/xmol(i,lt)
                zovl = min( 1._r8,zovl )
                psih = -5._r8 * zovl
                vds  = 2.e-3_r8*ustar(i,lt)
             end if
             dep_ra (i,lt,lchnk) = (vonkar - psih*cvar(i,lt))/(ustar(i,lt)*vonkar*cvar(i,lt))
             dep_rb (i,lt,lchnk) = (2._r8/(vonkar*ustar(i,lt))) * crb
             rds(i,lt) = 1._r8/vds
          end if
       end do
    end do

    !-------------------------------------------------------------------------------------
    ! surface resistance : depends on both land type and species
    ! land types are computed seperately, then resistance is computed as average of values
    ! following wesely rc=(1/(rs+rm) + 1/rlu +1/(rdc+rcl) + 1/(rac+rgs))**-1
    !
    ! compute rsmx = 1/(rs+rm) : multiply by 3 if surface is wet
    !-------------------------------------------------------------------------------------
    species_loop1 :  do ispec = 1,ndvel_gas
       if( has_dvel(ispec) ) then
          do lt = 1,n_land_type
             do i = 1,ncol
                if( fr_lnduse(i,lt) ) then
                   sndx = index_season(i,lt)
                   if( ispec == o3_ndx .or. ispec == so2_ndx ) then
                      rmx = 0._r8
                   else
                      rmx = 1._r8/(heff(i,ispec)/3000._r8 + 100._r8*foxd(ispec))
                   end if
                   cts = 1000._r8*exp( -max( 263.15_r8,sfc_temp(i) ) + 269.15_r8 )                 ! correction for frost
                   rgsx(i,lt,ispec) = cts + 1._r8/((heff(i,ispec)/(1.e5_r8*rgss(sndx,lt))) + (foxd(ispec)/rgso(sndx,lt)))
                   !-------------------------------------------------------------------------------------
                   ! special case for H2 and CO;; CH4 is set ot a fraction of dv(H2)
                   !-------------------------------------------------------------------------------------
                   if( ispec == h2_ndx .or. ispec == co_ndx .or. ispec == ch4_ndx ) then
                      if( ispec == co_ndx ) then
                         fact_h2 = 1.0_r8
                      elseif ( ispec == h2_ndx ) then
                         fact_h2 = 0.5_r8
                      elseif ( ispec == ch4_ndx ) then
                         fact_h2 = 50.0_r8
                      end if
                      !-------------------------------------------------------------------------------------
                      ! no deposition on snow, ice, desert, and water
                      !-------------------------------------------------------------------------------------
                      if( lt == 1 .or. lt == 7 .or. lt == 8 .or. sndx == 4 ) then
                         rgsx(i,lt,ispec) = large_value
                      else
                         var_soilw = max( .1_r8,min( soilw(i),.3_r8 ) )
                         if( lt == 3 ) then
                            var_soilw = log( var_soilw )
                         end if
                         dv_soil_h2 = h2_c(lt) + var_soilw*(h2_b(lt) + var_soilw*h2_a(lt))
                         if( dv_soil_h2 > 0._r8 ) then
                            rgsx(i,lt,ispec) = fact_h2/(dv_soil_h2*1.e-4_r8)
                         end if
                      end if
                   end if
                   if( lt == 7 ) then
                      rclx(i,lt,ispec) = large_value
                      rsmx(i,lt,ispec) = large_value
                      rlux(i,lt,ispec) = large_value
                   else
                      rs = ri(sndx,lt)*crs(i)
                      rsmx(i,lt,ispec) = (dewm(i)*rs*drat(ispec) + rmx)
                      !-------------------------------------------------------------------------------------
                      ! jfl : special case for PAN
                      !-------------------------------------------------------------------------------------
                      if( ispec == pan_ndx ) then
                         dv_pan =  c0_pan(lt) * (1._r8 - exp( -k_pan(lt)*(dewm(i)*rs*drat(ispec))*1.e-2_r8 ))
                         if( dv_pan > 0._r8 .and. sndx /= 4 ) then
                            rsmx(i,lt,ispec) = ( 1._r8/dv_pan )
                         end if
                      end if
                      rclx(i,lt,ispec) = cts + 1._r8/((heff(i,ispec)/(1.e5_r8*rcls(sndx,lt))) + (foxd(ispec)/rclo(sndx,lt)))
                      rlux(i,lt,ispec) = cts + rlu(sndx,lt)/(1.e-5_r8*heff(i,ispec) + foxd(ispec))
                   end if
                end if
             end do
          end do
       end if
    end do species_loop1


    species_loop2 : do ispec = 1,ndvel_gas
       if( has_dvel(ispec) ) then
          if( ispec /= o3_ndx .and. ispec /= so2_ndx ) then
             do lt = 1,n_land_type
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( dewm(i) >= 2.999_r8 ) then
!!$                               rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
!!$                                    + 1.e-7_r8*heff(i,ispec) + foxd(ispec)/rlux(i,lt,o3_ndx))
                                rlux(i,lt,ispec) = 1._r8/((1._r8/(3._r8*rlux(i,lt,ispec))) &
                                                 + 1.e-7_r8*heff(i,ispec) + foxd(ispec)/rlux_o3(i,lt))
                            end if
                         end if
                      end if
                   end do
                end if
             end do
          else if( ispec == so2_ndx ) then
             do lt = 1,n_land_type
                if( lt /= 7 ) then
                   do i = 1,ncol
                      if( fr_lnduse(i,lt) ) then
                         !-------------------------------------------------------------------------------------
                         ! no effect if sfc_temp < O C
                         !-------------------------------------------------------------------------------------
                         if( sfc_temp(i) > tmelt ) then
                            if( qs(i) <= spec_hum(i) ) then
                               rlux(i,lt,ispec) = 100._r8
                            end if
                            if( has_rain(i) ) then
                               !                               rlux(i,lt,ispec) = 1._r8/(2.e-4_r8 + (1._r8/(3._r8*rlu(index_season(i,lt),lt))))
                               rlux(i,lt,ispec) = 15._r8*rlu(index_season(i,lt),lt)/(5._r8 + 3.e-3_r8*rlu(index_season(i,lt),lt))
                            end if
                         end if
                      end if
                   end do
                end if
             end do
             do i = 1,ncol
                if( fr_lnduse(i,1) .and. dewm(i) >= 3._r8 ) then
                   rlux(i,1,ispec) = 50._r8
                end if
             end do
          end if
       end if
    end do species_loop2

    !-------------------------------------------------------------------------------------
    ! compute rc
    !-------------------------------------------------------------------------------------
    term(:ncol) = 1.e-2_r8 * pressure_10m(:ncol) / (r*tv(:ncol))

!!$    if( has_seasalt ) then
!!$       pmid(:ncol,1) = pressure_sfc(:ncol)
!!$       tfld(:ncol,1) = sfc_temp(:ncol)
!!$       call seasalt_sett_vel( pmid, tfld, rh, settling_velocity )
!!$    end if
    species_loop3 : do ispec = 1,ndvel_gas
       if( has_dvel(ispec) ) then

          wrk(:) = 0._r8
          do lt = 1,n_land_type

             do i = 1,ncol
                if (fr_lnduse(i,lt)) then


!                resc(i)=1._r8/(1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))
!#ifdef feil
                   resc(i) = 1._r8/(1._r8/rsmx(i,lt,ispec) + 1._r8/rlux(i,lt,ispec) &
                           + 1._r8/(rdc(i) + rclx(i,lt,ispec)) &
                           + 1._r8/(rac(index_season(i,lt),lt) + rgsx(i,lt,ispec)))
!#endif
                   resc(i) = max( 10._r8,resc(i) )
                   lnd_frc(i) = fraction_landuse(lonndx(i),latndx(i),lt)
                endif
             enddo
             !-------------------------------------------------------------------------------------
             ! 	... compute average deposition velocity
             !-------------------------------------------------------------------------------------
!os             select case( solsym(ispec) )
!             case( 'SO2' )

                if( lt == 7 ) then
                   do i=1,ncol
                      if( fr_lnduse(i,lt) ) then

!
                        wrk(i) = wrk(i) + lnd_frc(i)/(dep_ra(i,lt,lchnk) + dep_rb(i,lt,lchnk) + resc(i))

                      endif
                    end do
                else
                   do i = 1,ncol
                      if( index_season(i,lt) == 4 ) then
                         if( fr_lnduse(i,lt) ) then
                            wrk(i) = wrk(i) + lnd_frc(i) * 0.1e-2_r8
                         end if
                      else
                         if( fr_lnduse(i,lt) ) then
                            wrk(i) = wrk(i) + lnd_frc(i) * 0.6e-2_r8
                         end if
                      end if

                   end do
                end if
!             case( 'SO4' )
!                where( fr_lnduse(:,lt) )
!                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + rds(:,lt))
!                endwhere
!             case( 'NH4', 'NH4NO3' )
!                where( fr_lnduse(:,lt) )
!                   wrk(:) = wrk(:) + lnd_frc(:)/(dep_ra(:ncol,lt,lchnk) + 0.5_r8*rds(:,lt))
!                endwhere
!             case( 'SA1', 'SA2', 'SA3', 'SA4' )
! move this to progseasalt_intr module -- let dvel be zero here
!!$                if( ispec == sa1_ndx ) then
!!$                   m = 1
!!$                else if( ispec == sa2_ndx ) then
!!$                   m = 2
!!$                else if( ispec == sa3_ndx ) then
!!$                   m = 3
!!$                else if( ispec == sa4_ndx ) then
!!$                   m = 4
!!$                end if
!!$                where( fr_lnduse(:,lt) )
!!$                   wrk(:) = wrk(:) + lnd_frc(:) &
!!$                        *(1._r8/(dep_ra(:,lt)+dep_rb(:,lt)+dep_ra(:,lt)*dep_rb(:,lt)*settling_velocity(:,1,m)) &
!!$                        + settling_velocity(:,1,m)) ! here settling_velocity is defined for one lev - bottom lev
!!$                endwhere
                !-------------------------------------------------------------------------------------
                ! 	... special case for Pb (for consistency with offline code)
                !-------------------------------------------------------------------------------------
!             case( 'Pb' )
!                if( lt == 7 ) then
!                   where( fr_lnduse(:,lt) )
!                      wrk(:) = wrk(:) + lnd_frc(:) * 0.05e-2_r8
!                   endwhere
!                else
!                   where( fr_lnduse(:ncol,lt) )
!                      wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol) * 0.2e-2_r8
!                   endwhere
!                end if
!             case default
!                where( fr_lnduse(:ncol,lt) )
!                   wrk(:ncol) = wrk(:ncol) + lnd_frc(:ncol)/(dep_ra(:ncol,lt,lchnk) + dep_rb(:ncol,lt,lchnk) + resc(:ncol))
!                endwhere
!             end select
          end do
!          select case( trim( solsym(ispec) ) )
!          case( 'CB1', 'CB2', 'OC1', 'OC2' )
!             wrk(:ncol) = 0.10e-2_r8
!          end select

       end if
! OS Uses m/s in the flux routine so the scaling facctor is removed
          dvel(:ncol,ispec) = wrk(:ncol) 
!* scaling_to_cm_per_s
!          dflx(:ncol,ispec) = term(:ncol) * dvel(:ncol,ispec) * mmr(:ncol,plev,ispec)



    end do species_loop3

    !-------------------------------------------------------------------------------------
    ! 	... special adjustments
    !-------------------------------------------------------------------------------------
!    if( mpan_ndx > 0 ) then
!       if( has_dvel(mpan_ndx) ) then
!          dvel(:ncol,mpan_ndx) = dvel(:ncol,mpan_ndx)/3._r8
!          dflx(:ncol,mpan_ndx) = term(:ncol) * dvel(:ncol,mpan_ndx) * mmr(:ncol,plev,mpan_ndx)
!       end if
!    end if

  end subroutine dvel_so2

end module drydep_so2
