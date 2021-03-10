! File :         p_prep_obs.F90
!
!                Created: unknown
!
! Author:        unknown
!
! Purpose:       Read data from different data sources and convert it to
!                type(measurement).
!
! Description:   The code calls different subroutines for particular types of
!                input data, depending on the source, observation type and
!                format. The output is the array of type(measurement) that
!                contains, in particular, pivot points and bilinear
!                interpolation coefficients for each observation. This array
!                is written to "observations.uf" in binary format.
!
! Modifications: 30/01/2008 - Pavel Sakov gave a trim (formatted) -- sorry
!                  for that, could not stand -- and modified to allow in-situ
!                  argo data from ifremer.
!                02/09/2008 - Pavel Sakov added superobing for SST and SLA data
!                17/08/2010 PS - turned (3D) superobing on for Argo obs
!                02/03/2016 Yiguo WANG: available for EN4 dataset 
!                04/08/2016 Yiguo WANG modified m_read_CLS_SSH
!
program p_prep_obs
  use mod_measurement
  use mod_grid
  use m_read_micom_SST
  use m_read_micom_SSH
  use m_read_CLS_header
  use m_read_CLS_data
  use m_read_CLS_SST_grid
  use m_read_MET_SST_grid
  use m_read_CLS_TSLA_grid
  use m_read_CLS_SST
  use m_read_CLS_SSH
  use m_read_CLS_SLA
  use m_read_CLS_TSLA
  use m_read_MET_SST
  use m_read_HadI_SST
  use m_read_NOAA_SST
  use m_read_CERSAT_data
  use m_read_ifremer_argo
  use m_read_EN4_profile
  use m_read_FFI_glider
  use m_read_metno_icec
  use m_get_def_wet_point
  use m_write_wet_file
  use m_get_micom_grid
  use m_parse_blkdat
  use m_read_amsr_norsex
  use m_superobs
  use m_uobs
  use m_get_micom_dim
  implicit none

  integer, parameter :: STRLEN = 512

  type (measurement), allocatable :: data(:)
  type (measurement), allocatable :: obs(:)
  type (grid) :: gr

  integer :: nx, ny, nz
  real, allocatable, dimension(:,:) :: depths, modlat, modlon,mask
  integer, parameter :: maxobs = 5000000
  character(STRLEN) :: fname, fnamehdr, dataformat, producer
  character(len=3) :: form
  character(len=5) :: obstype

  integer :: nrobs
  integer :: grpoints, k
  real :: factor, var
  real :: rdummy, mindx, meandx

  logical :: data_eq_obs

  ! superobs
  logical :: dosuperob
  logical :: is3d
  integer :: nrsobs
  type(measurement), allocatable :: sobs(:)

  integer :: i
  integer :: nthisobs
  integer, allocatable, dimension(:) :: thisobs

  gr = default_grid
  data_eq_obs = .false.

  open(10, file = 'infile.data')
  read(10, '(a)') producer
  read(10, '(a)') obstype
  read(10, '(a)') fnamehdr
  read(10, '(a)') fname
  close(10)

  print *, 'Data producer: ', trim(producer)
  print *, 'Data to be processed for NorCPM: ', trim(obstype)
  print *, 'Observational error variance: ', trim(fnamehdr)
  print *, 'Filenames to be processed: ', trim(fname)
  print *, 'Result of processing is stored in temporary file "observation.uf"'

  ! Get grid dimensions
  call get_micom_dim(nx,ny,nz)
  !
  allocate(depths(nx, ny))
  allocate(mask(nx, ny))
  allocate(modlon(nx, ny))
  allocate(modlat(nx, ny))
  ! Read position and depth from model grid
  !
  call  get_micom_grid(modlon, modlat, depths, mindx, meandx, nx, ny)
 
  dosuperob = .false.
  is3d = .false.

  ! Fill the "data" array by calling subroutines specific for the producer
  ! and observation type
  !
  if (trim(producer) == 'Reynolds') then

     if (trim(obstype) == 'SST') then
        dosuperob = .true.
        call read_CLS_header(fnamehdr, gr, dataformat, form, factor, var)
        !grpoints = gr % nx * gr % ny
        grpoints = 64980
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CLS_data(fname, obstype, dataformat, gr, form, data, factor, var)
        print*, 'Reynolds- ', obstype, ' data has been scaled by a factor = ', factor  
     else
        stop 'ERROR: Reynolds only produce SST'
     endif

  else if (trim(Producer) == 'MET') then

     if (trim(obstype) == 'SST') then
        dosuperob = .true.
        call read_MET_SST_grid(fnamehdr, gr)
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_MET_SST(fname, gr, data)
     else
        stop 'ERROR: OSTIA (MET) only produces SST'
     endif
  else if (trim(Producer) == 'Had') then
     if (trim(obstype) == 'SST') then
        dosuperob = .false.
        !data has 1 degree resolution (300*180)
        !could have read it from the file but hardcoded from now
        gr%reg = .true.
        gr%order = 2
        gr%nx=360
        gr%ny=180
        gr%set=.true.
        grpoints = gr%nx*gr%ny
        allocate(data(grpoints))
        print * ,'Reading HadI SST data'
        call read_HadI_SST(fname,fnamehdr,data,modlon, modlat, depths,360,180, nrobs)
        allocate(obs(nrobs))
        data_eq_obs = .false.
        print *,'HadI Data read',nrobs
     else
        stop 'ERROR: (HadI) only produces SST'
     endif
  else if (trim(Producer) == 'NOAA') then
     if (trim(obstype) == 'SST') then
        dosuperob = .true.
        !data has 1 degree resolution (300*180)
        !could have read it from the file but hardcoded from now
        gr%reg = .true.
        gr%order = 2
        gr%nx=360
        gr%ny=180
        gr%set=.true.
        grpoints = gr%nx*gr%ny
        allocate(data(grpoints))
        print * ,'Reading NOAA SST data'
        call read_NOAA_SST(fname,fnamehdr,data,modlon, modlat, depths,360,180, nrobs)
        allocate(obs(nrobs))
        data_eq_obs = .false.
        print *,'NOAA Data read',nrobs
     else
        stop 'ERROR: (NOAA) only produces SST'
     endif
  else if (trim(Producer) == 'micom') then
     if (trim(obstype) == 'SST') then
        dosuperob = .false.
        mask(:,:)=0
        where(depths>0)mask=1;
        grpoints = sum(mask)
        print *,'Number of wet point',grpoints
        allocate(data(grpoints))
       !fnamehdr contains the month and fname the year file 
        call read_micom_SST(fname,fnamehdr,data,modlon, modlat, depths, nx, ny, nrobs)
        data_eq_obs = .true.
        print *,'Data read'
     elseif (trim(obstype) == 'SSH') then
        dosuperob = .false.
        mask(:,:)=0
        where(depths>0)mask=1;
        grpoints = sum(mask)
        print *,'Number of wet point',grpoints
        allocate(data(grpoints))
       !fnamehdr contains the month and fname the year file 
        call read_micom_SSH(fname,fnamehdr,data,modlon, modlat, depths, nx, ny)
        data_eq_obs = .true.
        print *,'Data read'
     else
        stop 'ERROR: micom only produces SST'
     endif
  else if (trim(Producer) == 'NSIDC-AMSR') then

     if (trim(obstype) == 'ICEC') then
        dosuperob = .true.
        call read_amsr_norsex(fname, gr, data, obstype)
        allocate (obs(maxobs))
     else
        print *, 'No ',obstype, ' data from:', Producer
        stop  'ERROR: p_prep_obs'
     endif

  else if (trim(Producer) == 'METNO') then

     if (trim(obstype) == 'ICEC') then
        dosuperob = .true.
        call read_metno_icec(fname, data, gr)
        allocate (obs(size(data)))
     else
        print *, 'There can be no ', obstype,' data from', Producer
        stop
     endif

  elseif (trim(producer) == 'CLS') then

     if (trim(obstype) == 'SLA') then
        dosuperob = .true.
        ! call read_CLS_SST_grid() here because SST data grid has the same
        ! structure
        call read_CLS_SST_grid(fnamehdr, gr)
        grpoints = gr % nx * gr % ny
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CLS_SLA(fname, gr, data)
  
     elseif (trim(obstype) == 'SSH') then
        dosuperob = .false.
        !data has 1 degree resolution (360*180)
        !could have read it from the file but hardcoded from now
        gr%reg = .true.
        gr%order = 2
        gr%nx=360
        gr%ny=180
        gr%set=.true.
        grpoints = gr%nx*gr%ny
        allocate(data(grpoints))
        print * ,'Reading CLS SSH data'
        call read_CLS_SSH(fname, data, modlon, modlat, depths, gr%nx, gr%ny, nrobs, nx, ny)
        data_eq_obs = .true.
        print *,'CLS Data read', nrobs

     elseif (trim(obstype) == 'SST') then
        dosuperob = .true.
        call read_CLS_SST_grid(fnamehdr, gr)
        grpoints = gr % nx * gr % ny
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CLS_SST(fname, gr, data)

     elseif (trim(obstype) == 'TSLA') then
        dosuperob = .true.
        call read_CLS_TSLA_grid(fnamehdr, gr)
        print *, 'read_CLS_TSLA_grid finished, total # of obs = ', gr % nx 
        grpoints = gr % nx 
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CLS_TSLA(fname,gr,data)
     else
        print *, 'data of type "', trim(obstype),'"  from producer "', producer, '" is not handled'
        stop  'ERROR: p_prep_obs'
     endif

  else if (trim(producer) == 'NSIDC') then
     if (trim(obstype) == 'ICEC') then
        dosuperob = .true.
        call read_CLS_header(fnamehdr, gr, dataformat, form, factor, var)
        grpoints = gr % nx * gr % ny
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CLS_data(fname, obstype, dataformat, gr, form, data, factor, var)
        print *, producer, obstype, 'data has been scaled by a factor = ', factor  
     else
        print *, 'no data of type "', trim(obstype),'"  from producer "', producer, '" is not handled'
        stop  'ERROR: p_prep_obs'
     endif

  else if (trim(producer) == 'CERSAT') then
     if (trim(obstype) == 'idrft') then
        dosuperob = .false.
        call read_CLS_header(fnamehdr, gr, dataformat, form, factor, var)
        grpoints = gr % nx ! NB - 2 vector components - irregular grid
        allocate(data(grpoints))
        allocate(obs(maxobs))
        call read_CERSAT_data(trim(fname), gr, data, grpoints, var)
        print *, producer, obstype, 'data has been scaled by a factor = ', factor  
     else
        print *, 'no data of type "', trim(obstype),'"  from producer "', producer, '" is not handled'
        stop 'ERROR: p_prep_obs'
     endif

  elseif (trim(producer) == 'IFREMER') then

     dosuperob = .true.
     is3d = .true.
     read(fnamehdr, *) var
     print *, 'variance =', var
     call read_ifremer_argo(fname, obstype, var, nx, ny, data)

     ! PS: This is a flag to denote that read_ifremer_argo() takes care of
     ! filling type(measurement) array "data" in a correct way, and it should
     ! not be re-processed by calling get_def_wet_point(). This may not match
     ! the ideology behind the workflow adopted in this module and may be
     ! changed in future.
     !
     data_eq_obs = .true.

  elseif (trim(producer) == 'FFI') then

     dosuperob = .true.
     is3d = .true.
     read(fnamehdr, *) var
     print *, 'variance =', var
     call read_FFI_glider(fname, obstype, var, nx, ny, data)
     data_eq_obs = .true.

  elseif (trim(producer) == 'EN4') then

     dosuperob = .true.
     is3d = .true.
     read(fnamehdr, *) var
     call read_EN4_profile(fname, obstype, var, nx, ny, data, nrobs)
     
     data_eq_obs = .true.
  else
     print *, 'unknown producer ', trim(producer), ' in "infile.data"'
     stop 'ERROR: p_prep_obs'
  endif


  if (.not. data_eq_obs) then
     ! Compute bilinear coefficients
     ! Extract the defined and wet data points
     ! Write locations to ijfile to be used in TECPLOT
     !
     call get_def_wet_point(obs, data, gr, depths, modlat, modlon, nrobs, nx, ny)
  else
     allocate(obs(nrobs))
     obs(:nrobs) = data

     print *,'Nb of data dumped is :',nrobs
  end if
  deallocate(data)

  print *,'Updated age ',minval(obs(:)%date),maxval(obs(:)%date)
  if (trim(obstype) == 'TSLA') then
     call set_re_TSLA(nrobs, obs, nx, ny, modlon, modlat)
  end if

  where (obs % d + 1.0 == obs % d)
     obs % status = .false.
  end where

  ! Superob dense data
  !
  if (dosuperob) then
     allocate(sobs(nrobs))
     print *,'Min age/max age',minval(obs(1:nrobs)%date),maxval(obs(1:nrobs)%date)
     call superob(obstype, nrobs, obs, nx, ny, modlon, modlat, nrsobs, sobs, is3d)
     
     deallocate(obs)
     allocate(obs(nrsobs))
     obs = sobs(1 : nrsobs)
     nrobs = nrsobs
     deallocate(sobs)
  end if

  if (nrobs .ge. maxobs) then
     print *, 'max No. of data reached, increase it!'
     stop 'ERROR: p_prep_obs'
  elseif (nrobs .le. 1) then
     print *, 'less than one observation in the whole dataset'
     !PS 4/9/2011 stop 'ERROR: p_prep_obs: Not worth the money'
  end if

  ! Write data to the binary file "observations.uf"
  !
  call write_wet_file(obs, nrobs)

  call uobs_get(obs(1 : nrobs) % id, nrobs, .true.)
  allocate(thisobs(nrobs))
  do i = 1, nuobs
     nthisobs = 0
     do k = 1, nrobs
        if (trim(unique_obs(i)) == trim(obs(k) % id)) then
           nthisobs = nthisobs + 1
           thisobs(nthisobs) = k
        end if
     end do

     if (nthisobs > 0) then
        call obs2nc(nthisobs, obs(thisobs(1 : nthisobs)))
     end if
  end do
  deallocate(thisobs)

  print *, 'Last observation:'
  print '(a)','   obs       var    id      lon   lat  depth   ipiv  jpiv   nsup'//&
         '  4-bilin-coeffs    active  orig (i,j)   dp    age orig_id'
  print '(2g10.2,a6,3f6.1,3i6,4f5.1,l5,2i7,f7.1,2i5)', obs(nrobs)

  deallocate(obs)
  deallocate(depths)
  deallocate(modlon)
  deallocate(modlat)

  print *, 'prep_obs: end of processing'
end program p_prep_obs


subroutine obs2nc(nobs, obs)
  use mod_measurement
  use nfw_mod
  implicit none

  integer, parameter :: STRLEN = 512

  integer, intent(in) :: nobs
  type(measurement), intent(in) :: obs(nobs)

  character(STRLEN) :: fname
  integer :: ncid, obsdimid(1), lon_id, lat_id, depth_id, d_id, var_id, age_id
  integer :: n_id, ipiv_id, jpiv_id
  integer :: n(nobs)

  ! Create netcdf file of observations
  !
  write(fname, '(a, a, a)') 'observations-', trim(obs(1) % id), '.nc'
  print *, 'dumping observations to "', trim(fname), '"'

  call nfw_create(fname, nf_clobber, ncid)

  call nfw_def_dim(fname, ncid, 'nobs', nobs, obsdimid(1))
  call nfw_def_var(fname, ncid, 'lon', nf_float, 1, obsdimid(1), lon_id)
  call nfw_def_var(fname, ncid,  'lat', nf_float, 1, obsdimid(1), lat_id)
  call nfw_def_var(fname, ncid, 'depth', nf_float, 1, obsdimid(1), depth_id)
  call nfw_def_var(fname, ncid, 'd', nf_float, 1, obsdimid(1), d_id)
  call nfw_def_var(fname, ncid, 'var', nf_float, 1, obsdimid(1), var_id)
  call nfw_def_var(fname, ncid, 'age', nf_int, 1, obsdimid(1), age_id)
  call nfw_def_var(fname, ncid, 'n', nf_int, 1, obsdimid(1), n_id)
  call nfw_def_var(fname, ncid, 'ipiv', nf_int, 1, obsdimid(1), ipiv_id)
  call nfw_def_var(fname, ncid, 'jpiv', nf_int, 1, obsdimid(1), jpiv_id)
  call nfw_enddef(fname, ncid)

  call nfw_put_var_double(fname, ncid, lon_id, obs(1:nobs) % lon)
  call nfw_put_var_double(fname, ncid, lat_id, obs(1:nobs) % lat)
  call nfw_put_var_double(fname, ncid, depth_id, obs(1:nobs) % depth)
  call nfw_put_var_double(fname, ncid, d_id, obs(1:nobs) % d)
  call nfw_put_var_double(fname, ncid, var_id, obs(1:nobs) % var)
  call nfw_put_var_int(fname, ncid, age_id, obs(1:nobs) % date)
  call nfw_put_var_int(fname, ncid, ipiv_id, obs(1:nobs) % ipiv)
  call nfw_put_var_int(fname, ncid, jpiv_id, obs(1:nobs) % jpiv)
  n = int(obs(1:nobs) % h)
  call nfw_put_var_int(fname, ncid, n_id, n)
  
  call nfw_close(fname, ncid)
end subroutine obs2nc


subroutine check_forland(data, depths, nrobs, ni, nj)
  use mod_measurement

  type (measurement), intent(inout), dimension(nrobs) :: data
  real, dimension(ni, nj), intent(in)  ::  depths
  integer, intent(in) :: nrobs
  integer, intent(in) :: ni, nj

  integer :: o, imin, jmin, imax, jmax, nmasked

  nmasked = 0
  do o = 1, nrobs
     imin = max(1, data(o) % ipiv - 1)
     jmin = max(1, data(o) % jpiv - 1)
     imax = min(ni, data(o) % ipiv + 2)
     jmax = min(nj, data(o) % jpiv + 2)
     if (any(depths(imin:imax,jmin:jmax) < 10.0 .or. depths(imin:imax,jmin:jmax) == depths(imin:imax,jmin:jmax) + 1.0)) then
        data(o) % status = .false.
        nmasked = nmasked + 1
     end if
  end do
  print *, "  check_forland(): ", nmasked, "obs close to land masked" 
end subroutine check_forland
