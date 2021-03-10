! File:          m_read_EN4_profile.F90
!
! Created:       17 July 2015
!
! Author:        Yiguo WANG (YW)
!                NERSC
!
! Purpose:       Read profile data from NetCDF files from EN4 into NorCPM
!                system.
!
! Description:   Data file(s) are defined by the string in the 4th line of
!                "infile.data". It should have the following format:
!                <BEGIN>
!                EN4
!                SAL | TEM 
!                <obs. error variance>
!                <File name(s) or a wildcard>
!                <END>
!                After that:
!                1. all profiles are read into two arrays,
!                   deph(1 : nlev, 1 : nprof) and v(1 : nlev, 1 : nprof), where
!                   nprof is the total number of profiles in all files, and
!                   nlev is the maximum number of horizontal levels for all
!                   profiles;
!                2. bad data with qc flags other than '1' is discarded;
!                3. dry or outside locations are discarded
!                4. if there close profiles (in the same grid cell), the best
!                   one (with most data or the most recent) is retained
!
!

module m_read_EN4_profile
  implicit none

  integer, parameter, private :: STRLEN = 512
  integer, parameter, private :: kdm = 42
  integer, parameter, private :: kdm1 = 42

  real, parameter, private :: DENS_DIFF_MIN = -0.02
  logical, parameter, private :: DISCARD_CLOSE = .false.

#ifdef ANOMALY
  ! anomaly
  real, parameter, private :: TEM_MIN = -6.0
  real, parameter, private :: TEM_MAX = 6.0
  real, parameter, private :: SAL_MIN = -3.0
  real, parameter, private :: SAL_MAX = 3.0
#else
  ! full field
  real, parameter, private :: TEM_MIN = -2.5
  real, parameter, private :: TEM_MAX = 35.0
  real, parameter, private :: SAL_MIN = 1.0
  real, parameter, private :: SAL_MAX = 50.0
#endif

  public read_EN4_profile

  private data_inquire
  private data_readfile
  private potential_density
  private get_pivot
  private data_variance
  private data_obsunc
contains

  subroutine read_EN4_profile(fnames, obstype, variance, nx, ny, data, nrobs)
    use mod_measurement
    use m_oldtonew
    use m_confmap
    use m_bilincoeff
    use m_pivotp
    use nfw_mod
    use ieee_arithmetic
    use m_get_micom_fld

    character(*), intent(in) :: fnames
    character(*), intent(in) :: obstype
    real(8), intent(in) :: variance
    integer, intent(in) :: nx, ny
    integer, intent(out) :: nrobs
    type(measurement), allocatable, intent(out) :: data(:)

    character(STRLEN) :: fname
    integer :: nfile, nprof, nlev

    real(8), allocatable :: juld(:)
    real(8), allocatable :: lat(:), lon(:)
    real(8), allocatable :: deph(:,:)
    real(8), allocatable :: temp(:,:), salt(:, :)
    
    character, allocatable :: pos_qc(:)
    character, allocatable :: temp_qc(:,:), salt_qc(:, :)
    character, allocatable :: prof_temp_qc(:),prof_saln_qc(:)

    integer, allocatable :: ipiv(:), jpiv(:)

    real(8), dimension(nx, ny) :: modlat, modlon
    real(8), dimension(nx, ny) :: depths, mxd1, mxd2, fice

    real(8), dimension(nx, ny, kdm1) :: obs_unc
    real(8), dimension(2, kdm1) :: d3z 

    real(8), dimension(360, 173, kdm) :: obs_mean
    real(8), dimension(360) :: dlon
    real(8), dimension(173) :: dlat
    real(8), dimension(kdm) :: ddepth
    real(8), allocatable, dimension(:) :: splinex, spliney, splined
    real(8), allocatable, dimension(:) :: splinexe, splineye
    integer :: lon_int, lat_int 

    integer :: f, l, p, np, k, ll
    integer :: ipp1,ipm1,jpp1,jpm1
    integer, allocatable :: mask(:)
    integer, allocatable :: mask2(:, :)
    integer, allocatable :: fid(:);
    integer, allocatable :: profid(:)
    integer, allocatable :: done(:)
    real(8) :: zmax, Q, Qbest, rho, rho_prev, rho_inc
    integer :: best
    integer :: p1
    
    integer ngood, ndata
    real(8) :: latnew, lonnew
    
    print *, 'BEGIN read_EN4_profile()'

    call data_inquire(fnames, nfile, nprof, nlev)
    print *, '  overall: nprof =', nprof, ', nlev =', nlev

    allocate(juld(nprof))
    allocate(lat(nprof))
    allocate(lon(nprof))
    allocate(temp(nlev, nprof))
    allocate(salt(nlev, nprof))
    allocate(deph(nlev, nprof))
    
    allocate(pos_qc(nprof))
    allocate(temp_qc(nlev, nprof))
    allocate(salt_qc(nlev, nprof))
    allocate(prof_temp_qc(nprof))
    allocate(prof_saln_qc(nprof))

    allocate(fid(nprof))
    allocate(profid(nprof))

    ! read pre-estimated obs uncertainties
    !
    fname = './obs_unc_'//trim(obstype)//'.nc'
    call data_obsunc(trim(fname), trim(obstype), d3z, obs_unc)
    
    ! read data
    !
    p = 1
    do f = 1, nfile
       call data_readfile(f, np, juld(p : nprof),&
            lat(p : nprof), lon(p : nprof), pos_qc(p : nprof), &
            temp(1 : nlev, p : nprof), temp_qc(1 : nlev, p : nprof), &
            salt(1 : nlev, p : nprof), salt_qc(1 : nlev, p : nprof), &
            deph(1 : nlev, p : nprof), prof_temp_qc(p : nprof), &
            prof_saln_qc(p : nprof))
       fid(p : p + np - 1) = f
       do l = 1, np
          profid(p + l - 1) = l
       end do
       p = p + np
    end do

    allocate(mask(nprof))
    mask(:) = 1
    allocate(mask2(nlev, nprof))
    mask2(:, :) = 1

#ifdef ANOMALY
    ! read climatology
    fname = 'mean_obs.nc'
    call data_obsmean(trim(fname), trim(obstype), ddepth, dlon, dlat, obs_mean)

    do p = 1, nprof
       if ((lat(p) .lt. minval(dlat)) .or. (lat(p) .gt. maxval(dlat))) then
          mask(p) = 0
          mask2(:, p) = 0
          cycle
       end if
       
       ! identify coordinate in obs_mean
       lon_int = modulo(nint(lon(p)), 360)
       if (lon_int .eq. 0) lon_int = 360
       lat_int = nint(lat(p)) + 84
       
       ! find available depth for obs_mean
       f = count(obs_mean(lon_int, lat_int, :) .gt. -10000.)
       if (f .lt. 2) then
          mask(p) = 0
          mask2(:, p) = 0
          cycle
       end if
       
       allocate(splinex(f), spliney(f),splined(f))
       l = 0
       do k = 1, kdm
          if (obs_mean(lon_int, lat_int, k) .gt. -10000.) then 
             l = l + 1
             splinex(l) = ddepth(k)
             spliney(l) = obs_mean(lon_int, lat_int, k)
          end if
       end do
       if (l .ne. f) print *, 'Error in spline: l /= f'
       call spline_pchip_set(f, splinex, spliney, splined)

       ! find available depths in obs
       f = count(deph(:, p) .lt. 10000.)
       if (f .lt. 1) then
          mask(p) = 0
          mask2(:, p) = 0
          deallocate(splinex, spliney, splined)
          cycle
       end if

       allocate(splinexe(f), splineye(f))       
       ll = 0
       do k = 1, nlev
          if (deph(k, p) .lt. 10000.) then
             ll = ll + 1
             splinexe(ll) = deph(k, p)
          end if
       end do
       if (ll .ne. f) print *, 'Error in spline: ll /= f'
       call spline_pchip_val(l, splinex, spliney, splined, &
            ll, splinexe, splineye)

       ! calcule anomaly
       ll = 0
       if (trim(obstype) == 'SAL') then
          do k = 1, nlev
             if (deph(k, p) .lt. maxval(splinex)) then
                ll = ll + 1
                salt(k, p) = salt(k, p) - splineye(ll)
             else if (deph(k, p) .lt. 10000.) then
                ll = ll + 1
                mask2(k, p) = 0
             else
                mask2(k, p) = 0
             end if
          end do
       else if (trim(obstype) == 'TEM') then
          do k = 1, nlev
             if (deph(k, p) .lt. maxval(splinex)) then
                ll = ll + 1
                temp(k, p) = temp(k, p) - splineye(ll)
             else if (deph(k, p) .lt. 10000.) then
                ll = ll + 1
                mask2(k, p) = 0
             else
                mask2(k, p) = 0
             end if
          end do
       else
          print *, 'ERROR: Anomaly: <obstype> should be <TEM> or <SAL>...'
          stop
       end if
       deallocate(splinex, spliney, splined, splinexe, splineye)
    end do
#endif

    ! mask <- pos_qc
    !
    where (pos_qc .ne. '1') mask = 0
    where ((lat .lt. -90) .or. (lat .gt. 89)) mask = 0
    where ((lon .lt. -180) .or. (lon .gt. 180)) mask = 0
    do p = 1, nprof
       if (mask(p) == 0) then
          mask2(:, p) = 0
       end if
    end do
    print *, '  after examining POS_QC:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    ! ipiv, jpiv
    !
    allocate(ipiv(nprof))
    allocate(jpiv(nprof))
    ipiv(:) = -999
    jpiv(:) = -999
    call get_pivot(nx, ny, nprof, lon, lat, ipiv, jpiv, depths, modlon, modlat)
    where (ipiv < 1 .or. jpiv < 1 .or. ipiv > nx - 1 .or. jpiv > ny - 1) mask = 0

    do p = 1, nprof
       if (ipiv(p) < 1 .or. jpiv(p) < 1) then 
          print *, 'negatif values for ipiv, jpiv'
          print *, 'ipiv=', ipiv
          print *, 'jpiv=', jpiv
          mask = 0
       end if
    end do

    do p = 1, nprof
       if (mask(p) == 0) then
          mask2(:, p) = 0
       end if
    end do
    print *, '  after calculating pivot points:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    ! Check land grid    
    !
    do p = 1, nprof
#ifdef MASK_LANDNEIGHBOUR
       ipm1=max(ipiv(p)-1,1)
       ipp1=min(ipiv(p)+1,nx)
       jpm1=max(jpiv(p)-1,1)
       jpp1=min(jpiv(p)+1,ny)
       if (any(depths(ipm1:ipp1, jpm1:jpp1) < 60) ) mask(p) = 0
#endif
       if (depths(ipiv(p), jpiv(p)) < 60) mask(p) = 0
    end do
    do p = 1, nprof
       if (mask(p) == 0) then
          mask2(:, p) = 0
       end if
    end do
    print *, '  after examinating model land points:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    ! Check for the observation being wet                                                             
    !                                           
    do p = 1, nprof
       if (mask(p) == 0) then
          cycle
       end if
       do l = 1, nlev
          if (mask2(l, p) == 0) then
             cycle
          end if
          if (deph(l, p) > depths(ipiv(p), jpiv(p)) .or.&
               deph(l, p) > depths(ipiv(p) + 1, jpiv(p)) .or.&
               deph(l, p) > depths(ipiv(p), jpiv(p) + 1) .or.&
               deph(l, p) > depths(ipiv(p) + 1, jpiv(p) + 1)) then
             mask2(l, p) = 0
          end if
       end do
       if (count(mask2(:, p) == 1) == 0) then
          mask(p) = 0
       end if
    end do
    print *, '  after examining for wet cells:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    ! Check for the observation above mixed layer
    !
    fname = 'forecast001'
    call get_micom_fld_new(trim(fname), mxd1, 'dp', 1, 1, nx, ny)
    call get_micom_fld_new(trim(fname), mxd2, 'dp', 2, 1, nx, ny)
    mxd1 = (mxd1 + mxd2) / 98060. ! Unit: meter
    do p = 1, nprof
       if (mask(p) == 0) then
          cycle
       end if
       do l = 1, nlev
          if (mask2(l, p) == 0) then
             cycle
          end if
          if (deph(l, p) < mxd1(ipiv(p), jpiv(p))) then
             mask2(l, p) = 0
          end if
       end do
       if (count(mask2(:, p) == 1) == 0) then
          mask(p) = 0
       end if
    end do
    print *, '  after examining for obs above the mixed layer:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    !
    ! Now examine 3D quality flags; set the mask for a profile to 0 if there
    ! are no good samples in this profile
    !

    ! <data>_qc
    !
    if (trim(obstype) == 'SAL') then
       where (prof_saln_qc .ne. '1') mask = 0
       do p = 1, nprof
          if (mask(p) == 0) then
             mask2(:, p) = 0
          end if
       end do

       do p = 1, nprof
          do l = 1, nlev
             if (salt_qc(l, p) /= '1') then
                mask2(l, p) = 0
             end if
             if ((salt(l, p) .lt. SAL_MIN) .or. (salt(l, p) .gt. SAL_MAX)) then
                mask2(l, p) = 0
             end if
          end do
          if (count(mask2(:, p) == 1) == 0) then
             mask(p) = 0
          end if
       end do
    else if (trim(obstype) == 'TEM') then
       where (prof_temp_qc .ne. '1') mask = 0
       do p = 1, nprof
          if (mask(p) == 0) then
             mask2(:, p) = 0
          end if
       end do

       do p = 1, nprof
          do l = 1, nlev
             if (temp_qc(l, p) /= '1') then
                mask2(l, p) = 0
             end if
             if ((temp(l, p) .lt. TEM_MIN) .or. (temp(l, p) .gt. TEM_MAX)) then
                mask2(l, p) = 0
             end if
          end do
          if (count(mask2(:, p) == 1) == 0) then
             mask(p) = 0
          end if
       end do
    end if
    print *, '  after examining prof_<data>_QC and <data>_QC:'
    print *, '    ', count(mask == 1), ' good profiles'
    print *, '    ', count(mask2 == 1), ' good obs'

    ! Finally, discard redundant observations
    ! This is a O(n^2) search, which can become a bit long when the number of
    ! examined profiles becomes really large (say, 10^4)
    !   
    if (DISCARD_CLOSE) then
       allocate(done(nprof))
       done = 0
       do p = 1, nprof
          if (mask(p) == 0 .or. done(p) == 1) then
             cycle
          end if
          np = 1
          profid(np) = p
          do p1 = p + 1, nprof
             if (ipiv(p1) == ipiv(p) .and. jpiv(p1) == jpiv(p)) then
                np = np + 1
                profid(np) = p1
                done(p1) = 1
             end if
          end do
          if (np > 1) then
             ! for each of close profiles find the depth range, number of points
             ! and the age
             Qbest = 0.0
             do p1 = 1, np
                zmax = 0.0
                ndata = 0
                do l = 1, nlev
                   if (mask2(l, profid(p1)) == 1) then
                      ndata = ndata + 1
                      if (deph(l, profid(p1)) > zmax) then
                         zmax =  deph(l, profid(p1))
                      end if
                   end if
                end do
                Q = min(zmax, 400.0) / 400.0 + min(ndata, 10) / 10
                if (Q > Qbest) then
                   best = p1
                end if
             end do
             do p1 = 1, np
                if (p1 == best) then
                   cycle
                end if
                mask(profid(p1)) = 0
                mask2(:, profid(p1)) = 0
             end do
          end if
       end do
       deallocate(done)
       print *, '  after discarding close profiles:'
       print *, '    ', count(mask == 1), ' good profiles'
       print *, '    ', count(mask2 == 1), ' good obs'
    end if ! DISCARD_CLOSE

    ! Read sea ice clim from model 
    fname = 'mean_mod'
    call get_micom_fld_new(trim(fname), fice, 'fice', 0, 1, nx, ny)

    ngood = count(mask2 == 1)
    nrobs = ngood
    allocate(data(ngood))
    ndata = 0
    do p = 1, nprof
       if (mask(p) == 0) then
          cycle
       end if
       do l = 1, nlev
          if (mask2(l, p) == 0) then
             cycle
          end if

          ndata = ndata + 1

          if (ndata > ngood) then
             print *, 'ERROR: read_EN4_profile(): programming error'
             print *, '  p =', p, ', l =', l
             print *, '  # data =', ndata, ', ngood =', ngood
             stop
          end if
       
          ! PS: I guess we should not bother about the cost of the
          ! comparisons below.
          !
          if (trim(obstype) == 'SAL') then
             data(ndata) % d = salt(l, p)
          else if (trim(obstype) == 'TEM') then
             data(ndata) % d = temp(l, p)
!             if (ipiv(p) == 18 .and. jpiv(p) == 137 .and. temp(l, p) < -5) print *, p
          else
             print *, 'ERROR: read_EN4_profile(): <obstype> should be <TEM> or <SAL>...'
             stop
          end if
          data(ndata) % id = obstype
          data(ndata) % lon = lon(p)
          data(ndata) % lat = lat(p)
          data(ndata) % depth = max(0.0, deph(l, p))
          if (variance > 0) then
             data(ndata) % var = variance
          else
             call data_variance(trim(obstype), data(ndata) % depth, data(ndata) % var)
          end if
          do k = 1,kdm1
             if (data(ndata) % depth >= d3z(1,k) .and. data(ndata) % depth < d3z(2,k)) then
                data(ndata) % var = max(data(ndata) % var, obs_unc(ipiv(p), jpiv(p), k))
                exit
             end if
          end do
          if (fice(ipiv(p), jpiv(p)) > 15.0) data(ndata) % var = 10.0 * data(ndata) % var 
          data(ndata) % ipiv = ipiv(p)
          data(ndata) % jpiv = jpiv(p)
          data(ndata) % ns = 0 ! for a point (not gridded) measurement
          data(ndata) % date = 0 ! assimilate synchronously

          data(ndata) % a1 = 1
          data(ndata) % a2 = 0
          data(ndata) % a3 = 0
          data(ndata) % a4 = 0
!          call bilincoeff1(modlon, modlat, nx, ny, lon(p), lat(p), ipiv(p),&
!               jpiv(p), data(ndata) % a1, data(ndata) % a2, data(ndata) % a3,&
!               data(ndata) % a4)

          data(ndata) % status = .true. ! (active)
          data(ndata) % i_orig_grid = p
          data(ndata) % j_orig_grid = l
          data(ndata) % orig_id = 0
          data(ndata) % h = 0
       end do
    end do

    if (ndata /= ngood) then
       print *, 'ERROR: read_EN4_profile(): programming error'
       print *, '  ndata =', ndata, ', ngood =', ngood
       stop
    end if

    deallocate(juld)
    deallocate(lat)
    deallocate(lon)
    deallocate(pos_qc)
    deallocate(fid)
    deallocate(profid)
    deallocate(temp)
    deallocate(salt)
    deallocate(temp_qc)
    deallocate(salt_qc)
    deallocate(deph)
    deallocate(prof_temp_qc)
    deallocate(prof_saln_qc)
    deallocate(mask)
    deallocate(mask2)
    deallocate(ipiv)
    deallocate(jpiv)

    print *, 'END read_EN4_profile()'

  end subroutine read_EN4_profile


  subroutine data_inquire(fnames, nfile, nprof, nlev)
    use nfw_mod

    character(*), intent(in) :: fnames
    integer, intent(inout) :: nfile, nprof, nlev

    character(STRLEN) :: command ! (there may be a limit of 80 on some systems)
    character(STRLEN) :: fname
    integer :: ios
    integer :: ncid
    integer :: id

    integer :: nprof_this, nlev_this

    nfile = 0
    nprof = 0
    nlev = 0

    command = 'ls '//trim(fnames)//' > infiles.txt'
    call system(command);

    nfile = 0
    open(10, file = 'infiles.txt')
    do while (.true.)
       read(10, fmt = '(a)', iostat = ios) fname
       if (ios /= 0) then
          exit
       end if

       nfile = nfile + 1
       print *, '  file #', nfile, ' = "', trim(fname), '"'

       call nfw_open(fname, nf_nowrite, ncid)

       ! nprof
       !
       call nfw_inq_dimid(fname, ncid, 'N_PROF', id)
       call nfw_inq_dimlen(fname, ncid, id, nprof_this)
       print *, '    nprof = ', nprof_this

       ! nlev
       !
       call nfw_inq_dimid(fname, ncid, 'N_LEVELS', id)
       call nfw_inq_dimlen(fname, ncid, id, nlev_this)
       print *, '    nlev = ', nlev_this
       
       nprof = nprof + nprof_this
       if (nlev_this > nlev) then
          nlev = nlev_this
       end if

       call nfw_close(fname, ncid)
    end do
    close(10)
  end subroutine data_inquire


  subroutine data_readfile(fid, nprof, juld_all, &
    lat_all, lon_all, pos_qc_all, temp_all, temp_qc_all, salt_all, salt_qc_all, &
    deph_all, prof_temp_qc_all, prof_saln_qc_all)
    use nfw_mod

    integer, intent(in) :: fid
    integer, intent(inout) :: nprof
    real(8), intent(inout), dimension(:) :: juld_all
    real(8), intent(inout), dimension(:) :: lat_all, lon_all
    character, intent(inout), dimension(:) :: pos_qc_all
    real(8), intent(inout), dimension(:,:) :: temp_all
    character, intent(inout), dimension(:,:) :: temp_qc_all
    real(8), intent(inout), dimension(:,:) :: salt_all
    character, intent(inout), dimension(:,:) :: salt_qc_all
    real(8), intent(inout), dimension(:,:) :: deph_all
    character, intent(inout), dimension(:) :: prof_temp_qc_all
    character, intent(inout), dimension(:) :: prof_saln_qc_all

    character(STRLEN) :: fname
    integer :: f
    integer :: ncid
    integer :: id
    integer :: nlev
    
    open(10, file = 'infiles.txt')
    do f = 1, fid
       read(10, fmt = '(a)') fname
    end do
    close(10)

    print *, '  reading "', trim(fname), '"'
    
    call nfw_open(fname, nf_nowrite, ncid)

    ! nprof
    !
    call nfw_inq_dimid(fname, ncid, 'N_PROF', id)
    call nfw_inq_dimlen(fname, ncid, id, nprof)

    ! nlev
    !
    call nfw_inq_dimid(fname, ncid, 'N_LEVELS', id)
    call nfw_inq_dimlen(fname, ncid, id, nlev)

    ! juld
    !
    call nfw_inq_varid(fname, ncid, 'JULD', id)
    call nfw_get_var_double(fname, ncid, id, juld_all(1 : nprof))

    ! lat
    !
    call nfw_inq_varid(fname, ncid, 'LATITUDE', id)
    call nfw_get_var_double(fname, ncid, id, lat_all(1 : nprof))

    ! lon
    !
    call nfw_inq_varid(fname, ncid, 'LONGITUDE', id)
    call nfw_get_var_double(fname, ncid, id, lon_all(1 : nprof))

    ! pos_qc
    !
    call nfw_inq_varid(fname, ncid, 'POSITION_QC', id)
    call nfw_get_var_text(fname, ncid, id, pos_qc_all(1 : nprof))

    ! temp
    !
    call nfw_inq_varid(fname, ncid, 'POTM_CORRECTED', id)
    call nfw_get_var_double(fname, ncid, id, temp_all(1 : nlev, 1 : nprof))

    ! temp_qc
    !
    call nfw_inq_varid(fname, ncid, 'POTM_CORRECTED_QC', id)
    call nfw_get_var_text(fname, ncid, id, temp_qc_all(1 : nlev, 1 : nprof))

    ! psal
    !
    call nfw_inq_varid(fname, ncid, 'PSAL_CORRECTED', id)
    call nfw_get_var_double(fname, ncid, id, salt_all(1 : nlev, 1 : nprof))

    ! psal_qc
    !
    call nfw_inq_varid(fname, ncid, 'PSAL_CORRECTED_QC', id)
    call nfw_get_var_text(fname, ncid, id, salt_qc_all(1 : nlev, 1 : nprof))
    
    ! deph 
    !
    call nfw_inq_varid(fname, ncid, 'DEPH_CORRECTED', id)
    call nfw_get_var_double(fname, ncid, id, deph_all(1 : nlev, 1 : nprof))

    ! profile_potm_qc
    !
    call nfw_inq_varid(fname, ncid, 'PROFILE_POTM_QC', id)
    call nfw_get_var_text(fname, ncid, id, prof_temp_qc_all(1 : nprof))

    ! profile_psal_qc
    !
    call nfw_inq_varid(fname, ncid, 'PROFILE_PSAL_QC', id)
    call nfw_get_var_text(fname, ncid, id, prof_saln_qc_all(1 : nprof))

    call nfw_close(fname, ncid)
  end subroutine data_readfile

  real(8) function potential_density(T, S)
    real(8), intent(in) :: T, S

    if (T < -2.0d0 .or. T > 40.0d0 .or. S < 0.0d0 .or. S > 42.0d0) then
       potential_density = -999.0d0
       return
    end if

    potential_density =&
         -9.20601d-2&
         + T * (5.10768d-2 + S * (- 3.01036d-3)&
         + T * (- 7.40849d-3 + T * 3.32367d-5 + S * 3.21931d-5))&
         + 8.05999d-1 * S
  end function potential_density

  ! Purpose: find pivot points from profile data 
  !
  subroutine get_pivot(nx, ny, nprof, lon, lat, ipiv, jpiv, depths, modlon, modlat)
    use m_get_micom_grid
    use m_get_micom_dim
    use m_pivotp_micom
    
    implicit none

    ! Grid dimensions
    integer, intent(in) :: nx, ny
    integer, intent(in) :: nprof 

    real   , intent(in)    :: lon(nprof), lat(nprof)
    integer, intent(inout) :: ipiv(nprof), jpiv(nprof)

    real, intent(inout), dimension(nx,ny) :: depths, modlon, modlat

    real, parameter :: onem=98060.
    
    character(len=80) :: filename
    logical          :: ex
    character(len=8) :: ctmp

    real :: meandx,mindx
    real, allocatable, dimension(:,:) :: min_r, max_r
    integer, allocatable, dimension(:,:) :: itw, &
         jtw, its, jts, itn, jtn, ite, jte
    integer :: p
    integer :: dimids(2)
    integer :: ncid, x_ID, y_ID, z_ID 
    integer :: vJPIV_ID, vIPIV_ID
    integer :: ncid2, jns_ID, ins_ID, inw_ID, jnw_ID,jnn_ID, inn_ID, ine_ID, jne_ID
    
    allocate(min_r(nx, ny))
    allocate(max_r(nx, ny))
    allocate(itw(nx, ny))
    allocate(jtw(nx, ny))
    allocate(its(nx, ny))
    allocate(jts(nx, ny))
    allocate(itn(nx, ny))
    allocate(jtn(nx, ny))
    allocate(ite(nx, ny))
    allocate(jte(nx, ny))

    ! Read position and depth from model grid
    !
    call get_micom_grid(modlon, modlat, depths, mindx, meandx, nx, ny)
    call ini_pivotp(modlon,modlat, nx, ny, min_r, max_r, itw, jtw, itn, jtn, &
         its, jts, ite, jte)
    do p=1,nprof
       ipiv(p) = 1
       jpiv(p) = 1
       if (lat(p) .ge. -90 .and. lat(p) .le. 89 .and. lon(p) .ge. -180 .and. lon(p) .le. 180) then
          call pivotp_micom_new(lon(p), lat(p), modlon, modlat, ipiv(p), jpiv(p), &
               nx, ny, min_r, max_r,itw, jtw, its, jts, itn, jtn, ite, jte)
       end if
    end do
  end subroutine get_pivot

  ! Define the observation error variance as in Xie and Zhu, 2010
  ! 
  subroutine data_variance(obstype, depth, var)

    implicit none
    
    character(*), intent(in) :: obstype
    real        , intent(in) :: depth

    real, intent(inout) :: var

    if (trim(obstype) == 'TEM') then
       var = 0.05 + 0.45 * exp(-0.002 * depth)
       var = var ** 2.0
    elseif(trim(obstype) == 'SAL') then
       var = 0.02 + 0.10 * exp(-0.008 * depth)
       var = var ** 2.0
    else
       print *, 'ERROR: data_variance(): the definition of variance is only available for <TEM> and <SAL>'
       print *, 'The inital variance in the file <infile.data> should be non-negative...'
       stop
    end if
  end subroutine data_variance

  ! Purpose: read obervation uncertainties (instrumental and 
  !          representativeness error) from the pre-estimation.
  ! Refer to Karspeck, A. (2016).
  !
  subroutine data_obsunc(fname, obstype, depth_bnd, field)
    use nfw_mod

    character(*), intent(in)            :: fname
    character(*), intent(in)                 :: obstype
    real(8), intent(inout), dimension(:,:)   :: depth_bnd
    real(8), intent(inout), dimension(:,:,:) :: field
    
    integer :: ncid
    integer :: id
    
    !print *, '  reading "', trim(fname), '"'
    call nfw_open(trim(fname), nf_nowrite, ncid)

    ! depth_bnd
    call nfw_inq_varid(trim(fname), ncid, 'depth_bnds', id)
    call nfw_get_var_double(trim(fname), ncid, id, depth_bnd)
    
    ! Obs_Unc
    call nfw_inq_varid(trim(fname), ncid, 'var_o', id)
    call nfw_get_var_double(trim(fname), ncid, id, field)

    call nfw_close(trim(fname), ncid)
  end subroutine data_obsunc

  ! Purpose: read obervation climatology
  !                                     
  subroutine data_obsmean(fname, obstype, depth, lon, lat, field)
    use nfw_mod

    character(*), intent(in)            :: fname
    character(*), intent(in)                 :: obstype
    real(8), intent(inout), dimension(:)     :: depth, lon, lat
    real(8), intent(inout), dimension(:,:,:) :: field

    integer :: ncid
    integer :: id
    !integer, dimension(4) :: ns, nc

!    print *, '  reading "', trim(fname), '"'
    call nfw_open(trim(fname), nf_nowrite, ncid)

    ! depth
    call nfw_inq_varid(trim(fname), ncid, 'depth', id)
    call nfw_get_var_double(trim(fname), ncid, id, depth)

    ! lon
    call nfw_inq_varid(trim(fname), ncid, 'lon', id)
    call nfw_get_var_double(trim(fname), ncid, id, lon)

    ! lat
    call nfw_inq_varid(trim(fname), ncid, 'lat', id)
    call nfw_get_var_double(trim(fname), ncid, id, lat)

    ! Obs_mean
    if (trim(obstype) == 'SAL') then
       call nfw_inq_varid(trim(fname), ncid, 'salinity', id)
       call nfw_get_var_double(trim(fname), ncid, id, field)
    elseif (trim(obstype) == 'TEM') then
       call nfw_inq_varid(trim(fname), ncid, 'temperature', id)
       call nfw_get_var_double(trim(fname), ncid, id, field(:,:,:))
       ! Convert from Kelvin to Celcius
       field = field - 273.15
    end if

    call nfw_close(trim(fname), ncid)
  end subroutine data_obsmean
    
end module  m_read_EN4_profile
