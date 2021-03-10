! File:          m_Generate_element_Si.F90
!     
! Created:       ???                  
!
! Last modified: 09/03/2016
!
! Purpose:       Calculation of HA_i ("S_i")
!
! Description:   Calculates HA_i for ensmeble memmber i for given data type.
!                                                      
! Modifications:
!                09/03/2016 Yiguo WANG: create get_S_spline() that performs
!                           the piecewise cubic Hermite interpolate for vertical
!                           interpolation
!                18/03/2016 Yiguo WANG: set SAL_MIN = 1.0, instead of SAL_MIN = 5.0
!                22/07/2016 Yiguo WANG: create get_climato_spline() for climotological
!                           model data 
!
!

module m_Generate_element_Si
  implicit none

  public Generate_element_Si
  public get_S
  public get_S_spline

  integer, parameter, private :: NONE = 0
  integer, parameter, private :: TEMPERATURE = 1
  integer, parameter, private :: SALINITY = 2

  real, parameter, private :: TEM_MIN = -2.5
  real, parameter, private :: TEM_MAX = 35.0
  real, parameter, private :: SAL_MIN = 1.0
  real, parameter, private :: SAL_MAX = 41.0

  logical, parameter, private :: VERT_INTERP_GRID = .true.

contains

  subroutine Generate_element_Si(S, obstype, fld, depths, nx, ny, nz, t)
    use mod_measurement
    use m_obs
    implicit none

    real, dimension(nobs), intent(inout) :: S ! input/output vector
    character(len=5), intent(in) :: obstype ! the model fld type in "fld"
    integer, intent(in) :: nx,ny,nz ! grid size
    real, intent(in) :: fld   (nx,ny) ! field to be placed in Si
    real, intent(in) :: depths(nx,ny) ! depth mask -- needed for support 
    integer, intent(in), optional :: t !time of fld

    integer :: iobs
    integer :: i, j, ip1, jp1
    integer :: ix, jy, imin, imax, jmin, jmax, cnt

    logical :: isprofile
    real :: depth
    integer :: ns

    real, parameter :: undef = 999.9 ! land points have value huge()

    ! TEM, GTEM, SAL and GSAL come from profiles
    isprofile = (trim(obstype) .eq. 'SAL' .or.&
         trim(obstype) .eq. 'GSAL' .or.&
         trim(obstype) .eq. 'TEM' .or.&
         trim(obstype) .eq. 'GTEM')

    do iobs = 1, nobs
       if (trim(obstype) == obs(iobs) % id) then
          if (trim(obstype) .ne. 'TSLA' .or. obs(iobs) % date == t) then
             ! Get model gridcell
             i = obs(iobs) % ipiv
             j = obs(iobs) % jpiv
             ip1 = min(i + 1, nx)
             jp1 = min(j + 1, ny)
             
             depth = obs(iobs) % depth
             
             !TODO: 1. check consistency for ns = 1 vs ns = 0
             !      2. check consistency of running from -ns to +ns (this can
             !         lead perhaps for averaginf over -1 0 1 = 3 x 3 instead
             !         of 2 x 2 grid cells if ns = 1
             if (depth .lt. 10.0 .and. .not. isprofile) then ! satellite data
                ns = obs(iobs) % ns
                if(ns .lt. 2) then ! point data : zero support
                   if (fld(i,j)>1000.) print *, 'WTF',i,j,depths(i,j)
                   S(iobs) = fld(i, j) 
                else ! data support assumed a square of 2ns * 2ns grid cells
                   imin = max( 1, i - ns)
                   imax = min(nx, i + ns)
                   jmin = max( 1, j - ns)
                   jmax = min(ny, j + ns)
                   cnt = 0
                   S(iobs) = 0.0
                   do jy = jmin, jmax
                      do ix = imin, imax
                         if (depths(ix, jy) > 1.0 .and. abs(fld(ix, jy)) < 10.0d3 .and. fld(ix, jy) /= 0.0d0 .and. fld(ix, jy) + 1.0d0 /= fld(ix, jy)) then 
                            S(iobs) = S(iobs) + fld(ix, jy)
                            cnt = cnt + 1
                         endif
                      enddo
                   enddo
                   
                   if (cnt == 0) then
                      print *, ' observation on land ', i, j, obs(iobs) % d
                      stop 'm_Generate_element_Sij: report bug to LB (laurentb@nersc.no)'
                   end if
                   S(iobs) = S(iobs) / real(cnt)
                endif

             elseif(isprofile) then      ! in-situ data (in depth)
                print *,'(m_Generate_element_Si does not handle profiles yet)'
                stop '(m_Generate_element_Si)'
             else
                stop 'Generate_element_Sij: not a profile but depth is deeper than 10m'
             endif
          end if ! obs and model are at similar time
       end if ! (trim(obstype) == obs(iobs) % id) then
    end do
  end subroutine Generate_element_Si

  ! Get S = HA for in-situ data. Linearly interpolate for obs positioned
  ! between the layer centres; otherwise use the layer value for the obs above
  ! the middle of the first layer or below the middle of the last layer.
  !
  ! Note - this procedure parses through all obs for each ensemble member
  ! to work out profiles. This indeed invlolves some redundancy because
  ! this work could be done only once. However, the penalty (I think) is
  ! quite small compared to the time required for reading the fields from
  ! files and does not worth modifying (and complicating) the code.
  !
  subroutine get_S(S, obstag, nobs, obs, iens)
    use mod_measurement
    use m_insitu
    use m_get_micom_dim
    use m_get_micom_fld
    !    use m_parse_blkdat
    use m_parameters
    implicit none

    real, dimension(nobs), intent(inout) :: S
    character(*), intent(in) :: obstag
    integer, intent(in) :: nobs
    type(measurement), dimension(nobs) :: obs
    integer, intent(in) :: iens

    real, parameter :: ONEMETER = 98060.0

    ! obs stuff
    !
    integer :: p, o
    integer, allocatable, dimension(:) :: ipiv, jpiv
    real, allocatable, dimension(:) :: a1, a2, a3, a4

    ! grid stuff
    !
    integer :: k
    integer :: ni, nj, nk
    real :: rdummy

    ! vertical stuff
    !
    real, allocatable, dimension(:) :: zgrid, zcentre, zgrid_prev, zcentre_prev
    real, allocatable, dimension(:) :: v, v_prev

    ! fields & I/O stuff
    !
    real, allocatable, dimension(:, :) :: dz2d, v2d, sstbias, mld, offset, z
    integer :: tlevel
    character(8) :: fieldtag
    character(3) :: cmem
    character(80) :: fname
    real, dimension(2, 2) :: dz_cell, v_cell
    real :: dz, depth, z0, z1, z01, delta
    integer :: field

    field = NONE

    if (nobs == 0) then
       return
    end if

    if (master .and. iens == 1) then
       if (VERT_INTERP_GRID) then
          print *, trim(obstag), ': vertical interpolation in grid space'
       else
          print *, trim(obstag), ': vertical interpolation in physical space'
       end if
    end if

    !
    ! 1. Identify profiles presented in "obs"
    !

    ! note that profiles are being used by each
    ! ensemble member...
    !
    call insitu_setprofiles(obstag, nobs, obs)

    allocate(ipiv(nprof))
    allocate(jpiv(nprof))
    allocate(a1(nprof))
    allocate(a2(nprof))
    allocate(a3(nprof))
    allocate(a4(nprof))
    allocate(zgrid(nprof))
    allocate(zgrid_prev(nprof))
    allocate(zcentre(nprof))
    allocate(zcentre_prev(nprof))
    allocate(v(nprof))
    allocate(v_prev(nprof))

    ipiv = obs(pstart(1 : nprof)) % ipiv
    jpiv = obs(pstart(1 : nprof)) % jpiv
    a1 = obs(pstart(1 : nprof)) % a1
    a2 = obs(pstart(1 : nprof)) % a2
    a3 = obs(pstart(1 : nprof)) % a3
    a4 = obs(pstart(1 : nprof)) % a4

    !
    ! 2. Map the observations for this ensemble member proceeding by layers
    !    to reduce I/O:
    !
    !    -cycle through layers
    !       -find the middle of this layer
    !       -cycle through profiles
    !          -for each obs between the middle of the prev layer and the
    !           middle of this layer
    !             -interpolate the field value
    !             -write to S
    !

    ! get grid dimensions
    !
!    call parse_blkdat('idm   ','integer', rdummy, ni)
!    call parse_blkdat('jdm   ','integer', rdummy, nj)
!    call parse_blkdat('kdm   ','integer', rdummy, nk)
    call get_micom_dim(ni, nj, nk)

    allocate(v2d(ni, nj))
    allocate(dz2d(ni, nj))

    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'GSAL') then
       fieldtag = 'saln    '
       field = SALINITY
    elseif (trim(obstag) == 'TEM' .or. trim(obstag) == 'GTEM') then
       fieldtag = 'temp    '
       field = TEMPERATURE
    else
       if (master) then
          print *, 'ERROR: get_S(): unknown observatioon tag "', trim(obstag), '"'
       end if
       stop
    end if
    write(cmem, '(i3.3)') iens
    fname = 'forecast'//cmem

    if (field == TEMPERATURE .and. prm_prmestexists('sstb')) then
       allocate(sstbias(ni, nj))
       allocate(mld(ni, nj))
       allocate(offset(ni, nj))
       allocate(z(ni, nj))
       z = 0.0d0

       tlevel = 1
       call get_micom_fld_new(trim(fname), sstbias, 'sstb ', 0, tlevel, ni, nj)
       if (tlevel == -1) then
          if (master) then
             print *, 'ERROR: get_micom_fld_new(): failed for "sstb"'
          end if
          stop
       end if
       call get_micom_fld_new(trim(fname), mld, 'dpmixl  ', 0, tlevel, ni, nj)
       if (tlevel == -1) then
          if (master) then
             print *, 'ERROR: get_micom_fld_new(): failed for "dpmixl"'
          end if
          stop
       end if
     end if

    ! cycle through layers
    !
    tlevel = 1
    do k = 1, nk + 1

       if (k == 1) then
          zgrid_prev = 0.0
          zcentre_prev = 0.0
       end if

       if (k <= nk) then

          ! read the depth and the requested field at this layer
          !
          call get_micom_fld_new(trim(fname), dz2d, 'dp      ', k, tlevel, ni, nj)
          if (tlevel == -1) then
             if (master) then
                print *, 'ERROR: get_micom_fld_new(): failed for "dp"'
             end if
             stop
          end if
          call get_micom_fld_new(trim(fname), v2d, fieldtag, k, tlevel, ni, nj)
          if (tlevel == -1) then
             if (master) then
                print *, 'ERROR: get_micom_fld_new(): failed for "', fieldtag, '"'
             end if
             stop
          end if
       end if

       ! calculate correction from SST bias at this depth
       !
       if (field == TEMPERATURE .and. prm_prmestexists('sstb')) then
          offset = 0.0d0
          z = z + dz2d / 2.0 ! at the middle of the layer
          where (mld > 0.0d0 .and. mld < 1.0d8) ! < 10000 m
             offset = sstbias * exp(-(z / mld) ** 2)
          end where
          v2d = v2d - offset
          z = z + dz2d / 2.0
       end if

       ! cycle through profiles
       !
       do p = 1, nprof
          if (k <= nk) then
             dz_cell(:, :) = dz2d(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1)
             dz = dz_cell(1, 1) * a1(p) + dz_cell(2, 1) * a2(p)&
                  + dz_cell(1, 2) * a3(p) + dz_cell(2, 2) * a4(p)
             dz = dz / ONEMETER
             zgrid(p) = zgrid_prev(p) + dz
             zcentre(p) = (zgrid_prev(p) + zgrid(p)) / 2.0
             v_cell(:, :) = v2d(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1)
             v(p) = v_cell(1, 1) * a1(p) + v_cell(2, 1) * a2(p)&
                  + v_cell(1, 2) * a3(p) + v_cell(2, 2) * a4(p)
          else
             ! for the lower half of the last layer -- just use the layer value
             ! (note that there was no reading in this case, so that
             ! v = v_prev)
             zcentre(p) = zgrid(p)
          end if

          if (k == 1) then
             v_prev(p) = v(p)
          end if

          ! cycle through the obs, pick the ones in between the middle of the
          ! previous layer and the middle of this layer, interpolate the
          ! ensemble field to their locations, and save the results in S
          !
          z0 = zcentre_prev(p)
          z1 = zcentre(p)
          z01 = zgrid_prev(p)
          if (z1 == z0) then
             cycle
          end if
          do while (pstart(p) <= pend(p))
             o = pstart(p)
             depth = obs(o) % depth

             ! check that this obs is within the current layer
             !
             if (depth > z1 .and. k <= nk) then
                exit ! next profile
             elseif (depth >= z0 .and. depth <= z1) then

                if (.not. VERT_INTERP_GRID) then
                   ! interpolate linearly in physical space
                   !
                   S(o) = (z1 - depth) / (z1 - z0) * v_prev(p) +&
                        (depth - z0) / (z1 - z0) * v(p)
                else
                   ! interpolate linearly in the grid space
                   !
                   if (depth < z01) then
                      delta = 0.5d0 * (depth - z0) / (z01 - z0)
                   else
                      delta = 0.5d0 + 0.5d0 * (depth - z01) / (z1 - z01)
                   end if
                   S(o) = (1.0d0 - delta) * v_prev(p) + delta * v(p)
                end if

                ! Here we check the range of interpolated ensemble values;
                ! the range of observed values is checked in insitu_QC().
                !
                if (field == SALINITY) then
                   if ((S(o) < SAL_MIN .or. S(o) > SAL_MAX) .and. master) then
                      print *, 'WARNING: get_S(): suspicious value (SAL): ',&
                           'iens =', iens, ', obs =', o, ', profile = ', p,&
                           'depth =', depth, ', S =', S(o)
                   end if
                else if (field == TEMPERATURE) then
                   if ((S(o) < TEM_MIN .or. S(o) > TEM_MAX) .and. master) then
                      print *, 'WARNING: get_S(): suspicious value (TEM): ',&
                           'iens =', iens, ', obs =', o, ', profile = ', p,&
                           'depth =', depth, ', S =', S(o)
                      print *, v_cell
                      print *, dz_cell
                      print *, delta, v_prev(p), v(p)
                      print *, dz
                      print *, a1(p), a2(p), a3(p), a4(p)
                      print *, ipiv(p), jpiv(p)
                      stop
                   end if
                end if
             else ! k == nk + 1
                S(o) = v(p)
             end if
                ! go to the next obs
                !
                pstart(p) = pstart(p) + 1
          end do ! o
       end do ! p
       zgrid_prev = zgrid
       zcentre_prev = zcentre
       v_prev = v
    end do ! k

    deallocate(dz2d)
    deallocate(v2d)
    deallocate(v_prev)
    deallocate(v)
    deallocate(zcentre_prev)
    deallocate(zcentre)
    deallocate(zgrid_prev)
    deallocate(zgrid)
    deallocate(a4)
    deallocate(a3)
    deallocate(a2)
    deallocate(a1)
    deallocate(jpiv)
    deallocate(ipiv)
    if (allocated(sstbias)) then
       deallocate(sstbias)
       deallocate(mld)
       deallocate(offset)
       deallocate(z)
    end if
  end subroutine get_S


  ! Get S = HA for in-situ data. Piecewise cubic Hermite interpolate for 
  ! obs positioned between the first and last layer centres; otherwise use 
  ! the layer value for the obs above the middle of the first layer or 
  ! below the middle of the last layer.
  !
  ! Note - this procedure parses through all obs for each ensemble member
  ! to work out profiles. This indeed invlolves some redundancy because
  ! this work could be done only once. However, the penalty (I think) is
  ! quite small compared to the time required for reading the fields from
  ! files and does not worth modifying (and complicating) the code.
  !
  subroutine get_S_spline(S, obstag, nobs, obs, iens)
    use mod_measurement
    use mod_eosfun
    use m_insitu
    use m_get_micom_dim
    use m_get_micom_fld
    use m_parameters

    implicit none

    real, dimension(nobs), intent(inout) :: S
    character(*), intent(in) :: obstag
    integer, intent(in) :: nobs
    type(measurement), dimension(nobs) :: obs
    integer, intent(in) :: iens

    real, parameter :: g = 980.6 ! cm/s^2

    ! obs stuff
    !
    integer :: p, o, o1, o2
    integer, allocatable, dimension(:) :: ipiv, jpiv
    real, allocatable, dimension(:) :: a1, a2, a3, a4

    ! grid stuff
    !
    integer :: k, m, n
    integer :: ni, nj, nk

    ! vertical stuff
    !
    real, allocatable, dimension(:) :: zgrid, zcentre
    real, allocatable, dimension(:) :: v, vd
    
    ! fields & I/O stuff
    !
    real, allocatable, dimension(:, :, :) :: dz3d, v3d
    real, allocatable, dimension(:, :, :) :: temp, saln
    integer :: tlevel
    character(8) :: fieldtag
    character(3) :: cmem
    character(80) :: fname
    real, dimension(2, 2) :: dz_cell, v_cell, th_cell, sa_cell
    real :: dz, th, sa, plo, dphi, alp1, alp2
    integer :: field

    field = NONE

    if (nobs == 0) then
       return
    end if

    if (master .and. iens == 1) then
       print *, trim(obstag), ': spline vertical interpolation in grid space'
    end if

    !
    ! Define coefficients for equation of state functions
    !
    call eosini

    !
    ! 1. Identify profiles presented in "obs"
    !

    ! note that profiles are being used by each 
    ! ensemble member...
    !
    call insitu_setprofiles(obstag, nobs, obs)

    allocate(ipiv(nprof))
    allocate(jpiv(nprof))
    allocate(a1(nprof))
    allocate(a2(nprof))
    allocate(a3(nprof))
    allocate(a4(nprof))

    ipiv = obs(pstart(1 : nprof)) % ipiv
    jpiv = obs(pstart(1 : nprof)) % jpiv
    a1 = obs(pstart(1 : nprof)) % a1
    a2 = obs(pstart(1 : nprof)) % a2
    a3 = obs(pstart(1 : nprof)) % a3
    a4 = obs(pstart(1 : nprof)) % a4

    !
    ! 2. Map the observations for this ensemble member proceeding by layers
    !    to reduce I/O:
    !

    ! get grid dimensions
    !
    call get_micom_dim(ni, nj, nk)

    allocate(v3d(ni, nj, nk))
    allocate(dz3d(ni, nj, nk))
    allocate(temp(ni, nj, nk))
    allocate(saln(ni, nj, nk))
    allocate(zgrid(nk))
    allocate(zcentre(nk))
    allocate(v(nk))
    allocate(vd(nk))

    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'GSAL') then
       fieldtag = 'saln    '
       field = SALINITY
    elseif (trim(obstag) == 'TEM' .or. trim(obstag) == 'GTEM') then
       fieldtag = 'temp    '
       field = TEMPERATURE
    else
       if (master) then
          print *, 'ERROR: get_S_spline(): unknown observatioon tag "', trim(obstag), '"'
       end if
       stop
    end if
    write(cmem, '(i3.3)') iens
    fname = 'forecast'//cmem

    tlevel = 1
    ! read the depth and the requested field at all layers                                                                            
    !                                                                                                       
    call get_micom_fld(trim(fname), dz3d, 'dp  ', tlevel, ni, nj, nk)
    call get_micom_fld(trim(fname), temp, 'temp', tlevel, ni, nj, nk)
    call get_micom_fld(trim(fname), saln, 'saln', tlevel, ni, nj, nk)

    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'GSAL') then
       v3d = saln
    elseif (trim(obstag) == 'TEM' .or. trim(obstag) == 'GTEM') then
       v3d = temp
    else
       if (master) then
          print *, 'ERROR: get_S_spline(): unknown observatioon tag "', trim(obstag), '"'
       end if
       stop
    end if
    !call get_micom_fld(trim(fname), v3d, fieldtag, tlevel, ni, nj, nk)

    ! cycle through profiles
    !
    do p = 1, nprof
       ! cycle through layers
       !
        do k = 1, nk
           if (k == 1) then
              m = 1
              zgrid   = 0.
              zcentre = 0.
              plo = 0.
           end if
           dz_cell(:, :) = dz3d(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1, k)
           dz = dz_cell(1, 1) * a1(p) + dz_cell(2, 1) * a2(p)&
                + dz_cell(1, 2) * a3(p) + dz_cell(2, 2) * a4(p)

           if (dz .lt. 1.) cycle ! strictly insrease

           th_cell = temp(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1, k)
           th = th_cell(1, 1) * a1(p) + th_cell(2, 1) * a2(p)&
                + th_cell(1, 2) * a3(p) + th_cell(2, 2) * a4(p)

           sa_cell = saln(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1, k)
           sa = sa_cell(1, 1) * a1(p) + sa_cell(2, 1) * a2(p)&
                + sa_cell(1, 2) * a3(p) + sa_cell(2, 2) * a4(p)

           plo = plo + dz
           call delphi(plo, plo - dz, th, sa, dphi, alp1, alp2)
           dz = dphi / g / 100 ! cm -> meter

           if (m == 1) then
              zcentre(m) = dz / 2.0  
              zgrid(m) = zcentre(m) + dz / 2.0  
           else
              zcentre(m) = zgrid(m-1) + dz / 2.0  
              zgrid(m) = zcentre(m) + dz / 2.0  
           end if
           v_cell(:, :) = v3d(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1, k)
           v(m) = v_cell(1, 1) * a1(p) + v_cell(2, 1) * a2(p)&
                + v_cell(1, 2) * a3(p) + v_cell(2, 2) * a4(p)
           m = m + 1
        end do ! k
        m = m - 1

        ! observation length in the profile
        n = pend(p) - pstart(p) + 1 

        ! first obs below the first model layer centre
        o1 = pstart(p)
        do while (o1 .le. pend(p) .and. obs(o1) % depth .lt. zcentre(1))
           o1 = o1 + 1
        end do

        ! last obs above the last model layer centre
        o2 = pend(p)
        do while (o2 .ge. pstart(p) .and. obs(o2) % depth .gt. zcentre(m))
           o2 = o2 - 1
        end do

        if (o1 .gt. o2+1) then
           if (master) then
              print *, 'ERROR: get_S_spline(): failed for first and last obs in the range'
              print *, 'pstart =', pstart(p)
              print *, 'pend =', pend(p)
              print *, 'o1 =', o1
              print *, 'o2 =', o2
              print *, 'zcentre =', zcentre
              print *, 'm =', m
              print *, 'obs % depth =', obs(pstart(p):pend(p)) % depth
              print *, a1(p), a2(p), a3(p), a4(p)
              print *, ipiv(p), jpiv(p)
           end if
           stop
        end if

        ! use the first value for the obs above the middle of the first layer        
        do o = pstart(p), o1-1
           S(o) = v(1)
        end do

        ! use the last value for the obs below the middle of the last layer
        do o = o2+1, pend(p)
           S(o) = v(m)
        end do

        if (o2 - o1 .ge. 0) then
           call spline_pchip_set(m, zcentre(1:m), v(1:m), vd(1:m))
           call spline_pchip_val(m, zcentre(1:m), v(1:m), vd(1:m), o2 - o1 + 1, &
                obs(o1:o2) % depth, S(o1:o2))
        end if

        ! Here we check the range of interpolated ensemble values;                                                                  
        ! the range of observed values is checked in insitu_QC().                                                                   
        !                                                                                                                           
        do o = pstart(p), pend(p)
           if (field == SALINITY) then
              if ((S(o) < SAL_MIN .or. S(o) > SAL_MAX)) then
!!$                 print *, 'WARNING: get_S(): suspicious value (SAL): ',&
!!$                      'iens =', iens, ', obs =', o, ', profile = ', p,&
!!$                      'depth =', obs(o) % depth, ', S =', S(o)
!!$                 print *, zcentre(1:m)
!!$                 print *, v(1:m)
!!$                 print *, a1(p), a2(p), a3(p), a4(p)
!!$                 print *, ipiv(p), jpiv(p)
                 if (S(o) < SAL_MIN) then
                    S(o) = SAL_MIN
                 else
                    S(o) = SAL_MAX
                 end if
              end if
           else if (field == TEMPERATURE) then
              if ((S(o) < TEM_MIN .or. S(o) > TEM_MAX)) then
!!$                 print *, 'WARNING: get_S(): suspicious value (TEM): ',&
!!$                      'iens =', iens, ', obs =', o, ', profile = ', p,&
!!$                      'depth =', obs(o) % depth, ', S =', S(o)
!!$                 print *, zcentre(1:m)
!!$                 print *, v(1:m)
!!$                 print *, a1(p), a2(p), a3(p), a4(p)
!!$                 print *, ipiv(p), jpiv(p)
                 if (S(o) < TEM_MIN) then
                    S(o) = TEM_MIN
                 else
                    S(o) = TEM_MAX
                 end if
              end if
           end if
        end do ! o
     end do ! p

     deallocate(dz3d, v3d, v, vd)
     deallocate(zcentre, zgrid)
     deallocate(a1, a2, a3, a4)
     deallocate(ipiv, jpiv)
   end subroutine get_S_spline

   ! Get climatological data. Piecewise cubic Hermite interpolate for                                                                           
   ! obs positioned between min and max avlaible depths; otherwise use
   ! the first value for the obs above the first model depth or
   ! the last value below the last model depth.
   !
   subroutine get_climato_spline(S, obstag, nobs, obs)
     use mod_measurement
     use m_insitu
     use m_get_micom_dim
     use m_get_micom_fld
     use m_parameters
     use ieee_arithmetic
     implicit none

     real, dimension(nobs), intent(inout) :: S
     character(*), intent(in) :: obstag
     integer, intent(in) :: nobs
     type(measurement), dimension(nobs) :: obs
     
     ! obs stuff
     ! 
     integer :: p, o, o1, o2
     integer, allocatable, dimension(:) :: ipiv, jpiv
     real, allocatable, dimension(:) :: a1, a2, a3, a4

     ! grid stuff
     !
     integer :: k, m, n
     integer :: ni, nj, nk

     ! vertical stuff
     !
     real, allocatable, dimension(:) :: zgrid, zcentre
     real, allocatable, dimension(:) :: v, vd
     
     ! fields & I/O stuff
     !
     real, allocatable, dimension(:) :: dz3d
     real, allocatable, dimension(:, :, :) :: v3d
     integer :: tlevel
     character(8) :: fieldtag

     character(80) :: fname
     real, dimension(2, 2) :: dz_cell, v_cell
     real :: dz
     integer :: field
     
     field = NONE
     
     if (nobs == 0) then
        return
     end if
     
     if (master) then
        print *, trim(obstag), ': climatology spline vertical interpolation in grid space'
     end if

     !
     ! 1. Identify profiles presented in "obs"
     !
     call insitu_setprofiles(obstag, nobs, obs)
     
     allocate(ipiv(nprof), jpiv(nprof))
     allocate(a1(nprof), a2(nprof), a3(nprof), a4(nprof))
     
     ipiv = obs(pstart(1 : nprof)) % ipiv
     jpiv = obs(pstart(1 : nprof)) % jpiv
     a1 = obs(pstart(1 : nprof)) % a1
     a2 = obs(pstart(1 : nprof)) % a2
     a3 = obs(pstart(1 : nprof)) % a3
     a4 = obs(pstart(1 : nprof)) % a4

    !
    ! 2. Map the observations for climatology proceeding by depth
    !

    ! get grid dimensions
    ! 
    call get_climato_dim(ni, nj, nk)
    
    allocate(v3d(ni, nj, nk))
    allocate(dz3d(nk))
    allocate(zgrid(nk))
    allocate(zcentre(nk))
    allocate(v(nk))
    allocate(vd(nk))

    if (trim(obstag) == 'SAL' .or. trim(obstag) == 'GSAL') then
       fieldtag = 'saln    '
       field = SALINITY
    elseif (trim(obstag) == 'TEM' .or. trim(obstag) == 'GTEM') then
       fieldtag = 'temp    '
       field = TEMPERATURE
    else
       if (master) then
          print *, 'ERROR: get_S(): unknown observatioon tag "', trim(obstag), '"'
       end if
       stop
    end if

    fname = 'mean_mod'
    tlevel = 1
    ! read the depth and the requested field at all layers
    !
    call get_micom_fld_1d(trim(fname), dz3d, 'depth', nk)
    call get_micom_fld(trim(fname), v3d, trim(fieldtag)//'lvl', tlevel, ni, nj, nk)

    ! cycle through profiles
    do p = 1, nprof
       ! cycle through layers
       m = 1
       do k = 1, nk
          dz = dz3d(k)
          zcentre(m) = dz

          v_cell(:, :) = v3d(ipiv(p) : ipiv(p) + 1, jpiv(p) : jpiv(p) + 1, k)
          v(m) = v_cell(1, 1) * a1(p) + v_cell(2, 1) * a2(p)&
               + v_cell(1, 2) * a3(p) + v_cell(2, 2) * a4(p)
          if (count(v_cell < -32766.) .gt. 0) cycle
          m = m + 1
       end do ! k
       m = m - 1

       ! observation length in the profile
       n = pend(p) - pstart(p) + 1
       
       ! first obs below the first model depth
       o1 = pstart(p)
       do while (o1 .le. pend(p) .and. obs(o1) % depth .lt. zcentre(1))
          o1 = o1 + 1
       end do

       ! first obs above the last model depth
       o2 = pend(p)
       do while (o2 .ge. pstart(p) .and. obs(o2) % depth .gt. zcentre(m))
          o2 = o2 - 1
       end do

       if (o1 .gt. o2+1) then
          if (master) then
             print *, 'ERROR: get_climato_spline(): failed for first and last obs in the range'
             print *, 'pstart =', pstart(p)
             print *, 'pend =', pend(p)
             print *, 'o1 =', o1
             print *, 'o2 =', o2
             print *, 'zcentre =', zcentre
             print *, 'value = ', v
             print *, 'm =', m
             print *, 'obs % depth =', obs(pstart(p):pend(p)) % depth
             print *, a1(p), a2(p), a3(p), a4(p)
             print *, ipiv(p), jpiv(p)
          end if
          stop
        end if

        ! use the first value for the obs above the first depth
        do o = pstart(p), o1-1
           S(o) = v(1)
        end do

        ! use the last value for the obs below the last depth
        do o = o2+1, pend(p)
           S(o) = v(m)
        end do

        if (o2 - o1 .ge. 0) then
           call spline_pchip_set(m, zcentre(1:m), v(1:m), vd(1:m))
           call spline_pchip_val(m, zcentre(1:m), v(1:m), vd(1:m), o2 - o1 + 1, &
                obs(o1:o2) % depth, S(o1:o2))
        end if

        ! Here we check the range of interpolated values;
        ! the range of observed values is checked in insitu_QC().
        !
        do o = pstart(p), pend(p)
           if (ieee_is_nan(S(o))) then
              print *, 'WARNING: get_Climato_S(): find nan: ',&
                   'obs =', o, ', profile = ', p,&
                   'depth =', obs(o) % depth, ', S =', S(o)
              print *, zcentre(1:m)
              print *, v(1:m)
              print *, a1(p), a2(p), a3(p), a4(p)
              print *, ipiv(p), jpiv(p)
           end if
           if (field == SALINITY) then
              if ((S(o) < SAL_MIN .or. S(o) > SAL_MAX)) then
                 if (S(o) < SAL_MIN) then
                    S(o) = SAL_MIN
                 else
                    S(o) = SAL_MAX
                 end if
              end if
           else if (field == TEMPERATURE) then
              if ((S(o) < TEM_MIN .or. S(o) > TEM_MAX)) then
                 if (S(o) < TEM_MIN) then
                    S(o) = TEM_MIN
                 else
                    S(o) = TEM_MAX
                 end if
              end if
           end if
        end do ! o
     end do ! p
     
     if (master) then
        print *, trim(obstag), ': End of climatology splining vertical interpolation in grid space'
     end if

     deallocate(dz3d, v3d, v, vd)
     deallocate(zcentre, zgrid)
     deallocate(a1, a2, a3, a4)
     deallocate(ipiv, jpiv)
   end subroutine get_climato_spline

end module m_Generate_element_Si
