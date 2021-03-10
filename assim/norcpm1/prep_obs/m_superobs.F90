! File:          m_superobs.F90
!
! Created:       02 Sep 2008
!
! Author:        Pavel Sakov
!                NERSC
!
! Purpose:       Superobing observations to model grid (one observation/profile only per grid)
!
! Description:   Adapted from algorithm of Pavel SAKOV
!                Conducts the following operations:
!                  - determine the number of observations with this tag
!                  - sort observations according to the pivot values
!                  - calculate superobs
!
! Modifications: 14.10.2009 PS: added cycle over the data age, so that it only
!                               superobs the data of the same age. Do not set 
!                               age if you assiilate data of different age in
!                               one go.
!                15.11.2009 PS: fixed a defect at l.102: it should be
!                               "thisob = obs_now(sorted(o))", not
!                               "thisob = obs_now(o)"
!                17.11.2009 PS: extended to handle the 3D case 
!                02.03.2016 Yiguo Wang: modify the definition  of kpiv using 
!                               the intrinsic function INT instead of z2k  
!
module m_superobs
  use mod_measurement
  use m_bilincoeff
  implicit none

  integer, parameter, private :: STRLEN = 512
  logical, parameter, private :: TEST = .false.

  contains

  subroutine superob(obstag, nobs, obs, ni, nj, modlon, modlat, nnewobs, newobs, is3d)
    character(*), intent(in) :: obstag
    integer, intent(in) :: nobs
    type(measurement), intent(inout), dimension(:) :: obs
    integer, intent(in) :: ni, nj
    real, dimension(:,:), intent(in) :: modlon, modlat
    integer, intent(inout) :: nnewobs
    type(measurement), intent(inout), dimension(:) :: newobs
    logical, intent(in), optional :: is3d

    integer :: age_min, age_max, nobs_total, nobs_now, age
    integer :: o, iprev, jprev, kprev, ii, ii_now
    logical, dimension(nobs) :: mask
    integer, dimension(nobs) :: sorted
    type(measurement), dimension(nobs) :: obs_now
    real(8), dimension(1) :: nobs_real
    type(measurement) :: thisob
    real :: n, nmax, valsum, valsqsum, varinvsum, lonsum, latsum, depthsum, valmax, valmin
    real :: a1sum, a2sum, a3sum, a4sum
    integer :: nlon_pos, nlon_neg
    real :: lonsum_abs
    integer, dimension(nobs) :: kpiv ! vertical index for 3D case
    integer, dimension(nobs) :: ids ! ids of obs contributing to this superob
    integer :: fid

    print *, 'BEGIN superob()'

    ! find the range of the data age
    !
    age_min = minval(obs (1:nobs) % date)
    age_max = maxval(obs (1:nobs) % date)
    print *, '  min age =', age_min
    print *, '  max age =', age_max

    ! get the total number of observations to process
    !
    mask = .false.
    do o = 1, nobs
       if (trim(obs(o) % id) == trim(obstag)) then
          mask(o) = .true.
       end if
       obs(o) % orig_id = o
    end do
    nobs_total = count(mask)
    print *, '  total # of obs of all types =', nobs
    print *, '  total # of obs of type "', trim(obstag), '" =', nobs_total

    if (TEST) then
       open(101, file = 'superobs.txt', access = 'sequential', status = 'replace')
    end if

    ii = 0
    do age = age_min, age_max
       ! trim() prevents vectorising below
       mask = .false.
       do o = 1, nobs
          if (trim(obs(o) % id) == trim(obstag) .and. obs(o) % date == age .and. obs(o) % status) then
             mask(o) = .true.
          end if
       end do

       nobs_now = count(mask)
       print *, '  age =', age
       print *, '    nobs =', nobs_now

       if (nobs_now == 0) then
          cycle
       end if

       obs_now(1 : nobs_now) = pack(obs(1 : nobs), mask)

       nobs_real(1) = nobs_now
       if (.not. present(is3d) .or. .not. is3d) then
          kpiv = 0
          call sortgriddedobs(nobs_real, obs_now % ipiv, obs_now % jpiv, sorted)
       else
          kpiv = z2k(obs_now % depth) !int(obs_now % depth)
          call sortgriddedobs3d(nobs_real, obs_now % ipiv, obs_now % jpiv,&
               kpiv, sorted)
       end if
 
       iprev = 0
       jprev = 0
       kprev = 0
       nmax = 0
       ii_now = 0
       do o = 1, nobs_now + 1
          if (o <= nobs_now) then
             thisob = obs_now(sorted(o))
          else
             thisob % ipiv = -1 ! to force write of the previous measurement
          end if
          if (thisob % ipiv /= iprev .or. thisob % jpiv /= jprev .or. kpiv(sorted(o)) /= kprev) then
             if (ii_now > 0) then ! write the previous measurement
                newobs(ii) % d = valsum / n
                if (is3d) then
                   newobs(ii) % var = varinvsum / n
                else
                   newobs(ii) % var = 1.0d0 / varinvsum
                end if
                newobs(ii) % id = obstag
                if (nlon_pos == 0 .or. nlon_neg == 0 .or. lonsum_abs / n < 90.0d0) then
                   newobs(ii) % lon = lonsum / n
                else
                   lonsum = lonsum + real(nlon_neg) * 360.0d0;
                   newobs(ii) % lon = lonsum / n
                   if (newobs(ii) % lon > 180.0d0) then
                      newobs(ii) % lon = newobs(ii) % lon - 360.0d0
                   end if
                end if
                newobs(ii) % lat = latsum / n
                newobs(ii) % depth = depthsum / n
                newobs(ii) % ipiv = iprev
                newobs(ii) % jpiv = jprev
                newobs(ii) % ns = 0 ! not 100% sure
                newobs(ii) % a1 = a1sum / n
                newobs(ii) % a2 = a2sum / n
                newobs(ii) % a3 = a3sum / n
                newobs(ii) % a4 = a4sum / n
                newobs(ii) % status = .true.
                newobs(ii) % i_orig_grid = -1
                newobs(ii) % j_orig_grid = -1
                newobs(ii) % h = n
                newobs(ii) % date = age
                newobs(ii) % orig_id = ids(1) ! ID of the first ob
                nmax = max(n, nmax)
                if (TEST) then
                   write(101, '(a, g10.3)') 'total # of obs = ', n
                   write(101, '(a, i6)') '  index = ', ii
                   write(101, '(a, g10.3)') '  d = ', newobs(ii) % d
                   write(101, '(a, g10.3)') '  var = ', newobs(ii) % var
                   write(101, '(a, g10.3)') '  lon = ', newobs(ii) % lon
                   write(101, '(a, g10.3)') '  lat = ', newobs(ii) % lat
                   write(101, '(a, i4)') '  ipiv = ', newobs(ii) % ipiv
                   write(101, '(a, i4)') '  jpiv = ', newobs(ii) % jpiv
                   write(101, '(a, g10.3)') '  depth = ', newobs(ii) % depth
                   write(101, '(a, g10.3)') '  a1 = ', newobs(ii) % a1
                   write(101, '(a, g10.3)') '  a2 = ', newobs(ii) % a2
                   write(101, '(a, g10.3)') '  a3 = ', newobs(ii) % a3
                   write(101, '(a, g10.3)') '  a4 = ', newobs(ii) % a4
                   write(101, '(a)') '---'
                   call superobs_dump(trim(obstag), ii, ids, int(n))
                end if
             end if
             if (o > nobs_now) then
                exit
             end if
             ii = ii + 1
             ii_now = ii_now + 1
             if (TEST) then
                write(101, '(a, i6)') 'new superob, index = ', ii
             end if
             n = 0.0
             valsum = 0.0d0
             valsqsum = 0.0d0
             varinvsum = 0.0d0
             lonsum = 0.0d0
             latsum = 0.0d0
             depthsum = 0.0
             a1sum = 0.0d0
             a2sum = 0.0d0
             a3sum = 0.0d0
             a4sum = 0.0d0
             valmax = -1.0d+20
             valmin = 1.0d+20
             iprev = thisob % ipiv
             jprev = thisob % jpiv
             kprev = kpiv(sorted(o))
             nlon_pos = 0
             nlon_neg = 0
             lonsum_abs = 0.0d0
          end if
          n = n + 1.0
          valsum = valsum + thisob % d
          valsqsum = valsqsum + (thisob % d) ** 2
          if (is3d) then
             varinvsum = varinvsum + thisob % var
          else
             varinvsum = varinvsum + 1.0 / thisob % var
          end if
          lonsum = lonsum + thisob % lon
          lonsum_abs = lonsum_abs + abs(thisob % lon)
          if (thisob % lon >= 0.0) then
             nlon_pos = nlon_pos + 1
          else
             nlon_neg = nlon_neg + 1
          end if
          latsum = latsum + thisob % lat
          depthsum = depthsum + thisob % depth
          a1sum = a1sum + thisob % a1
          a2sum = a2sum + thisob % a2
          a3sum = a3sum + thisob % a3
          a4sum = a4sum + thisob % a4
          valmin = min(valmin, thisob % d)
          valmax = max(valmax, thisob % d)
          ids(int(n)) = thisob % orig_id;
          if (TEST) then
             write(101, '(a, i6)') '  obs index = ', sorted(o)
             write(101, '(a, g10.3)') '    d = ', thisob % d
             write(101, '(a, g10.3)') '    var = ', thisob % var
             write(101, '(a, g10.3)') '    lon = ', thisob % lon
             write(101, '(a, g10.3)') '    lat = ', thisob % lat
             write(101, '(a, i4)') '    ipiv = ', thisob % ipiv
             write(101, '(a, i4)') '    jpiv = ', thisob % jpiv
             write(101, '(a, g10.3)') '    depth = ', thisob % depth
             write(101, '(a, g10.3)') '    a1 = ', thisob % a1
             write(101, '(a, g10.3)') '    a2 = ', thisob % a2
             write(101, '(a, g10.3)') '    a3 = ', thisob % a3
             write(101, '(a, g10.3)') '    a4 = ', thisob % a4
          end if
       end do ! obs for this age
       print *, '    nsuperobs =', ii_now
    end do ! age
    if (TEST) then
       close(101)
    end if

    nnewobs = ii
    print *, '  superobing("', trim(obstag), '"):'
    print *, '  ', nobs, 'observations ->', nnewobs, 'observations'
    print *, '  ', 'max # of obs found in a grid cell  =', int(nmax)
    print *, 'END superob()'
  end subroutine superob

  function z2k(z)
    real, intent(in), dimension(:) :: z
    integer, dimension(size(z)) :: z2k

    real, dimension(2, 42) :: depth_bnds
    integer :: i, k, nz
    

    nz = size(z)
    depth_bnds = reshape((/0.0, 10.0475, 10.0475, 20.1158, 20.1158, 30.214, 30.214, 40.3553, 40.3553, 50.5586, &
         50.5586, 60.8509, 60.8509, 71.2711, 71.2711, 81.8753, 81.8753, 92.7437, 92.7437, 103.9913, &
         103.9913, 115.7825, 115.7825, 128.3522, 128.3522, 142.0345, 142.0345, 157.3019, 157.3019, 174.8189,&
         174.8189, 195.5097, 195.5097, 220.6418, 220.6418, 251.9186, 251.9186, 291.5638, 291.5638, 342.3651,&
         342.3651, 407.6244, 407.6244, 490.9494, 490.9494, 595.8501, 595.8501, 725.1709, 725.1709, 880.5102,&
         880.5102, 1061.846, 1061.846, 1267.542, 1267.542, 1494.729, 1494.729, 1739.887, 1739.887, 1999.413,&
         1999.413, 2270.017, 2270.017, 2548.924, 2548.924, 2833.926, 2833.926, 3123.331, 3123.331, 3415.883,&
         3415.883, 3710.664, 3710.664, 4007.015, 4007.015, 4304.469, 4304.469, 4602.695, 4602.695, 4901.459,&
         4901.459, 5200.599, 5200.599, 5500./), shape(depth_bnds))
    
    do i = 1, nz
       do k = 1, 42
          if (depth_bnds(1, k) <= z(i) .and. z(i) < depth_bnds(2, k)) then
             z2k(i) = k
             exit
          end if
       end do
    end do
  end function z2k


  subroutine superobs_dump(tag, id, ids, n)
    use nfw_mod

    character(*) :: tag
    integer, intent(in) :: id
    integer, intent(in) :: ids(n)
    integer, intent(in) :: n

    character(STRLEN) :: fname
    character(64) :: dname
    character(64) :: vname
    integer :: ncid, did(1), vid

    if (id > NF_MAX_DIMS) then
       return
    end if

    write(fname, '(a, a, a)') 'superobs-', trim(tag), '.nc'
    if (id == 1) then
       print *, 'dumping obs ids for each superob to "', trim(fname), '"'
       call nfw_create(fname, nf_clobber, ncid)
    else
       call nfw_open(fname, nf_write, ncid)
       call nfw_redef(fname, ncid)
    end if

    write(dname, '(a,i0)') 'd', id
    call nfw_def_dim(fname, ncid, trim(dname), n, did(1))
    write(vname, '(a,i0)') 'v', id
    call nfw_def_var(fname, ncid, trim(vname), nf_int, 1, did(1), vid)
    call nfw_enddef(fname, ncid)

    call nfw_put_var_int(fname, ncid, vid, ids)

    call nfw_close(fname, ncid)
  end subroutine superobs_dump

end module m_superobs
