module m_read_CLS_SSH
  
contains

  subroutine read_CLS_SSH(fname,data,modlon,modlat,depths,dlon,dlat,nrobs,nx,ny)
  use mod_measurement
  use m_pivotp_micom
  use mod_grid
  use nfw_mod

  implicit none

  integer, intent(in) :: nx, ny
  integer, intent(in) :: dlon,dlat
  integer, intent(out) :: nrobs
  type (measurement), intent(inout)  :: data(:)
  real, dimension(nx,ny), intent(in) :: depths,modlon,modlat
  character(len=80), intent(in) :: fname

  integer :: vLON_ID,vLAT_ID
  integer :: ncid,vSSH_ID,i,j,k
  integer :: irec
  integer :: ipiv, jpiv
  integer :: ipp1,ipm1,jpp1,jpm1
  integer, dimension(nx,ny) :: itw, jtw, its, jts, itn, jtn, ite, jte

  logical :: ex, wet

  real :: lon, lat,ssh,ssh_sq, wetsill
  real, dimension(nx,ny) :: min_r, max_r, obs_unc
  real(4), dimension(1) :: scalefac, addoffset, fillvalue, fillvalue2
  real(4), allocatable :: vssh(:,:,:,:), vssh2(:,:,:,:)
  real   , allocatable :: vlongitude(:), vlatitude(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call ini_pivotp(modlon,modlat, nx, ny, min_r, max_r, itw, jtw, itn, jtn, its, jts, ite, jte)
  ipiv=1
  jpiv=1

  ! Read  observation file
  allocate(vssh(dlon,dlat,1,1))
  allocate(vlongitude(dlon))
  allocate(vlatitude(dlat))
  
  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'Data file ', fname, ' not found.'
     stop
  end if
  call nfw_open(fname, nf_nowrite, ncid)
  call nfw_inq_varid(fname, ncid,'height', vSSH_ID)
  call nfw_inq_varid(fname, ncid,'lon', vLON_ID)
  call nfw_inq_varid(fname, ncid,'lat', vLAT_ID)
  call nfw_get_var_real(fname, ncid, vSSH_ID, vssh)
  call nfw_get_var_double(fname, ncid, vLON_ID, vlongitude)
  call nfw_get_var_double(fname, ncid, vLAT_ID, vlatitude)
  call nfw_get_att_real(fname, ncid, vSSH_ID, 'add_offset', addoffset)
  call nfw_get_att_real(fname, ncid, vSSH_ID, 'scale_factor', scalefac)
  call nfw_get_att_real(fname, ncid, vSSH_ID, '_FillValue', fillvalue)
  call nfw_close(fname, ncid)

  where (vssh(:,:,1,1) .ne. fillvalue(1))
     vssh(:,:,1,1) =  vssh(:,:,1,1) * scalefac(1) + addoffset(1)
  end where

  !Read the monthly mean 
  allocate(vssh2(dlon,dlat,1,1))

  print *, 'Start  reading anom'
  call nfw_open('mean_obs.nc', nf_nowrite, ncid)
  print *, 'openning ID'
  call nfw_inq_varid('mean_obs.nc', ncid, 'height', vSSH_ID)
  print *, 'reading ID'
  call nfw_get_var_real('mean_obs.nc', ncid, vSSH_ID, vssh2)
  call nfw_get_att_real('mean_obs.nc', ncid, vSSH_ID, '_FillValue', fillvalue2)
  print *, 'closing'
  call nfw_close('mean_obs.nc', ncid)
  print *, 'Finished  reading anom'

  where (vssh(:,:,1,1) .ne. fillvalue(1) .and. vssh2(:,:,1,1) .ne. fillvalue2(1))
     vssh(:,:,1,1) = vssh(:,:,1,1) - vssh2(:,:,1,1)
  elsewhere
     vssh(:,:,1,1) = fillvalue(1)
  end where

  ! read pre-estimated obs errors
  call nfw_open('./obs_unc_SSH.nc', nf_nowrite, ncid)
  call nfw_inq_varid('./obs_unc_SSH.nc', ncid, 'var_o', vSSH_ID)
  call nfw_get_var_double('./obs_unc_SSH.nc', ncid, vSSH_ID, obs_unc)
  call nfw_close('./obs_unc_SSH.nc', ncid)

  print *,'Nb obs mem'
  nrobs=1
  wetsill = 200.  ! Discarding data in shallow waters 
  do j = 1, dlat
     do i = 1, dlon
        call pivotp_micom(vlongitude(i), vlatitude(j), modlon, modlat, ipiv, jpiv, &
             nx, ny, min_r, max_r,itw, jtw, itn, jtn, its, jts, ite, jte)
#ifdef MASK_LANDNEIGHBOUR
        ipm1=max(ipiv-1,1)
        ipp1=min(ipiv+1,nx)
        jpm1=max(jpiv-1,1)
        jpp1=min(jpiv+1,ny)
        if (any(depths(ipm1:ipp1, jpm1:jpp1) < wetsill) ) cycle
#endif
        if (depths(ipiv, jpiv) < wetsill ) cycle
        if (vssh(i,j,1,1) .ne. fillvalue(1)) then
           data(nrobs)%d = vssh(i,j,1,1)
           data(nrobs)%ipiv = ipiv
           data(nrobs)%jpiv = jpiv
           !regular grid [-179.5 -> 179.5] & [89.5 -> -89.5]
           data(nrobs)%lon = vlongitude(i) 
           data(nrobs)%lat = vlatitude(j)
           data(nrobs)%a1 = 1
           data(nrobs)%a2 = 0
           data(nrobs)%a3 = 0
           data(nrobs)%a4 = 0
           data(nrobs)%ns = 0
           data(nrobs)%depth = 0
           data(nrobs)%date = 0
           data(nrobs)%id ='SSH'
           data(nrobs)%orig_id =0
           data(nrobs)%i_orig_grid = -1
           data(nrobs)%j_orig_grid = -1
           data(nrobs)%h = 1
           data(nrobs)%status = .true.
           data(nrobs)%var = max(0.0025, obs_unc(ipiv, jpiv))
           nrobs=nrobs+1
        endif
     enddo   ! dlat
  enddo    ! dlon
  nrobs=nrobs-1
  print *,'Max,min obs',maxval(data(:)%d),minval(data(:)%d),maxval(data(:)%lon),minval(data(:)%lat)
  print *,'Max,min age',maxval(data(:)%date),minval(data(:)%date)
 ! print *,'Nb of obs',nrobs
 ! print *,'Max,min lon',maxval(data(:)%lon),minval(data(:)%lon)
 ! print *,'Max,min lat',maxval(data(:)%lat),minval(data(:)%lat)
end subroutine read_CLS_SSH
end module m_read_CLS_SSH
