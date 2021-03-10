 module m_read_micom_SSH

contains

  subroutine read_micom_SSH(fname,cmonth,data,modlon, modlat, depths,nx,ny)
  use mod_measurement
  use mod_grid
  use nfw_mod

  implicit none

  integer, intent(in) :: nx,ny
  type (measurement), intent(inout)  :: data(:)
  real, dimension(nx,ny), intent(in) :: modlon,modlat,depths
  character(len=80), intent(in) :: fname,cmonth
  real(4) :: vssh(nx,ny,12)
  integer :: ncid,vSSH_ID,i,j,k,imonth
  logical :: ex, found, fleeting
  real :: lon, lat
  real(4), dimension(1) :: scalefac, addoffset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read  observation file
  read(cmonth,'(i3.3)') imonth
  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'Data file ', fname, ' not found.'
     stop
  end if
  call nfw_open(fname, nf_nowrite, ncid)
  call nfw_inq_varid(fname, ncid,'sealv', vSSH_ID)
  call nfw_get_var_real(fname, ncid, vSSH_ID, vssh)
  call nfw_get_att_real(fname, ncid, vSSH_ID, 'add_offset', addoffset)
  call nfw_get_att_real(fname, ncid, vSSH_ID, 'scale_factor', scalefac)
  k=1
  do j = 1, ny
     do i = 1, nx
       if (depths(i,j)>0 .and. vssh(i,j,imonth)*scalefac(1)+addoffset(1)>-10 ) then
        data(k)%d = vssh(i,j,imonth)*scalefac(1)+addoffset(1)
        data(k)%ipiv = i
        data(k)%jpiv = j
        data(k)%lon = modlon(i,j)
        data(k)%lat = modlat(i,j)
        data(k)%a1 = 1
        data(k)%a2 = 0
        data(k)%a3 = 0
        data(k)%a4 = 0
        data(k)%ns = 0
        data(k)%var = 0.0009
        data(k)%depth = 0
        data(k)%date = 0
        data(k)%status = .true.
        data(k)%id = 'SSH'
        data(k)%orig_id =0
        data(k)%i_orig_grid = -1
        data(k)%j_orig_grid = -1
        data(k)%h = 1
        k=k+1
      endif
     enddo   ! ny
  enddo    ! nx
  !print *,'number of obs',k-1
  !print *,'Max,min obs',maxval(data(:)%d),minval(data(:)%d),maxval(depths(:,:)),minval(depths(:,:))
end subroutine read_micom_SSH
end module m_read_micom_SSH
