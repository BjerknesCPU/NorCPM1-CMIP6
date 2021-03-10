 module m_read_micom_SST

contains

  subroutine read_micom_SST(fname,cmonth,data,modlon, modlat, depths,nx,ny,nrobs)
  use mod_measurement
  use mod_grid
  use nfw_mod

  implicit none

  integer, intent(in) :: nx,ny
  integer, intent(out) :: nrobs
  type (measurement), intent(inout)  :: data(:)
  real, dimension(nx,ny), intent(in) :: modlon,modlat,depths
  character(len=80), intent(in) :: fname,cmonth
  real(4) :: vsst(nx,ny,1)
  real(4) :: vfice(nx,ny)
  integer :: ncid,vSST_ID,i,j,k,imonth
  integer :: vFICE_ID
  logical :: ex, found, fleeting
  real :: lon, lat
  real(4), dimension(1) :: scalefac, addoffset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read  observation file
  vfice(:,:)=0
  inquire (file='mask.nc', exist=ex)
  if (ex) then
     print *,'Masking sst under cice'
     call nfw_open('mask.nc', nf_nowrite, ncid)
     call nfw_inq_varid('mask.nc', ncid,'ficem', vFICE_ID)
     call nfw_get_var_real('mask.nc', ncid, vFICE_ID, vfice)
     call nfw_close('mask.nc', ncid)
     !no offset or scaling factor
  endif

  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'Data file ', fname, ' not found.'
     stop
  end if
  call nfw_open(fname, nf_nowrite, ncid)
  call nfw_inq_varid(fname, ncid,'sst', vSST_ID)
  call nfw_get_var_real(fname, ncid, vSST_ID, vsst)
  !call nfw_get_att_real(fname, ncid, vSST_ID, 'add_offset', addoffset)
  !call nfw_get_att_real(fname, ncid, vSST_ID, 'scale_factor', scalefac)
  k=1
  do j = 1, ny
     do i = 1, nx
       if (depths(i,j)>0 .and. vsst(i,j,1)>-1.81 .and. vfice(i,j)==0) then
        data(k)%d = vsst(i,j,1)
        data(k)%ipiv = i
        data(k)%jpiv = j
        data(k)%lon = modlon(i,j)
        data(k)%lat = modlat(i,j)
        data(k)%a1 = 1
        data(k)%a2 = 0
        data(k)%a3 = 0
        data(k)%a4 = 0
        data(k)%ns = 0
        data(k)%var = 0.01
        data(k)%depth = 0
        data(k)%date = 0
        data(k)%status = .true.
        data(k)%id ='SST'
        data(k)%orig_id =0
        data(k)%i_orig_grid = -1
        data(k)%j_orig_grid = -1
        data(k)%h = 1
        k=k+1
       endif
     enddo   ! ny
  enddo    ! nx
  nrobs=k-1
  !print *,'Max,min obs',maxval(data(:)%d),minval(data(:)%d),maxval(depths(:,:)),minval(depths(:,:))
end subroutine read_micom_SST
end module m_read_micom_SST
