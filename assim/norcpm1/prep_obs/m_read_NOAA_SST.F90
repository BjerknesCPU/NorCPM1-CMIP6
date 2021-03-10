 module m_read_NOAA_SST

contains

  subroutine read_NOAA_SST(fname,cens,data,modlon,modlat,depths,dlon,dlat,nrobs)
  use mod_measurement
  use mod_grid
  use nfw_mod

  implicit none

  integer, intent(in) :: dlon,dlat
  integer, intent(out) :: nrobs
  type (measurement), intent(inout)  :: data(:)
  real, dimension(dlon,dlat), intent(in) :: depths,modlon,modlat
  character(len=80), intent(in) :: fname
  character(len=3), intent(in) :: cens
  real(4) ,allocatable :: vsst(:,:,:), vsic(:,:,:) ,vsst2(:,:,:),verr(:,:,:)
  real(4) ,allocatable :: vlongitude(:), vlatitude(:)
  integer :: vLON_ID,vLAT_ID,vERR_ID
  integer :: ncid,vSST_ID,vSST2_ID,i,j,k,imonth
  integer :: vsic_ID,irec,nens
  integer, allocatable :: ns(:), nc(:)
  logical :: ex, ice_status
  real :: lon, lat,sst,sst_sq
  real(4), dimension(1) :: scalefac_sst, addoffset_sst, scalefac_sic, addoffset_sic
  integer :: dimids(2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read  observation file
  read(cens,*) nens
  allocate(vsst(dlon,dlat,1))
  allocate(vsic(dlon,dlat,1))
  allocate(verr(dlon,dlat,1))
  allocate(vlongitude(dlon))
  allocate(vlatitude(dlat))
  allocate(ns(3))
  allocate(nc(3))

  ns(1)=1
  ns(2)=1
  ns(3)=1
  nc(1)=dlon
  nc(2)=dlat
  nc(3)=1

  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'Data file ', fname, ' not found.'
     stop
  end if
  call nfw_open(fname, nf_nowrite, ncid)
  call nfw_inq_varid(fname, ncid,'sst', vSST_ID)
  call nfw_inq_varid(fname, ncid,'lon', vLON_ID)
  call nfw_inq_varid(fname, ncid,'icec', vSIC_ID)
  call nfw_inq_varid(fname, ncid,'err', vERR_ID)
  call nfw_inq_varid(fname, ncid,'lat', vLAT_ID)
  call nfw_get_vara_real(fname, ncid, vSST_ID, ns, nc, vsst)
  call nfw_get_att_real(fname, ncid, vSST_ID, 'add_offset', addoffset_sst)
  call nfw_get_att_real(fname, ncid, vSST_ID, 'scale_factor', scalefac_sst)
  call nfw_get_vara_real(fname, ncid, vSIC_ID, ns, nc, vsic)
  call nfw_get_att_real(fname, ncid, vSIC_ID, 'add_offset', addoffset_sic)
  call nfw_get_att_real(fname, ncid, vSIC_ID, 'scale_factor', scalefac_sic)
  call nfw_get_vara_real(fname, ncid, vERR_ID, ns, nc, verr)
  call nfw_get_var_real(fname, ncid, vLON_ID, vlongitude)
  call nfw_get_var_real(fname, ncid, vLAT_ID, vlatitude)
  call nfw_close(fname, ncid)
  do j = 1, dlat
     do i = 1, dlon
        vsst(i,j,1)=vsst(i,j,1)*scalefac_sst(1)+addoffset_sst(1)
        vsic(i,j,1)=vsic(i,j,1)*scalefac_sic(1)+addoffset_sic(1)
     enddo   
  enddo    
#ifdef ANOMALY
!Read the monthly mean 
  allocate(vsst2(dlon,dlat,1))
  print *, 'Reading anom'
  call nfw_open('mean_obs.nc', nf_nowrite, ncid)
  call nfw_inq_varid('mean_obs.nc', ncid,'sst', vSST2_ID)
  call nfw_get_vara_real('mean_obs.nc', ncid, vSST2_ID, ns, nc, vsst2)
  call nfw_close('mean_obs.nc', ncid)
  vsst(:,:,1)=vsst(:,:,1)-vsst2(:,:,1)
#endif
  nrobs=0
  do j = 1, dlat
     do i = 1, dlon
          if (vsic(i,j,1).eq.0.) then
           nrobs=nrobs+1
           data(nrobs)%d = vsst(i,j,1)
           data(nrobs)%ipiv = 0
           data(nrobs)%jpiv = 0
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
           data(nrobs)%id ='SST'
           data(nrobs)%orig_id =0
           data(nrobs)%i_orig_grid = -1
           data(nrobs)%j_orig_grid = -1
           data(nrobs)%h = 1
           data(nrobs)%date = 0
           data(nrobs)%status = .true.
           data(nrobs)%var = max(real(verr(i,j,1)),0.01)
          endif
     enddo   ! dlat
  enddo    ! dlon
  print *,'Max,min obs',maxval(data(:)%d),minval(data(:)%d),maxval(data(:)%lon),minval(data(:)%lat)
  print *,'Max,min age',maxval(data(:)%date),minval(data(:)%date)
  print *,'Nb of obs',nrobs
 ! print *,'Max,min lon',maxval(data(:)%lon),minval(data(:)%lon)
 ! print *,'Max,min lat',maxval(data(:)%lat),minval(data(:)%lat)
end subroutine read_NOAA_SST
end module m_read_NOAA_SST
