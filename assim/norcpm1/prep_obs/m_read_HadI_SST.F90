 module m_read_HadI_SST

contains

  subroutine read_HadI_SST(fname,cens,data,modlon,modlat,depths,dlon,dlat,nrobs)
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
  real(4) ,allocatable :: vsst(:,:,:,:), vsic(:,:,:,:) ,vsst2(:,:,:,:)
  real(4) ,allocatable :: vlongitude(:), vlatitude(:)
  integer :: vLON_ID,vLAT_ID
  integer :: ncid,vSST_ID,i,j,k,imonth
  integer :: vsic_ID,irec,nens
  integer, allocatable :: ns(:), nc(:)
  logical :: ex, ice_status
  real :: lon, lat,sst,sst_sq
  real(4), dimension(1) :: scalefac, addoffset, undef
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read  observation file
  read(cens,*) nens
  nens=10
  allocate(vsst(dlon,dlat,1,nens))
  allocate(vsic(dlon,dlat,1,nens))
  allocate(vlongitude(dlon))
  allocate(vlatitude(dlat))
  allocate(ns(4))
  allocate(nc(4))

  ns(1)=1
  ns(2)=1
  ns(3)=1
  ns(4)=1
  nc(1)=dlon
  nc(2)=dlat
  nc(3)=1
  nc(4)=nens

  inquire (file=fname, exist=ex)
  if (.not. ex) then
     print *, 'Data file ', fname, ' not found.'
     stop
  end if
  call nfw_open(fname, nf_nowrite, ncid)
  call nfw_inq_varid(fname, ncid,'sst', vSST_ID)
  call nfw_inq_varid(fname, ncid,'longitude', vLON_ID)
  call nfw_inq_varid(fname, ncid,'sic', vSIC_ID)
  call nfw_inq_varid(fname, ncid,'latitude', vLAT_ID)
  call nfw_get_att_real(fname, ncid, vSST_ID, '_FillValue', undef)
  call nfw_get_vara_real(fname, ncid, vSST_ID, ns, nc, vsst)
  call nfw_get_vara_real(fname, ncid, vSIC_ID, ns, nc, vsic)
  nc(4)=1
  call nfw_get_var_real(fname, ncid, vLON_ID, vlongitude)
  call nfw_get_var_real(fname, ncid, vLAT_ID, vlatitude)
  !call nfw_get_vara_real(fname, ncid, vLON_ID, 1, dlon, vlongitude)
  !call nfw_get_vara_real(fname, ncid, vLAT_ID, 1, dlat, vlatitude)
  call nfw_close(fname, ncid)
  !Convert from Kelvin to Celcius
  where (vsst.ne.undef(1)) vsst=vsst-273.15 
#ifdef ANOMALY
!Read the monthly mean 
  allocate(vsst2(dlon,dlat,1,1))
  nc(4)=1
  print *, 'Start  reading anom'
  call nfw_open('mean_obs.nc', nf_nowrite, ncid)
  print *, 'openning ID'
  call nfw_inq_varid('mean_obs.nc', ncid,'sst', vSST_ID)
  print *, 'reading ID'
  call nfw_get_vara_real('mean_obs.nc', ncid, vSST_ID, ns, nc, vsst2)
  print *, 'closing'
  call nfw_close('mean_obs.nc', ncid)
  print *, 'Finished  reading anom'
  !Convert from Kelvin to Celcius
  where (vsst2.ne.undef(1)) vsst2=vsst2-273.15
  do k = 1, nens
     where (vsst2(:,:,1,1).ne.undef(1)) vsst(:,:,1,k)=vsst(:,:,1,k)-vsst2(:,:,1,1)
  enddo    ! dlon
#endif
  print *,'Nb obs mem'
  nrobs=1
  do j = 1, dlat
     do i = 1, dlon
        sst=0.
        sst_sq=0.
        ice_status=.true.
        do k = 1, nens
          if (ice_status .and. vsic(i,j,1,k).eq.0. .and. vsst(i,j,1,k).ne.undef(1)) then
           !convert from Kelvin to Celcius
           sst=sst+vsst(i,j,1,k)
           sst_sq=sst_sq+vsst(i,j,1,k)**2
           !Only fill with realistic value for the last member
           data(nrobs)%d = 9999.
           data(nrobs)%ipiv = i
           data(nrobs)%jpiv = j
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
           data(nrobs)%status = .false.
           if (k.eq.nens) then
              data(nrobs)%status = .true.
              ! calculate the variance of the ensemble of obs
              !if nens=1 -> division by 0
              sst=sst/real(nens)
              sst_sq=sst_sq/real(nens)
              data(nrobs)%d = sst
              !Add a min erro because some places obs spread is null
             ! data(nrobs)%var = max(sst_sq-sst**2,0.01) -> min error of 0.1
             ! degree
              data(nrobs)%var = max(sst_sq-sst**2,0.01)
              nrobs=nrobs+1
           endif
          else
           ice_status=.false.
          endif
        enddo   ! dlon
     enddo   ! dlat
  enddo    ! dlon
  nrobs=nrobs-1
  print *,'Max,min obs',maxval(data(:)%d),minval(data(:)%d),maxval(data(:)%lon),minval(data(:)%lat)
  print *,'Max,min age',maxval(data(:)%date),minval(data(:)%date)
  print *,'Nb of obs',nrobs
 ! print *,'Max,min lon',maxval(data(:)%lon),minval(data(:)%lon)
 ! print *,'Max,min lat',maxval(data(:)%lat),minval(data(:)%lat)
end subroutine read_HadI_SST
end module m_read_HadI_SST
