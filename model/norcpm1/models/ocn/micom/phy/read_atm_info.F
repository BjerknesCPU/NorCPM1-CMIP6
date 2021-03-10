      subroutine read_atm_info(filename,londim,latdim,rlon,rlat,land)
c
c --- read lon,lat and land from netcdf file
c
      implicit none
c
      character filename*120
      integer londim,latdim
      real rlon(londim),rlat(latdim),land(londim,latdim)
c
#include "netcdf.inc"
c
      integer status,ncid,lonid,latid,fieldid,i,j,
     .        start(3),count(3),stride(3)
      real*4 lon(londim),lat(latdim)
      integer*2 ifield(londim,latdim)
c
c --- open netcdf file
      status=nf_open(filename,nf_nowrite,ncid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_open: ',filename,',',nf_strerror(status)
        stop
      endif
c
c --- get id of variables
      status=nf_inq_varid(ncid,'lon',lonid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_inq_varid: lon:',nf_strerror(status)
        stop
      endif
      status=nf_inq_varid(ncid,'lat',latid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_inq_varid: lat:',nf_strerror(status)
        stop
      endif
      status=nf_inq_varid(ncid,'land',fieldid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_inq_varid: land:',nf_strerror(status)
        stop
      endif
c
c --- read lat/lon variables
      status=nf_get_var_real(ncid,lonid,lon)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_var_real: lon:',nf_strerror(status)
        stop
      endif
      status=nf_get_var_real(ncid,latid,lat)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_var_real: lat:',nf_strerror(status)
        stop
      endif
c
      do i=1,londim
        rlon(i)=lon(i)
      enddo
      do j=1,latdim
        rlat(j)=lat(j)
      enddo
c
c --- read field
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=londim
      count(2)=latdim
      count(3)=1
      stride(1)=1
      stride(2)=1
      stride(3)=1
      status=nf_get_vars_int2(ncid,fieldid,start,count,stride,ifield)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_vars_int2: land:',nf_strerror(status)
        stop
      endif
c
c --- close netcdf file
      status=nf_close(ncid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_close:',nf_strerror(status)
        stop
      endif
c
c --- convert field to real land mask
c$OMP PARALLEL DO
      do j=1,latdim
        do i=1,londim
          land(i,j)=ifield(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end
