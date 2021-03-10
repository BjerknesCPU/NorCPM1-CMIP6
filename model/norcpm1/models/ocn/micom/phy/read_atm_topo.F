      subroutine read_atm_topo(filename,londim,latdim,topo)
c
c --- read topography from netcdf file
c
      implicit none
c
      character filename*120
      integer londim,latdim
      real topo(londim,latdim)
c
#include "netcdf.inc"
c
#if   defined(NCEP)
      integer start(3),count(3),stride(3)
#elif defined(ERA)
      integer start(4),count(4),stride(4)
#endif
      integer status,ncid,fieldid,i,j
      real*4 field(londim,latdim)
c
c --- open netcdf file
      status=nf_open(filename,nf_nowrite,ncid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_open: ',filename,',',nf_strerror(status)
        stop
      endif
c
c --- get id of topography variable
      status=nf_inq_varid(ncid,'hgtsfc',fieldid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_inq_varid: hgtsfc:',nf_strerror(status)
        stop
      endif
c
c --- read field
      start(1)=1
      start(2)=1
      start(3)=1
#ifdef ERA
      start(4)=1
#endif
      count(1)=londim
      count(2)=latdim
      count(3)=1
#ifdef ERA
      count(4)=1
#endif
      stride(1)=1
      stride(2)=1
      stride(3)=1
#ifdef ERA
      stride(4)=1
#endif
      status=nf_get_vars_real(ncid,fieldid,start,count,stride,field)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_vars_real: hgtsfc:',nf_strerror(status)
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
c --- convert field to real topography field
c$OMP PARALLEL DO
      do j=1,latdim
        do i=1,londim
          topo(i,j)=field(i,latdim+1-j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
      end
