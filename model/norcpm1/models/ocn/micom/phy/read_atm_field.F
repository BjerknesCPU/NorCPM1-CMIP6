      subroutine read_atm_field(filename,fieldname,londim,latdim,
     .                          timelevel,missing_value,rfield)
c
c --- Read field from NetCDF file at selected time level
c
      use mod_xc
c
      implicit none
c
      character filename*120,fieldname*5
      integer londim,latdim,timelevel
      real rfield(londim,latdim),missing_value
c
#include "netcdf.inc"
c
      integer status,ncid,timeid,fieldid,start(3),count(3),stride(3),i,j
      integer*2 ifield(londim,latdim),i2_mval
      real*8 time
      real*4 offset,scale_factor
c
c --- open netcdf file
      status=nf_open(filename,nf_nowrite,ncid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_open: ',filename,',',nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- get id of time variable
      status=nf_inq_varid(ncid,'time',timeid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_inq_varid: time:',nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- read time variable at selected time level
      status=nf_get_var1_double(ncid,timeid,timelevel,time)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_var1_double: time:',nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
      if (time.eq.0.) then
        write(*,'(2a)') 'read_atm_field: ',filename
        write(*,'(a)') '  Time variable is zero at selected time level!'
        write(*,'(a)') '  Probable cause is a corrupted file.'
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- get id of field variable
      status=nf_inq_varid(ncid,fieldname,fieldid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_inq_varid: ',fieldname,':',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- read field
      start(1)=1
      start(2)=1
      start(3)=timelevel
      count(1)=londim
      count(2)=latdim
      count(3)=1
      stride(1)=1
      stride(2)=1
      stride(3)=1
#if   defined(NCEP)
      status=nf_get_vars_int2(ncid,fieldid,start,count,stride,ifield)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_get_vars_int2: ',fieldname,':',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- read offset, scale factor, and value of no data
      status=nf_get_att_real(ncid,fieldid,'add_offset',offset)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: add_offset:',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
      status=nf_get_att_real(ncid,fieldid,'scale_factor',scale_factor)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: scale_factor:',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
      status=nf_get_att_int2(ncid,fieldid,'missing_value',i2_mval)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: missing_value:',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
c --- scale and add offset to field
c$OMP PARALLEL DO
      do j=1,latdim
        do i=1,londim
          if (ifield(i,j).eq.i2_mval) then
            rfield(i,j)=missing_value
          else
            rfield(i,j)=ifield(i,j)*scale_factor+offset
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
#elif defined(ERA)
      status=nf_get_vars_double(ncid,fieldid,start,count,stride,rfield)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_get_vars_double: ',fieldname,':',
     .                  nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
#endif
c
c --- close netcdf file
      status=nf_close(ncid)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_close:',nf_strerror(status)
        call xchalt('(read_atm_field)')
               stop '(read_atm_field)'
      endif
c
      return
      end
