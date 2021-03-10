      subroutine rdotr3(ncid,varname,kd,fld_ow,fld_oe,fld_os,fld_on)
c
      implicit none
c
#include "dimensions.h"
#include "param_nest.h"
#include "netcdf.inc"
c
      integer ncid,kd
      character varname*(*)
      real fld_ow(nbz,jdm,kd),fld_oe(nbz,jdm,kd),
     .     fld_os(idm,nbz,kd),fld_on(idm,nbz,kd)
c
      real mval,fval
      parameter (mval=-1.e12,fval=-1.e13)
c
      integer mask_o(idm_o,jdm_o),start(4),count(4),
     .        status,varid,i,j,k,ibz
      integer*2 i2fld(idm_o,jdm_o),i2mval
      real*4 offset,sclfac
      real fld_o(idm_o,jdm_o),fld(idm,jdm),untfac
c
c --- define unit conversion factor
      if     (varname.eq.'utot'.or.varname.eq.'vtot'.or.
     .        varname.eq.'sealv') then
        untfac=1.e2
      elseif (varname.eq.'fice') then
        untfac=1.e-2
      else
        untfac=1.
      endif
c
c --- dummy mask for coarse outer grid
c$OMP PARALLEL DO
      do j=1,jdm_o
        do i=1,idm_o
          mask_o(i,j)=1
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- get id of field variable
      status=nf_inq_varid(ncid,varname,varid)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_inq_varid: ',varname,': ',
     .                  nf_strerror(status)
        stop
      endif
c
c --- read offset, scale factor, and value of no data
      status=nf_get_att_real(ncid,varid,'add_offset',offset)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: add_offset: ',
     .                  nf_strerror(status)
        stop
      endif
      status=nf_get_att_real(ncid,varid,'scale_factor',sclfac)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: scale_factor: ',
     .                  nf_strerror(status)
        stop
      endif
      status=nf_get_att_int2(ncid,varid,'_FillValue',i2mval)
      if (status.ne.nf_noerr) then
        write(*,'(2a)') 'nf_get_att_real: _FillValue: ',
     .                  nf_strerror(status)
        stop
      endif
c
      do k=1,kd
c
c --- - read field
        start(1)=1
        start(2)=1
        start(3)=k
        start(4)=1
        count(1)=idm_o
        count(2)=jdm_o
        count(3)=1
        count(4)=1
        status=nf_get_vara_int2(ncid,varid,start,count,i2fld)
        if (status.ne.nf_noerr) then
          write(*,'(4a)') 'nf_get_vara_int2: ',varname,': ',
     .                    nf_strerror(status)
          stop
        endif
c
c --- - scale and add offset to field
c$OMP PARALLEL DO
        do j=1,jdm_o
          do i=1,idm_o
            if (i2fld(i,j).eq.i2mval) then
              fld_o(i,j)=mval
            else
              fld_o(i,j)=(i2fld(i,j)*sclfac+offset)*untfac
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
        if     (varname.eq.'utot') then
          call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                     mval,fval,mask_o,iu,fld_o,fld)
        elseif (varname.eq.'vtot') then
          call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                     mval,fval,mask_o,iv,fld_o,fld)
        else
          call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                     mval,fval,mask_o,ip,fld_o,fld)
        endif
c
c --- - extract boundary values
        if (varname.eq.'utot') then
c$OMP PARALLEL DO
          do j=1,jj
            do ibz=1,nbz
              fld_ow(ibz,j,k)=fld(ibz+1,j)
              fld_oe(ibz,j,k)=fld(ii-ibz,j)
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO
          do j=1,jj
            do ibz=1,nbz
              fld_ow(ibz,j,k)=fld(ibz,j)
              fld_oe(ibz,j,k)=fld(ii-ibz,j)
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        if (varname.eq.'vtot') then
c$OMP PARALLEL DO
          do ibz=1,nbz
            do i=1,ii
              fld_os(i,ibz,k)=fld(i,ibz+1)
              fld_on(i,ibz,k)=fld(i,jj-ibz)
            enddo
          enddo
c$OMP END PARALLEL DO
        else
c$OMP PARALLEL DO
          do ibz=1,nbz
            do i=1,ii
              fld_os(i,ibz,k)=fld(i,ibz)
              fld_on(i,ibz,k)=fld(i,jj-ibz)
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
c
      enddo
c
      return
      end
