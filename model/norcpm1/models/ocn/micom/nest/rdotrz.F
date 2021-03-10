      subroutine rdotrz(ncid,depths_o,depths,z_ow,z_oe,z_os,z_on)
c
      implicit none
c
#include "dimensions.h"
#include "param_nest.h"
#include "netcdf.inc"
c
      integer ncid
      real depths_o(idm_o,jdm_o),depths(idm,jdm),
     .     z_ow(nbz,jdm,kdm+1),z_oe(nbz,jdm,kdm+1),
     .     z_os(idm,nbz,kdm+1),z_on(idm,nbz,kdm+1)
c
      real mval,fval
      parameter (mval=-1.e12,fval=-1.e13)
c
      integer mask_o(idm_o,jdm_o),start(4),count(4),
     .        status,varid,i,j,k,ibz
      integer*2 i2fld(idm_o,jdm_o),i2mval
      real*4 offset,sclfac
      real z_o(idm_o,jdm_o),z(idm,jdm),untfac,q
      character varname*2
c
c --- define unit conversion factor
      untfac=1.e2
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
      varname='dz'
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
c$OMP PARALLEL DO
      do j=1,jdm_o
        do i=1,idm_o
          z_o(i,j)=mval
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
        do ibz=1,nbz
          z_ow(ibz,j,kk+1)=-depths(ibz,j)*1.e2
          z_oe(ibz,j,kk+1)=-depths(ii-ibz,j)*1.e2
        enddo
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do ibz=1,nbz
        do i=1,ii
          z_os(i,ibz,kk+1)=-depths(i,ibz)*1.e2
          z_on(i,ibz,kk+1)=-depths(i,jj-ibz)*1.e2
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=kk,1,-1
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
c$OMP PARALLEL DO PRIVATE(q)
        do j=1,jdm_o
          do i=1,idm_o
            if (i2fld(i,j).eq.i2mval) then
              z_o(i,j)=mval
            else
              q=(i2fld(i,j)*sclfac+offset)*untfac
              if (q.gt.1.e-6) then
                if (z_o(i,j).eq.mval) then
                  z_o(i,j)=-depths_o(i,j)*1.e2+q
                else
                  z_o(i,j)=z_o(i,j)+q
                endif
              endif
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                   mval,fval,mask_o,ip,z_o,z)
c
c --- - extract boundary values
c$OMP PARALLEL DO
        do j=1,jj
          do ibz=1,nbz
            if (z(ibz,j).eq.mval) then
              z_ow(ibz,j,k)=z_ow(ibz,j,k+1)
            else
              z_ow(ibz,j,k)=z(ibz,j)
            endif
            if (z(ii-ibz,j).eq.mval) then
              z_oe(ibz,j,k)=z_oe(ibz,j,k+1)
            else
              z_oe(ibz,j,k)=z(ii-ibz,j)
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do ibz=1,nbz
          do i=1,ii
            if (z(i,ibz).eq.mval) then
              z_os(i,ibz,k)=z_os(i,ibz,k+1)
            else
              z_os(i,ibz,k)=z(i,ibz)
            endif
            if (z(i,jj-ibz).eq.mval) then
              z_on(i,ibz,k)=z_on(i,ibz,k+1)
            else
              z_on(i,ibz,k)=z(i,jj-ibz)
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
      enddo
c
      do k=2,kk
c$OMP PARALLEL DO
        do j=1,jj
          do ibz=1,nbz
            z_ow(ibz,j,k)=max(min(z_ow(ibz,j,k),z_ow(ibz,j,k-1)),
     .                        z_ow(ibz,j,kk+1))
            z_oe(ibz,j,k)=max(min(z_oe(ibz,j,k),z_oe(ibz,j,k-1)),
     .                        z_oe(ibz,j,kk+1))
          enddo
        enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
        do ibz=1,nbz
          do i=1,ii
            z_os(i,ibz,k)=max(min(z_os(i,ibz,k),z_os(i,ibz,k-1)),
     .                          z_os(i,ibz,kk+1))
            z_on(i,ibz,k)=max(min(z_on(i,ibz,k),z_on(i,ibz,k-1)),
     .                          z_on(i,ibz,kk+1))
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      return
      end
