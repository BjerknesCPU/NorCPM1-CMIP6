      subroutine rdotru(ncid,ub_ow,ub_oe,ub_os,ub_on,
     .                       vb_ow,vb_oe,vb_os,vb_on)
c
      implicit none
c
#include "dimensions.h"
#include "param_nest.h"
#include "netcdf.inc"
c
      integer ncid
      real ub_ow(jdm),ub_oe(jdm),ub_os(idm),ub_on(idm),
     .     vb_ow(jdm),vb_oe(jdm),vb_os(idm),vb_on(idm)
c
      real mval,fval
      parameter (mval=-1.e12,fval=-1.e13)
c
      integer mask_o(idm_o,jdm_o),start(3),count(3),
     .        nfac,status,varid,i,j,k
      integer*2 i2fld(idm_o,jdm_o),i2mval
      real*4 offset,sclfac
      real fld_o(idm_o,jdm_o),fld(idm,jdm),untfac,q
      character varname*5
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
c --- grid scale factor between outer and inner grid
      nfac=(idm-3)/(i2o-i1o+1)
c
c --- get id of field variable
      varname='ubaro'
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
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=idm_o
      count(2)=jdm_o
      count(3)=1
      status=nf_get_vara_int2(ncid,varid,start,count,i2fld)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_get_vara_int2: ',varname,': ',
     .                  nf_strerror(status)
        stop
      endif
c
c --- scale and add offset to field
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
      call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                 mval,fval,mask_o,iu,fld_o,fld)
c
c --- extract boundary values
c$OMP PARALLEL DO
      do j=1,jj
        ub_ow(j)=fld(2,j)
        ub_oe(j)=fld(ii-1,j)
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do i=1,ii
        ub_os(i)=fld(i,1)
        ub_on(i)=fld(i,jj-1)
      enddo
c$OMP END PARALLEL DO
c
c --- ensure "flux conservation" trough boundary
c$OMP PARALLEL DO PRIVATE(q)
      do j=j1o,j2o
c
        if (fld_o(i1o,j).ne.mval) then
          q=fld_o(i1o,j)*nfac
          do k=1,nfac
            q=q-ub_ow((j-j1o)*nfac+k+1)
          enddo
          q=q/nfac
          do k=1,nfac
            ub_ow((j-j1o)*nfac+k+1)=ub_ow((j-j1o)*nfac+k+1)+q
          enddo
        endif
c
        if (fld_o(i2o+1,j).ne.mval) then
          q=fld_o(i2o+1,j)*nfac
          do k=1,nfac
            q=q-ub_oe((j-j1o)*nfac+k+1)
          enddo
          q=q/nfac
          do k=1,nfac
            ub_oe((j-j1o)*nfac+k+1)=ub_oe((j-j1o)*nfac+k+1)+q
          enddo
        endif
c
      enddo
c$OMP END PARALLEL DO
c
c --- get id of field variable
      varname='vbaro'
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
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=idm_o
      count(2)=jdm_o
      count(3)=1
      status=nf_get_vara_int2(ncid,varid,start,count,i2fld)
      if (status.ne.nf_noerr) then
        write(*,'(4a)') 'nf_get_vara_int2: ',varname,': ',
     .                  nf_strerror(status)
        stop
      endif
c
c --- scale and add offset to field
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
      call intpol_nest(idm_o,jdm_o,idm,jdm,i1o,i2o,j1o,j2o,
     .                 mval,fval,mask_o,iv,fld_o,fld)
c
c --- extract boundary values
c$OMP PARALLEL DO
      do j=1,jj
        vb_ow(j)=fld(1,j)
        vb_oe(j)=fld(ii-1,j)
      enddo
c$OMP END PARALLEL DO
c$OMP PARALLEL DO
      do i=1,ii
        vb_os(i)=fld(i,2)
        vb_on(i)=fld(i,jj-1)
      enddo
c$OMP END PARALLEL DO
c
c --- ensure "flux conservation" trough boundary
c$OMP PARALLEL DO PRIVATE(q)
      do i=i1o,i2o
c
        if (fld_o(i,j1o).ne.mval) then
          q=fld_o(i,j1o)*nfac
          do k=1,nfac
            q=q-vb_os((i-i1o)*nfac+k+1)
          enddo
          q=q/nfac
          do k=1,nfac
            vb_os((i-i1o)*nfac+k+1)=vb_os((i-i1o)*nfac+k+1)+q
          enddo
        endif
c
        if (fld_o(i,j2o+1).ne.mval) then
          q=fld_o(i,j2o+1)*nfac
          do k=1,nfac
            q=q-vb_on((i-i1o)*nfac+k+1)
          enddo
          q=q/nfac
          do k=1,nfac
            vb_on((i-i1o)*nfac+k+1)=vb_on((i-i1o)*nfac+k+1)+q
          enddo
        endif
c
      enddo
c$OMP END PARALLEL DO
c
      return
      end
