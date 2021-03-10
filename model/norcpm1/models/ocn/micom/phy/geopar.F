      subroutine geopar
c
c --- set up model parameters related to geography (NERSC version)
c
c --- version 2.8 -- cyclic and noncyclic b.c. combined
c
      use mod_xc
      use mod_nctools
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
#include "common_forc.h"
#include "common_ice.h"
#include "common_asf.h"
#include "netcdf.inc"
c
      real mval,fval
      parameter (mval=-1.e12,fval=-1.e13)
c
      real, dimension(itdm,jtdm) :: tmpg
      real rnwp,rmxnbp,rtnbp,rnbp,dx2,dy2,btdtmx,
     .     umaxmin,vmaxmin,umaxmax,vmaxmax
      integer, dimension(3) :: start,count
      integer i,j,k,l,status,ncid,dimid,varid,nwp,kmax
c
#ifdef TKE
#  include "common_tke.h"
#endif
#ifdef TRC
#  include "param_trc.h"
#  include "common_trc.h"
      integer nt
#endif
c
c --- ------------------------------------------------------------------
c --- read grid information from grid.nc
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading grid information from ',
     .                    path(1:path_len)//'grid.nc'
        call flush(lp)
c
c --- - open netcdf file
        status=nf_open(path(1:path_len)//'grid.nc',nf_nowrite,
     .                 ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open: grid.nc: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
c --- - check dimensions
        status=nf_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: x: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimlen(ncid,dimid,i)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: x: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: y: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimlen(ncid,dimid,j)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: y: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in grid.nc'
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
c --- - read bathymetry
        status=nf_inq_varid(ncid,'pdepth',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: pdepth: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: pdepth: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
c --- - count number of wet points for subsequent xcsum testing
        nwp=0
        do j=1,jtdm
          do i=1,itdm
            if (tmpg(i,j).gt.0.) nwp=nwp+1
          enddo
        enddo
      endif
c
      call xcaput(tmpg,depths,1)
c
c --- read grid coordinates
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qlon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qlon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: qlon: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,qlon,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qlat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qlat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: qlat: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,qlat,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'plon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: plon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: plon: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,plon,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'plat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: plat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: plat: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,plat,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'ulon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: ulon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: ulon: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,ulon,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'ulat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: ulat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: ulat: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,ulat,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vlon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vlon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: vlon: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,vlon,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vlat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vlat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: vlat: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,vlat,1)
c
      start(1)=1
      start(2)=1
      count(1)=itdm
      count(2)=jtdm
      count(3)=1
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qclon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qclon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: qclon: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,qclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qclat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qclat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: qclat: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,qclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'pclon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: pclon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: pclon: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,pclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'pclat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: pclat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: pclat: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,pclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'uclon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: uclon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: uclon: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,uclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'uclat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: uclat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: uclat: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,uclat(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vclon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vclon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vclon: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,vclon(1-nbdy,1-nbdy,k),1)
      enddo
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vclat',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vclat: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      do k=1,4
        if (mnproc.eq.1) then
          start(3)=k
          status=nf_get_vara_double(ncid,varid,start,count,tmpg)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_get_vara_double: vclat: ',
     .                       nf_strerror(status)
            call xchalt('(geopar)')
                   stop '(geopar)'
          endif
        endif
        call xcaput(tmpg,vclat(1-nbdy,1-nbdy,k),1)
      enddo
c
c --- read scale factors
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qdx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qdx: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: qdx: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scqx,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qdy',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qdy: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: qdy: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scqy,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'pdx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: pdx: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: pdx: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scpx,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'pdy',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: pdy: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: pdy: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scpy,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'udx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: udx: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: udx: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scux,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'udy',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: udy: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: udy: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scuy,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vdx',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vdx: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: vdx: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scvx,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'vdy',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: vdy: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: vdy: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scvy,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'qarea',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: qarea: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: qarea: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scq2,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'parea',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: parea: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: parea: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scp2,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'uarea',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: uarea: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: uarea: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scu2,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'varea',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: varea: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: varea: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,scv2,1)
c
#ifdef CCSMCOUPLED
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'angle',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: angle: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: angle: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,angle,1)
c
      if (mnproc.eq.1) then
        status=nf_inq_varid(ncid,'cplmask',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: cplmask: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      if (mnproc.eq.1) then
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: cplmask: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
      call xcaput(tmpg,util3,1)
#endif
c
c --- close grid information file
c
      if (mnproc.eq.1) then
        status=nf_close(ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_close: grid.nc: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
c
c --- ------------------------------------------------------------------
c --- define grid point coordinates, scale factors, and Coriolis
c --- parameter
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do i=1,ii
c
          scqx(i,j)=scqx(i,j)*1.e2
          scqy(i,j)=scqy(i,j)*1.e2
          scpx(i,j)=scpx(i,j)*1.e2
          scpy(i,j)=scpy(i,j)*1.e2
          scux(i,j)=scux(i,j)*1.e2
          scuy(i,j)=scuy(i,j)*1.e2
          scvx(i,j)=scvx(i,j)*1.e2
          scvy(i,j)=scvy(i,j)*1.e2
          scq2(i,j)=scq2(i,j)*1.e4
          scp2(i,j)=scp2(i,j)*1.e4
          scu2(i,j)=scu2(i,j)*1.e4
          scv2(i,j)=scv2(i,j)*1.e4
c
          scq2i(i,j)=1./max(1.,scq2(i,j))
          scp2i(i,j)=1./max(1.,scp2(i,j))
          scuxi(i,j)=1./max(1.,scux(i,j))
          scvyi(i,j)=1./max(1.,scvy(i,j))
          scuyi(i,j)=1./max(1.,scuy(i,j))
          scvxi(i,j)=1./max(1.,scvx(i,j))
c
          corioq(i,j)=sin(qlat(i,j)/radian)*4.*pi/86164.
          coriop(i,j)=sin(plat(i,j)/radian)*4.*pi/86164.
          betafp(i,j)=cos(plat(i,j)/radian)*4.*pi/(86164.*rearth)
c
        enddo
      enddo
c
c --- ------------------------------------------------------------------
c --- determine do-loop limits for u,v,p,q points
c --- ------------------------------------------------------------------
c
      call bigrid(depths)
c
      call xctilr(qlat, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(qlon, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(plat, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(plon, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(ulat, 1,1, nbdy,nbdy, halo_us)
      call xctilr(ulon, 1,1, nbdy,nbdy, halo_us)
      call xctilr(vlat, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(vlon, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scqx, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scqy, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scpx, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scpy, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scux, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scuy, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvx, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scvy, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scq2, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scp2, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scu2, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scv2, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scq2i, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(scp2i, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(scuxi, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvyi, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(scuyi, 1,1, nbdy,nbdy, halo_us)
      call xctilr(scvxi, 1,1, nbdy,nbdy, halo_vs)
      call xctilr(corioq, 1,1, nbdy,nbdy, halo_qs)
      call xctilr(coriop, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(betafp, 1,1, nbdy,nbdy, halo_ps)
      call xctilr(qclat, 1,4, nbdy,nbdy, halo_qs)
      call xctilr(qclon, 1,4, nbdy,nbdy, halo_qs)
      call xctilr(pclat, 1,4, nbdy,nbdy, halo_ps)
      call xctilr(pclon, 1,4, nbdy,nbdy, halo_ps)
      call xctilr(uclat, 1,4, nbdy,nbdy, halo_us)
      call xctilr(uclon, 1,4, nbdy,nbdy, halo_us)
      call xctilr(vclat, 1,4, nbdy,nbdy, halo_vs)
      call xctilr(vclon, 1,4, nbdy,nbdy, halo_vs)
#ifdef CCSMCOUPLED
      call xctilr(util3, 1,1, nbdy,nbdy, halo_ps)
c
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          cplmsk(i,j)=nint(util3(i,j))
        enddo
      enddo
#endif
c
c --- ------------------------------------------------------------------
c --- set mask used for global sums
c --- ------------------------------------------------------------------
c
      if (nreg.eq.2.and.nproc.eq.jpr) then
        do j=1-nbdy,jj-1
          do i=1-nbdy,ii+nbdy
            ips(i,j)=ip(i,j)
          enddo
        enddo
        do j=jj,jj+nbdy
          do i=1-nbdy,ii+nbdy
            ips(i,j)=0
          enddo
        enddo
      else
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            ips(i,j)=ip(i,j)
          enddo
        enddo
      endif
c
      do j=1,jj
        do i=1,ii
          util1(i,j)=ip(i,j)
        enddo
      enddo
      call xcsum(rnwp,util1,ip)
      do j=1,jj
        do i=1,ii
          util1(i,j)=depths(i,j)*scp2(i,j)
          util2(i,j)=scp2(i,j)
        enddo
      enddo
      call xcsum(avgbot,util1,ips)
      call xcsum(area,  util2,ips)
      avgbot=avgbot/area
      if     (mnproc.eq.1) then
        if (nwp.ne.nint(rnwp)) then
          write (lp,'(a)') ' xcsum test failed!'
          write (lp,'(a,i7)') ' number of wet points:',nwp
          write (lp,'(a,i7)') ' xcsum on ocean mask: ',nint(rnwp)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        write (lp,100) avgbot,area
        call flush(lp)
 100    format(' mean basin depth (m) and area (10^6 km^2):',f9.1,
     .         -16p,f9.1)
      endif
c
c --- ------------------------------------------------------------------
c --- set mask for grid cells connected to the world ocean
c --- ------------------------------------------------------------------
c
      do j=1,jj
        do i=1,ii
          util1(i,j)=ips(i,j)
        enddo
      enddo
      call xcsum(rnwp,util1,ips)
      if (mnproc.eq.1) then
        write (lp,*) 'Number of wet points',nint(rnwp)
        call flush(lp)
      endif
c
      do j=1,jj
        do i=1,ii
          if (ips(i,j).eq.1) then
            util1(i,j)=fval
          else
            util1(i,j)=mval
          endif
        enddo
      enddo
c
      k=0
      rmxnbp=0.
      rtnbp=0.
      do
        k=k+1
        call xcaget(tmpg,util1,1)
        if (mnproc.eq.1) then
          do l=1,itdm*jtdm
            j=(l-1)/itdm+1
            i=l-(j-1)*itdm
            if (tmpg(i,j).eq.fval) then
              tmpg(i,j)=k
              exit
            endif
          enddo
        endif
        call xcaput(tmpg,util1,1)
        call fill_global(mval,fval,halo_ps,util1)
        do j=1,jj
          do i=1,ii
            if (util1(i,j).eq.mval.or.util1(i,j).eq.fval) then
              util2(i,j)=0.
            else
              if (nint(util1(i,j)).eq.k) then
                util2(i,j)=1.
              else
                util2(i,j)=0.
              endif
            endif
          enddo
        enddo
        call xcsum(rnbp,util2,ips)
        if (mnproc.eq.1) then
          write (lp,*) 'Number of basin points',nint(rnbp)
          call flush(lp)
        endif
        if (rnbp.gt.rmxnbp) then
          rmxnbp=rnbp
          kmax=k
        endif
        rtnbp=rtnbp+rnbp
        if (nint(rtnbp-rnwp).eq.0) exit
      enddo
c
      do j=1,jj
        do i=1,ii
          if (util1(i,j).eq.mval) then
            util1(i,j)=0.
          else
            if (nint(util1(i,j)).eq.kmax) then
              util1(i,j)=1.
            else
              util1(i,j)=0.
            endif
          endif
        enddo
      enddo
      call xctilr(util1, 1,1, nbdy,nbdy, halo_ps)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          ipwocn(i,j)=nint(util1(i,j))
        enddo
      enddo
c      
      call ncfopn('ipwocn.nc','w','c',1)
      call ncdims('x',itdm)
      call ncdims('y',jtdm)
      call ncwrti('ipwocn','x y',ipwocn,ip,1)
      call ncfcls
c
c --- ------------------------------------------------------------------
c --- read tidal wave energy dissipation divided by bottom buoyancy
c --- frequency
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,'(2a)') ' reading tidal dissipation data from ',
     .                    path(1:path_len)//'tidal_dissipation.nc'
        call flush(lp)
c
c --- - open netcdf file
        status=nf_open(path(1:path_len)//'tidal_dissipation.nc',
     .                 nf_nowrite,ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_open: tidal_dissipation.nc: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
c --- - check dimensions
        status=nf_inq_dimid(ncid,'x',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: x: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimlen(ncid,dimid,i)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: x: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimid(ncid,'y',dimid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimid: y: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_inq_dimlen(ncid,dimid,j)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_dimlen: y: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        if (i.ne.itdm.or.j.ne.jtdm) then
          write (lp,*) 'wrong dimensions in tidal_dissipation.nc'
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
c --- - read tidal dissipation data
        status=nf_inq_varid(ncid,'twedon',varid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_inq_varid: twedon: ',nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
        status=nf_get_var_double(ncid,varid,tmpg)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_get_var_double: twedon: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
c
        status=nf_close(ncid)
        if (status.ne.nf_noerr) then
          write(lp,'(2a)') 'nf_close: tidal_dissipation.nc: ',
     .                     nf_strerror(status)
          call xchalt('(geopar)')
                 stop '(geopar)'
        endif
      endif
c
      call xcaput(tmpg,twedon,1)
c
c --- change units from W*s/m**2=kg/s**2 to g/s**2
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          twedon(i,j)=twedon(i,j)*1.e3
        enddo
        enddo
      enddo
c
      call xctilr(twedon, 1,1, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- Determine upper bound of lateral diffusivity based on numerical
c --- stability concerns
c --- ------------------------------------------------------------------
c
      do 170 j=1-nbdy,jj+nbdy
      do 170 i=1-nbdy,ii+nbdy
      dx2=scpx(i,j)*scpx(i,j)
      dy2=scpy(i,j)*scpy(i,j)
      difmxp(i,j)=.9*.5*dx2*dy2/max(1.,(dx2+dy2)*(baclin+baclin))
      dx2=scqx(i,j)*scqx(i,j)
      dy2=scqy(i,j)*scqy(i,j)
 170  difmxq(i,j)=.9*.5*dx2*dy2/max(1.,(dx2+dy2)*(baclin+baclin))
c
c --- ------------------------------------------------------------------
c --- estimate maximum barotropic time step
c --- ------------------------------------------------------------------
c
      btdtmx=86400.
      do 180 j=1,jj
      do 180 l=1,isp(j)
      do 180 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
 180  btdtmx=min(btdtmx,
     .           scpx(i,j)*scpy(i,j)
     .           /sqrt(g*depths(i,j)*100.*(scpx(i,j)*scpx(i,j)
     .                                    +scpy(i,j)*scpy(i,j))))
      call xcminr(btdtmx)
      if (mnproc.eq.1) then
        write (lp,*) 'estimated max. barotropic time step:',
     .               btdtmx/sqrt(2.)
        call flush(lp)
      endif
c
c --- ------------------------------------------------------------------
c --- set maximum velocities allowable ensuring stability of the upwind
c --- scheme
c --- ------------------------------------------------------------------
c
      umaxmin=huge
      vmaxmin=huge
      umaxmax=0.
      vmaxmax=0.
      do 190 j=1,jj
      do 191 l=1,isu(j)
      do 191 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      umax(i,j)=.9*.125*min(scp2(i-1,j),scp2(i,j))/(scuy(i,j)*baclin)
      umaxmin=min(umaxmin,umax(i,j))
 191  umaxmax=max(umaxmax,umax(i,j))
      do 190 l=1,isv(j)
      do 190 i=max(1,ifv(j,l)),min(ii,ilv(j,l))
      vmax(i,j)=.9*.125*min(scp2(i,j-1 ),scp2(i,j))/(scvx(i,j)*baclin)
      vmaxmin=min(vmaxmin,vmax(i,j))
 190  vmaxmax=max(vmaxmax,vmax(i,j))
c
      call xctilr(umax, 1,1, nbdy,nbdy, halo_us)
      call xctilr(vmax, 1,1, nbdy,nbdy, halo_vs)
c
      call xcminr(umaxmin)
      call xcmaxr(umaxmax)
      call xcminr(vmaxmin)
      call xcmaxr(vmaxmax)
c
      if (mnproc.eq.1) then
        write (lp,*) 'min/max umax:',umaxmin,umaxmax
        write (lp,*) 'min/max vmax:',vmaxmin,vmaxmax
        call flush(lp)
      endif
c
c --- ------------------------------------------------------------------
c --- initialize some arrays
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO
      do 209 j=1-nbdy,jj+nbdy
      do 209 i=1-nbdy,ii+nbdy
      p(i,j,1)=huge
      pu(i,j,1)=huge
      pv(i,j,1)=huge
      phi(i,j,1)=huge
      do k=1,3
      ubflxs(i,j,k)=huge
      vbflxs(i,j,k)=huge
      enddo
      do k=1,2
      pb(i,j,k)=huge
      ubflx(i,j,k)=huge
      vbflx(i,j,k)=huge
      pb_mn(i,j,k)=huge
      ubflx_mn(i,j,k)=huge
      vbflx_mn(i,j,k)=huge
      ub(i,j,k)=huge
      vb(i,j,k)=huge
      ubflxs_p(i,j,k)=huge
      vbflxs_p(i,j,k)=huge
      pbu(i,j,k)=huge
      pbv(i,j,k)=huge
      pvtrop(i,j,k)=huge
      enddo
      pb_p(i,j)=huge
      pbu_p(i,j)=huge
      pbv_p(i,j)=huge
      ubcors_p(i,j)=huge
      vbcors_p(i,j)=huge
      defor1(i,j)=huge
      defor2(i,j)=huge
      pvtrop_o(i,j)=huge
      utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfxm(i,j,1)=huge
      pgfxm(i,j,2)=huge
      pgfxm_o(i,j)=huge
      xixp(i,j,1)=huge
      xixp(i,j,2)=huge
      xixp_o(i,j)=huge
      xixm(i,j,1)=huge
      xixm(i,j,2)=huge
      xixm_o(i,j)=huge
      pgfym(i,j,1)=huge
      pgfym(i,j,2)=huge
      pgfym_o(i,j)=huge
      xiyp(i,j,1)=huge
      xiyp(i,j,2)=huge
      xiyp_o(i,j)=huge
      xiym(i,j,1)=huge
      xiym(i,j,2)=huge
      xiym_o(i,j)=huge
      difwgt(i,j)=huge
      sealv(i,j)=huge
      taux(i,j)=huge
      tauy(i,j)=huge
      swa(i,j)=huge
      nsf(i,j)=huge
      hmltfz(i,j)=huge
      hmlt(i,j)=huge
      dfl(i,j)=huge
      lip(i,j)=huge
      sop(i,j)=huge
      eva(i,j)=huge
      rnf(i,j)=huge
      rfi(i,j)=huge
      fmltfz(i,j)=huge
      sfl(i,j)=huge
      ztx(i,j)=huge
      mty(i,j)=huge
      ustarw(i,j)=huge
      tsi(i,j)=huge
      slp(i,j)=huge
      abswnd(i,j)=huge
      albw(i,j)=huge
      frzpot(i,j)=huge
      mltpot(i,j)=huge
      atmco2(i,j)=huge
      flxco2(i,j)=huge
      tsi_tda(i,j)=huge
      tml_tda(i,j)=huge
      sml_tda(i,j)=huge
      alb_tda(i,j)=huge
      fice_tda(i,j)=huge
      ssu_tda(i,j)=huge
      ssv_tda(i,j)=huge
      alb(i,j)=huge
      rnfres(i,j)=huge
      rnfflx(i,j)=huge
      rfiflx(i,j)=huge
      eiacc(i,j)=huge
      pracc(i,j)=huge
      ficem(i,j)=huge
      hicem(i,j)=huge
      tsrfm(i,j)=huge
      hsnwm(i,j)=huge
      ticem(i,j)=huge
      ustari(i,j)=0.
      tauxice(i,j)=0.
      tauyice(i,j)=0.
      uicem(i,j)=0.
      vicem(i,j)=0.
      iagem(i,j)=huge
#if defined(CLIM) || defined(SYN)
      cd_d(i,j)=huge
      ch_d(i,j)=huge
      ce_d(i,j)=huge
      wg2_d(i,j)=huge
      cd_m(i,j)=huge
      ch_m(i,j)=huge
      ce_m(i,j)=huge
      wg2_m(i,j)=huge
      rhoa(i,j)=huge
#endif
      do 209 k=1,kk
      u  (i,j,k   )=huge
      u  (i,j,k+kk)=huge
      v  (i,j,k   )=huge
      v  (i,j,k+kk)=huge
      uflx(i,j,k   )=huge
      uflx(i,j,k+kk)=huge
      utflx(i,j,k   )=huge
      utflx(i,j,k+kk)=huge
      usflx(i,j,k   )=huge
      usflx(i,j,k+kk)=huge
      vflx(i,j,k   )=huge
      vflx(i,j,k+kk)=huge
      vtflx(i,j,k   )=huge
      vtflx(i,j,k+kk)=huge
      vsflx(i,j,k   )=huge
      vsflx(i,j,k+kk)=huge
      umfltd(i,j,k   )=huge
      umfltd(i,j,k+kk)=huge
      vmfltd(i,j,k   )=huge
      vmfltd(i,j,k+kk)=huge
      dp (i,j,k   )=huge
      dp (i,j,k+kk)=huge
      dpu(i,j,k   )=huge
      dpu(i,j,k+kk)=huge
      dpuold(i,j,k)=huge
      dpv(i,j,k   )=huge
      dpv(i,j,k+kk)=huge
      dpvold(i,j,k)=huge
      p (i,j,k+1)=huge
      pu(i,j,k+1)=huge
      pv(i,j,k+1)=huge
      phi(i,j,k+1)=huge
c
      temp(i,j,k)=huge
      saln(i,j,k)=huge
      sigma(i,j,k   )=huge
      sigma(i,j,k+kk)=huge
      pgfx(i,j,k   )=huge
      pgfx(i,j,k+kk)=huge
      pgfy(i,j,k   )=huge
      pgfy(i,j,k+kk)=huge
      pgfxo(i,j,k)=huge
      pgfyo(i,j,k)=huge
      diaflx(i,j,k)=0.
#ifdef TKE
      tke(i,j,k)=huge
      tke(i,j,k+kk)=huge
      gls_psi(i,j,k)=huge
      gls_psi(i,j,k+kk)=huge
#endif
#ifdef TRC
      do nt=1,ntr
      trc(i,j,k,nt)=huge
      enddo
#endif
 209  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 210 j=1,jj+1
      do 210 l=1,isq(j)
      do 210 i=max(1,ifq(j,l)),min(ii+1,ilq(j,l))
      do k=1,2
      pb(i  ,j  ,k)=0.
      pb(i-1,j  ,k)=0.
      pb(i  ,j-1,k)=0.
      pb(i-1,j-1,k)=0.
      pb_mn(i  ,j  ,k)=0.
      pb_mn(i-1,j  ,k)=0.
      pb_mn(i  ,j-1,k)=0.
      pb_mn(i-1,j-1,k)=0.
      enddo
      pb_p(i  ,j  )=0.
      pb_p(i-1,j  )=0.
      pb_p(i  ,j-1)=0.
      pb_p(i-1,j-1)=0.
      p(i  ,j  ,1)=0.
      p(i-1,j  ,1)=0.
      p(i  ,j-1,1)=0.
      p(i-1,j-1,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,k   )=0.
      dp(i  ,j  ,k+kk)=0.
      dp(i-1,j  ,k   )=0.
      dp(i-1,j  ,k+kk)=0.
      dp(i  ,j-1,k   )=0.
      dp(i  ,j-1,k+kk)=0.
      dp(i-1,j-1,k   )=0.
 210  dp(i-1,j-1,k+kk)=0.
c$OMP END PARALLEL DO
c
      call xctilr(pb,   1,   2, nbdy,nbdy, halo_ps)
      call xctilr(pb_mn,1,   2, nbdy,nbdy, halo_ps)
      call xctilr(pb_p, 1,   1, nbdy,nbdy, halo_ps)
      call xctilr(p,    1,   1, nbdy,nbdy, halo_ps)
      call xctilr(dp,   1,2*kk, nbdy,nbdy, halo_ps)
c
c$OMP PARALLEL DO
      do 211 j=1-nbdy,jj+nbdy
      do 211 l=1,isp(j)
      do 211 i=ifp(j,l),ilp(j,l)
 211  p(i,j,1)=0.
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 212 j=1,jj
      do 212 l=1,isp(j)
      do 212 i=max(1,ifp(j,l)),min(ii,ilp(j,l))
      defor2(i  ,j  )=0.
      defor2(i+1,j  )=0.
      defor2(i  ,j+1)=0.
 212  defor2(i+1,j+1)=0.
c$OMP END PARALLEL DO
c
c --- initialize  u,ub,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  pbu,dpu,utotn  upstream and downstream of p
c --- points as well as at lateral neighbors of interior u points.
c
c$OMP PARALLEL DO
      do 156 j=0,jj+1
      do 156 l=1,isu(j)
      do 156 i=max(1,ifu(j,l)),min(ii,ilu(j,l))
      pu(i,j,1)=0.
c
      utotn(i,j-1)=0.
      utotn(i,j+1)=0.
c
      do k=1,2
      pbu(i,j-1,k)=0.
      pbu(i,j+1,k)=0.
      enddo
      pbu_p(i,j-1)=0.
      pbu_p(i,j+1)=0.
c
      do 156 k=1,kk
      dpu(i,j-1,k   )=0.
      dpu(i,j-1,k+kk)=0.
      dpuold(i,j-1,k)=0.
      dpu(i,j+1,k   )=0.
      dpu(i,j+1,k+kk)=0.
      dpuold(i,j+1,k)=0.
 156  continue
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do 158 j=1,jj
      do 158 l=1,isp(j)
      do 158 i=max(1,ifp(j,l)),min(ii,ilp(j,l)+1)
      do k=1,3
      ubflxs(i,j,k)=0.
      enddo
      do k=1,2
      ubflx(i,j,k)=0.
      ubflx_mn(i,j,k)=0.
      ub(i,j,k)=0.
      ubflxs_p(i,j,k)=0.
      pbu(i,j,k)=0.
      enddo
      pbu_p(i,j)=0.
      ubcors_p(i,j)=0.
      utotn (i,j)=0.
      utotm(i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      do 158 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
      dpuold(i,j,k)=0.
      uflx(i,j,k   )=0.
      uflx(i,j,k+kk)=0.
      utflx(i,j,k   )=0.
      utflx(i,j,k+kk)=0.
      usflx(i,j,k   )=0.
      usflx(i,j,k+kk)=0.
      umfltd(i,j,k   )=0.
      umfltd(i,j,k+kk)=0.
      u(i,j,k   )=0.
 158  u(i,j,k+kk)=0.
c$OMP END PARALLEL DO
c
      call xctilr(pu,     1,   1, nbdy,nbdy, halo_us)
      call xctilr(utotn,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(pbu,    1,   2, nbdy,nbdy, halo_us)
      call xctilr(pbu_p,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(dpu,    1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(dpuold, 1,  kk, nbdy,nbdy, halo_us)
      call xctilr(ub,     1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflx,  1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflx_mn, 1,   2, nbdy,nbdy, halo_us)
      call xctilr(ubflxs, 1,   3, nbdy,nbdy, halo_us)
      call xctilr(ubflxs_p, 1, 2, nbdy,nbdy, halo_us)
      call xctilr(ubcors_p, 1, 1, nbdy,nbdy, halo_us)
      call xctilr(utotm,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux,  1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux2, 1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflux3, 1,   1, nbdy,nbdy, halo_us)
      call xctilr(uja,    1,   1, nbdy,nbdy, halo_us)
      call xctilr(ujb,    1,   1, nbdy,nbdy, halo_us)
      call xctilr(uflx,   1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(utflx,  1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(usflx,  1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(umfltd, 1,2*kk, nbdy,nbdy, halo_us)
      call xctilr(u   ,   1,2*kk, nbdy,nbdy, halo_us)
#ifdef TKE
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          uflxtke(i,j)=uflux(i,j)
#  ifdef GLS
          uflxgls_psi(i,j)=uflux(i,j)
#  endif
        enddo
      enddo
#endif
#ifdef TRC
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          do nt=1,ntr
            uflxtr(nt,i,j)=uflux(i,j)
          enddo
        enddo
      enddo
#endif
c
c --- initialize  v,vb,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  pby,dpv,vtotn  upstream and downstream of p
c --- points as well as at lateral neighbors of interior v points.
c
      do 166 i=0,ii+1
      do 166 l=1,jsv(i)
      do 166 j=max(1,jfv(i,l)),min(jj,jlv(i,l))
      pv(i,j,1)=0.
c
      vtotn (i-1,j)=0.
      vtotn (i+1,j)=0.
c
      do k=1,2
      pbv(i-1,j,k)=0.
      pbv(i+1,j,k)=0.
      enddo
      pbv_p(i-1,j)=0.
      pbv_p(i+1,j)=0.
c
      do 166 k=1,kk
      dpv(i-1,j,k   )=0.
      dpv(i-1,j,k+kk)=0.
      dpvold(i-1,j,k)=0.
c
      dpv(i+1,j,k   )=0.
      dpv(i+1,j,k+kk)=0.
      dpvold(i+1,j,k)=0.
 166  continue
c
      do 168 i=1,ii
      do 168 l=1,jsp(i)
      do 168 j=max(1,jfp(i,l)),min(jj,jlp(i,l)+1)
      do k=1,3
      vbflxs(i,j,k)=0.
      enddo
      do k=1,2
      vbflx(i,j,k)=0.
      vbflx_mn(i,j,k)=0.
      vb(i,j,k)=0.
      vbflxs_p(i,j,k)=0.
      pbv(i,j,k)=0.
      enddo
      pbv_p(i,j)=0.
      vbcors_p(i,j)=0.
      vtotn (i,j)=0.
      vtotm(i,j)=0.
      vflux (i,j)=0.
      vflux2(i,j)=0.
      vflux3(i,j)=0.
      via(i,j)=0.
      vib(i,j)=0.
c
      do 168 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
      dpvold(i,j,k)=0.
      vflx(i,j,k   )=0.
      vflx(i,j,k+kk)=0.
      vtflx(i,j,k   )=0.
      vtflx(i,j,k+kk)=0.
      vsflx(i,j,k   )=0.
      vsflx(i,j,k+kk)=0.
      vmfltd(i,j,k   )=0.
      vmfltd(i,j,k+kk)=0.
      v(i,j,k   )=0.
 168  v(i,j,k+kk)=0.
c
      call xctilr(pv,     1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vtotn,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(pbv,    1,   2, nbdy,nbdy, halo_vs)
      call xctilr(pbv_p,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(dpv,    1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(dpvold, 1,  kk, nbdy,nbdy, halo_vs)
      call xctilr(vb,     1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflx,  1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflx_mn, 1,   2, nbdy,nbdy, halo_vs)
      call xctilr(vbflxs, 1,   3, nbdy,nbdy, halo_vs)
      call xctilr(vbflxs_p, 1, 2, nbdy,nbdy, halo_vs)
      call xctilr(vbcors_p, 1, 1, nbdy,nbdy, halo_vs)
      call xctilr(vtotm,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux,  1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux2, 1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflux3, 1,   1, nbdy,nbdy, halo_vs)
      call xctilr(via,    1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vib,    1,   1, nbdy,nbdy, halo_vs)
      call xctilr(vflx,   1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vtflx,  1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vsflx,  1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(vmfltd, 1,2*kk, nbdy,nbdy, halo_vs)
      call xctilr(v   ,   1,2*kk, nbdy,nbdy, halo_vs)
#ifdef TKE
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          vflxtke(i,j)=vflux(i,j)
#  ifdef GLS
          vflxgls_psi(i,j)=vflux(i,j)
#  endif
        enddo
      enddo
#endif
#ifdef TRC
      do j=1-nbdy,jdm+nbdy
        do i=1-nbdy,idm+nbdy
          do nt=1,ntr
            vflxtr(nt,i,j)=vflux(i,j)
          enddo
        enddo
      enddo
#endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'geopar:'
        endif
        call chksummsk(depths,ip,1,'depths')
        call chksummsk(plat,ip,1,'plat')
        call chksummsk(plon,ip,1,'plon')
        call chksummsk(pclat,ip,4,'pclat')
        call chksummsk(pclon,ip,4,'pclon')
        call chksummsk(corioq,iq,1,'corioq')
        call chksummsk(coriop,ip,1,'coriop')
        call chksummsk(betafp,ip,1,'betafp')
        call chksummsk(scqx,iq,1,'scqx')
        call chksummsk(scqy,iq,1,'scqy')
        call chksummsk(scpx,ip,1,'scpx')
        call chksummsk(scpy,ip,1,'scpy')
        call chksummsk(scux,iu,1,'scux')
        call chksummsk(scuy,iu,1,'scuy')
        call chksummsk(scvx,iv,1,'scvx')
        call chksummsk(scvy,iv,1,'scvy')
        call chksummsk(scq2,iq,1,'scq2')
        call chksummsk(scp2,ip,1,'scp2')
        call chksummsk(scu2,iu,1,'scu2')
        call chksummsk(scv2,iv,1,'scv2')
        call chksummsk(scp2i,ip,1,'scp2i')
        call chksummsk(scq2i,iq,1,'scq2i')
        call chksummsk(scuxi,iu,1,'scuxi')
        call chksummsk(scvyi,iv,1,'scvyi')
        call chksummsk(scuyi,iu,1,'scuyi')
        call chksummsk(scvxi,iv,1,'scvxi')
        call chksummsk(difmxp,ip,1,'difmxp')
        call chksummsk(difmxq,iq,1,'difmxq')
        call chksummsk(umax,iu,1,'umax')
        call chksummsk(vmax,iv,1,'vmax')
c       call chksummsk(p,ip,kk+1,'p')
c       call chksummsk(pu,iu,kk+1,'pu')
c       call chksummsk(pv,iv,kk+1,'pv')
c       call chksummsk(pb,ip,2,'pb')
c       call chksummsk(pb_mn,ip,2,'pb')
c       call chksummsk(pb_p,ip,1,'pb_p')
c       call chksummsk(ub,iu,2,'ub')
c       call chksummsk(vb,iv,2,'vb')
c       call chksummsk(pbu,iu,2,'pbu')
c       call chksummsk(pbu_p,iu,1,'pbu_p')
c       call chksummsk(pbv,iv,2,'pbv')
c       call chksummsk(pbv_p,iv,1,'pbv_p')
c       call chksummsk(ubflx,iu,2,'ubflx')
c       call chksummsk(ubflx_mn,iu,2,'ubflx')
c       call chksummsk(vbflx,iv,2,'vbflx')
c       call chksummsk(vbflx_mn,iv,2,'vbflx')
c       call chksummsk(ubflxs,iu,3,'ubflxs')
c       call chksummsk(ubflxs_p,iu,2,'ubflxs_p')
c       call chksummsk(vbflxs,iv,3,'vbflxs')
c       call chksummsk(vbflxs_p,iv,2,'vbflxs_p')
c       call chksummsk(ubcors_p,iu,1,'ubcors_p')
c       call chksummsk(vbcors_p,iv,1,'vbcors_p')
c       call chksummsk(utotm,iu,1,'utotm')
c       call chksummsk(vtotm,iv,1,'vtotm')
c       call chksummsk(utotn,iu,1,'utotn')
c       call chksummsk(vtotn,iv,1,'vtotn')
c       call chksummsk(uflux,iu,1,'uflux')
c       call chksummsk(vflux,iv,1,'vflux')
c       call chksummsk(uflux2,iu,1,'uflux2')
c       call chksummsk(vflux2,iv,1,'vflux2')
c       call chksummsk(uflux3,iu,1,'uflux3')
c       call chksummsk(vflux3,iv,1,'vflux3')
c       call chksummsk(uja,iu,1,'uja')
c       call chksummsk(ujb,iu,1,'ujb')
c       call chksummsk(via,iv,1,'via')
c       call chksummsk(vib,iv,1,'vib')
c       call chksummsk(u,iu,2*kk,'u')
c       call chksummsk(v,iv,2*kk,'v')
c       call chksummsk(uflx,iu,2*kk,'uflx')
c       call chksummsk(vflx,iv,2*kk,'vflx')
c       call chksummsk(dp,ip,2*kk,'dp')
c       call chksummsk(dpu,iu,2*kk,'dpu')
c       call chksummsk(dpuold,iu,kk,'dpuold')
c       call chksummsk(dpv,iv,2*kk,'dpv')
c       call chksummsk(dpvold,iv,kk,'dpvold')
      endif
c
      return
      end
c
c
c> Revision history:
c>
c> May  1997 - extended list of variables set to 'huge' on land
c> May  2000 - changed loop 56 from i/j to j/i to improve memory layout
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
