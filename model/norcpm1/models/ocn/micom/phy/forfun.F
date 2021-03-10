      subroutine forfun
c
c --- ------------------------------------------------------------------
c --- Generate monthly climatological fields and prepare SST/SSS
c --- relaxation
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_atm.h"
#include "common_atm.h"
#include "common_forc.h"
#include "common_geo.h"
#include "netcdf.inc"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field
      real, dimension(atm_idm) :: rlon
      real, dimension(atm_jdm) :: rlat
      real dx2,dy2
      real*4, dimension(atm_idm,atm_jdm) :: tmp_r4
      real, dimension(itdm,jtdm) :: tmp2da,tmp2db
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
      integer, dimension(3) :: start,count
      integer i,j,k,l,atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t,mbdy,
     .        status,ncid,dimid,varid
      character filename*120
c
c --- ------------------------------------------------------------------
c --- Prepare interpolation of atmospheric surface fields
c --- ------------------------------------------------------------------
c
c --- get global coordinates of p-points
      call xcaget(tmp2da,plat,1)
      call xcaget(tmp2db,plon,1)
c
      if (mnproc.eq.1) then
c
c --- - read atmospheric lat/lon and land mask information
        filename=path(1:path_len)//'land.sfc.gauss.nc'
        call read_atm_info(filename,atm_idm,atm_jdm,rlon,rlat,atm_field)
c$OMP PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_lon(i,j)=rlon(i)
            atm_lat(i,j)=rlat(j)
            if (atm_field(i,j).lt..5) then
              atm_msk(i,j)=1 ! ocean
            else
              atm_msk(i,j)=0 ! land
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - read atmospheric model topography
        filename=path(1:path_len)//'hgt.sfc.nc'
        call read_atm_topo(filename,atm_idm,atm_jdm,atm_topo)
c
c --- - read interpolation weights if they exist, otherwise compute them
        filename=path(1:path_len)//'atm_intwgt.uf'
        open (unit=nfu,file=filename,form='unformatted',status='old',
     .        err=10)
        read (nfu) atm_idm_t,atm_jdm_t,idm_t,jdm_t,atm_nwgt_t
        if (atm_idm_t.ne.atm_idm.or.atm_jdm_t.ne.atm_jdm.or.
     .      idm_t.ne.itdm.or.jdm_t.ne.jtdm.or.
     .      atm_nwgt_t.ne.atm_nwgt) then
          close (unit=nfu)
        else
          read (nfu) atm_wgt,atm_iwgt,atm_jwgt
          close (unit=nfu)
          goto 20
        endif
c 
 10     write (lp,*)
     .    'forfun: computing atmospheric interpolation weights...'
c
c --- - compute weights
        call compweights(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                   atm_lat,atm_lon,tmp2da,tmp2db,
     .                   atm_wgt,atm_iwgt,atm_jwgt)
        filename=path(1:path_len)//'atm_intwgt.uf'
        open (unit=nfu,file=filename,form='unformatted')
        write (nfu) atm_idm,atm_jdm,itdm,jtdm,atm_nwgt
        write (nfu) atm_wgt,atm_iwgt,atm_jwgt
        close (unit=nfu)
c
 20     continue
c
      endif
c
c --- get mask for the full domain
      do j=1,jj
        do i=1,ii
          util1(i,j)=ip(i,j)
        enddo
      enddo
      call xcaget(tmp2da,util1,1)
      if (mnproc.eq.1) then
        do j=1,jtdm
          do i=1,itdm
            itp(i,j)=tmp2da(i,j)
          enddo
        enddo
      endif
c
c --- Compute smooting weight atm_ice_swgt. For stability
c --- atm_ice_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
      atm_ice_swgt=huge
c$OMP PARALLEL DO PRIVATE(dx2,dy2) REDUCTION(min:atm_ice_swgt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          atm_ice_swgt=min(atm_ice_swgt,.5*dx2*dy2/(dx2+dy2))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcminr(atm_ice_swgt)
      atm_ice_swgt=.9*atm_ice_swgt
c
c --- number of smoothing iterations is choosen to get a geographical
c --- extent of the smoothing independent of the grid resolution
      atm_ice_nsmt=nint(atm_ice_csmt/atm_ice_swgt)
      if (mnproc.eq.1) then
        write (lp,*) 'atm_ice_swgt',atm_ice_swgt,
     .               'atm_ice_nsmt',atm_ice_nsmt
      endif
c
c --- ------------------------------------------------------------------
c --- Read climatological ice concentration and surface temperature from
c --- atmospheric model
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*) 'reading atm. climatological ice concentration...'
        open (unit=nfu,file=path(1:path_len)//'icec_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read(nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_icec(i,j,k)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_icec(1,1,k),tmp2da)
        endif
        call xcaput(tmp2da,ricclm(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   ricclm(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(ricclm, 1,12, nbdy,nbdy, halo_ps)
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological surface temperature...'
        open (unit=nfu,file=path(1:path_len)//'skt_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_sktclm(i,j,k)=tmp_r4(i,j)+.0065*atm_topo(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_sktclm(1,1,k),tmp2da)
        endif
        call xcaput(tmp2da,sstclm(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   sstclm(1-nbdy,1-nbdy,k))
c
c --- - create smoothing mask - smooth where ice conc. is above 0.5
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            if (ricclm(i,j,k).gt..5.and.ip(i,j).eq.1) then
              smtmsk(i,j)=1
            else
              smtmsk(i,j)=0
            endif
          enddo
        enddo
c
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(sstclm(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                sstclm(1-nbdy,1-nbdy,k),smtmsk,
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(sstclm, 1,12, nbdy,nbdy, halo_ps)
c
c
      if (srxday.gt.epsil) then
c
c --- ------------------------------------------------------------------
c --- - read monthly climatological SSS (psu)
c --- ------------------------------------------------------------------
c
        if (mnproc.eq.1) then
          if     (srxsrc.eq.'PHC3.0') then
            filename='sss_clim_phc3.0.nc'
          elseif (srxsrc.eq.'CORE') then
            filename='sss_clim_core.nc'
          endif
          write (lp,'(2a)') ' reading monthly climatological SSS from ',
     .                      path(1:path_len)//trim(filename)
          call flush(lp)
c
c --- --- open netcdf file
          status=nf_open(path(1:path_len)//trim(filename),
     .                   nf_nowrite,ncid)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_open: '//trim(filename),
     .                       nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
c
c --- --- check dimensions
          status=nf_inq_dimid(ncid,'x',dimid)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_inq_dimid: x: ',nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
          status=nf_inq_dimlen(ncid,dimid,i)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_inq_dimlen: x: ',nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
          status=nf_inq_dimid(ncid,'y',dimid)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_inq_dimid: y: ',nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
          status=nf_inq_dimlen(ncid,dimid,j)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_inq_dimlen: y: ',nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,*) 'wrong dimensions in '//trim(filename)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
c
          status=nf_inq_varid(ncid,'sss',varid)
          if (status.ne.nf_noerr) then
            write(lp,'(2a)') 'nf_inq_varid: sss: ',nf_strerror(status)
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
c
          start(1)=1
          start(2)=1
          count(1)=itdm
          count(2)=jtdm
          count(3)=1
c
        endif
c
        do k=1,12
          if (mnproc.eq.1) then
            start(3)=k
            status=nf_get_vara_double(ncid,varid,start,count,tmp2da)
            if (status.ne.nf_noerr) then
              write(lp,'(2a)') 'nf_get_vara_double: sss: ',
     .                         nf_strerror(status)
              call xchalt('(forfun)')
                     stop '(forfun)'
            endif
          endif
          call xcaput(tmp2da,sssclm(1-nbdy,1-nbdy,k),1)
        enddo
c
      else
        do k=1,12
          do j=1,jj
            do i=1,ii
              sssclm(i,j,k)=huge
            enddo
          enddo
        enddo
      endif
c
      call xctilr(sssclm, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- initialize diagnosing/application of relaxation fluxes
c --- ------------------------------------------------------------------
c
      if (aptflx) then
        if (mnproc.eq.1) then
          open (unit=nfu,file=path(1:path_len)//'tflxdi.uf',
     .          form='unformatted')
          read (nfu) i,j
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,*) 'wrong dimensions in tflxdi.uf'
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
        endif
        do k=1,48
          if (mnproc.eq.1) then
            read (nfu) tmp2da
          endif
          call xcaput(tmp2da,tflxap(1-nbdy,1-nbdy,k),1)
        enddo
        if (mnproc.eq.1) then
          close (unit=nfu)
        endif
        call xctilr(tflxap, 1,48, nbdy,nbdy, halo_ps)
      endif
      if (apsflx) then
        if (mnproc.eq.1) then
          open (unit=nfu,file=path(1:path_len)//'sflxdi.uf',
     .          form='unformatted')
          read (nfu) i,j
          if (i.ne.itdm.or.j.ne.jtdm) then
            write (lp,*) 'wrong dimensions in sflxdi.uf'
            call xchalt('(forfun)')
                   stop '(forfun)'
          endif
        endif
        do k=1,48
          if (mnproc.eq.1) then
            read (nfu) tmp2da
          endif
          call xcaput(tmp2da,sflxap(1-nbdy,1-nbdy,k),1)
        enddo
        if (mnproc.eq.1) then
          close (unit=nfu)
        endif
        call xctilr(sflxap, 1,48, nbdy,nbdy, halo_ps)
      endif
c
      if (ditflx.or.disflx) then
        do k=1,48
          nflxdi(k)=0
        enddo
        if (ditflx) then
c$OMP PARALLEL DO
          do k=1,48
            do j=1-nbdy,jj+nbdy
              do i=1-nbdy,ii+nbdy
                tflxdi(i,j,k)=0.
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
        if (disflx) then
c$OMP PARALLEL DO
          do k=1,48
            do j=1-nbdy,jj+nbdy
              do i=1-nbdy,ii+nbdy
                sflxdi(i,j,k)=0.
              enddo
            enddo
          enddo
c$OMP END PARALLEL DO
        endif
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'forfun:'
        endif
        call chksummsk(ricclm,ip,12,'ricclm')
        call chksummsk(sstclm,ip,12,'sstclm')
        call chksummsk(sssclm,ip,12,'sssclm')
      endif
c
      return
      end
