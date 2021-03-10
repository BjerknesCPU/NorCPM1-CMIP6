      subroutine forfun_clim
c
c --- ------------------------------------------------------------------------
c --- load atmospheric monthly climatological forcing fields
c --- ------------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#include "param_atm.h"
#include "common_atm.h"
#include "common_clim.h"
#include "common_geo.h"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field
      real dx2,dy2
      real*4, dimension(atm_abdm,atm_idm,atm_jdm) :: rnf_weight
      real*4, dimension(atm_idm,atm_jdm) :: tmp_r4
      real*4 rw4
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,12) :: smtmsk
      real, dimension(itdm,jtdm) :: tmp2d
      integer i,j,k,l,mbdy,il,jl
      integer*2, dimension(atm_abdm,atm_idm,atm_jdm) ::
     .  rnf_ocdpi,rnf_ocdpj
      integer*2 rn2,ri2,rj2
c
c --- ------------------------------------------------------------------
c --- create smoothing mask - smooth where ice conc. is above 0.5
c --- ------------------------------------------------------------------
c
      do k=1,12
c$OMP PARALLEL DO
        do j=1-nbdy,jj+nbdy
          do i=1-nbdy,ii+nbdy
            if (ricclm(i,j,k).gt..5.and.ip(i,j).eq.1) then
              smtmsk(i,j,k)=1
            else
              smtmsk(i,j,k)=0
            endif
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c --- ------------------------------------------------------------------
c --- read short-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological short-wave radiation flux...'
        open (unit=nfu,file=path(1:path_len)//'dswrf_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   dswrfl(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(dswrfl(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                dswrfl(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(dswrfl, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read net long-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write (lp,*)
     .    'reading atm. climatological long-wave radiation flux...'
        open (unit=nfu,file=path(1:path_len)//'nlwrs_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   nlwrfs(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(nlwrfs(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                nlwrfs(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(nlwrfs, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read total cloud cover [0-100%]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological total cloud cover...'
        open (unit=nfu,file=path(1:path_len)//'tcdc_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   clouds(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(clouds(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                clouds(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
c
c --- - convert range of cloudiness from 0-100 to 0-1
c$OMP PARALLEL DO
          do j=1,jj
            do l=1,isp(j)
            do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
              clouds(i,j,k)=clouds(i,j,k)*1.e-2
            enddo
            enddo
          enddo
c$OMP END PARALLEL DO
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(clouds, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read precipitation rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological precipitation rate...'
        open (unit=nfu,file=path(1:path_len)//'prate_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,precip(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   precip(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(precip(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                precip(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(precip, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read latent heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological latent heat flux...'
        open (unit=nfu,file=path(1:path_len)//'lhtfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   lhtflx(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(lhtflx(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                lhtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(lhtflx, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read sensible heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological sensible heat flux...'
        open (unit=nfu,file=path(1:path_len)//'shtfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   shtflx(1-nbdy,1-nbdy,k))
        mbdy=1
        do i=1,atm_ice_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(shtflx(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_ice_swgt,
     .                shtflx(1-nbdy,1-nbdy,k),smtmsk(1-nbdy,1-nbdy,k),
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(shtflx, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read surface pressure [Pa]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological surface pressure...'
        open (unit=nfu,file=path(1:path_len)//'pres_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=
     .          tmp_r4(i,j)*exp(9.81*atm_topo(i,j)
     .          /(287.*(atm_sktclm(i,j,k)-.00325*atm_topo(i,j))))
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   slpres(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(slpres, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read momentum flux [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological momentum flux...'
        open (unit=nfu,file=path(1:path_len)//'momfl_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,taud(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_ps,
     .                   taud(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(taud, 1,12, nbdy,nbdy, halo_ps)
c
c --- ------------------------------------------------------------------
c --- read momentum flux components [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        write(lp,*)
     .    'reading atm. climatological u-component of momentum flux...'
        open (unit=nfu,file=path(1:path_len)//'uflx_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_pv,
     .                   tauxd(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      if (mnproc.eq.1) then
        write(lp,*)
     .    'reading atm. climatological v-component of momentum flux...'
        open (unit=nfu,file=path(1:path_len)//'vflx_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
          call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                  atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                  atm_msk,itp,atm_field,tmp2d)
        endif
        call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,k),1)
        call fill_global(atm_mval,atm_fval,halo_pv,
     .                   tauyd(1-nbdy,1-nbdy,k))
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      do k=1,12
c
c --- - change sign of momentum flux
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            tauxd(i,j,k)=-tauxd(i,j,k)
            tauyd(i,j,k)=-tauyd(i,j,k)
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
c
c --- - rotate the vector components
        call xctilr(tauxd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
        call xctilr(tauyd(1-nbdy,1-nbdy,k), 1,1, 1,1, halo_pv)
        call uvrotr2g(atm_mval,atm_fval,plat,plon,iu,iv,
     .                tauxd(1-nbdy,1-nbdy,k),tauyd(1-nbdy,1-nbdy,k))
c
      enddo
c
      call xctilr(tauxd, 1,12, nbdy,nbdy, halo_uv)
      call xctilr(tauyd, 1,12, nbdy,nbdy, halo_vv)
c
c --- ------------------------------------------------------------------
c --- read runoff rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
c --- Compute smooting weight atm_rnf_swgt. For stability
c --- atm_rnf_swgt < .5*dx^2*dy^2/(dx^2+dy^2).
      atm_rnf_swgt=huge
c$OMP PARALLEL DO PRIVATE(dx2,dy2) REDUCTION(min:atm_rnf_swgt)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dx2=scpx(i,j)*scpx(i,j)
          dy2=scpy(i,j)*scpy(i,j)
          atm_rnf_swgt=min(atm_rnf_swgt,.5*dx2*dy2/(dx2+dy2))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
      call xcminr(atm_rnf_swgt)
      atm_rnf_swgt=.9*atm_rnf_swgt
c
c --- number of smoothing iterations is choosen to get a geographical
c --- extent of the smoothing independent of the grid resolution
      atm_rnf_nsmt=nint(atm_rnf_csmt/atm_rnf_swgt)
      if (mnproc.eq.1) then
        write (lp,*) 'atm_rnf_swgt',atm_rnf_swgt,
     .               'atm_rnf_nsmt',atm_rnf_nsmt
      endif
c
      if (mnproc.eq.1) then
        write(lp,*) 'reading atm. climatological runoff rate...'
c
c --- - read file containing ocean discarge points/weights for land
c --- - areas
        open (unit=nfu,file=path(1:path_len)//'runoffweights.uf',
     .        form='unformatted',status='old',action='read')
        do j=1,atm_jdm
          do i=1,atm_idm
            read (nfu) rn2
            do k=1,rn2
              read (nfu) rw4,ri2,rj2
              rnf_weight(k,i,j)=rw4
              rnf_ocdpi(k,i,j)=ri2
              rnf_ocdpj(k,i,j)=rj2
            enddo
            do k=rn2+1,atm_abdm
              rnf_weight(k,i,j)=0.
              rnf_ocdpi(k,i,j)=0
              rnf_ocdpj(k,i,j)=0
            enddo
          enddo
        enddo
        close (unit=nfu)
c
        open (unit=nfu,file=path(1:path_len)//'runof_1968-1996.uf',
     .        form='unformatted')
      endif
c
      do k=1,12
        if (mnproc.eq.1) then
          read (nfu) tmp_r4
c$OMP PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              atm_field(i,j)=tmp_r4(i,j)
            enddo
          enddo
c$OMP END PARALLEL DO
c
c --- --- place runoff at ocean discharge points
c$OMP PARALLEL DO
          do j=1,jtdm
            do i=1,itdm
              tmp2d(i,j)=0.
            enddo
          enddo
c$OMP END PARALLEL DO
          do j=1,atm_jdm
            do i=1,atm_idm
              do l=1,atm_abdm
                if (rnf_weight(l,i,j).gt.0.) then
                  il=rnf_ocdpi(l,i,j)
                  jl=rnf_ocdpj(l,i,j)
                  tmp2d(il,jl)=tmp2d(il,jl)
     .              +atm_field(i,j)*rnf_weight(l,i,j)
                endif
              enddo
            enddo
          enddo
        endif
        call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,k),1)
c
c --- - smooth the runoff field
        mbdy=1
        do i=1,atm_rnf_nsmt
          if (mbdy.eq.1) then
            mbdy=nbdy
            call xctilr(runoff(1-nbdy,1-nbdy,k), 1,1, nbdy,nbdy,halo_ps)
          else
            mbdy=mbdy-1
          endif
          call smooth(mbdy,atm_rnf_swgt,
     .                runoff(1-nbdy,1-nbdy,k),ip,
     .                scuy,scvx,scuxi,scvyi,scp2i)
        enddo
c
c --- - convert unit of runoff from kg/m^2/day to kg/m^2/s
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k)=runoff(i,j,k)/86400.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
      if (mnproc.eq.1) then
        close (unit=nfu)
      endif
c
      call xctilr(runoff, 1,12, nbdy,nbdy, halo_ps)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'forfun_clim:'
        endif
        call chksummsk(dswrfl,ip,12,'dswrfl')
        call chksummsk(nlwrfs,ip,12,'nlwrfs')
        call chksummsk(clouds,ip,12,'clouds')
        call chksummsk(precip,ip,12,'precip')
        call chksummsk(lhtflx,ip,12,'lhtflx')
        call chksummsk(shtflx,ip,12,'shtflx')
        call chksummsk(slpres,ip,12,'slpres')
        call chksummsk(taud,ip,12,'taud')
        call chksummsk(tauxd,iu,12,'tauxd')
        call chksummsk(tauyd,iv,12,'tauyd')
        call chksummsk(runoff,ip,12,'runoff')
      endif
c
      return
      end
