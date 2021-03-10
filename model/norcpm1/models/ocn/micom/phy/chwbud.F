      subroutine chwbud
c
c --- ------------------------------------------------------------------
c --- Compute global surface budgets of climatological heat and
c --- freshwater
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_geo.h"
#include "common_clim.h"
#include "common_forc.h"
#include "common_ice.h"
c
      real prc_sum,eva_sum,rnf_sum,swa_sum,lwa_sum,lht_sum,sht_sum,
     .     fwf_fac,dangle,garea,le,albedo,fac,
     .     swa_ave,lwa_ave,lht_ave,sht_ave,crnf,cswa,q
      integer i,j,k,l
c
      fwf_fac=1.e-3 ! conversion factor kg/m^2/s -> m/s
c
c --- sum up freshwater fluxes
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
        enddo
        enddo
      enddo
c
      do k=1,12
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea=scp2(i,j)*1.e-4 ! [m^2]
c
c --- ----- freshwater fluxes [m/s]
            util1(i,j)=util1(i,j)+precip(i,j,k)*fwf_fac*garea
            le=(2.501-.00237*(sstclm(i,j,k)-273.15))*1.e6
            util2(i,j)=util2(i,j)-(lhtflx(i,j,k)/le)*1.e-3*garea
            util3(i,j)=util3(i,j)+runoff(i,j,k)*fwf_fac*garea
          enddo
          enddo
        enddo
      enddo
      call xcsum(prc_sum,util1,ip)
      call xcsum(eva_sum,util2,ip)
      call xcsum(rnf_sum,util3,ip)
c
      fac=1.e-6/12.
      prc_sum=prc_sum*fac
      eva_sum=eva_sum*fac
      rnf_sum=rnf_sum*fac
c
c --- sum up heat fluxes
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          util1(i,j)=0.
          util2(i,j)=0.
          util3(i,j)=0.
          util4(i,j)=0.
        enddo
        enddo
      enddo

      do k=1,12
c
c --- - compute albedo
        dangle=8.*atan(1.)*real(30*k-15)/360.
        call albw_eval(dangle,plat,clouds(1-nbdy,1-nbdy,k),albw_d,albw)
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            garea=scp2(i,j)*1.e-4 ! [m^2]
c
c --- ----- heat fluxes
            albedo=albs_f*ricclm(i,j,k)+albw(i,j)*(1.-ricclm(i,j,k))
            util1(i,j)=util1(i,j)+dswrfl(i,j,k)*(1.-albedo)*garea
            util2(i,j)=util2(i,j)-nlwrfs(i,j,k)*garea
            util3(i,j)=util3(i,j)-lhtflx(i,j,k)*garea
            util4(i,j)=util4(i,j)-shtflx(i,j,k)*garea
          enddo
          enddo
        enddo
      enddo
c
      call xcsum(swa_sum,util1,ip)
      call xcsum(lwa_sum,util2,ip)
      call xcsum(lht_sum,util3,ip)
      call xcsum(sht_sum,util4,ip)
c
      fac=1.e4/(12.*area)
      swa_ave=swa_sum*fac
      lwa_ave=lwa_sum*fac
      lht_ave=lht_sum*fac
      sht_ave=sht_sum*fac
c
      if (mnproc.eq.1) then
      write (lp,*)
      write (lp,*) 'Global precipitation:         ',prc_sum,'Sv'
      write (lp,*) 'Global evaporation:           ',eva_sum,'Sv'
      write (lp,*) 'Global runoff:                ',rnf_sum,'Sv'
      write (lp,*) 'Global balance of freshwater: ',
     .            prc_sum+eva_sum+rnf_sum,'Sv'
      write (lp,*)
      write (lp,*) 'Global mean short-wave radiation: ',swa_ave,'W/m^2'
      write (lp,*) 'Global mean long-wave radiation:  ',lwa_ave,'W/m^2'
      write (lp,*) 'Global mean latent heat flux:     ',lht_ave,'W/m^2'
      write (lp,*) 'Global mean sensible heat flux:   ',sht_ave,'W/m^2'
      write (lp,*) 'Global balance of mean heat-flux: ',
     .           swa_ave+lwa_ave+lht_ave+sht_ave,'W/m^2'
      write (lp,*)
      endif
c
      crnf=-(prc_sum+eva_sum)/rnf_sum
      cswa=-(lwa_ave+lht_ave+sht_ave)/swa_ave
      do k=1,12
        do j=1-nbdy,jj+nbdy
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,k)=runoff(i,j,k)*crnf
            dswrfl(i,j,k)=dswrfl(i,j,k)*cswa
          enddo
          enddo
        enddo
      enddo
      if (mnproc.eq.1) then
      write (lp,*) 'Runoff has been adjusted by a factor',crnf
      write (lp,*) 'Short-wave radiation has been adjusted by a factor',
     .              cswa
      write (lp,*)
      endif
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'chwbud:'
        endif
        call chksummsk(runoff,ip,12,'runoff')
        call chksummsk(dswrfl,ip,12,'dswrfl')
      endif
c
      return
      end
