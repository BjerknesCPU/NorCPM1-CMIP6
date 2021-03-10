      subroutine asflux
c
c --- ------------------------------------------------------------------
c --- Compute air-sea fluxes. Same routine is used both for
c --- climatological and synoptic forcing
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_forc.h"
#ifdef CLIM
#  include "common_clim.h"
#else
#  include "common_syn.h"
#endif
#include "common_asf.h"
#include "common_geo.h"
#include "common_clndr.h"
c
c --- Parameters:
c ---   dtmax - maximum near surface temperature gradient [K]
c ---   dqmax - maximum near surface specific humidity gradient [kg/kg]
c
      real dtmax,dqmax
      parameter (dtmax=30.,dqmax=.05)
c
      integer i,j,l,niter
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: taufac,cc
      real rntda,tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice,
     .     tml,sml,albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa,
     .     tsrf_m,qsrf_m,dqsrf_m,dangle
c
      real intpol,swtfrz,qsati,qsatw,rhoair,dqsati,dqsatw
      external intpol,swtfrz,qsati,qsatw,rhoair,dqsati,dqsatw
c
      rntda=1./real(ntda)
      ntda=0
c
c$OMP PARALLEL DO PRIVATE(
c$OMP+  tau_d,prcp,dswrf,nlwrs,shtfl,lhtfl,slpr,tsrf_d,rice,tml,sml,
c$OMP+  albedo,fice,tice_f,tml_d,tsi_d,qsrf_d,le,ua,sa,ta,qa,niter,
c$OMP+  tsrf_m,qsrf_m,dqsrf_m)
c$OMP+  SHARED(xx,zu,zt,zq)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- ------------------------------------------------------------------
c --- --- Interpolate the climatological/synoptic fields
c --- ------------------------------------------------------------------
c
          tau_d =intpol(taud  (i,j,ll1),taud  (i,j,ll2),taud  (i,j,ll3),
     .                  taud  (i,j,ll4),taud  (i,j,ll5),xx)
          prcp  =intpol(precip(i,j,ll1),precip(i,j,ll2),precip(i,j,ll3),
     .                  precip(i,j,ll4),precip(i,j,ll5),xx)
          dswrf =intpol(dswrfl(i,j,ll1),dswrfl(i,j,ll2),dswrfl(i,j,ll3),
     .                  dswrfl(i,j,ll4),dswrfl(i,j,ll5),xx)
          nlwrs =intpol(nlwrfs(i,j,ll1),nlwrfs(i,j,ll2),nlwrfs(i,j,ll3),
     .                  nlwrfs(i,j,ll4),nlwrfs(i,j,ll5),xx)
          shtfl =intpol(shtflx(i,j,ll1),shtflx(i,j,ll2),shtflx(i,j,ll3),
     .                  shtflx(i,j,ll4),shtflx(i,j,ll5),xx)
          lhtfl =intpol(lhtflx(i,j,ll1),lhtflx(i,j,ll2),lhtflx(i,j,ll3),
     .                  lhtflx(i,j,ll4),lhtflx(i,j,ll5),xx)
          slpr  =intpol(slpres(i,j,ll1),slpres(i,j,ll2),slpres(i,j,ll3),
     .                  slpres(i,j,ll4),slpres(i,j,ll5),xx)
#ifdef CLIM
          tsrf_d=intpol(sstclm(i,j,ll1),sstclm(i,j,ll2),sstclm(i,j,ll3),
     .                  sstclm(i,j,ll4),sstclm(i,j,ll5),xx)
          rice  =intpol(ricclm(i,j,ll1),ricclm(i,j,ll2),ricclm(i,j,ll3),
     .                  ricclm(i,j,ll4),ricclm(i,j,ll5),xx)
#endif
#ifdef SYN
          tsrf_d=intpol(tmpsfc(i,j,ll1),tmpsfc(i,j,ll2),tmpsfc(i,j,ll3),
     .                  tmpsfc(i,j,ll4),tmpsfc(i,j,ll5),xx)
          rice  =intpol(ricec (i,j,ll1),ricec (i,j,ll2),ricec (i,j,ll3),
     .                  ricec (i,j,ll4),ricec (i,j,ll5),xx)
#endif
          rnf(i,j)
     .          =intpol(runoff(i,j,ll1),runoff(i,j,ll2),runoff(i,j,ll3),
     .                  runoff(i,j,ll4),runoff(i,j,ll5),xx)
c
          prcp =max(0.,prcp)
          dswrf=max(0.,dswrf)
          rice =max(0.,min(1.,rice))
c
c --- ------------------------------------------------------------------
c --- --- Get averaged quantities obtained using the previous surface
c --- --- fluxes
c --- ------------------------------------------------------------------
c
          tsi(i,j)=tsi_tda(i,j)*rntda
          tml=tml_tda(i,j)*rntda
          sml=sml_tda(i,j)*rntda
          albedo=alb_tda(i,j)*rntda
          fice=fice_tda(i,j)*rntda
c
          tsi_tda(i,j)=0.
          tml_tda(i,j)=0.
          sml_tda(i,j)=0.
          alb_tda(i,j)=0.
          fice_tda(i,j)=0.
c
c --- --- Freezing temperature of sea water
          tice_f=swtfrz(sml)+t0deg
c
c --- ------------------------------------------------------------------
c --- --- Compute the atmospheric state by using the prescribed momentum
c --- --- and heat fluxes and the prescribed sea surface state
c --- ------------------------------------------------------------------
c
          tml_d=max(tsrf_d,tice_f)
          tsi_d=max(200.,(tsrf_d-(1.-rice)*tml_d)/max(rice,1.e-6))
          qsrf_d=rice*qsati(tsi_d,slpr)+(1.-rice)*qsatw(tml_d,slpr)
          le=(2.501-0.00237*(tsrf_d-273.15))*1.e6
c
c --- --- Make sure wind stress is not too small compared to latent and
c --- --- sensible heat fluxes
          sa=max(abs(shtfl)/(rhoa(i,j)*cpair*ch_d(i,j)*dtmax),
     .           abs(lhtfl)/(rhoa(i,j)*le*ce_d(i,j)*dqmax))
          tau_d=max(tau_d,rhoa(i,j)*cd_d(i,j)*sa*sa)
c
c --- --- First guess on the atmospheric state by using the transfer
c --- --- coefficients and density from the previous time step
          ua=sqrt(.5*(-wg2_d(i,j)
     .                +sqrt(wg2_d(i,j)*wg2_d(i,j)
     .                     +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
          sa=sqrt(ua*ua+wg2_d(i,j))
          ta=tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
          qa=qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
          rhoa(i,j)=rhoair(ta,qa,slpr)
c
c --- --- Iteration loop for estimating transfer coefficients and
c --- --- atmospheric state
          do niter=1,tciter
c
c --- ----- update the transfer coefficients and gustiness
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_d,qsrf_d,rice,
     .                  cd_d(i,j),ch_d(i,j),ce_d(i,j),wg2_d(i,j))
c
c --- ----- update the atmospheric state
            ua=sqrt(.5*(-wg2_d(i,j)
     .                  +sqrt(wg2_d(i,j)*wg2_d(i,j)
     .                       +4.*(tau_d/(rhoa(i,j)*cd_d(i,j)))**2)))
            sa=sqrt(ua*ua+wg2_d(i,j))
            ta=tsrf_d-.0098*zt-shtfl/(rhoa(i,j)*cpair*ch_d(i,j)*sa)
            qa=qsrf_d-lhtfl/(rhoa(i,j)*le*ce_d(i,j)*sa)
            rhoa(i,j)=rhoair(ta,qa,slpr)
c
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Update transfer coefficients and gustiness with the computed
c --- --- atmospheric state the models ocean state
c --- ------------------------------------------------------------------
c
          tsrf_m=fice*tsi(i,j)+(1.-fice)*tml
          qsrf_m=fice*qsati(tsi(i,j),slpr)+(1.-fice)*qsatw(tml,slpr)
c
          do niter=1,tciter
            call bulktf(ua,zu,ta,zt,qa,zq,tsrf_m,qsrf_m,fice,
     .                  cd_m(i,j),ch_m(i,j),ce_m(i,j),wg2_m(i,j))
          enddo
c
c --- ------------------------------------------------------------------
c --- --- Compute correction of the wind stress on the surface and wind
c --- --- generated tke [m/s]
c --- ------------------------------------------------------------------
c
          sa=sqrt(ua*ua+wg2_m(i,j))
          taufac(i,j)=rhoa(i,j)*cd_m(i,j)*sa*ua/tau_d
c
c --- --- Wind generated TKE
          ustarw(i,j)=sqrt(cd_m(i,j)*sa*ua*rhoa(i,j)/rhowat)
c
c --- ------------------------------------------------------------------
c --- --- Compute heat fluxes (positive downward) [W/m^2]
c --- ------------------------------------------------------------------
c
          swa(i,j)=dswrf*(1.-albedo)
          le=(2.501-.00237*(tsrf_m-273.15))*1.e6
          nsf(i,j)=rhoa(i,j)*cpair*ch_m(i,j)*sa*(ta+0.0098*zt-tsrf_m)
     .            +rhoa(i,j)*ce_m(i,j)*le*sa*(qa-qsrf_m)
     .            -nlwrs-4.*emiss*stefanb*ta**3*(tsrf_m-tsrf_d)
c
c --- ------------------------------------------------------------------
c --- --- Compute evaporation (positive downward) [kg/m^2/s]
c --- ------------------------------------------------------------------
c
          eva(i,j)=rhoa(i,j)*ce_m(i,j)*sa*(qa-qsrf_m)
c
c --- ------------------------------------------------------------------
c --- --- Compute derivative of non-solar flux by surface temperature
c --- --- [W/m^2/K]
c --- ------------------------------------------------------------------
c
          dqsrf_m=fice*dqsati(tsi(i,j),slpr)+(1.-fice)*dqsatw(tml,slpr)
          dfl(i,j)=-rhoa(i,j)*cpair*ch_m(i,j)*sa
     .             -rhoa(i,j)*ce_m(i,j)*le*sa*dqsrf_m
     .             -4.*emiss*stefanb*ta**3
c
c --- ------------------------------------------------------------------
c --- --- Split solid and liquid precipitation (positive downward)
c --- --- [kg/m^2/s]
c --- ------------------------------------------------------------------
c
          if (ta.lt.t0deg) then
            lip(i,j)=0.
            sop(i,j)=prcp
          else
            lip(i,j)=prcp
            sop(i,j)=0.
          endif
c
c --- ------------------------------------------------------------------
c --- --- Sea level pressure [Pa] and wind speed at measurement height
c --- --- [m/s]
c --- ------------------------------------------------------------------
c
          slp(i,j)=slpr
          abswnd(i,j)=sa
c
c --- ------------------------------------------------------------------
c --- --- Frozen runoff is not computed.
c --- ------------------------------------------------------------------
c
          rfi(i,j)=0.
c
c --- ------------------------------------------------------------------
c --- --- If requested, apply correction to precipitation and runoff to
c --- --- balance the fresh water budget
c --- ------------------------------------------------------------------
c
          if (sprfac) then
            lip(i,j)=lip(i,j)*prfac
            sop(i,j)=sop(i,j)*prfac
            rnf(i,j)=rnf(i,j)*prfac
            rfi(i,j)=rfi(i,j)*prfac
          endif
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute the surface wind stress [N/m^2]
c --- ------------------------------------------------------------------
c
      call xctilr(taufac,  1,1, 1,1, halo_ps)
c
c$OMP PARALLEL DO SHARED(xx)
      do j=1,jj
        do l=1,isu(j)
        do i=max(1,ifu(j,l)),min(ii,ilu(j,l))
          ztx(i,j)=.5*(taufac(i,j)+taufac(i-1,j))
     .             *intpol(tauxd(i,j,ll1),tauxd(i,j,ll2),tauxd(i,j,ll3),
     .                     tauxd(i,j,ll4),tauxd(i,j,ll5),xx)
        enddo
        enddo
        do l=1,isv(j)
        do i=max(1,ifv(j,l)),min(ii,ilv(j,l))
          mty(i,j)=.5*(taufac(i,j)+taufac(i,j-1))
     .             *intpol(tauyd(i,j,ll1),tauyd(i,j,ll2),tauyd(i,j,ll3),
     .                     tauyd(i,j,ll4),tauyd(i,j,ll5),xx)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Compute open water albedo
c --- ------------------------------------------------------------------
c 
c --- Time interpolation of cloud cover
c$OMP PARALLEL DO SHARED(x)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          cc(i,j)=
     .      intpol(clouds(i,j,ll1),clouds(i,j,ll2),clouds(i,j,ll3),
     .             clouds(i,j,ll4),clouds(i,j,ll5),xx)
          cc(i,j)=max(0.,min(1.,cc(i,j)))
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Obtain albedo for open water
      dangle=8.*atan(1.)*real(nday_of_year-1)/real(nday_in_year)
      call albw_eval(dangle,plat,cc,albw_d,albw)
c
      if (csdiag) then
        if (mnproc.eq.1) then
          write (lp,*) 'asflux:'
        endif
        call chksummsk(swa,ip,1,'swa')
        call chksummsk(nsf,ip,1,'nsf')
        call chksummsk(dfl,ip,1,'dfl')
        call chksummsk(lip,ip,1,'lip')
        call chksummsk(sop,ip,1,'sop')
        call chksummsk(eva,ip,1,'eva')
        call chksummsk(ztx,iu,1,'ztx')
        call chksummsk(mty,iv,1,'mty')
        call chksummsk(rnf,ip,1,'rnf')
        call chksummsk(ustarw,ip,1,'ustarw')
        call chksummsk(tsi,ip,1,'tsi')
        call chksummsk(slp,ip,1,'slp')
        call chksummsk(abswnd,ip,1,'abswnd')
        call chksummsk(albw,ip,1,'albw')
      endif
c
      return
      end
