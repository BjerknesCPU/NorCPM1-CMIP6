      subroutine rdatm_syn
c
c --- ------------------------------------------------------------------
c --- Read and interpolate daily atmospheric forcing fields. The
c --- interpolated fields of heat fluxes, cloud cover, precipitation,
c --- and surface temperature are smoothed over ice covered regions.
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_atm.h"
#include "common_atm.h"
#include "common_atmsyn.h"
#include "common_syn.h"
#include "common_geo.h"
c
      real, dimension(atm_idm,atm_jdm) :: atm_field,atm_skt,atm_tau
      real, dimension(itdm,jtdm) :: tmp2d
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: smtmsk
      integer i,j,l,nday_of_year1,nyear1,mbdy,iii,jjj
      character cyear*4,filename*120,fieldname*5
c
c --- The actual fields to be read are 2 days ahead
      nday_of_year1=nday_of_year+2
      if (nday_of_year1.gt.nday_in_year) then
        nday_of_year1=nday_of_year1-nday_in_year
        nyear1=nyear+1
      else
        nyear1=nyear
      endif
c
c --- Do not go beyond the start of the atm. time series
#if   defined(NCEP)
      if (nyear1.lt.1948) then
        nday_of_year1=1
        nyear1=1948
      endif
#elif defined(ERA)
      if (nyear1.lt.1958) then
        nday_of_year1=1
        nyear1=1958
      endif
#endif
c
      write(cyear,'(i4)') nyear1
c
c --- Rearrange the time level indexes
      i=ll1
      ll1=ll2
      ll2=ll3
      ll3=ll4
      ll4=ll5
      ll5=i
c
      if (mnproc.eq.1) then
        write (lp,'(a,i3,a,i4,a)')
     .    'Reading and interpolating atm. forcing fields for day ',
     .    nday_of_year1,' of year ',nyear1,'...'
        write (lp,'(a,5i2)') 'Time level indexes:',ll1,ll2,ll3,ll4,ll5
      endif
c
c --- ------------------------------------------------------------------
c --- Read/interpolate ice concentration [0-1]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'ICECsfc/icec.sfc.gauss.'//cyear//'.nc'
        fieldname='icec'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,ricec(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 ricec(1-nbdy,1-nbdy,ll5))
c
c --- create smoothing mask - smooth where ice conc. is above 0.5
c$OMP PARALLEL DO
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          if (ricec(i,j,ll5).gt..5.and.ip(i,j).eq.1) then
            smtmsk(i,j)=1
          else
            smtmsk(i,j)=0
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate short-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'DSWRFsfc/dswrf.sfc.gauss.'//cyear//'.nc'
        fieldname='dswrf'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,dswrfl(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 dswrfl(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(dswrfl(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              dswrfl(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- Adjust short-wave radiation field
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          dswrfl(i,j,ll5)=dswrfl(i,j,ll5)*atm_cswa
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate net long-wave radiation flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'NLWRSsfc/nlwrs.sfc.gauss.'//cyear//'.nc'
        fieldname='nlwrs'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,nlwrfs(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 nlwrfs(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(nlwrfs(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              nlwrfs(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- ------------------------------------------------------------------
c --- Read/interpolate total cloud cover [0-100%]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'TCDCeatm/tcdc.eatm.gauss.'//cyear//'.nc'
        fieldname='tcdc'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,clouds(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 clouds(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(clouds(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              clouds(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- Convert range of cloudiness from 0-100 to 0-1
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          clouds(i,j,ll5)=clouds(i,j,ll5)*1.e-2
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- Read/interpolate precipitation rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'PRATEsfc/prate.sfc.gauss.'//cyear//'.nc'
        fieldname='prate'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,precip(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 precip(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(precip(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              precip(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- ------------------------------------------------------------------
c --- Read/interpolate latent heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'LHTFLsfc/lhtfl.sfc.gauss.'//cyear//'.nc'
        fieldname='lhtfl'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,lhtflx(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 lhtflx(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(lhtflx(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              lhtflx(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- ------------------------------------------------------------------
c --- Read/interpolate sensible heat net flux [W/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'SHTFLsfc/shtfl.sfc.gauss.'//cyear//'.nc'
        fieldname='shtfl'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,shtflx(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 shtflx(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(shtflx(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              shtflx(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- ------------------------------------------------------------------
c --- Read/interpolate sea surface temperature [K]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'SKTsfc/skt.sfc.gauss.'//cyear//'.nc'
        fieldname='skt'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_skt)
c
c --- - Compute sea level temperature
c$OMP PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_skt(i,j)=atm_skt(i,j)+.0065*atm_topo(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_skt,tmp2d)
      endif
      call xcaput(tmp2d,tmpsfc(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 tmpsfc(1-nbdy,1-nbdy,ll5))
      mbdy=1
      do i=1,atm_ice_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(tmpsfc(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_ice_swgt,
     .              tmpsfc(1-nbdy,1-nbdy,ll5),smtmsk,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- ------------------------------------------------------------------
c --- Read/interpolate surface pressure [Pa]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'PRESsfc/pres.sfc.gauss.'//cyear//'.nc'
        fieldname='pres'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
c
c --- Compute sea level pressure
c$OMP PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_field(i,j)=atm_field(i,j)
     .                      *exp(9.81*atm_topo(i,j)
     .                      /(287.*(atm_skt(i,j)-.00325*atm_topo(i,j))))
          enddo
        enddo
c$OMP END PARALLEL DO
c
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
      endif
      call xcaput(tmp2d,slpres(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 slpres(1-nbdy,1-nbdy,ll5))
c
c --- ------------------------------------------------------------------
c --- Read/interpolate momentum flux [N/m^2]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'UFLXsfc/uflx.sfc.gauss.'//cyear//'.nc'
        fieldname='uflx'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
c$OMP PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_tau(i,j)=atm_field(i,j)*atm_field(i,j)
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcaput(tmp2d,tauxd(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_pv,
     .                 tauxd(1-nbdy,1-nbdy,ll5))
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'VFLXsfc/vflx.sfc.gauss.'//cyear//'.nc'
        fieldname='vflx'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_field,tmp2d)
c$OMP PARALLEL DO
        do j=1,atm_jdm
          do i=1,atm_idm
            atm_tau(i,j)=sqrt(atm_tau(i,j)
     .                        +atm_field(i,j)*atm_field(i,j))
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
      call xcaput(tmp2d,tauyd(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_pv,
     .                 tauyd(1-nbdy,1-nbdy,ll5))
      if (mnproc.eq.1) then
        call intpol2d(atm_idm,atm_jdm,itdm,jtdm,atm_nwgt,
     .                atm_wgt,atm_iwgt,atm_jwgt,atm_mval,atm_fval,
     .                atm_msk,itp,atm_tau,tmp2d)
      endif
      call xcaput(tmp2d,taud(1-nbdy,1-nbdy,ll5),1)
      call fill_global(atm_mval,atm_fval,halo_ps,
     .                 taud(1-nbdy,1-nbdy,ll5))
c
c --- Change sign of momentum flux
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          tauxd(i,j,ll5)=-tauxd(i,j,ll5)
          tauyd(i,j,ll5)=-tauyd(i,j,ll5)
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Rotate the vector components
      call xctilr(tauxd(1-nbdy,1-nbdy,ll5), 1,1, 1,1, halo_pv)
      call xctilr(tauyd(1-nbdy,1-nbdy,ll5), 1,1, 1,1, halo_pv)
      call uvrotr2g(atm_mval,atm_fval,plat,plon,iu,iv,
     .              tauxd(1-nbdy,1-nbdy,ll5),tauyd(1-nbdy,1-nbdy,ll5))
c
c --- ------------------------------------------------------------------
c --- Read/distribute runoff rate [kg/m^2/s]
c --- ------------------------------------------------------------------
c
      if (mnproc.eq.1) then
        filename=atm_path(1:atm_path_len)
     .           //'RUNOFsfc/runof.sfc.gauss.'//cyear//'.nc'
        fieldname='runof'
        call read_atm_field(filename,fieldname,atm_idm,atm_jdm,
     .                      nday_of_year1,atm_mval,atm_field)
c
c --- - Place runoff at ocean discharge points
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
                iii=rnf_ocdpi(l,i,j)
                jjj=rnf_ocdpj(l,i,j)
                tmp2d(iii,jjj)=tmp2d(iii,jjj)
     .            +atm_field(i,j)*rnf_weight(l,i,j)
              endif
            enddo
          enddo
        enddo
      endif
      call xcaput(tmp2d,runoff(1-nbdy,1-nbdy,ll5),1)
c
c --- Multiply runoff by 2 for fields prior to May 1957 
      if ( nyear1.lt.1957.or.
     .    (nyear1.eq.1957.and.nday_of_year1.le.120)) then
c$OMP PARALLEL DO
        do j=1,jj
          do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            runoff(i,j,ll5)=runoff(i,j,ll5)*2.
          enddo
          enddo
        enddo
c$OMP END PARALLEL DO
      endif
c
c --- Adjust runoff field
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,ll5)=runoff(i,j,ll5)*atm_crnf
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- Smooth the runoff field
      mbdy=1
      do i=1,atm_rnf_nsmt
        if (mbdy.eq.1) then
          mbdy=nbdy
          call xctilr(runoff(1-nbdy,1-nbdy,ll5), 1,1, nbdy,nbdy,halo_ps)
        else
          mbdy=mbdy-1
        endif
        call smooth(mbdy,atm_rnf_swgt,
     .              runoff(1-nbdy,1-nbdy,ll5),ip,
     .              scuy,scvx,scuxi,scvyi,scp2i)
      enddo
c
c --- Convert unit of runoff from kg/m^2/day to kg/m^2/s
c$OMP PARALLEL DO
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
          runoff(i,j,ll5)=runoff(i,j,ll5)/86400.
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      if (csdiag) then
        call chksummsk(ricec(1-nbdy,1-nbdy,ll5),ip,1,'ricec')
        call chksummsk(dswrfl(1-nbdy,1-nbdy,ll5),ip,1,'dswrfl')
        call chksummsk(nlwrfs(1-nbdy,1-nbdy,ll5),ip,1,'nlwrfs')
        call chksummsk(clouds(1-nbdy,1-nbdy,ll5),ip,1,'clouds')
        call chksummsk(precip(1-nbdy,1-nbdy,ll5),ip,1,'precip')
        call chksummsk(lhtflx(1-nbdy,1-nbdy,ll5),ip,1,'lhtflx')
        call chksummsk(shtflx(1-nbdy,1-nbdy,ll5),ip,1,'shtflx')
        call chksummsk(tmpsfc(1-nbdy,1-nbdy,ll5),ip,1,'tmpsfc')
        call chksummsk(slpres(1-nbdy,1-nbdy,ll5),ip,1,'slpres')
        call chksummsk(tauxd(1-nbdy,1-nbdy,ll5),ip,1,'tauxd')
        call chksummsk(tauyd(1-nbdy,1-nbdy,ll5),ip,1,'tauyd')
        call chksummsk(taud(1-nbdy,1-nbdy,ll5),ip,1,'taud')
        call chksummsk(runoff(1-nbdy,1-nbdy,ll5),ip,1,'runoff')
      endif
c
      return
      end
