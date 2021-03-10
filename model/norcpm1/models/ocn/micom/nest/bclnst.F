      subroutine bclnst(m,n,mm,nn,k1m,k1n)
c
      implicit none
c
      integer m,n,mm,nn,k1m,k1n
c
#include "dimensions.h"
#include "common_blocks.h"
#include "common_ice.h"
#include "param_nest.h"
#include "common_nest.h"
c
      integer i,j,k,l,kn,ibz
      real pbfac(idm),wt1,wt2,ws1,ws2,dpw1,dpw21,dpw22,q
c
      real sig
      external sig
c
      wt2=real(nstep-nst31)/(nst32-nst31)
      wt1=1.-wt2
c
c --- ------------------------------------------------------------------
c --- "western" and "eastern" boundary
c --- ------------------------------------------------------------------
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(ws2,ws1,i,dpw1,dpw21,dpw22,q)
        do j=1,jj
          do ibz=1,nbz
            ws2=wnst(ibz)
            ws1=1.-ws2
c
            i=ibz
            if (ip(i,j).eq.1) then
              dpw1=max(onemm,dp(i,j,kn)*ws1)
              dpw21=min(dp(i,j,kn),dp_ow1(ibz,j,k))*ws2*wt1
              dpw22=min(dp(i,j,kn),dp_ow2(ibz,j,k))*ws2*wt2
              q=1./(dpw1+dpw21+dpw22)
              temp(i,j,kn)=(dpw1*temp(i,j,kn)+dpw21*temp_ow1(ibz,j,k)
     .                     +dpw22*temp_ow2(ibz,j,k))*q
              saln(i,j,kn)=(dpw1*saln(i,j,kn)+dpw21*saln_ow1(ibz,j,k)
     .                     +dpw22*saln_ow2(ibz,j,k))*q
c             temp(i,j,kn)=ws1*temp(i,j,kn)
c    .          +ws2*(wt1*temp_ow1(ibz,j,k)+wt2*temp_ow2(ibz,j,k))
c             saln(i,j,kn)=ws1*saln(i,j,kn)
c    .          +ws2*(wt1*saln_ow1(ibz,j,k)+wt2*saln_ow2(ibz,j,k))
              dp(i,j,kn)=ws1*dp(i,j,kn)
     .          +ws2*(wt1*dp_ow1(ibz,j,k)+wt2*dp_ow2(ibz,j,k))
              sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
            i=ii-ibz
            if (ip(i,j).eq.1) then
              dpw1=max(onemm,dp(i,j,kn)*ws1)
              dpw21=min(dp(i,j,kn),dp_oe1(ibz,j,k))*ws2*wt1
              dpw22=min(dp(i,j,kn),dp_oe2(ibz,j,k))*ws2*wt2
              q=1./(dpw1+dpw21+dpw22)
              temp(i,j,kn)=(dpw1*temp(i,j,kn)+dpw21*temp_oe1(ibz,j,k)
     .                     +dpw22*temp_oe2(ibz,j,k))*q
              saln(i,j,kn)=(dpw1*saln(i,j,kn)+dpw21*saln_oe1(ibz,j,k)
     .                     +dpw22*saln_oe2(ibz,j,k))*q
c             temp(i,j,kn)=ws1*temp(i,j,kn)
c    .          +ws2*(wt1*temp_oe1(ibz,j,k)+wt2*temp_oe2(ibz,j,k))
c             saln(i,j,kn)=ws1*saln(i,j,kn)
c    .          +ws2*(wt1*saln_oe1(ibz,j,k)+wt2*saln_oe2(ibz,j,k))
              dp(i,j,kn)=ws1*dp(i,j,kn)
     .          +ws2*(wt1*dp_oe1(ibz,j,k)+wt2*dp_oe2(ibz,j,k))
              sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(i,pbfac,kn)
      do j=1,jj
c
        do ibz=1,nbz
c
          i=ibz
          if (ip(i,j).eq.1) then
            pbfac(i)=pb(i,j,n)/p(i,j,kk+1)
          endif
c
          i=ii-ibz
          if (ip(i,j).eq.1) then
            pbfac(i)=pb(i,j,n)/p(i,j,kk+1)
          endif
c
        enddo
c
        do k=1,kk
          kn=k+nn
          do ibz=1,nbz
c
            i=ibz
            if (ip(i,j).eq.1) then
              dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
            i=ii-ibz
            if (ip(i,j).eq.1) then
              dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
          enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do ibz=1,nbz
c
          i=ibz+1
          if (iu(i,j).eq.1) then
            util1(i,j)=0.
          endif
          i=ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=0.
          endif
c
          i=ii-ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=0.
          endif
          i=ii-ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=0.
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(ws2,ws1,i)
        do j=1,jj
          do ibz=1,nbz
            ws2=wnst(ibz)
            ws1=1.-ws2
c
            i=ibz+1
            if (iu(i,j).eq.1) then
              u(i,j,kn)=ws1*u(i,j,kn)
     .          +ws2*(wt1*utot_ow1(ibz,j,k)+wt2*utot_ow2(ibz,j,k))
              dpu(i,j,kn)=
     .          max(0.,min(pbu(i,j,n),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .                -min(pbu(i,j,n),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
              util1(i,j)=util1(i,j)+u(i,j,kn)*dpu(i,j,kn)
            endif
            i=ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=ws1*v(i,j,kn)
     .          +ws2*(wt1*vtot_ow1(ibz,j,k)+wt2*vtot_ow2(ibz,j,k))
              dpv(i,j,kn)=
     .          max(0.,min(pbv(i,j,n),.5*(p(i,j  ,k+1)+p(i,j-1,k+1))) 
     .                -min(pbv(i,j,n),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
              util2(i,j)=util2(i,j)+v(i,j,kn)*dpv(i,j,kn)
            endif
c
            i=ii-ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=ws1*u(i,j,kn)
     .          +ws2*(wt1*utot_oe1(ibz,j,k)+wt2*utot_oe2(ibz,j,k))
              dpu(i,j,kn)=
     .          max(0.,min(pbu(i,j,n),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .                -min(pbu(i,j,n),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
              util1(i,j)=util1(i,j)+u(i,j,kn)*dpu(i,j,kn)
            endif
            i=ii-ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=ws1*v(i,j,kn)
     .          +ws2*(wt1*vtot_oe1(ibz,j,k)+wt2*vtot_oe2(ibz,j,k))
              dpv(i,j,kn)=
     .          max(0.,min(pbv(i,j,n),.5*(p(i,j  ,k+1)+p(i,j-1,k+1))) 
     .                -min(pbv(i,j,n),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
              util2(i,j)=util2(i,j)+v(i,j,kn)*dpv(i,j,kn)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(i)
      do j=1,jj
        do ibz=1,nbz
c
          i=ibz+1
          if (iu(i,j).eq.1) then
            util1(i,j)=util1(i,j)/pbu(i,j,n)
          endif
          i=ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=util2(i,j)/pbv(i,j,n)
          endif
c
          i=ii-ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=util1(i,j)/pbu(i,j,n)
          endif
          i=ii-ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=util2(i,j)/pbv(i,j,n)
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(i)
        do j=1,jj
          do ibz=1,nbz
c
            i=ibz+1
            if (iu(i,j).eq.1) then
              u(i,j,kn)=u(i,j,kn)-util1(i,j)
            endif
            i=ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=v(i,j,kn)-util2(i,j)
            endif
c
            i=ii-ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=u(i,j,kn)-util1(i,j)
            endif
            i=ii-ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=v(i,j,kn)-util2(i,j)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(ws2,ws1,i)
      do j=1,jj
        do ibz=1,nbz
          ws2=wnst(ibz)
          ws1=1.-ws2
c
          i=ibz
          if (ip(i,j).eq.1) then
            ficem(i,j)=ws1*ficem(i,j)
     .        +ws2*(wt1*ficem_ow1(ibz,j)+wt2*ficem_ow2(ibz,j))
            hicem(i,j)=ws1*hicem(i,j)
     .        +ws2*(wt1*hicem_ow1(ibz,j)+wt2*hicem_ow2(ibz,j))
            hsnwm(i,j)=ws1*hsnwm(i,j)
     .        +ws2*(wt1*hsnwm_ow1(ibz,j)+wt2*hsnwm_ow2(ibz,j))
            iagem(i,j)=ws1*iagem(i,j)
     .        +ws2*(wt1*iagem_ow1(ibz,j)+wt2*iagem_ow2(ibz,j))
          endif
c
          i=ii-ibz
          if (ip(i,j).eq.1) then
            ficem(i,j)=ws1*ficem(i,j)
     .        +ws2*(wt1*ficem_oe1(ibz,j)+wt2*ficem_oe2(ibz,j))
            hicem(i,j)=ws1*hicem(i,j)
     .        +ws2*(wt1*hicem_oe1(ibz,j)+wt2*hicem_oe2(ibz,j))
            hsnwm(i,j)=ws1*hsnwm(i,j)
     .        +ws2*(wt1*hsnwm_oe1(ibz,j)+wt2*hsnwm_oe2(ibz,j))
            iagem(i,j)=ws1*iagem(i,j)
     .        +ws2*(wt1*iagem_oe1(ibz,j)+wt2*iagem_oe2(ibz,j))
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- "southern" and "northern" boundary
c --- ------------------------------------------------------------------
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(ws2,ws1,i,dpw1,dpw21,dpw22,q)
        do ibz=1,nbz
          ws2=wnst(ibz)
          ws1=1.-ws2
          do i=1,ii
c
            j=ibz
            if (ip(i,j).eq.1) then
              dpw1=max(onemm,dp(i,j,kn)*ws1)
              dpw21=min(dp(i,j,kn),dp_os1(i,ibz,k))*ws2*wt1
              dpw22=min(dp(i,j,kn),dp_os2(i,ibz,k))*ws2*wt2
              q=1./(dpw1+dpw21+dpw22)
              temp(i,j,kn)=(dpw1*temp(i,j,kn)+dpw21*temp_os1(i,ibz,k)
     .                     +dpw22*temp_os2(i,ibz,k))*q
              saln(i,j,kn)=(dpw1*saln(i,j,kn)+dpw21*saln_os1(i,ibz,k)
     .                     +dpw22*saln_os2(i,ibz,k))*q
c             temp(i,j,kn)=ws1*temp(i,j,kn)
c    .          +ws2*(wt1*temp_os1(i,ibz,k)+wt2*temp_os2(i,ibz,k))
c             saln(i,j,kn)=ws1*saln(i,j,kn)
c    .          +ws2*(wt1*saln_os1(i,ibz,k)+wt2*saln_os2(i,ibz,k))
              dp(i,j,kn)=ws1*dp(i,j,kn)
     .          +ws2*(wt1*dp_os1(i,ibz,k)+wt2*dp_os2(i,ibz,k))
              sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
            j=jj-ibz
            if (ip(i,j).eq.1) then
              dpw1=max(onemm,dp(i,j,kn)*ws1)
              dpw21=min(dp(i,j,kn),dp_on1(i,ibz,k))*ws2*wt1
              dpw22=min(dp(i,j,kn),dp_on2(i,ibz,k))*ws2*wt2
              q=1./(dpw1+dpw21+dpw22)
              temp(i,j,kn)=(dpw1*temp(i,j,kn)+dpw21*temp_on1(i,ibz,k)
     .                     +dpw22*temp_on2(i,ibz,k))*q
              saln(i,j,kn)=(dpw1*saln(i,j,kn)+dpw21*saln_on1(i,ibz,k)
     .                     +dpw22*saln_on2(i,ibz,k))*q
c             temp(i,j,kn)=ws1*temp(i,j,kn)
c    .          +ws2*(wt1*temp_on1(i,ibz,k)+wt2*temp_on2(i,ibz,k))
c             saln(i,j,kn)=ws1*saln(i,j,kn)
c    .          +ws2*(wt1*saln_on1(i,ibz,k)+wt2*saln_on2(i,ibz,k))
              dp(i,j,kn)=ws1*dp(i,j,kn)
     .          +ws2*(wt1*dp_on1(i,ibz,k)+wt2*dp_on2(i,ibz,k))
              sigma(i,j,kn)=sig(temp(i,j,kn),saln(i,j,kn))
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(j,pbfac,kn)
      do ibz=1,nbz
c
        do i=1,ii
c
          j=ibz
          if (ip(i,j).eq.1) then
            pbfac(i)=pb(i,j,n)/p(i,j,kk+1)
          endif
c
          j=jj-ibz
          if (ip(i,j).eq.1) then
            pbfac(i)=pb(i,j,n)/p(i,j,kk+1)
          endif
c
        enddo
c
        do k=1,kk
          kn=k+nn
          do i=1,ii
c
            j=ibz
            if (ip(i,j).eq.1) then
              dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
            j=jj-ibz
            if (ip(i,j).eq.1) then
              dp(i,j,kn)=dp(i,j,kn)*pbfac(i)
              p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
            endif
c
          enddo
        enddo
c
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO PRIVATE(j)
      do ibz=1,nbz
        do i=1,ii
c
          j=ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=0.
          endif
          j=ibz+1
          if (iv(i,j).eq.1) then
            util2(i,j)=0.
          endif
c
          j=jj-ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=0.
          endif
          j=jj-ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=0.
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(ws2,ws1,j)
        do ibz=1,nbz
          ws2=wnst(ibz)
          ws1=1.-ws2
          do i=1,ii
c
            j=ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=ws1*u(i,j,kn)
     .          +ws2*(wt1*utot_os1(i,ibz,k)+wt2*utot_os2(i,ibz,k))
              dpu(i,j,kn)=
     .          max(0.,min(pbu(i,j,n),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .                -min(pbu(i,j,n),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
              util1(i,j)=util1(i,j)+u(i,j,kn)*dpu(i,j,kn)
            endif
            j=ibz+1
            if (iv(i,j).eq.1) then
              v(i,j,kn)=ws1*v(i,j,kn)
     .          +ws2*(wt1*vtot_os1(i,ibz,k)+wt2*vtot_os2(i,ibz,k))
              dpv(i,j,kn)=
     .          max(0.,min(pbv(i,j,n),.5*(p(i,j  ,k+1)+p(i,j-1,k+1))) 
     .                -min(pbv(i,j,n),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
              util2(i,j)=util2(i,j)+v(i,j,kn)*dpv(i,j,kn)
            endif
c
            j=jj-ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=ws1*u(i,j,kn)
     .          +ws2*(wt1*utot_on1(i,ibz,k)+wt2*utot_on2(i,ibz,k))
              dpu(i,j,kn)=
     .          max(0.,min(pbu(i,j,n),.5*(p(i  ,j,k+1)+p(i-1,j,k+1))) 
     .                -min(pbu(i,j,n),.5*(p(i  ,j,k  )+p(i-1,j,k  ))))
              pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
              util1(i,j)=util1(i,j)+u(i,j,kn)*dpu(i,j,kn)
            endif
            j=jj-ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=ws1*v(i,j,kn)
     .          +ws2*(wt1*vtot_on1(i,ibz,k)+wt2*vtot_on2(i,ibz,k))
              dpv(i,j,kn)=
     .          max(0.,min(pbv(i,j,n),.5*(p(i,j  ,k+1)+p(i,j-1,k+1))) 
     .                -min(pbv(i,j,n),.5*(p(i,j  ,k  )+p(i,j-1,k  ))))
              pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
              util2(i,j)=util2(i,j)+v(i,j,kn)*dpv(i,j,kn)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(j)
      do ibz=1,nbz
        do i=1,ii
c
          j=ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=util1(i,j)/pbu(i,j,n)
          endif
          j=ibz+1
          if (iv(i,j).eq.1) then
            util2(i,j)=util2(i,j)/pbv(i,j,n)
          endif
c
          j=jj-ibz
          if (iu(i,j).eq.1) then
            util1(i,j)=util1(i,j)/pbu(i,j,n)
          endif
          j=jj-ibz
          if (iv(i,j).eq.1) then
            util2(i,j)=util2(i,j)/pbv(i,j,n)
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
      do k=1,kk
        kn=k+nn
c$OMP PARALLEL DO PRIVATE(j)
        do ibz=1,nbz
          do i=1,ii
c
            j=ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=u(i,j,kn)-util1(i,j)
            endif
            j=ibz+1
            if (iv(i,j).eq.1) then
              v(i,j,kn)=v(i,j,kn)-util2(i,j)
            endif
c
            j=jj-ibz
            if (iu(i,j).eq.1) then
              u(i,j,kn)=u(i,j,kn)-util1(i,j)
            endif
            j=jj-ibz
            if (iv(i,j).eq.1) then
              v(i,j,kn)=v(i,j,kn)-util2(i,j)
            endif
c
          enddo
        enddo
c$OMP END PARALLEL DO
      enddo
c
c$OMP PARALLEL DO PRIVATE(ws2,ws1,j)
      do ibz=1,nbz
        ws2=wnst(ibz)
        ws1=1.-ws2
        do i=1,ii
c
          j=ibz
          if (ip(i,j).eq.1) then
            ficem(i,j)=ws1*ficem(i,j)
     .        +ws2*(wt1*ficem_os1(i,ibz)+wt2*ficem_os2(i,ibz))
            hicem(i,j)=ws1*hicem(i,j)
     .        +ws2*(wt1*hicem_os1(i,ibz)+wt2*hicem_os2(i,ibz))
            hsnwm(i,j)=ws1*hsnwm(i,j)
     .        +ws2*(wt1*hsnwm_os1(i,ibz)+wt2*hsnwm_os2(i,ibz))
            iagem(i,j)=ws1*iagem(i,j)
     .        +ws2*(wt1*iagem_os1(i,ibz)+wt2*iagem_os2(i,ibz))
          endif
c
          j=jj-ibz
          if (ip(i,j).eq.1) then
            ficem(i,j)=ws1*ficem(i,j)
     .        +ws2*(wt1*ficem_on1(i,ibz)+wt2*ficem_on2(i,ibz))
            hicem(i,j)=ws1*hicem(i,j)
     .        +ws2*(wt1*hicem_on1(i,ibz)+wt2*hicem_on2(i,ibz))
            hsnwm(i,j)=ws1*hsnwm(i,j)
     .        +ws2*(wt1*hsnwm_on1(i,ibz)+wt2*hsnwm_on2(i,ibz))
            iagem(i,j)=ws1*iagem(i,j)
     .        +ws2*(wt1*iagem_on1(i,ibz)+wt2*iagem_on2(i,ibz))
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
