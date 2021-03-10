      subroutine forfun_syn
c
c --- ------------------------------------------------------------------------
c --- Prepare utilization of daily atmospheric forcing
c --- ------------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "param_atm.h"
#include "common_atm.h"
#include "common_atmsyn.h"
#include "common_syn.h"
#include "common_clndr.h"
c
      integer i,j,k,l
      integer*2 rn2,ri2,rj2
      real dx2,dy2
      real*4 rw4
c
      if (mnproc.eq.1) then
c
c --- - Read file containing ocean discarge points for land areas
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
      endif
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
c
c --- Initial indexes of the 5 time levels
      ll1=1
      ll2=2
      ll3=3
      ll4=4
      ll5=5
c
c --- Go back 4 days so that all necessary fields are read
c
      call pstdat(nyear,nmonth,nday,4)
      call inicln
      do i=1,4
        call rdatm_syn
        call updcln(0)
      enddo
c
      return
      end
