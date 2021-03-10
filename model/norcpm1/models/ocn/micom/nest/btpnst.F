      subroutine btpnst(m,n,nl,lll,wn)
c
      implicit none
c
      integer m,n,nl,lll
      real wn
c
#include "dimensions.h"
#include "common_blocks.h"
#include "param_nest.h"
#include "common_nest.h"
c
      integer i,j,l
      real wt1,wt2,pbo,ubo,vbo,xi,pbi
c
      wt2=(nstep-1+real(lll)/lstep-nst31)/(nst32-nst31)
      wt1=1.-wt2
c     if (lll.le.lstep) write (lp,*) 'btpnst',wt2
c
c --- "Western" boundary
      i=2
      do l=1,jsu(i)
      do j=jfu(i,l),jlu(i,l)
        pb(i-1,j,nl)=wt1*pb_ow1(1,j)+wt2*pb_ow2(1,j)
        pbo=wt1*pb_ow1(2,j)+wt2*pb_ow2(2,j)
        ubo=wt1*ub_ow1(j)+wt2*ub_ow2(j)
        xi=(1.-wn)*xixp(i,j,m)+wn*xixp(i,j,n)
        pbi=pb(i,j,nl)
        ub(i,j,nl)=ubo+sqrt(xi/pbi)*(pbo-pbi)
      enddo
      enddo
      i=1
      do l=1,jsv(i)
      do j=jfv(i,l),jlv(i,l)
        vb(i,j,nl)=wt1*vb_ow1(j)+wt2*vb_ow2(j)
      enddo
      enddo
c
c --- "Eastern" boundary
      i=ii-1
      do l=1,jsu(i)
      do j=jfu(i,l),jlu(i,l)
        pb(i,j,nl)=wt1*pb_oe1(1,j)+wt2*pb_oe2(1,j)
        pbo=wt1*pb_oe1(2,j)+wt2*pb_oe2(2,j)
        ubo=wt1*ub_oe1(j)+wt2*ub_oe2(j)
        xi=(1.-wn)*xixm(i,j,m)+wn*xixm(i,j,n)
        pbi=pb(i-1,j,nl)
        ub(i,j,nl)=ubo+sqrt(xi/pbi)*(pbi-pbo)
      enddo
      enddo
      i=ii-1
      do l=1,jsv(i)
      do j=jfv(i,l),jlv(i,l)
        vb(i,j,nl)=wt1*vb_oe1(j)+wt2*vb_oe2(j)
      enddo
      enddo
c
c --- "Southern" boundary
      j=1
      do l=1,isu(j)
      do i=ifu(j,l),ilu(j,l)
        ub(i,j,nl)=wt1*ub_os1(i)+wt2*ub_os2(i)
      enddo
      enddo
      j=2
      do l=1,isv(j)
      do i=ifv(j,l),ilv(j,l)
        pb(i,j-1,nl)=wt1*pb_os1(i,1)+wt2*pb_os2(i,1)
        pbo=wt1*pb_os1(i,2)+wt2*pb_os2(i,2)
        vbo=wt1*vb_os1(i)+wt2*vb_os2(i)
        xi=(1.-wn)*xiyp(i,j,m)+wn*xiyp(i,j,n)
        pbi=pb(i,j,nl)
        vb(i,j,nl)=vbo+sqrt(xi/pbi)*(pbo-pbi)
      enddo
      enddo
c
c --- "Northern" boundary
      j=jj-1
      do l=1,isu(j)
      do i=ifu(j,l),ilu(j,l)
        ub(i,j,nl)=wt1*ub_on1(i)+wt2*ub_on2(i)
      enddo
      enddo
      j=jj-1
      do l=1,isv(j)
      do i=ifv(j,l),ilv(j,l)
        pb(i,j,nl)=wt1*pb_on1(i,1)+wt2*pb_on2(i,1)
        pbo=wt1*pb_on1(i,2)+wt2*pb_on2(i,2)
        vbo=wt1*vb_on1(i)+wt2*vb_on2(i)
        xi=(1.-wn)*xiym(i,j,m)+wn*xiym(i,j,n)
        pbi=pb(i,j-1,nl)
        vb(i,j,nl)=vbo+sqrt(xi/pbi)*(pbi-pbo)
      enddo
      enddo
c
      return
      end
