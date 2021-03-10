      module mod_nudge_ts
c
      implicit none 
#ifdef NUDGE_TS
c
c --- namelist parameters
      logical, save :: nudge_t=.false.,nudge_s=.false.
      real, save :: nudge_nlim=90.,nudge_slim=-90.,nudge_taper=0., 
     .  nudge_tau=1.
      character(len=500), save :: nudge_path
      namelist /NUDGETS/ nudge_path,nudge_t,nudge_s, 
     .                   nudge_nlim,nudge_slim,nudge_taper,nudge_tau 

      contains 

      subroutine nudge_ts(m,n,mm,nn,k1m,k1n)     
c----------------------------------------------------------------------
c     Description: SST and SSS nudging 
c       This routine modfies the upper temperature and salinity layer 
c       according to y_new = (1-dt/tau) X y_old + dt/tau X y_clim.  
c----------------------------------------------------------------------
      use mod_xc
      use dimensions
      use mod_dia
      use mod_nctools
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "common_geo.h"
c     
      integer m,n,mm,nn,k1m,k1n
c
      integer :: jday,m1,m2,k,k2,l,i,j,klast
      integer, save :: m1_old=0,m2_old=0
      real    :: a1,a1tmp,a2,zbot,ztop,w,wacc
      real, save :: r, fillvalue=1e+20
      integer, parameter :: mons(12)=(/1,32,60,91,121,152,182,213,244,
     .  274,305,335/),nlev=33
      real, parameter :: monc(12)=(/16.0,45.5,75.0,105.5,136.0,166.5,
     .  197.0,228.0,258.5,289.0,319.5,350.0/)
      real,parameter :: eps=1e-20 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy), save :: mask=1 
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,nlev), save :: 
     .  tlev,tlev1,tlev2,slev,slev1,slev2
      real, dimension (1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,kdm) :: 
     .  tlay,slay,mask3d,util 
      real, save :: levtop(nlev),levbot(nlev), lev(nlev)= (/0, 10, 20, 
     .  30, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 600, 700, 
     .  800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1750, 2000, 2500, 
     .  3000, 3500, 4000, 4500, 5000, 5500/)
      logical, save :: first=.true.
c     
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: ind1,ind2
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy,ddm) :: wghts,
     .  wghtsflx
c
c --- check if nudging should be done
      if (.not.nudge_t.and..not.nudge_s) return
c
c --- compute tapering mask etc 
      if (first) then 
        first=.false.
        r=baclin/(nudge_tau*3600*24)
        where(plat.lt.nudge_slim.or.plat.gt.nudge_nlim) mask=0 
        where(plat.ge.nudge_slim.and.plat.le.nudge_slim+nudge_taper) 
     .    mask=(plat-nudge_slim)/nudge_taper 
        where(plat.le.nudge_nlim.and.plat.ge.nudge_nlim-nudge_taper) 
     .    mask=1-(plat-nudge_nlim+nudge_taper)/nudge_taper 
c --- - compute depth bounds
        levbot(1)=0
        do k=1,nlev
          levtop(k)=levbot(max(k-1,1))
          levbot(k)=0.5*(lev(k)+lev(min(k+1,nlev)))
        end do 
        levbot(nlev)=8000.
      endif 
c
c --- compute interpolation weigths
      write(lp,*) nyear,'-',nmonth,'-',nday
      jday=MONS(nmonth)+nday-1
      if (jday.eq.MONC(nmonth)) then
        m1=nmonth
        m2=nmonth
        a1=0.5
        a2=0.5
      elseif (jday.lt.MONC(nmonth)) then
        m1=mod(nmonth-2+12,12)+1
        m2=nmonth
        a1=(MONC(m2)-jday)
        if (jday-MONC(m1).gt.0) then
          a2=jday-MONC(m1)
        else
          a2=jday-MONC(m1)+365
        endif
      elseif (jday.gt.MONC(nmonth)) then
        m1=nmonth
        m2=mod(nmonth,12)+1
        a2=(jday-MONC(m1))
        if (MONC(m2)-jday.gt.0) then
          a1=(MONC(m2)-jday)
        else
          a1=(MONC(m2)-jday+365)
        endif
      endif
      a1tmp=a1/(a1+a2)
      a2=a2/(a1+a2)
      a1=a1tmp
      if (abs(a1+a2-1).gt.1e-6) stop 'nudging error'
      write(lp,*) 'm1=',m1,' a1=',a1,'m2=',m2,' a2=',a2
c
      if (m1.ne.m1_old.or.m2.ne.m2_old) then
c 
c --- - read nudging data             
        if (nudge_t) then 
          write(lp,*) 'T nudging read month 1:',m1
          call ncfopn(trim(nudge_path),'r',' ',m1)
          call ncread('Temp',tlev1,ip,0,fillvalue)   
          call ncfcls     
c 
          write(lp,*) 'T nudging read month 2:',m2
          call ncfopn(trim(nudge_path),'r',' ',m2)
          call ncread('Temp',tlev2,ip,0,fillvalue)        
          call ncfcls 
        else 
          tlev1=0.
          tlev2=0.
        end if 
c
        if (nudge_s) then
          write(lp,*) 'S nudging read month 1:',m1
          call ncfopn(trim(nudge_path),'r',' ',m1)
          call ncread('Salt',slev1,ip,0,fillvalue)
          call ncfcls
c 
          write(lp,*) 'S nudging read month 2:',m2
          call ncfopn(trim(nudge_path),'r',' ',m2)
          call ncread('Salt',slev2,ip,0,fillvalue)
          call ncfcls
        else 
          slev1=0.
          slev2=0.
        end if
      end if
c
c --- interpolate to model layers
      tlev=a1*tlev1+a2*tlev2
      slev=a1*slev1+a2*slev2
      do j=1,jj
        do l=1,isp(j)
          do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
            zbot=0
            klast=1
            do k=1,kk
              ztop=zbot
              zbot=zbot+(phi(i,j,k)-phi(i,j,k+1))*1e-2/g
              tlay(i,j,k)=0.
              slay(i,j,k)=0.
              mask3d(i,j,k)=0.
              wacc=0.
              do k2=klast,nlev 
c                write(lp,*) 'DEBUG k,k2,tlev,slev=',k,k2
c                write(lp,*) 'DEBUG tlev=',tlev(i,j,k2)
c                write(lp,*) 'DEBUG slev=',slev(i,j,k2)
                if (levbot(k2).lt.ztop) then 
                  klast=k2
                  cycle 
                else if (levtop(k2).gt.zbot) then 
                  exit 
                else if (tlev(i,j,k2).ne.fillvalue.and.
     .            tlev(i,j,k2).gt.-2..and.tlev(i,j,k2).lt.50..and.
     .            slev(i,j,k2).gt.eps.and.slev(i,j,k2).lt.50..and.
     .            slev(i,j,k2).ne.fillvalue) then  
c                 write(lp,*) 'DEBUG k,k2,ztop,zbot,levtop,levbot=',k,
c    .              k2,ztop,zbot,levtop(k2),levbot(k2)
c                 write(lp,*) 'DEBUG k,k2,tlev,slev=',k,k2,
c    .              tlev(i,j,k2),slev(i,j,k2)
                  w=(min(zbot,levbot(k2))-max(ztop,levtop(k2)))
     .              /max(zbot-ztop,eps)
                  tlay(i,j,k)=tlay(i,j,k)+w*tlev(i,j,k2)
                  slay(i,j,k)=slay(i,j,k)+w*slev(i,j,k2) 
                  wacc=wacc+w
                end if                        
              enddo 
              if (wacc.gt.eps) then  
                tlay(i,j,k)=tlay(i,j,k)/wacc
                slay(i,j,k)=slay(i,j,k)/wacc
                mask3d(i,j,k)=mask(i,j)
              endif 
            enddo 
          enddo
        enddo
      enddo
c
c --- apply nudging 
      if (nudge_t) then      
        util=temp(:,:,1+mm:kk+mm)    
c       write(lp,*) 'DEBUG fillvalue=',fillvalue
c       write(lp,*) 'DEBUG mm,k1m=',mm,k1m
c       write(lp,*) 'DEBUG levbot=',levbot
c       write(lp,*) 'DEBUG levtop=',levtop
c       write(lp,*) 'DEBUG zbot=',zbot
c       write(lp,*) 'DEBUG mask=',mask(1,1) 
c       write(lp,*) 'DEBUG mask3d=',mask3d(1,1,:) 
c       write(lp,*) 'DEBUG temp=',util(1,1,:) 
c       write(lp,*) 'DEBUG tlev=',tlev(1,1,:)
c       write(lp,*) 'DEBUG tlay=',tlay(1,1,:)
c       do j=1,jj
c         do l=1,isp(j)
c           do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c             do k=1,kk
c               if (mask3d(i,j,k).gt.0.5) temp(i,j,k+mm)=
c    .            (1-r)*temp(i,j,k+mm)+r*tlay(i,j,k)
c             enddo
c           enddo
c         enddo
c       enddo
        where (mask3d.gt.0.5) temp(:,:,1+mm:kk+mm)=(1-r)*util+r*tlay
        util=(temp(:,:,1+mm:kk+mm)-util)/baclin
        call acclyr(ACC_TNUDGE,util,util,0,'p')
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
          call acclvl(ACC_TNUDGELVL,util,'p',k,ind1,ind2,wghts)
        enddo
      end if  
c
      if (nudge_s) then
        util=saln(:,:,1+mm:kk+mm)
c       write(lp,*) 'DEBUG saln=',util(1,1,:)
c       write(lp,*) 'DEBUG slev=',slev(1,1,:)
c       write(lp,*) 'DEBUG slay=',slay(1,1,:)
c       do j=1,jj
c         do l=1,isp(j)
c           do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c             do k=1,kk
c               if (mask3d(i,j,k).gt.0.5) saln(i,j,k+mm)=
c    .            (1-r)*saln(i,j,k+mm)+r*slay(i,j,k)
c               if (slay(i,j,k).lt.20.and.mask3d(i,j,k).gt.0.5) then 
c                 write(*,*) 'DEBUG i,j,k=',i,j,k
c                 write(*,*) 'DEBUG slev=',slev 
c                 write(*,*) 'DEBUG slay=',slay 
c                 write(*,*) 'DEBUG smod=',util(i,j,:)
c                 stop 
c               endif 
c             enddo
c           enddo 
c         enddo 
c       enddo
        where (mask3d.gt.0.5) saln(:,:,1+mm:kk+mm)=(1-r)*util+r*slay
        util=(saln(:,:,1+mm:k1m)-util)/baclin
        call acclyr(ACC_SNUDGE,util,util,0,'p')
        do k=1,kk
          call diazlv('p',k,mm,nn,ind1,ind2,wghts,wghtsflx)
          call acclvl(ACC_SNUDGELVL,util,'p',k,ind1,ind2,wghts)
        enddo
      end if
c
c --- remember m1 and m2
      m1_old=m1
      m2_old=m2
c
      end subroutine nudge_ts
#endif
      end module mod_nudge_ts
