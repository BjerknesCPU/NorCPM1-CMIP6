      subroutine intccn1to3 (ncol, nlons, ind, lev, kcomp, xsupin, xctin, &
                             Nnat, xfccn, cxs)

!    Optimized (CCM-version) by Egil Stoeren/NoSerC 2002-05-02 (comments:ces)
!    Updated for use of new CCN-tables by Alf Kirkevaag 2003-05-02, and again 
!    with OC included 2003-02-09. Rewritten for CAM by Trude Storelvmo and 
!    Alf Kirkevaag, November 2004. Modified for new aerosol schemes by Alf 
!    Kirkevaag in January 2006.

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use ppgrid
      use opttab
      use const

      implicit none

      real(r8), intent(in) :: Nnat(pcols)        ! Modal number concentration

      real(r8), intent(in) :: xsupin(pcols)      ! relative humidity (>1)
      real(r8), intent(in) :: xctin(pcols)	 ! total internally mixed conc. (ug/m3)	
      real(r8), intent(out) :: xfccn(pcols)      ! CCN/Nnat
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff
      real(r8) xsup(pcols), xct(pcols)

!ces: integer arrays isup1, isup2, ict1, ict2, 
!    substituted with scalar variables with the same name.

      integer lon, long

      integer i, isup, isupx, ictot, isup1, isup2, ict1, ict2

      real(r8) f1, f2, f11, f12, f21, f22, a , b 

!ces: New local variables introduced by (or inspired by) Egil Stoeren:
      real(r8) esssf10, ess



!     Initialize excess mass cxs, wrt. maximum allowed internal mixing
      do lon=1,ncol
        cxs(lon,kcomp) = 0.0_r8
        xct(lon)  = 0.0_r8
      enddo

!ces: All loops "do long=1,nlons" combined to one loop:

!      do lon=1,ncol
      do long=1,nlons
       lon=ind(long)
       xfccn(lon) = 0.0_r8
       if(Nnat(lon).gt.0.0_r8) then    ! diffrent cdnc if removed!?

!  write(*,*) 'xsup, xct =', xsup(lon), xct(lon)

	xsup(lon) = min(max(xsupin(lon),S(1)),S(9))
	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))

        camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

        cxs(lon,kcomp)  = max(0.0_r8,camdiff)

!        if(lon.eq.1) then
!          write(*,*) 'ind: k, kcomp, cxs =', lev, kcomp, cxs(lon,kcomp)
!        endif

	isup=1
	ess = xsup(lon)
	do while (isup.lt.8.and.(ess.lt.S(isup).or.ess.gt.S(isup+1)))
	 isup=isup+1
	enddo
        isup1=isup
        isup2=isup+1

        ictot=1
	ess = xct(lon)
	do while (ictot.lt.15.and.(ess.lt.cate(kcomp,ictot).or. &
                ess.gt.cate(kcomp,ictot+1)))
	 ictot=ictot+1
	enddo
	ict1=ictot
	ict2=ictot+1

      esssf10= 1.0_r8/(cate(kcomp,ict2)-cate(kcomp,ict1))

!  interpolated (in 2 dimensions) ccn-fraction for each kcomp:

      f11=fff1to3(kcomp,isup1,ict1)
      f12=fff1to3(kcomp,isup1,ict2)
      f21=fff1to3(kcomp,isup2,ict1)
      f22=fff1to3(kcomp,isup2,ict2)

      f1=((cate(kcomp,ict2)-xct(lon))*f11+ &
          (xct(lon)-cate(kcomp,ict1))*f12)*esssf10
      f2=((cate(kcomp,ict2)-xct(lon))*f21+ &
          (xct(lon)-cate(kcomp,ict1))*f22)*esssf10


!     old, linear interpolation: 
!t
      xfccn(lon)=((S(isup2)-xsup(lon))*f1+ &
          (xsup(lon)-S(isup1))*f2)         &
          /(S(isup2)-S(isup1))    
!t
      
!     new, average of linear and exponential interpolation:
!t       a=(log(f2)-log(f1))/(S(isup2)-S(isup1))
!t       b=(S(isup2)*log(f1)-S(isup1)*log(f2))/(S(isup2)-S(isup1))
!t       xfccn(lon)=0.5_r8*(e**(a*xsup(lon)+b) & 
!t          +((S(isup2)-xsup(lon))*f1+(xsup(lon)-S(isup1))*f2) &
!t          /(S(isup2)-S(isup1)))

!         if(xfccn(lon).le.0._r8) then
!          write(*,*)'INTCCNB:xsup,xct=',xsup(lon),xct(lon)
!          write(*,*)'INTCCNB:S1,S2=',S(isup1),S(isup2)
!          write(*,*)'INTCCNB:lon,kcomp,xfccn=',lon,kcomp,xfccn(lon)
!          write(*,*)'INTCCNB:Nnat(lon)=',Nnat(lon)
!         endif

!         write(*,*)'INTCCNB:lon,kcomp,xfccn=',lon,kcomp,xfccn(lon)

       endif    ! Nnatk.gt.0
      end do   ! lon

      return
      end

