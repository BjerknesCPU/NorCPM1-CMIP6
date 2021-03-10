      subroutine intccn4 (ncol, nlons, ind, lev, kcomp, xsupin, xctin, Nnat, &
                          xfacin, xfaqin, xfccn, cxs)

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

      real(r8), intent(in) :: Nnat(pcols)           ! Modal number concentration

      real(r8), intent(in) :: xsupin(pcols)         ! relative humidity (>1)
      real(r8), intent(in) :: xctin(pcols)	    ! total internally mixed (cond+aq) SO4conc. (ug/m3)	
      real(r8), intent(in) :: xfacin(pcols)         ! Cbc/(Cbc+Coc) for the background mode
      real(r8), intent(in) :: xfaqin(pcols)         ! = Cso4a2/(Cso4a1+Cso4a2)
      real(r8), intent(out) :: xfccn(pcols)         ! CCN/Nnat
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff
      real(r8) xsup(pcols), xct(pcols), xfac(pcols), xfaq(pcols)
!ces: integer arrays isup1, isup2, ict1, ict2, ifaq1 and ifaq2
!     substituted with scalar variables with the same name.

      integer lon, long

      integer i, isup, isupx, ictot, ifac, ifaq, &
        isup1, isup2, ict1, ict2, ifac1, ifac2, ifaq1, ifaq2

      real(r8) f1, f2, f11, f12, f21, f22, f111, f121, f112, f211, &
        f122, f212, f221, f222, f1111, f1112, f1121, f1122, f1211, & 
        f1212, f1221, f1222, f2111, f2112, f2121, f2122, f2211, & 
        f2212, f2221, f2222, a , b

!ces: New local variables introduced by (or inspired by) Egil Stoeren:
      real(r8) esssf1, esssf2, esssf3, esssf7, esssf8, esssf9, esssf10, ess


!     Initialize excess mass cxs, wrt. maximum allowed internal mixing
      do lon=1,ncol
        cxs(lon,kcomp) = 0.0_r8
        xct(lon)  = 0.0_r8
        xfac(lon) = 0.0_r8
        xfaq(lon) = 0.0_r8
      enddo

!ces: All loops "do long=1,nlons" combined to one loop:

!      do lon=1,ncol
      do long=1,nlons
       lon=ind(long)
       xfccn(lon) = 0.0_r8
       if(Nnat(lon).gt.0.0_r8) then    ! diffrent cdnc if removed!?

!  write(*,*) 'xsup, xct =', xsup(lon), xct(lon)
!  write(*,*) 'xfac, xfaq =', xfac(lon), xfaq(lon)

	xsup(lon) = min(max(xsupin(lon),S(1)),S(9))
	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))
	xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
	xfaq(lon) = min(max(xfaqin(lon),faq(1)),faq(6))

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

	ifac=1
	ess = xfac(lon)
	do while (ifac.lt.5.and.(ess.lt.fac(ifac).or. &
                 ess.gt.fac(ifac+1)))
	 ifac=ifac+1
	enddo
	ifac1=ifac
	ifac2=ifac+1

	ifaq=1
	ess = xfaq(lon)
	do while (ifaq.lt.5.and.(ess.lt.faq(ifaq) &
                 .or.ess.gt.faq(ifaq+1)))
	 ifaq=ifaq+1
        enddo
	ifaq1=ifaq
	ifaq2=ifaq+1

      esssf1 = (faq(ifaq2)-xfaq(lon))
      esssf2 = (xfaq(lon)-faq(ifaq1))
      esssf3 = 1.0_r8/(faq(ifaq2)-faq(ifaq1))

      esssf7 = (fac(ifac2)-xfac(lon))
      esssf8 = (xfac(lon)-fac(ifac1))
      esssf9 = 1.0_r8/(fac(ifac2)-fac(ifac1))
      esssf10= 1.0_r8/(cate(kcomp,ict2)-cate(kcomp,ict1))


!  interpolated (in 4 dimensions) ccn-fraction:

      f1111=fff4(isup1,ict1,ifac1,ifaq1)
      f1112=fff4(isup1,ict1,ifac1,ifaq2)
      f1121=fff4(isup1,ict1,ifac2,ifaq1)
      f1122=fff4(isup1,ict1,ifac2,ifaq2)
      f1211=fff4(isup1,ict2,ifac1,ifaq1)
      f1212=fff4(isup1,ict2,ifac1,ifaq2)
      f1221=fff4(isup1,ict2,ifac2,ifaq1)
      f1222=fff4(isup1,ict2,ifac2,ifaq2)
      f2111=fff4(isup2,ict1,ifac1,ifaq1)
      f2112=fff4(isup2,ict1,ifac1,ifaq2)
      f2121=fff4(isup2,ict1,ifac2,ifaq1)
      f2122=fff4(isup2,ict1,ifac2,ifaq2)
      f2211=fff4(isup2,ict2,ifac1,ifaq1)
      f2212=fff4(isup2,ict2,ifac1,ifaq2)
      f2221=fff4(isup2,ict2,ifac2,ifaq1)
      f2222=fff4(isup2,ict2,ifac2,ifaq2)

      f111=esssf1*f1111+esssf2*f1112
      f112=esssf1*f1121+esssf2*f1122
      f121=esssf1*f1211+esssf2*f1212
      f122=esssf1*f1221+esssf2*f1222
      f211=esssf1*f2111+esssf2*f2112
      f212=esssf1*f2121+esssf2*f2122
      f221=esssf1*f2211+esssf2*f2212
      f222=esssf1*f2221+esssf2*f2222

      f11 =esssf7*f111+esssf8*f112
      f12 =esssf7*f121+esssf8*f122
      f21 =esssf7*f211+esssf8*f212
      f22 =esssf7*f221+esssf8*f222

      f1=((cate(kcomp,ict2)-xct(lon))*f11+ &
          (xct(lon)-cate(kcomp,ict1))*f12)*esssf3*esssf9*esssf10
      f2=((cate(kcomp,ict2)-xct(lon))*f21+ &
          (xct(lon)-cate(kcomp,ict1))*f22)*esssf3*esssf9*esssf10

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
end subroutine intccn4

