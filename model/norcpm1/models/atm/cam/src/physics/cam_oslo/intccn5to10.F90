      subroutine intccn5to10 (ncol, nlons, ind, lev, kcomp, xsupin, xctin, Nnat, &
                         xfacin, xfbcin, xfaqin, xfccn, cxs)

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

      real(r8), intent(in) :: xsupin(pcols)        ! relative humidity (>1)
      real(r8), intent(in) :: xctin(pcols)	   ! total internally mixed conc. (ug/m3)	
      real(r8), intent(in) :: xfacin(pcols)        ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
      real(r8), intent(in) :: xfbcin(pcols)        ! = Cbc/(Cbc+Coc)
      real(r8), intent(in) :: xfaqin(pcols)        ! = Cso4a2/(Cso4a1+Cso4a2)
      real(r8), intent(out) :: xfccn(pcols)      ! CCN/Nnat
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff
      real(r8) xsup(pcols), xct(pcols), xfac(pcols),  xfbc(pcols), xfaq(pcols)  

!ces: integer arrays isup1, isup2, ict1, ict2, ifbc1, ifbc2, ifaq1 and
!    ifaq2 substituted with scalar variables with the same name.

      integer lon, long

      integer i, isup, isupx, ictot, ifac, ifbc, ifaq, &
        isup1, isup2, ict1, ict2, ifac1, ifac2, &
        ifbc1, ifbc2, ifaq1, ifaq2

      real(r8) f1, f2, f11, f12, f21, f22, f111, f121, f112, f211, &
        f122, f212, f221, f222, f1111, f1112, f1121, f1122, f1211, & 
        f1212, f1221, f1222, f2111, f2112, f2121, f2122, f2211, & 
        f2212, f2221, f2222, f11111, f11112, f11121, f11122, f11211, & 
        f11212, f11221, f11222, f12111, f12112, f12121, f12122, & 
        f12211, f12212, f12221, f12222, f21111, f21112, f21121, & 
        f21122, f21211, f21212, f21221, f21222, f22111, f22112, & 
        f22121, f22122, f22211, f22212, f22221, f22222, a , b, c, d

!ces: New local variables introduced by (or inspired by) Egil Stoeren:
      real(r8) esssf1, esssf2, esssf3, esssf4, esssf5, esssf6, esssf7, &
       esssf8, esssf9, esssf10, esssf3x6x9x10, ess



!     Initialize excess mass cxs, wrt. maximum allowed internal mixing
      do lon=1,ncol
        cxs(lon,kcomp) = 0.0_r8
        xct(lon)  = 0.0_r8
        xfac(lon) = 0.0_r8
        xfbc(lon) = 0.0_r8
        xfaq(lon) = 0.0_r8
!        xfccn(lon) = 0.0_r8
      enddo

!ces: All loops "do long=1,nlons" combined to one loop:

!      do lon=1,ncol
      do long=1,nlons
       lon=ind(long)
       xfccn(lon) = 0.0_r8

       if(Nnat(lon).gt.0.0_r8) then    ! diffrent cdnc if removed!?

!  write(*,*) 'xsup, xct =', xsup(lon), xct(lon)
!  write(*,*) 'xfac, xfbc, xfaq =', xfac(lon), xfbc(lon), xfaq(lon)

	xsup(lon) = min(max(xsupin(lon),S(1)),S(9))
	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cat(kcomp,1)),cat(kcomp,6))
	xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
	xfbc(lon) = min(max(xfbcin(lon),fbc(1)),fbc(6))
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
	do while (ictot.lt.5.and.(ess.lt.cat(kcomp,ictot).or. &
                ess.gt.cat(kcomp,ictot+1)))
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

	ifbc=1
	ess = xfbc(lon)
	do while (ifbc.lt.5.and.(ess.lt.fbc(ifbc).or. &
                 ess.gt.fbc(ifbc+1)))
	 ifbc=ifbc+1
	enddo
	ifbc1=ifbc
	ifbc2=ifbc+1

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
      esssf4 = (fbc(ifbc2)-xfbc(lon))
      esssf5 = (xfbc(lon)-fbc(ifbc1))
      esssf6 = 1.0_r8/(fbc(ifbc2)-fbc(ifbc1))
      esssf7 = (fac(ifac2)-xfac(lon))
      esssf8 = (xfac(lon)-fac(ifac1))
      esssf9 = 1.0_r8/(fac(ifac2)-fac(ifac1))
      esssf10= 1.0_r8/(cat(kcomp,ict2)-cat(kcomp,ict1))
      esssf3x6x9x10 = esssf3*esssf6*esssf9*esssf10

!  interpolated (in 5 dimensions) ccn-fraction for each kcomp:

      f11111=fff(kcomp,isup1,ict1,ifac1,ifbc1,ifaq1)
      f11112=fff(kcomp,isup1,ict1,ifac1,ifbc1,ifaq2)
      f11121=fff(kcomp,isup1,ict1,ifac1,ifbc2,ifaq1)
      f11122=fff(kcomp,isup1,ict1,ifac1,ifbc2,ifaq2)
      f11211=fff(kcomp,isup1,ict1,ifac2,ifbc1,ifaq1)
      f11212=fff(kcomp,isup1,ict1,ifac2,ifbc1,ifaq2)
      f11221=fff(kcomp,isup1,ict1,ifac2,ifbc2,ifaq1)
      f11222=fff(kcomp,isup1,ict1,ifac2,ifbc2,ifaq2)
      f12111=fff(kcomp,isup1,ict2,ifac1,ifbc1,ifaq1)
      f12112=fff(kcomp,isup1,ict2,ifac1,ifbc1,ifaq2)
      f12121=fff(kcomp,isup1,ict2,ifac1,ifbc2,ifaq1)
      f12122=fff(kcomp,isup1,ict2,ifac1,ifbc2,ifaq2)
      f12211=fff(kcomp,isup1,ict2,ifac2,ifbc1,ifaq1)
      f12212=fff(kcomp,isup1,ict2,ifac2,ifbc1,ifaq2)
      f12221=fff(kcomp,isup1,ict2,ifac2,ifbc2,ifaq1)
      f12222=fff(kcomp,isup1,ict2,ifac2,ifbc2,ifaq2)
      f21111=fff(kcomp,isup2,ict1,ifac1,ifbc1,ifaq1)
      f21112=fff(kcomp,isup2,ict1,ifac1,ifbc1,ifaq2)
      f21121=fff(kcomp,isup2,ict1,ifac1,ifbc2,ifaq1)
      f21122=fff(kcomp,isup2,ict1,ifac1,ifbc2,ifaq2)
      f21211=fff(kcomp,isup2,ict1,ifac2,ifbc1,ifaq1)
      f21212=fff(kcomp,isup2,ict1,ifac2,ifbc1,ifaq2)
      f21221=fff(kcomp,isup2,ict1,ifac2,ifbc2,ifaq1)
      f21222=fff(kcomp,isup2,ict1,ifac2,ifbc2,ifaq2)
      f22111=fff(kcomp,isup2,ict2,ifac1,ifbc1,ifaq1)
      f22112=fff(kcomp,isup2,ict2,ifac1,ifbc1,ifaq2)
      f22121=fff(kcomp,isup2,ict2,ifac1,ifbc2,ifaq1)
      f22122=fff(kcomp,isup2,ict2,ifac1,ifbc2,ifaq2)
      f22211=fff(kcomp,isup2,ict2,ifac2,ifbc1,ifaq1)
      f22212=fff(kcomp,isup2,ict2,ifac2,ifbc1,ifaq2)
      f22221=fff(kcomp,isup2,ict2,ifac2,ifbc2,ifaq1)
      f22222=fff(kcomp,isup2,ict2,ifac2,ifbc2,ifaq2)

      f1111=esssf1*f11111+esssf2*f11112
      f1112=esssf1*f11121+esssf2*f11122
      f1121=esssf1*f11211+esssf2*f11212
      f1122=esssf1*f11221+esssf2*f11222
      f1211=esssf1*f12111+esssf2*f12112
      f1212=esssf1*f12121+esssf2*f12122
      f1221=esssf1*f12211+esssf2*f12212
      f1222=esssf1*f12221+esssf2*f12222
      f2111=esssf1*f21111+esssf2*f21112
      f2112=esssf1*f21121+esssf2*f21122
      f2121=esssf1*f21211+esssf2*f21212
      f2122=esssf1*f21221+esssf2*f21222
      f2211=esssf1*f22111+esssf2*f22112
      f2212=esssf1*f22121+esssf2*f22122
      f2221=esssf1*f22211+esssf2*f22212
      f2222=esssf1*f22221+esssf2*f22222

      f111=esssf4*f1111+esssf5*f1112
      f112=esssf4*f1121+esssf5*f1122
      f121=esssf4*f1211+esssf5*f1212
      f122=esssf4*f1221+esssf5*f1222
      f211=esssf4*f2111+esssf5*f2112
      f212=esssf4*f2121+esssf5*f2122
      f221=esssf4*f2211+esssf5*f2212
      f222=esssf4*f2221+esssf5*f2222

      f11 =esssf7*f111+esssf8*f112
      f12 =esssf7*f121+esssf8*f122
      f21 =esssf7*f211+esssf8*f212
      f22 =esssf7*f221+esssf8*f222

      f1=((cat(kcomp,ict2)-xct(lon))*f11+ &
          (xct(lon)-cat(kcomp,ict1))*f12) &
            *esssf3x6x9x10
      f2=((cat(kcomp,ict2)-xct(lon))*f21+ &
          (xct(lon)-cat(kcomp,ict1))*f22) &
            *esssf3x6x9x10

!     old, linear interpolation: 
!t
!      xfccn(lon)=((S(isup2)-xsup(lon))*f1+ &
!          (xsup(lon)-S(isup1))*f2)         &
!          /(S(isup2)-S(isup1))    
!t      
!     new, average of linear and exponential interpolation:
       a=(log(f2)-log(f1))/(S(isup2)-S(isup1))
       b=(S(isup2)*log(f1)-S(isup1)*log(f2))/(S(isup2)-S(isup1))

       xfccn(lon)=0.5_r8*(e**(a*xsup(lon)+b) & 
          +((S(isup2)-xsup(lon))*f1+(xsup(lon)-S(isup1))*f2) &
          /(S(isup2)-S(isup1)))

!       write(6,*) 'xfccn ',lon,a,b,e,xsup(lon)

       endif    ! Nnatk.gt.0

      end do   ! lon

      return
end subroutine intccn5to10

