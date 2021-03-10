      subroutine intlog5to10 (ncol, nlons, ind, lev, kcomp, xctin, Nnat, &
                         xfacin, xfbcin, xfaqin, cxs, xstdv, xrk)

!Created by Trude Storelvmo, fall 2007, based on method of A. Kirkevag. 
!This subroutine gives as output the "new" modal radius and standard deviation 
!for a given aerosol mode, kcomp 1-5. These parameters are calculated for a 
!best lognormal fit approximation of the aerosol size distribution. 
! This because the aerosol activation routine (developed by Abdul-Razzak & Ghan
!, 2000) requires the size distribution to be described by lognormal modes.
!

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use ppgrid
      use const
      use opttab,   only: nbmodes, nmodes, nbmp1,cat,fbc,fac,faq

      implicit none

      real(r8), intent(in) :: Nnat(pcols)          ! Modal number concentration
      real(r8), intent(in) :: xctin(pcols)	   ! total internally mixed conc. (ug/m3)	
      real(r8), intent(in) :: xfacin(pcols)        ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
      real(r8), intent(in) :: xfbcin(pcols)        ! = Cbc/(Cbc+Coc)
      real(r8), intent(in) :: xfaqin(pcols)        ! = Cso4a2/(Cso4a1+Cso4a2)
      real(r8), intent(out) :: xstdv(pcols)        ! log10 of standard deviation of lognormal fit
      real(r8), intent(out) :: xrk(pcols)          ! Modal radius of lognormal fit
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff, eps
      real(r8) xct(pcols), xfac(pcols),  xfbc(pcols), xfaq(pcols)  

      integer lon, long

      integer i, ictot, ifac, ifbc, ifaq, &
        ict1, ict2, ifac1, ifac2, &
        ifbc1, ifbc2, ifaq1, ifaq2

      real(r8) r1, r2, r11, r12, r21, r22, r111, r121, r112, r211, &
        r122, r212, r221, r222, r1111, r1112, r1121, r1122, r1211, & 
        r1212, r1221, r1222, r2111, r2112, r2121, r2122, r2211, & 
        r2212, r2221, r2222, e, s1, s2, s11, s12, s21, s22, s111, s121, s112, s211, &
        s122, s212, s221, s222, s1111, s1112, s1121, s1122, s1211, & 
        s1212, s1221, s1222, s2111, s2112, s2121, s2122, s2211, & 
        s2212, s2221, s2222

      real esssf1, esssf2, esssf3, esssf4, esssf5, esssf6, esssf7, &
       esssf8, esssf9, esssf10, esssf3x6x9x10, ess

       !parameter (e=2.718281828_r8, eps=1.0e-60)_r8
       parameter (e=2.718281828_r8, eps=1.0e-10_r8)


!     Initialize excess mass cxs, wrt. maximum allowed internal mixing
      do lon=1,ncol
        cxs(lon,kcomp) = 0.0_r8
        xct(lon)  = 0.0_r8
        xfac(lon) = 0.0_r8
        xfbc(lon) = 0.0_r8
        xfaq(lon) = 0.0_r8
      enddo

!ces: All loops "do long=1,nlons" combined to one loop:

!      do lon=1,ncol
      do long=1,nlons
       lon=ind(long)
       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8
       if(Nnat(lon).gt.0.0_r8) then    ! diffrent cdnc if removed!?

	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cat(kcomp,1)),cat(kcomp,6))
	xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
	xfbc(lon) = min(max(xfbcin(lon),fbc(1)),fbc(6))
	xfaq(lon) = min(max(xfaqin(lon),faq(1)),faq(6))

        camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

        cxs(lon,kcomp)  = max(0.0_r8,camdiff)

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

!  interpolated (in 4 dimensions) radius and standard deviation for each kcomp:

      r1111=rrr(kcomp,ict1,ifac1,ifbc1,ifaq1)
      r1112=rrr(kcomp,ict1,ifac1,ifbc1,ifaq2)
      r1121=rrr(kcomp,ict1,ifac1,ifbc2,ifaq1)
      r1122=rrr(kcomp,ict1,ifac1,ifbc2,ifaq2)
      r1211=rrr(kcomp,ict1,ifac2,ifbc1,ifaq1)
      r1212=rrr(kcomp,ict1,ifac2,ifbc1,ifaq2)
      r1221=rrr(kcomp,ict1,ifac2,ifbc2,ifaq1)
      r1222=rrr(kcomp,ict1,ifac2,ifbc2,ifaq2)
      r2111=rrr(kcomp,ict2,ifac1,ifbc1,ifaq1)
      r2112=rrr(kcomp,ict2,ifac1,ifbc1,ifaq2)
      r2121=rrr(kcomp,ict2,ifac1,ifbc2,ifaq1)
      r2122=rrr(kcomp,ict2,ifac1,ifbc2,ifaq2)
      r2211=rrr(kcomp,ict2,ifac2,ifbc1,ifaq1)
      r2212=rrr(kcomp,ict2,ifac2,ifbc1,ifaq2)
      r2221=rrr(kcomp,ict2,ifac2,ifbc2,ifaq1)
      r2222=rrr(kcomp,ict2,ifac2,ifbc2,ifaq2)
      s1111=sss(kcomp,ict1,ifac1,ifbc1,ifaq1)
      s1112=sss(kcomp,ict1,ifac1,ifbc1,ifaq2)
      s1121=sss(kcomp,ict1,ifac1,ifbc2,ifaq1)
      s1122=sss(kcomp,ict1,ifac1,ifbc2,ifaq2)
      s1211=sss(kcomp,ict1,ifac2,ifbc1,ifaq1)
      s1212=sss(kcomp,ict1,ifac2,ifbc1,ifaq2)
      s1221=sss(kcomp,ict1,ifac2,ifbc2,ifaq1)
      s1222=sss(kcomp,ict1,ifac2,ifbc2,ifaq2)
      s2111=sss(kcomp,ict2,ifac1,ifbc1,ifaq1)
      s2112=sss(kcomp,ict2,ifac1,ifbc1,ifaq2)
      s2121=sss(kcomp,ict2,ifac1,ifbc2,ifaq1)
      s2122=sss(kcomp,ict2,ifac1,ifbc2,ifaq2)
      s2211=sss(kcomp,ict2,ifac2,ifbc1,ifaq1)
      s2212=sss(kcomp,ict2,ifac2,ifbc1,ifaq2)
      s2221=sss(kcomp,ict2,ifac2,ifbc2,ifaq1)
      s2222=sss(kcomp,ict2,ifac2,ifbc2,ifaq2)

      r111=esssf1*r1111+esssf2*r1112
      r112=esssf1*r1121+esssf2*r1122
      r121=esssf1*r1211+esssf2*r1212
      r122=esssf1*r1221+esssf2*r1222
      r211=esssf1*r2111+esssf2*r2112
      r212=esssf1*r2121+esssf2*r2122
      r221=esssf1*r2211+esssf2*r2212
      r222=esssf1*r2221+esssf2*r2222
      s111=esssf1*s1111+esssf2*s1112
      s112=esssf1*s1121+esssf2*s1122
      s121=esssf1*s1211+esssf2*s1212
      s122=esssf1*s1221+esssf2*s1222
      s211=esssf1*s2111+esssf2*s2112
      s212=esssf1*s2121+esssf2*s2122
      s221=esssf1*s2211+esssf2*s2212
      s222=esssf1*s2221+esssf2*s2222

      r11=esssf4*r111+esssf5*r112
      r12=esssf4*r121+esssf5*r122
      r21=esssf4*r211+esssf5*r212
      r22=esssf4*r221+esssf5*r222
      s11=esssf4*s111+esssf5*s112
      s12=esssf4*s121+esssf5*s122
      s21=esssf4*s211+esssf5*s212
      s22=esssf4*s221+esssf5*s222

      r1 =esssf7*r11+esssf8*r12
      r2 =esssf7*r21+esssf8*r22
      s1 =esssf7*s11+esssf8*s12
      s2 =esssf7*s21+esssf8*s22

      xrk(lon)=((cat(kcomp,ict2)-xct(lon))*r1+ &
          (xct(lon)-cat(kcomp,ict1))*r2) &
            *esssf3x6x9x10
      xstdv(lon)=((cat(kcomp,ict2)-xct(lon))*s1+ &
          (xct(lon)-cat(kcomp,ict1))*s2) &
            *esssf3x6x9x10

       endif    ! Nnatk.gt.0
      end do   ! lon

      return
end subroutine intlog5to10

