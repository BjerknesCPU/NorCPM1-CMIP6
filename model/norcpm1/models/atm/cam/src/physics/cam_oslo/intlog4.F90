      subroutine intlog4 (ncol, nlons, ind, lev, kcomp, xctin, Nnat, &
                          xfacin, xfaqin, cxs, xstdv, xrk)

!   Created by Trude Storelvmo, fall 2007. This subroutine gives as output   
!   the "new" modal radius and standard deviation for aerosol mode kcomp=4.
!   These parameters are calculated for a best lognormal fit approximation of
!   the aerosol size distribution. This because the aerosol activation routine
!   (developed by Abdul-Razzak & Ghan, 2000) requiers the size distribution 
!   to be described 
! 	by lognormal modes.

      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use ppgrid
      use const
      use opttab,   only: nbmodes, nmodes, nbmp1,cate,fac,faq

      implicit none

      real(r8), intent(in) :: Nnat(pcols)           ! Modal number concentration
      real(r8), intent(in) :: xctin(pcols)	    ! total internally mixed conc. (ug/m3)	
      real(r8), intent(in) :: xfacin(pcols)         ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
      real(r8), intent(in) :: xfaqin(pcols)         ! = Cso4a2/(Cso4a1+Cso4a2)
      real(r8), intent(out) :: xstdv(pcols)         ! log10 of standard deviation for lognormal fit
      real(r8), intent(out) :: xrk(pcols)           ! Modal radius for lognormal fit
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff, eps
      real(r8) xct(pcols), xfac(pcols), xfaq(pcols)
!ces: integer arrays ict1, ict2, ifaq1 and ifaq2
!     substituted with scalar variables with the same name.

      integer lon, long

      integer i, ictot, ifac, ifaq, &
        ict1, ict2, ifac1, ifac2, ifaq1, ifaq2

      real(r8) r1, r2, r11, r12, r21, r22, r111, r121, r112, r211, &
        r122, r212, r221, r222 ,e,  s1, s2, s11, s12, s21, s22, s111, &
	s121, s112, s211, s122, s212, s221, s222

!ces: New local variables introduced by (or inspired by) Egil Stoeren:
      real esssf1, esssf2, esssf3, esssf7, esssf8, esssf9, esssf10, ess

       !parameter (e=2.718281828_r8, eps=1.0e-60_r8)
       parameter (e=2.718281828_r8, eps=1.0e-60_r8)


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
!ak+corr31mar08
       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8
!ak-corr31mar08
       if(Nnat(lon).gt.0.0_r8) then    ! diffrent cdnc if removed!?

	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))
	xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
	xfaq(lon) = min(max(xfaqin(lon),faq(1)),faq(6))

        camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

        cxs(lon,kcomp)  = max(0.0_r8,camdiff)

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


!  interpolated (in 3 dimensions) modal radius and standard deviation:

      r111=rrr4(ict1,ifac1,ifaq1)
      r112=rrr4(ict1,ifac1,ifaq2)
      r121=rrr4(ict1,ifac2,ifaq1)
      r122=rrr4(ict1,ifac2,ifaq2)
      r211=rrr4(ict2,ifac1,ifaq1)
      r212=rrr4(ict2,ifac1,ifaq2)
      r221=rrr4(ict2,ifac2,ifaq1)
      r222=rrr4(ict2,ifac2,ifaq2)
      s111=sss4(ict1,ifac1,ifaq1)
      s112=sss4(ict1,ifac1,ifaq2)
      s121=sss4(ict1,ifac2,ifaq1)
      s122=sss4(ict1,ifac2,ifaq2)
      s211=sss4(ict2,ifac1,ifaq1)
      s212=sss4(ict2,ifac1,ifaq2)
      s221=sss4(ict2,ifac2,ifaq1)
      s222=sss4(ict2,ifac2,ifaq2)

      r11=esssf1*r111+esssf2*r112
      r12=esssf1*r121+esssf2*r122
      r21=esssf1*r211+esssf2*r212
      r22=esssf1*r221+esssf2*r222
      s11=esssf1*s111+esssf2*s112
      s12=esssf1*s121+esssf2*s122
      s21=esssf1*s211+esssf2*s212
      s22=esssf1*s221+esssf2*s222

      r1 =esssf7*r11+esssf8*r12
      r2 =esssf7*r21+esssf8*r22
      s1 =esssf7*s11+esssf8*s12
      s2 =esssf7*s21+esssf8*s22

      xrk(lon)=((cate(kcomp,ict2)-xct(lon))*r1+ &
          (xct(lon)-cate(kcomp,ict1))*r2)*esssf3*esssf9*esssf10
      xstdv(lon)=((cate(kcomp,ict2)-xct(lon))*s1+ &
          (xct(lon)-cate(kcomp,ict1))*s2)*esssf3*esssf9*esssf10


       endif    ! Nnatk.gt.0
      end do   ! lon

      return
end subroutine intlog4

