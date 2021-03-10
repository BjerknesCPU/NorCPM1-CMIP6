      subroutine intlog1to3 (ncol, nlons, ind, lev, kcomp, xctin, &
                             Nnat, cxs, xstdv, xrk)

!	Created by Trude Storelvmo, fall 2007. This subroutine gives as output   
!	the "new" modal radius and standard deviation for a given aerosol mode, kcomp 
!       1-3. These parameters are calculated for a best lognormal fit approximation of
!       the aerosol size distribution. This because the aerosol activation routine
!       (developed by Abdul-Razzak & Ghan, 2000) requiers the size distribution to be described 
! 	by lognormal modes.
!


      use shr_kind_mod, only: r8 => shr_kind_r8
      use pmgrid
      use ppgrid
      use const
      use opttab,   only: nbmodes, nmodes, cate

      implicit none

      real(r8), intent(in) :: Nnat(pcols)        ! Modal number concentration
      real(r8), intent(in) :: xctin(pcols)	 ! total internally mixed conc. (ug/m3)	
      real(r8), intent(out) :: xstdv(pcols)      ! log10 of standard deviation for lognormal fit
      real(r8), intent(out) :: xrk(pcols)        ! Modal radius for lognormal fit
      real(r8), intent(inout) :: cxs(pcols,nbmodes) ! excess (modal) internally mixed conc.
      integer, intent(in) :: ncol
      integer, intent(in) :: nlons
      integer, intent(in) :: ind(pcols)
      integer, intent(in) :: lev
      integer, intent(in) :: kcomp

      real(r8) xctsave, camdiff, eps
      real(r8) xct(pcols)

      integer lon, long

      integer i, ictot, ict1, ict2

      real(r8) r1, r2, s1, s2, a , b ,e

      real esssf10, ess

       !parameter (e=2.718281828_r8, eps=1.0e-60)
       parameter (e=2.718281828_r8, eps=1.0e-10_r8)

!     Initialize excess mass cxs, wrt. maximum allowed internal mixing
      do lon=1,ncol
        cxs(lon,kcomp) = 0.0_r8
        xct(lon)  = 0.0_r8
      enddo

      do long=1,nlons
       lon=ind(long)
       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8
       if(Nnat(lon).gt.0.0_r8) then   

	xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))

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

!	Interpolation:

      esssf10= 1.0_r8/(cate(kcomp,ict2)-cate(kcomp,ict1))
	
      r1=rrr1to3(kcomp,ict1)
      r2=rrr1to3(kcomp,ict2)
      s1=sss1to3(kcomp,ict1)
      s2=sss1to3(kcomp,ict2)

	xstdv(lon) = ((cate(kcomp,ict2)-xct(lon))*s1 + &
	(xct(lon)-cate(kcomp,ict1))*s2)*esssf10
	xrk(lon) = ((cate(kcomp,ict2)-xct(lon))*r1 + &
	(xct(lon)-cate(kcomp,ict1))*r2)*esssf10
       endif    ! Nnatk.gt.0
      end do   ! lon

      return
      end

