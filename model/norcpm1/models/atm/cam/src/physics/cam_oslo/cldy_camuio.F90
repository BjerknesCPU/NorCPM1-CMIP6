      subroutine cldy_camuio(lchnk,   ncol,    ncnst,          &
                   so2,       t,       pmid,    cldfrc,      &
                   qc,        conicw, khet,  &
                   lo3,       lh2o2,  frach2) 

! Calculate the pH, solubility constants of SO2 and H2O2, and the 
! aqueous-phase reaction rate. Based on a routine in CAM3. ØS Jan 2011

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use physconst, only : mwdry
use aqrat,     only : calchet

implicit none
!include 'chemrates.h'   

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: ncnst                ! number of tracers 
   real(r8), intent(in) :: so2(pcols,pver)      ! SO2 
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction
   real(r8), intent(in) :: qc(pcols,pver)       ! Liquid cloud water
   real(r8), intent(in) :: conicw (pcols,pver)  ! in cloud water mixing ratio for convection
!   real(r8), intent(in) :: icwmr2 (pcols,pver)  ! in cloud water mixing ration for hack  scheme
   real(r8), intent(out) :: khet(pcols,pver)  ! Oxidation rate
   real(r8), intent(in) :: lo3(pcols,pver) ! O3 concentration
   real(r8), intent(in) :: lh2o2(pcols,pver) ! H2O2 concentration
   real(r8), intent(out) :: frach2(pcols,pver) ! H2O2 loss rate

! Local variables
      real(r8) :: ekh2o2(pcols,pver)     !H2O2 Henry's Law coefficient
      real(r8) :: ekso2(pcols,pver)      !SO2  effective Henry's Law coefficient
!      real(r8) :: frach2(pcols,pver)
      integer :: ind(pcols,pver)     !indices for cloudy grid points
      integer :: ncpts(pver)         !number of cloudy grid points
      real(r8) :: o3(pcols,pver),h2o2(pcols,pver)
!      real(r8) :: cldr(pcols,pver)
!c Local variables:

      integer :: i, k
!      integer loni(pcols)
      real(r8) :: cwat(pcols,pver)       !cloud water
!      real rwat(pcols,pver)       !rain
      real(r8) :: twat(pcols,pver)       !total liquid water = cwat 
!+ rwat
      real(r8) :: molso2 (pcols)  !so2 in mol/mol units At the moment in the model

      integer :: i2                  !indice for cloudy grid points
      integer :: ilw                 !number of cloudy grid points

      real(r8) :: tso2(pcols)              !so2 at cloudy grid points
      real(r8) :: ttwat(pcols)              !twat at cloudy grid points
      real(r8) :: temp(pcols)               !t at cloudy grid points
      real(r8) :: tprs(pcols)               !pres at cloudy grid points
!      real :: tph(pcols)                !ph at cloudy grid points
!      real :: thion(pcols)              !hion at cloudy grid points
!      real :: theff(pcols)              !ekso2 at cloudy grid points
!os
      real(r8) :: hso2(pcols),hs1(pcols),hs2(pcols),ho3(pcols),hh2o2(pcols)
      real(r8) :: xkh2o2(pcols),xk1o3(pcols),xk2o3(pcols),xk3o3(pcols)
      real(r8) :: xso2(pcols),xh2o2(pcols),xo3(pcols),xcldv(pcols)
      real(r8) :: xkhet(pcols),frh2o2(pcols)
!      real rform(pcols,pver)

!os
      real(r8) :: weight
      real(r8) :: tc                     !temperature in deg C


! Function:
      real (r8) :: e298, dhr, tt, ek
      ek(e298, dhr, tt) = e298*exp(dhr*(1._r8/tt - 1._r8/298._r8))
!      kh2o2(:,:)=0._r8
      frach2(:,:)=0._r8

      do k=1,pver
       do i=1,ncol
!c         ekh2o2(i,k) = ek(7.4e4_r8, 6621._r8, tair(i,k))
!c         ekso2(i,k)  = ek(1.23_r8,  3120._r8, tair(i,k)) 

          ekh2o2(i,k) = ek(7.4e4_r8, 7302._r8, t(i,k))
          ekso2(i,k)  = ek(1.23_r8,  3147._r8, t(i,k)) 
!c Determine cloud water and rain water mixing ratios
        tc     = t(i,k) - 273.16_r8
        weight = max(0.,min(-tc*0.04_r8,1.0_r8)) ! fraction of condensate that is ice
        cwat(i,k) = (qc(i,k)/max(cldfrc(i,k), 0.00001_r8) +  &
!c add suspended water from convection and do aqchem in only the liquid phase
             conicw(i,k))*(1._r8-weight)


!        if(tair(i,k) .gt. 273.16_r8) then
!         weight = 0._r8                    ! fraction of condensate that is ice
!        else
!         weight = 1._r8
!        endif
!        rwat(i,k) = qr(i,k)/max(cldv(i,k), 0.00001_r8)
!     $                *(1._r8-weight) ! aqchem in only the liquid phase

!c Sum cwat and rwat 
        twat(i,k) = cwat(i,k) 
!c+ rwat(i,k)
!	rform(i,k)=cmfdqr(i,k)+nrain(i,k)
       end do
      end do

!c Initialize arrays
      do k=1,pver
       do i=1,ncol
!        ph(i,k) = 99._r8
!        hion(i,k) = 0._r8
        ind(i,k) = 0
        khet(i,k)=0._r8
!        cldr(i,k)=0._r8
       end do
      end do


!c Find which grid points have liquid water
      do k=1,pver
       ncpts(k) = 0
       ilw = 0
       do i=1,ncol
!c        if (rform(i,k).gt.0._r8 .and. cldf(i,k).ge.0.02_r8 .and. 
!c     &    cwat(i,k).ge.1.e-12_r8) then
!c           cldr(i,k)=cldf(i,k)
!c	   end if
!c	if cmfdqr(i,k)
        if(cldfrc(i,k) .ge. 0.02_r8 .and. twat(i,k) .ge. 1.e-12_r8) then
          ilw = ilw + 1
          ind(ilw,k) = i
        endif
       end do




       ncpts(k) = ilw       !assign number of cloudy grid points to array

       do i2=1,ilw
!        loni(i2) = ind(i2,k)
        tso2(i2) = so2(ind(i2,k),k)
!c        tso4(i2) = so4(ind(i2,k),k)
        ttwat(i2) = twat(ind(i2,k),k)
        temp(i2) = t(ind(i2,k),k)
        tprs(i2) = pmid(ind(i2,k),k)
        xh2o2(i2) = lh2o2(ind(i2,k),k) 
        xo3(i2) = lo3(ind(i2,k),k)
        xcldv(i2) = cldfrc(ind(i2,k),k) 
!c Set pH as a function of HSO3- and SO4=
!c  Say NH4+ = 1.5 SO4= and NO3- = 0.5 SO4=
!c  Thus, H+ = HSO3- + SO4=  as done in ECHAM
!c
!COS tso2 multiplied by since it is given as amount of S
! OS 15/7 2003 Not multiplied at this time due to mol/mol concenntrations

!        molso2 = 2._r8*tso2(i2)  * 28.966_r8/64._r8
!        xso2(i2)=molso2

         xso2(i2)=2._r8*tso2(i2) * mwdry/64.1_r8
!c        molarso4 = tso4(i2)/(max(ttwat(i2), 1.e-30_r8)*96.e-3_r8) 
!c        ahso3 = ek(1.23_r8,3120._r8,temp(i2)) * ek(1.3e-2_r8,2015._r8,temp(i2)) *
!c     _            molso2 * tprs(i2)/101325._r8
!c        thion(i2) = (molarso4 + sqrt(max(molarso4*molarso4 + 4._r8*ahso3,
!c     _             1.e-60)))*0.5_r8
! 	thion(i2)=1.e-4_r8
!        tph(i2) = -alog10(thion(i2))



!c        write(6,*) 'ph ',thion(i2),tph(i2)

!c        theff(i2) = ek(1.23_r8, 3120._r8, temp(i2)) * (1._r8 + 
!c     _              ek(1.3e-2_r8, 2015._r8, temp(i2))/thion(i2) * (1._r8 +
!c     _              ek(6.31e-8_r8, 1505._r8, temp(i2))/thion(i2)))
!        theff(i2) = ek(1.23_r8, 3120._r8, temp(i2)) * (1._r8 + 
!     _              ek(1.3e-2_r8, 1960._r8, temp(i2))/thion(i2) * (1._r8 +
!     _              ek(6.6e-8_r8, 1500._r8, temp(i2))/thion(i2)))
!c	if (tso4(i2).gt.1.e-12_r8) then
!c	write(6,*) i2,k,temp(i2),thion(i2),tph(i2),ahso3,
!c     - theff(i2),tso4(i2),ekso2(ind(i2,k),k)
!c	end if
       end do
       do i2=1,ilw
	 hso2(i2)=ekso2(ind(i2,k),k)
	 hs1(i2)=ek(1.3e-2_r8,1960._r8,temp(i2))
	 hs2(i2)=ek(6.6e-8_r8,1500._r8,temp(i2))
         ho3(i2)=ek(1.13e-2_r8,2538._r8,temp(i2))
         hh2o2(i2)=ekh2o2(ind(i2,k),k)
         xkh2o2(i2)=ek(7.5e7_r8,-4430._r8,temp(i2))
         xk1o3(i2)=2.4e4_r8
         xk2o3(i2)=ek(3.7e5_r8,-5530._r8,temp(i2))
         xk3o3(i2)=ek(1.5e9_r8,-5280._r8,temp(i2))
       end do

       call calchet(hso2,hs1,hs2,ho3,hh2o2,  &
          xkh2o2,xk1o3,xk2o3,xk3o3,xkhet,xcldv,  &
          temp,ttwat,xso2,xh2o2,xo3,ilw,frh2o2)	

       do i2=1,ilw
!         hion(ind(i2,k),k) = thion(i2)  
!         ph(ind(i2,k),k) = tph(i2)  
!         ekso2(ind(i2,k),k) = theff(i2)
!        if (xso2(i2).gt.1.e-10_r8) &
         khet(ind(i2,k),k) = xkhet(i2)
         frach2(ind(i2,k),k) = frh2o2(i2)
!         if (xkhet(i2).lt.1.e-6_r8) &
!    write(6,*) 'cldy ',i,k,temp(i2),ttwat(i2),xkhet(i2)
       end do

      end do



      return
      end


