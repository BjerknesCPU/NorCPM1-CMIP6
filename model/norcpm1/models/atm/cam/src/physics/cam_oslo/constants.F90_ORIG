
subroutine constants
!
! A number of constants used in the emission and size-calculation in CAM-Oslo
! ØS Jan 2011
!
!
        use shr_kind_mod, only: r8 => shr_kind_r8
        use physconst,    only: pi
	use const
        use constituents, only: pcnst
        use aerosoldef

	implicit none


	integer kcomp,i
	real e, rhorbc, ntotint
!	real(r8),dimension(10) :: rk
!	real(r8),dimension(10) :: logsk
!	real(r8),dimension(10) :: rhob
	real(r8),dimension(0:10) :: rk
	real(r8),dimension(0:10) :: logsk
	real(r8),dimension(0:10) :: rhob
        e = 2.718281828_r8
        dlogr = 0.1_r8
        sq2pi=1._r8/sqrt(2.0_r8*pi)
!       modal parameters (for local use, and in koagsub)
        rk(0)    = effsize(l_bc_ax)*1.e6_r8
        logsk(0) = log10(sgpart(l_bc_ax))
        rhob(0)  = rhopart(l_bc_ax)       ! mostly not in use (rhorbc in stead)
        rk(1)    = effsize(l_so4_n)*1.e6_r8
        logsk(1) = log10(sgpart(l_so4_n))
        rhob(1)  = rhopart(l_so4_n)
        rk(2)    = effsize(l_bc_n)*1.e6_r8
        logsk(2) = log10(sgpart(l_bc_n))
        rhob(2)  = rhopart(l_bc_n)
        rk(3)    = effsize(l_om_ni)*1.e6_r8
        logsk(3) = log10(sgpart(l_om_ni))
        rhob(3)  = rhopart(l_om_ni)
        rk(4)    = effsize(l_om_ni)*1.e6_r8
        logsk(4) = log10(sgpart(l_om_ni))
        rhob(4)  = rhopart(l_om_ni)
        rk(5)    = effsize(l_so4_pr)*1.e6_r8
        logsk(5) = log10(sgpart(l_so4_pr))
        rhob(5)  = 1841.0_r8                   
        rk(6)    = effsize(l_dst_a2)*1.e6_r8 
        logsk(6) = log10(sgpart(l_dst_a2)) 
        rhob(6)  = rhopart(l_dst_a2) 
        rk(7)    = effsize(l_dst_a3)*1.e6_r8 
        logsk(7) = log10(sgpart(l_dst_a3))
        rhob(7)  = rhopart(l_dst_a3)
        rk(8)    = effsize(l_ss_a1)*1.e6_r8
        logsk(8) = log10(sgpart(l_ss_a1)) 
        rhob(8)  = rhopart(l_ss_a1) 
        rk(9)    = effsize(l_ss_a2)*1.e6_r8
        logsk(9) = log10(sgpart(l_ss_a2))
        rhob(9)  = rhopart(l_ss_a2) 
        rk(10)   = effsize(l_ss_a3)*1.e6_r8
        logsk(10)= log10(sgpart(l_ss_a3))
        rhob(10) = rhopart(l_ss_a3) 

!      Value of imax, which with dlogr=log(r(i+1))-log(r(i))  
!      gives rmin=0.001 and rmax=20 in microns: 
      imax=1+nint(4.3_r8/dlogr) ! rmax=20

      bcint=0.0_r8
      rc = 0.05_r8 
      irc=0
      do i=0,imax
        r(i)=10.0_r8**(dlogr*(i-1.0_r8)-3.0_r8)
        rp(i)=r(i)*10.0_r8**(dlogr/2.0_r8)
!        new by A. Burud:
        l10rd(i) = log10(rc/r(i))/dlogr
        l10rp(i) = log10(rp(i)/rc)
        if(irc==0.and.rp(i)>=rc) irc=i
!
        if(r(i)<=rk(2)) then
          rhorbc=rhob(2)
        else
          rhorbc=rhob(2)*(rk(2)/r(i))**0.5_r8
        endif
        bcint=bcint+dlogr*rhorbc*r(i)**3 &
!             *exp(-0.5_r8*(log10(r(i)/rk(5))/logsk(5))**2)
             *exp(-0.5_r8*(log10(r(i)/rk(0))/logsk(0))**2)
      enddo
!     write(*,*) 'irc=', irc

!     molal weights for ammonium sulphate, sulphate and
!     sulfuric acid, respectively:

      Ms=132.14_r8
      Mso4=96.06_r8
      Msv=98.08_r8

!     diffusion coefficient, termal velocity and mean free path
!     for sufuric acid, H2SO4:

      diff=1.1e7_r8                 ! um^2/s
      th=2.43e8_r8                  ! um/s                  
      mfv=3.5e-2_r8                 ! um
!      do kcomp = 1,10
      do kcomp = 0,10

!       other diffusion og coagulation parameters:
!        call koagsub(kcomp,rhob,rk) 
        if(kcomp>0) call koagsub(kcomp,rhob,rk) 
!       and finally, "normalized" size distribution for each background aerosol mode:
        do i=0,imax
          nk(kcomp,i) = sq2pi*(1.0_r8/logsk(kcomp)) &
                       *exp(-0.5_r8*(log10(r(i)/rk(kcomp))/logsk(kcomp))**2)
!cos+
           normnk(kcomp,i) =dlogr*nk(kcomp,i)
!cos-
        enddo 

      enddo  ! kcomp

!       conversion constants for use in convaer.F90 and pmxsub.F90
!
        efact_ss1 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_ss_a1)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_ss_a1))**3*rhopart(l_ss_a1))   
        efact_ss2 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_ss_a2)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_ss_a2))**3*rhopart(l_ss_a2))   
        efact_ss3 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_ss_a3)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_ss_a3))**3*rhopart(l_ss_a3))   
!NA        efact_dst1 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_dst_a1)))**2) &
!NA                   *(4.0_r8/3.0_r8)*pi*(effsize(l_dst_a1))**3*rhopart(l_dst_a1))   
        efact_dst1 = 0.0_r8 ! the dst_1 mode does not exist in the aerocomB emissions
        efact_dst2 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_dst_a2)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_dst_a2))**3*rhopart(l_dst_a2))   
        efact_dst3 = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_dst_a3)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_dst_a3))**3*rhopart(l_dst_a3))   
        efact_so4n = (Msv/Mso4)*1e-15_r8/(e**(4.5_r8*(log(sgpart(l_so4_n)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_so4_n))**3*rhopart(l_so4_n))
        efact_bcn  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_bc_n)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_bc_n))**3*rhopart(l_bc_n)) 
        efact_bcni = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_bc_ni)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_bc_ni))**3*rhopart(l_bc_ni))   
        efact_bcax = 3e3*log10(sgpart(l_bc_ax))*sq2pi/(2.0_r8*bcint)
        efact_omn  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_ni)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ni))**3*rhopart(l_om_ni))
        efact_omni = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_ni)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ni))**3*rhopart(l_om_ni))
        efact_so4na=efact_so4n
!os
!       Note The values of bc_n is used for calculating the number of bca modes.
!       This is due to the fact that bc_a aitken mode has supposedly grown from
!       bc_n and thus kept the number of bc_n particles. Similar for so4_na.

        efact_bca  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_bc_n)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_bc_n))**3*rhopart(l_bc_n))

!        efact_oma  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_a)))**2) &
!                     *(4.0_r8/3.0_r8)*pi*(effsize(l_om_a))**3*rhopart(l_om_a))

        efact_bcai  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_bc_ai)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_bc_ai))**3*rhopart(l_bc_ai))

        efact_omai  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_ai)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ai))**3*rhopart(l_om_ai))   
!os
        efact_so4pr  = 1e-15_r8/(e**(4.5_r8*(log(sgpart(l_so4_pr)))**2) &
                     *(4.0_r8/3.0_r8)*pi*(effsize(l_so4_pr))**3*rhopart(l_so4_pr))



!#ifdef PROGSSLT
! Fitted coefficients for MÃ¥rtensson sea salt emission parameterization
! defined for SS modes Rp = 0.022um (1.59), 0.130 (1.59), 0.740 um (2.0)
! Include a conversion from number flux to mass flux kg / m2 s
	ssa1 = -3.3551e-09_r8 / efact_ss1
	ssa2 = 1.1768e-10_r8 / efact_ss2
	ssa3 = -1.6675e-09_r8 / efact_ss3
	ssb1 = 1.0554e-06_r8 / efact_ss1
	ssb2 = -1.1369e-08_r8 / efact_ss2
	ssb3 = 2.2879e-07_r8 / efact_ss3
	ssc3 = 3.0608e-12_r8 / efact_ss3

!#endif

	return
      end subroutine constants




















