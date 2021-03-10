module const

!-----------------------------------------------------------------------------
!Module containing subroutines constants, koagsub and parmix and declaration
!of the variables required by them. Updated with one internally mixed and one
!externally mixed OC mode November 2004.
!-----------------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
!
implicit none
!
public
save
!ts++ From the old common block const.h
	real(r8) Ms
	real(r8) Mso4
	real(r8) Msv
	real(r8) bcint
	real(r8) r(0:100)
	real(r8) rp(0:100)
	real(r8) diff
	real(r8) th
	real(r8) mfv
	real(r8) Dmp(0:100)
        real(r8) Dmp_dst(0:100)
        real(r8) Dmp_bcn(0:100)
        real(r8) Dmp_omn(0:100)
        real(r8) Dmp_ni(0:100)
	real(r8) Dmpa(10,0:100)
	real(r8) Kp12(10,0:101)
	real(r8) Kp12oc(10,0:101)
	real(r8) Kp12s4(10,0:101)
	real(r8) Kp12ax(10,0:101)
	real(r8) nk(0:10,0:100)
        real(r8) rc
        integer irc 
!ts--
!ts++ from the old common block mc.h
      real(r8) cldmin
      real(r8) dlogr
 !     real, parameter :: pi=3.141592654_r8
!ts--
!ts++ From the old common block ccntab.h
        real(r8) :: ssup(9)
        real(r8) :: sfccn(9)

        real(r8) :: sup  (10,11664)

        real(r8) :: fff1to3 (3,9,16)
        real(r8) :: rrr1to3 (3,16)	!TS: Modal radius array, mode 1 - 3	 
        real(r8) :: sss1to3 (3,16)	!TS: Standard deviation array, Mode 1 -3
        real(r8) :: catot1to3(3,144)
        real(r8) :: calog1to3(3,16)	!TS: Array for reading catot from file
        real(r8) :: fccn1to3 (3,144)
	real(r8) :: rk1to3 (3,16)	!TS: Array for reading modal radius from file
	real(r8) :: stdv1to3 (3,16)	!TS: Array for reading std. dev. from file

!4        real(r8) :: fff4 (9,16,6)
!4        real(r8) :: catot4(864)
!4        real(r8) :: fccn4 (864)
!4        real(r8) :: frac4 (864)
        real(r8) :: fff4 (9,16,6,6)
        real(r8) :: rrr4 (16,6,6)	!TS: Modal radius array, mode 4	 	
        real(r8) :: sss4 (16,6,6)	!TS: Modal radius array, mode 4	 
        real(r8) :: catot4(5184)
        real(r8) :: calog4(576)	!TS: Same as catot4, but for initlogn.F90
        real(r8) :: fccn4 (5184)
        real(r8) :: rk4 (576)		!TS: Array for reading modal radius from file
        real(r8) :: stdv4 (576)	!TS: Array for reading std. dev. from file
        real(r8) :: frac4 (5184)
        real(r8) :: fraq4 (5184)
        real(r8) :: fraclog4 (576)	!TS: Same as frac4, but for initlogn.F90 
        real(r8) :: fraqlog4 (576)	!TS: Same as fraq4, but for initlogn.F90

        real(r8) :: fff (5:10,9,6,6,6,6)
        real(r8) :: sss (5:10,6,6,6,6)	!TS: Standard deviation array, mode 5 - 10
        real(r8) :: rrr (5:10,6,6,6,6)	!TS: Modal radius array, mode 5 - 10
        real(r8) :: catot (5:10,11664)
        real(r8) :: calog (5:10,1296)  !TS: Same as catot, but for initlogn.F90
        real(r8) :: fccn5to10 (5:10,11664)
        real(r8) :: rk5to10 (5:10,1296) !TS: Array for reading modal radius from file
        real(r8) :: stdv5to10 (5:10,1296) !TS: Array for reading std. dev. from file
        real(r8) :: frac5to10 (5:10,11664)
        real(r8) :: fraclog5to10 (5:10,1296) !TS: Same as frac5to10, but for initlogn.F90 
        real(r8) :: fabc5to10 (5:10,11664)
        real(r8) :: fabclog5to10 (5:10,1296) !TS: Same as fabc5to10, but for initlogn.F90 
        real(r8) :: fraq5to10 (5:10,11664)
        real(r8) :: fraqlog5to10 (5:10,1296) !TS: Same as fraq5to10, but for initlogn.F90 

        real(r8) :: ffs(9)
        real(r8) :: ffoc(9)
!ts--

	integer imax

        real(r8) sq2pi
        real(r8) efact_so4n, efact_so4na,efact_bcn, efact_bcax, efact_omn
        real(r8) efact_bcni,efact_omni
        real(r8) efact_so4pr
        real(r8) efact_ss1, efact_ss2, efact_ss3
        real(r8) efact_dst1, efact_dst2, efact_dst3 
        real(r8) efact_bca,efact_oma,efact_bcai,efact_omai

!      new by A. Burud :
        real(r8) l10rd(0:100)
        real(r8) l10rp(0:100)
        real(r8) normnk(0:10,0:100)
        real(r8) dmpxr(0:100)

! Fitted coefficients for Mårtensson sea salt emission parameterization
! defined for SS modes Rp = 0.022um (1.59), 0.130 (1.59), 0.740 um (2.0)
	real(r8) :: ssa1, ssa2, ssa3
	real(r8) :: ssb1, ssb2, ssb3
        real(r8) :: ssc3



end module const




