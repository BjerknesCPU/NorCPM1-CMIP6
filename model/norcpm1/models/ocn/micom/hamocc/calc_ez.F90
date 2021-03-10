      SUBROUTINE CALC_EZ(kpie,kpje,kpke,pddpo,ptiestw)
!
!**********************************************************************
!
!**** *CALC_EZ* - .
!
!     Karen Assmann          *BCCR*           02.12.05
! 
!     find lowest mass containing layer in the euphotic zone
!

      USE mo_biomod
      USE mo_control_bgc

      implicit none

      REAL :: pddpo(kpie,kpje,kpke),ptiestw(kpie,kpje,kpke+1)
      INTEGER :: kpie,kpje,kpke,i,j,k

! ******************************************************************

!$OMP PARALLEL DO
      DO 1321 j=1,kpje
      DO 1321 i=1,kpie
         kwrbioz(i,j)=1
      DO 1321 k=2,kpke
         IF(pddpo(i,j,k).GT.1.e-20.and.ptiestw(i,j,k).lt.99.) THEN
!         IF(pddpo(i,j,k).GT.1.e-6.and.ptiestw(i,j,k+1).le.100.) THEN
            kwrbioz(i,j)=k
         ENDIF
1321  CONTINUE
!$OMP END PARALLEL DO

      RETURN
      END
