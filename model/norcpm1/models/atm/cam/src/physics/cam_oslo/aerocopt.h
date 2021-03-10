!    For subroutines initaeropt and intaeropt1to3,4,6to10:

!old common /aerocopt1/ bep1to3, bep4, bep0, bep5to10
!old       real(r8) bep1to3(28,10,16,3)
!old       real(r8) bep4(28,10,16,6,6)
!old       real(r8) bep0(28)
!old       real(r8) bep5to10(28,10,6,6,6,6,5:10)

       common /aerocopt1/ bep1to3, bep4, bep5to10

!       real(r8) bep1to3(54,10,16,3)
!       real(r8) bep4(54,10,16,6,6)
!       real(r8) bep5to10(54,10,6,6,6,6,5:10)
       real(r8) bep1to3(38,10,16,3)
       real(r8) bep4(38,10,16,6,6)
       real(r8) bep5to10(38,10,6,6,6,6,5:10)
