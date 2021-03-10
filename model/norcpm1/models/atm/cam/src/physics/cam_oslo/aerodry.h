!    For subroutine initdryp and intdrypar:

       common /dryarr1/                               &
       a5to10cintbg, a5to10cintbg05, a5to10cintbg125, &
       a5to10cintbc, a5to10cintbc05, a5to10cintbc125, & 
       a5to10cintoc, a5to10cintoc05, a5to10cintoc125, &
       a5to10cintsc, a5to10cintsc05, a5to10cintsc125, &        
       a5to10cintsa, a5to10cintsa05, a5to10cintsa125, & 
  a5to10aaeros, a5to10aaerol, a5to10vaeros, a5to10vaerol, &
       a1to3cintbg, a1to3cintbg05, a1to3cintbg125,    &
       a1to3cintsc, a1to3cintsc05, a1to3cintsc125,    &        
  a1to3aaeros, a1to3aaerol, a1to3vaeros, a1to3vaerol, &
       a4cintbg, a4cintbg05, a4cintbg125,             &
       a4cintbc, a4cintbc05, a4cintbc125,             &
       a4cintsc, a4cintsc05, a4cintsc125,             &        
       a4cintsa, a4cintsa05, a4cintsa125,             &        
  a4aaeros, a4aaerol, a4vaeros, a4vaerol,             &
      a0cintbg, a0cintbg05, a0cintbg125,             &
  a0aaeros, a0aaerol, a0vaeros, a0vaerol

  real(r8)  a5to10cintbg(6,6,6,6,5:10),    &
            a5to10cintbg05(6,6,6,6,5:10),  &
            a5to10cintbg125(6,6,6,6,5:10), &
            a5to10cintbc(6,6,6,6,5:10),    &
            a5to10cintbc05(6,6,6,6,5:10),  &
            a5to10cintbc125(6,6,6,6,5:10), &
            a5to10cintoc(6,6,6,6,5:10),    &
            a5to10cintoc05(6,6,6,6,5:10),  &
            a5to10cintoc125(6,6,6,6,5:10), &
            a5to10cintsc(6,6,6,6,5:10),    &
            a5to10cintsc05(6,6,6,6,5:10),  &
            a5to10cintsc125(6,6,6,6,5:10), &
            a5to10cintsa(6,6,6,6,5:10),    &
            a5to10cintsa05(6,6,6,6,5:10),  &
            a5to10cintsa125(6,6,6,6,5:10), &
            a5to10aaeros(6,6,6,6,5:10),    &
            a5to10aaerol(6,6,6,6,5:10),    &
            a5to10vaeros(6,6,6,6,5:10),    &
            a5to10vaerol(6,6,6,6,5:10)

  real(r8)  a1to3cintbg(16,3),    &
            a1to3cintbg05(16,3),  &
            a1to3cintbg125(16,3), &
            a1to3cintsc(16,3),    &
            a1to3cintsc05(16,3),  &
            a1to3cintsc125(16,3), &
            a1to3aaeros(16,3),    &
            a1to3aaerol(16,3),    &
            a1to3vaeros(16,3),    &
            a1to3vaerol(16,3)

!4  real(r8)  a4cintbg(16,6),    &
!4            a4cintbg05(16,6),  &
!4            a4cintbg125(16,6), &
!4            a4cintbc(16,6),    &
!4            a4cintbc05(16,6),  &
!4            a4cintbc125(16,6), &
!4            a4cintsc(16,6),    &
!4            a4cintsc05(16,6),  &
!4            a4cintsc125(16,6), &
!4            a4aaeros(16,6),    &
!4            a4aaerol(16,6),    &
!4            a4vaeros(16,6),    &
!4            a4vaerol(16,6)

  real(r8)  a4cintbg(16,6,6),    &
            a4cintbg05(16,6,6),  &
            a4cintbg125(16,6,6), &
            a4cintbc(16,6,6),    &
            a4cintbc05(16,6,6),  &
            a4cintbc125(16,6,6), &
            a4cintsc(16,6,6),    &
            a4cintsc05(16,6,6),  &
            a4cintsc125(16,6,6), &
            a4cintsa(16,6,6),    &
            a4cintsa05(16,6,6),  &
            a4cintsa125(16,6,6), &
            a4aaeros(16,6,6),    &
            a4aaerol(16,6,6),    &
            a4vaeros(16,6,6),    &
            a4vaerol(16,6,6)

  real(r8)  a0cintbg, a0cintbg05, a0cintbg125, &
            a0aaeros, a0aaerol, a0vaeros, a0vaerol

