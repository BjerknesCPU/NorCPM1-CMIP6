subroutine initaeropt

!     Tabulating the 'aerocomk'-files to save computing time.

   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer
   use opttab,   only: cate, cat, fac, faq, fbc, rh, nbmodes, nmodes 

   implicit none

#include <aerocopt.h>
#include <aerocopt2.h>

      integer kcomp, irelh, ictot, ifac, ifbc, ifaq
      integer ic, ifil, lin
      real(r8) catot, relh, frac, fabc, fraq,                  & 
               bext440, babs440, bext500, babs500, babs550,    &
               bext670, babs670, bext870, babs870,             &
               bebg440, babg440, bebg500, babg500, babg550,    &
               bebg670, babg670, bebg870, babg870,             &
               bebc440, babc440, bebc500, babc500, babc550,    &
               bebc670, babc670, bebc870, babc870,             &
               beoc440, baoc440, beoc500, baoc500, baoc550,    &
               beoc670, baoc670, beoc870, baoc870,             &
               besu440, basu440, besu500, basu500, basu550,    &
               besu670, basu670, besu870, basu870,             &
               bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
               beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
               backscat550

      real(r8) rh2(10)
      real(r8) :: eps2 = 1.e-2_r8
      real(r8) :: eps4 = 1.e-4_r8
      real(r8) :: eps6 = 1.e-6_r8
      real(r8) :: eps7 = 1.e-7_r8

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      open(10,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk1.out'  &
             ,form='formatted',status='old')
      open(11,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk2.out'  &
             ,form='formatted',status='old')
      open(12,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk3.out'  &
             ,form='formatted',status='old')
      open(13,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk4.out' &
             ,form='formatted',status='old')
      open(14,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk5.out'  &
             ,form='formatted',status='old')
      open(15,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk6.out'  &
             ,form='formatted',status='old')
      open(16,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk7.out'  &
             ,form='formatted',status='old')
      open(17,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk8.out'  &
             ,form='formatted',status='old')
      open(18,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk9.out'  &
             ,form='formatted',status='old')
      open(19,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk10.out' &
             ,form='formatted',status='old')
      open(20,file='inputdata/atm/cam/camoslo/aerocomapr09/aerocomk0.out' &
             ,form='formatted',status='old')


!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do ifil = 5,10
        do lin = 1,12960     ! 10x6x6x6x6

          read(9+ifil,993) kcomp, relh, catot, frac, fabc, fraq, &
!               bext440, babs440, bext500, babs500, babs550,    &
!               bext670, babs670, bext870, babs870,             &
!               bebg440, babg440, bebg500, babg500,             &
!               bebg670, babg670, bebg870, babg870,             &
!               bebc440, babc440, bebc500, babc500,             &
!               bebc670, babc670, bebc870, babc870,             &
!               beoc440, baoc440, beoc500, baoc500,             &
!               beoc670, baoc670, beoc870, baoc870,             &
!               besu440, basu440, besu500, basu500,             &
!               besu670, basu670, besu870, basu870,             &
!               bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
!               beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
!               backscat550, babg550, babc550, baoc550, basu550
               bext440, bext500, bext670, bext870,             &
               bebg440, bebg500, bebg670, bebg870,             &
               bebc440, bebc500, bebc670, bebc870,             &
               beoc440, beoc500, beoc670, beoc870,             &
               besu440, besu500, besu670, besu870,             &
               babs440, babs500, babs550, babs670, babs870,    &
               bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
               beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
               backscat550, babg550, babc550, baoc550, basu550

!pga feil i antall siffer i tabellene:
!old          if(abs(relh-0.99_r8)<eps6) then          
!old            relh=0.995_r8
!old          endif    

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 11
	   endif
	  end do
   11     continue

 	  do ic=1,6
	   if(abs(catot-cat(kcomp,ic))<eps6) then
	    ictot=ic
	    goto 21
	   endif
	  end do
   21     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 31
	   endif
	  end do
   31     continue

 	  do ic=1,6
	   if(abs(fabc-fbc(ic))<eps4) then
	    ifbc=ic
	    goto 41
	   endif
	  end do
   41     continue

	  do ic=1,6
	   if(abs(fraq-faq(ic))<eps4) then
	    ifaq=ic
	    goto 51
	   endif
	  end do
   51     continue

           bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bext440 ! unit km^-1
           bep5to10(2,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bext500
           bep5to10(3,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bext670 
           bep5to10(4,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bext870
           bep5to10(5,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg440
           bep5to10(6,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg500
           bep5to10(7,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg670
           bep5to10(8,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg870
           bep5to10(9,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc440
           bep5to10(10,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc500
           bep5to10(11,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc670
           bep5to10(12,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc870
           bep5to10(13,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc440
           bep5to10(14,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc500
           bep5to10(15,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc670
           bep5to10(16,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc870
           bep5to10(17,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu440
           bep5to10(18,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu500
           bep5to10(19,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu670
           bep5to10(20,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu870
           bep5to10(21,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babs440
           bep5to10(22,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babs500
           bep5to10(23,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babs550
           bep5to10(24,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babs670
           bep5to10(25,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babs870
           bep5to10(26,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg550lt1
           bep5to10(27,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebg550gt1
           bep5to10(28,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc550lt1
           bep5to10(29,irelh,ictot,ifac,ifbc,ifaq,kcomp)=bebc550gt1
           bep5to10(30,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc550lt1
           bep5to10(31,irelh,ictot,ifac,ifbc,ifaq,kcomp)=beoc550gt1
           bep5to10(32,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu550lt1
           bep5to10(33,irelh,ictot,ifac,ifbc,ifaq,kcomp)=besu550gt1
           bep5to10(34,irelh,ictot,ifac,ifbc,ifaq,kcomp)=backscat550
           bep5to10(35,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babg550
           bep5to10(36,irelh,ictot,ifac,ifbc,ifaq,kcomp)=babc550
           bep5to10(37,irelh,ictot,ifac,ifbc,ifaq,kcomp)=baoc550
           bep5to10(38,irelh,ictot,ifac,ifbc,ifaq,kcomp)=basu550

        end do
      end do

    do kcomp=5,10
    do irelh=1,10
    do ictot=1,6
    do ifac=1,6
    do ifaq=1,6
     if(bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
      write(*,*) 'bep5to10 =', kcomp, irelh, ictot, ifac, ifbc, ifaq, &
                               bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)
      write(*,*) 'Error in initialization of bep5to10'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo

        write(*,*)'mode 5-10 ok'


!       Modes 1 to 3 (SO4/BC/OC + condesate from H2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do ifil = 1,3
        do lin = 1,160     ! 10x16

          read(9+ifil,994) kcomp, relh, catot, &
!           bext440, babs440, bext500, babs500, babs550,    &
!           bext670, babs670, bext870, babs870,             &
!           bebg440, babg440, bebg500, babg500,             &
!           bebg670, babg670, bebg870, babg870,             &
!           bebc440, babc440, bebc500, babc500,             &
!           bebc670, babc670, bebc870, babc870,             &
!           beoc440, baoc440, beoc500, baoc500,             &
!           beoc670, baoc670, beoc870, baoc870,             &
!           besu440, basu440, besu500, basu500,             &
!           besu670, basu670, besu870, basu870,             &
!           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
!           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
!           backscat550, babg550, babc550, baoc550, basu550
           bext440, bext500, bext670, bext870,             &
           bebg440, bebg500, bebg670, bebg870,             &
           bebc440, bebc500, bebc670, bebc870,             &
           beoc440, beoc500, beoc670, beoc870,             &
           besu440, besu500, besu670, besu870,             &
           babs440, babs500, babs550, babs670, babs870,    &
           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
           backscat550, babg550, babc550, baoc550, basu550


!pga feil i antall siffer i tabellene:
!old          if(abs(relh-0.99_r8)<eps6) then          
!old            relh=0.995_r8
!old          endif    

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 61
	   endif
	  end do
   61     continue

 	  do ic=1,16
	   if(abs(catot-cate(kcomp,ic))<eps7) then
	    ictot=ic
	    goto 71
	   endif
	  end do
   71     continue

           bep1to3(1,irelh,ictot,kcomp)=bext440 ! unit km^-1
           bep1to3(2,irelh,ictot,kcomp)=bext500
           bep1to3(3,irelh,ictot,kcomp)=bext670 
           bep1to3(4,irelh,ictot,kcomp)=bext870
           bep1to3(5,irelh,ictot,kcomp)=bebg440
           bep1to3(6,irelh,ictot,kcomp)=bebg500
           bep1to3(7,irelh,ictot,kcomp)=bebg670
           bep1to3(8,irelh,ictot,kcomp)=bebg870
           bep1to3(9,irelh,ictot,kcomp)=bebc440  !=0
           bep1to3(10,irelh,ictot,kcomp)=bebc500 !=0
           bep1to3(11,irelh,ictot,kcomp)=bebc670 !=0
           bep1to3(12,irelh,ictot,kcomp)=bebc870 !=0
           bep1to3(13,irelh,ictot,kcomp)=beoc440 !=0
           bep1to3(14,irelh,ictot,kcomp)=beoc500 !=0
           bep1to3(15,irelh,ictot,kcomp)=beoc670 !=0
           bep1to3(16,irelh,ictot,kcomp)=beoc870 !=0
           bep1to3(17,irelh,ictot,kcomp)=besu440
           bep1to3(18,irelh,ictot,kcomp)=besu500
           bep1to3(19,irelh,ictot,kcomp)=besu670
           bep1to3(20,irelh,ictot,kcomp)=besu870
           bep1to3(21,irelh,ictot,kcomp)=babs440
           bep1to3(22,irelh,ictot,kcomp)=babs500
           bep1to3(23,irelh,ictot,kcomp)=babs550
           bep1to3(24,irelh,ictot,kcomp)=babs670
           bep1to3(25,irelh,ictot,kcomp)=babs870
           bep1to3(26,irelh,ictot,kcomp)=bebg550lt1
           bep1to3(27,irelh,ictot,kcomp)=bebg550gt1
           bep1to3(28,irelh,ictot,kcomp)=bebc550lt1 !=0
           bep1to3(29,irelh,ictot,kcomp)=bebc550gt1 !=0
           bep1to3(30,irelh,ictot,kcomp)=beoc550lt1 !=0
           bep1to3(31,irelh,ictot,kcomp)=beoc550gt1 !=0
           bep1to3(32,irelh,ictot,kcomp)=besu550lt1
           bep1to3(33,irelh,ictot,kcomp)=besu550gt1
           bep1to3(34,irelh,ictot,kcomp)=backscat550
           bep1to3(35,irelh,ictot,kcomp)=babg550
           bep1to3(36,irelh,ictot,kcomp)=babc550 !=0
           bep1to3(37,irelh,ictot,kcomp)=baoc550 !=0
           bep1to3(38,irelh,ictot,kcomp)=basu550

        end do
      end do

    do kcomp=1,3
    do irelh=1,10
    do ictot=1,16
     if(bep1to3(1,irelh,ictot,kcomp)<=0.0_r8) then
      write(*,*) 'bep1to3 =', irelh, ictot, bep1to3(1,irelh,ictot,kcomp)
      write(*,*) 'Error in initialization of bep1to3'
      stop
     endif
    enddo
    enddo
    enddo

        write(*,*)'mode 1-3 ok' 


!       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 4
!4        do lin = 1,960     ! 10x16x6
        do lin = 1,5760     ! 10x16x6x6

!4          read(9+ifil,995) kcomp, relh, catot, frac, &
          read(9+ifil,995) kcomp, relh, catot, frac, fraq, &
!           bext440, babs440, bext500, babs500, babs550,    &
!           bext670, babs670, bext870, babs870,             &
!           bebg440, babg440, bebg500, babg500,             &
!           bebg670, babg670, bebg870, babg870,             &
!           bebc440, babc440, bebc500, babc500,             &
!           bebc670, babc670, bebc870, babc870,             &
!           beoc440, baoc440, beoc500, baoc500,             &
!           beoc670, baoc670, beoc870, baoc870,             &
!           besu440, basu440, besu500, basu500,             &
!           besu670, basu670, besu870, basu870,             &
!           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
!           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
!           backscat550, babg550, babc550, baoc550, basu550
           bext440, bext500, bext670, bext870,             &
           bebg440, bebg500, bebg670, bebg870,             &
           bebc440, bebc500, bebc670, bebc870,             &
           beoc440, beoc500, beoc670, beoc870,             &
           besu440, besu500, besu670, besu870,             &
           babs440, babs500, babs550, babs670, babs870,    &
           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
           backscat550, babg550, babc550, baoc550, basu550


!pga feil i antall siffer i tabellene:
!old          if(abs(relh-0.99_r8)<eps6) then          
!old            relh=0.995_r8
!old          endif    

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 81
	   endif
	  end do
   81     continue

 	  do ic=1,16
!skummelt	   if(abs(catot-cate(kcomp,ic))<eps7) then
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 91
	   endif
	  end do
   91     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 101
	   endif
	  end do
  101     continue

	  do ic=1,6
	   if(abs(fraq-faq(ic))<eps4) then
	    ifaq=ic
	    goto 111
	   endif
	  end do
  111     continue

           bep4(1,irelh,ictot,ifac,ifaq)=bext440 ! unit km^-1
           bep4(2,irelh,ictot,ifac,ifaq)=bext500
           bep4(3,irelh,ictot,ifac,ifaq)=bext670 
           bep4(4,irelh,ictot,ifac,ifaq)=bext870
           bep4(5,irelh,ictot,ifac,ifaq)=bebg440
           bep4(6,irelh,ictot,ifac,ifaq)=bebg500
           bep4(7,irelh,ictot,ifac,ifaq)=bebg670
           bep4(8,irelh,ictot,ifac,ifaq)=bebg870
           bep4(9,irelh,ictot,ifac,ifaq)=bebc440
           bep4(10,irelh,ictot,ifac,ifaq)=bebc500
           bep4(11,irelh,ictot,ifac,ifaq)=bebc670
           bep4(12,irelh,ictot,ifac,ifaq)=bebc870
           bep4(13,irelh,ictot,ifac,ifaq)=beoc440 !=0
           bep4(14,irelh,ictot,ifac,ifaq)=beoc500 !=0
           bep4(15,irelh,ictot,ifac,ifaq)=beoc670 !=0
           bep4(16,irelh,ictot,ifac,ifaq)=beoc870 !=0
           bep4(17,irelh,ictot,ifac,ifaq)=besu440
           bep4(18,irelh,ictot,ifac,ifaq)=besu500
           bep4(19,irelh,ictot,ifac,ifaq)=besu670
           bep4(20,irelh,ictot,ifac,ifaq)=besu870
           bep4(21,irelh,ictot,ifac,ifaq)=babs440
           bep4(22,irelh,ictot,ifac,ifaq)=babs500
           bep4(23,irelh,ictot,ifac,ifaq)=babs550
           bep4(24,irelh,ictot,ifac,ifaq)=babs670
           bep4(25,irelh,ictot,ifac,ifaq)=babs870
           bep4(26,irelh,ictot,ifac,ifaq)=bebg550lt1
           bep4(27,irelh,ictot,ifac,ifaq)=bebg550gt1
           bep4(28,irelh,ictot,ifac,ifaq)=bebc550lt1
           bep4(29,irelh,ictot,ifac,ifaq)=bebc550gt1
           bep4(30,irelh,ictot,ifac,ifaq)=beoc550lt1 !=0
           bep4(31,irelh,ictot,ifac,ifaq)=beoc550gt1 !=0
           bep4(32,irelh,ictot,ifac,ifaq)=besu550lt1
           bep4(33,irelh,ictot,ifac,ifaq)=besu550gt1
           bep4(34,irelh,ictot,ifac,ifaq)=backscat550
           bep4(35,irelh,ictot,ifac,ifaq)=babg550
           bep4(36,irelh,ictot,ifac,ifaq)=babc550 !=0
           bep4(37,irelh,ictot,ifac,ifaq)=baoc550 !=0
           bep4(38,irelh,ictot,ifac,ifaq)=basu550

        end do

    do irelh=1,10
    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(bep4(1,irelh,ictot,ifac,ifaq)<=0.0_r8) then
      write(*,*) 'bep4 =', irelh, ictot, ifac, ifaq, bep4(1,irelh,ictot,ifac,ifaq)
      write(*,*) 'Error in initialization of bep4'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

        write(*,*)'mode 4 ok'


!       Mode 0, BC(ax)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

          ifil = 11

          read(9+ifil,996) kcomp, relh,  &
           bex440, bax440, bex500, bax500, bax550, bex670, bax670, &
           bex870, bax870, bex550lt1, bex550gt1, backscx550

     if(bex440<=0.0_r8) then
      write(*,*) 'bex440 =', bex440
      write(*,*) 'Error in initialization of bex1'
      stop
     endif

        write(*,*)'mode 0 ok'

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!  993 format(I2,f6.3,3e10.3,f5.2,54e10.3)
!  994 format(I2,f6.3,e10.3,54e10.3)
!  995 format(I2,f6.3,3e10.3,54e10.3)
  993 format(I2,f6.3,3e10.3,f5.2,38e10.3)
  994 format(I2,f6.3,e10.3,38e10.3)
  995 format(I2,f6.3,3e10.3,38e10.3)
  996 format(I2,f6.3,12e11.4)

      do ifil=10,20
        close (ifil)
      end do 

      return
end subroutine initaeropt

