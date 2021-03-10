subroutine initdryp

   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer
   use opttab,   only: cate, cat, fac, faq, fbc, rh, nbmodes, nmodes 

   implicit none

!     Tabulating the 'aerodryk'-files to save computing time. Routine
!     originally made by  Alf Kirkevaag, and modified for new aerosol 
!     schemes in January 2006.

#include <aerodry.h>

      integer kcomp, irelh, ictot, ifac, ifbc, ifaq
      integer ic, ifil, lin
      real(r8) catot, relh, frac, fabc, fraq, &
          cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,   &  
          cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,   &      
          cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol 
      real(r8) :: eps2 = 1.e-2_r8
      real(r8) :: eps4 = 1.e-4_r8
      real(r8) :: eps6 = 1.e-6_r8
      real(r8) :: eps7 = 1.e-7_r8

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      open(10,file='inputdata/atm/cam/camoslo/aertab/aerodryk1.out' &
             ,form='formatted',status='old')
      open(11,file='inputdata/atm/cam/camoslo/aertab/aerodryk2.out' &
             ,form='formatted',status='old')
      open(12,file='inputdata/atm/cam/camoslo/aertab/aerodryk3.out' &
            ,form='formatted',status='old')
!2oct07      open(13,file='inputdata/atm/cam/camoslo/aertab/aerodryk4.out' &
      open(13,file='inputdata/atm/cam/camoslo/aertab/aerodryk4.out' &
             ,form='formatted',status='old')
      open(14,file='inputdata/atm/cam/camoslo/aertab/aerodryk5.out' &
             ,form='formatted',status='old')
      open(15,file='inputdata/atm/cam/camoslo/aertab/aerodryk6.out' &
             ,form='formatted',status='old')
      open(16,file='inputdata/atm/cam/camoslo/aertab/aerodryk7.out' &
             ,form='formatted',status='old')
      open(17,file='inputdata/atm/cam/camoslo/aertab/aerodryk8.out' &
             ,form='formatted',status='old')
      open(18,file='inputdata/atm/cam/camoslo/aertab/aerodryk9.out' &
             ,form='formatted',status='old')
      open(19,file='inputdata/atm/cam/camoslo/aertab/aerodryk10.out' & 
             ,form='formatted',status='old')
      open(20,file='inputdata/atm/cam/camoslo/aertab/aerodryk0.out' & 
             ,form='formatted',status='old')


!       Modes 5 to 10 (mineral and seasalt-modes + cond./coag./aq.)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do ifil = 5,10
        do lin = 1,1296     ! 6x6x6x6

          read(9+ifil,993) kcomp, catot, frac, fabc, fraq,            &
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &      
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

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

           a5to10cintbg(ictot,ifac,ifbc,ifaq,kcomp)=cintbg
           a5to10cintbg05(ictot,ifac,ifbc,ifaq,kcomp)=cintbg05
           a5to10cintbg125(ictot,ifac,ifbc,ifaq,kcomp)=cintbg125

           a5to10cintbc(ictot,ifac,ifbc,ifaq,kcomp)=cintbc 
           a5to10cintbc05(ictot,ifac,ifbc,ifaq,kcomp)=cintbc05
           a5to10cintbc125(ictot,ifac,ifbc,ifaq,kcomp)=cintbc125

           a5to10cintoc(ictot,ifac,ifbc,ifaq,kcomp)=cintoc
           a5to10cintoc05(ictot,ifac,ifbc,ifaq,kcomp)=cintoc05
           a5to10cintoc125(ictot,ifac,ifbc,ifaq,kcomp)=cintoc125

           a5to10cintsc(ictot,ifac,ifbc,ifaq,kcomp)=cintsc
           a5to10cintsc05(ictot,ifac,ifbc,ifaq,kcomp)=cintsc05
           a5to10cintsc125(ictot,ifac,ifbc,ifaq,kcomp)=cintsc125

           a5to10cintsa(ictot,ifac,ifbc,ifaq,kcomp)=cintsa
           a5to10cintsa05(ictot,ifac,ifbc,ifaq,kcomp)=cintsa05
           a5to10cintsa125(ictot,ifac,ifbc,ifaq,kcomp)=cintsa125

           a5to10aaeros(ictot,ifac,ifbc,ifaq,kcomp)=aaeros
           a5to10aaerol(ictot,ifac,ifbc,ifaq,kcomp)=aaerol
           a5to10vaeros(ictot,ifac,ifbc,ifaq,kcomp)=vaeros
           a5to10vaerol(ictot,ifac,ifbc,ifaq,kcomp)=vaerol

        end do  ! lin
      end do    ! ifil

    do kcomp=5,10
    do ictot=1,6
    do ifac=1,6
    do ifaq=1,6
     if(a5to10cintbg(ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
      write(*,*) 'a5to10cintbg =', kcomp, ictot, ifac, ifbc, ifaq, a5to10cintbg(ictot,ifac,ifbc,ifaq,kcomp)
      write(*,*) 'Error in initialization of a5to10cintbg'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

!        write(*,*)'mode 5-10 ok'


!       Modes 1 to 3 (SO4/BC/OC + condesate from H2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do ifil = 1,3
        do lin = 1,16     ! 16

          read(9+ifil,994) kcomp, catot,                              &
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &      
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

 	  do ic=1,16
	   if(abs(catot-cate(kcomp,ic))<eps7) then
	    ictot=ic
	    goto 61
	   endif
	  end do
   61     continue

!         no ifac-, ifbc- or ifaq-dependency for these modes,
!         since all catot comes from condensed sulfate

           a1to3cintbg(ictot,kcomp)=cintbg
           a1to3cintbg05(ictot,kcomp)=cintbg05
           a1to3cintbg125(ictot,kcomp)=cintbg125

           a1to3cintsc(ictot,kcomp)=cintsc
           a1to3cintsc05(ictot,kcomp)=cintsc05
           a1to3cintsc125(ictot,kcomp)=cintsc125

           a1to3aaeros(ictot,kcomp)=aaeros
           a1to3aaerol(ictot,kcomp)=aaerol
           a1to3vaeros(ictot,kcomp)=vaeros
           a1to3vaerol(ictot,kcomp)=vaerol

        end do  ! lin
      end do    ! ifil

    do kcomp=1,3
    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(a1to3cintbg(ictot,kcomp)<=0.0_r8) then
      write(*,*) 'a1to3cintbg =', kcomp, ictot, a1to3cintbg(ictot,kcomp)
      write(*,*) 'Error in initialization of a1to3cintbg'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

!        write(*,*)'mode 1-3 ok'


!       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 4
!4        do lin = 1,96     ! 16x6
        do lin = 1,576     ! 16x6x6

!4          read(9+ifil,995) kcomp, catot, frac,                        &
          read(9+ifil,995) kcomp, catot, frac, fraq,                  &
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &      
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

 	  do ic=1,16
!skummelt	   if(abs(catot-cate(kcomp,ic))<eps7) then
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 71
	   endif
	  end do
   71     continue

 	  do ic=1,6
	   if(abs(frac-fac(ic))<eps4) then
	    ifac=ic
	    goto 81
	   endif
	  end do
   81     continue

	  do ic=1,6
	   if(abs(fraq-faq(ic))<eps4) then
	    ifaq=ic
	    goto 91
	   endif
	  end do
   91     continue

!         no ifbc-dependency for this mode, since all catot 
!         comes from condensed or wet-phase sulfate 

!4           a4cintbg(ictot,ifac)=cintbg
!4           a4cintbg05(ictot,ifac)=cintbg05
!4           a4cintbg125(ictot,ifac)=cintbg125
!4           a4cintbc(ictot,ifac)=cintbc 
!4           a4cintbc05(ictot,ifac)=cintbc05
!4           a4cintbc125(ictot,ifac)=cintbc125
!4           a4cintsc(ictot,ifac)=cintsc
!4           a4cintsc05(ictot,ifac)=cintsc05
!4           a4cintsc125(ictot,ifac)=cintsc125
!4           a4aaeros(ictot,ifac)=aaeros
!4           a4aaerol(ictot,ifac)=aaerol
!4           a4vaeros(ictot,ifac)=vaeros
!4           a4vaerol(ictot,ifac)=vaerol
!4           a4aaeros(ictot,ifac)=aaeros
!4           a4aaerol(ictot,ifac)=aaerol
!4           a4vaeros(ictot,ifac)=vaeros
!4           a4vaerol(ictot,ifac)=vaerol

           a4cintbg(ictot,ifac,ifaq)=cintbg
           a4cintbg05(ictot,ifac,ifaq)=cintbg05
           a4cintbg125(ictot,ifac,ifaq)=cintbg125

           a4cintbc(ictot,ifac,ifaq)=cintbc 
           a4cintbc05(ictot,ifac,ifaq)=cintbc05
           a4cintbc125(ictot,ifac,ifaq)=cintbc125

           a4cintsc(ictot,ifac,ifaq)=cintsc
           a4cintsc05(ictot,ifac,ifaq)=cintsc05
           a4cintsc125(ictot,ifac,ifaq)=cintsc125

           a4cintsa(ictot,ifac,ifaq)=cintsa
           a4cintsa05(ictot,ifac,ifaq)=cintsa05
           a4cintsa125(ictot,ifac,ifaq)=cintsa125

     if(cintsa<cintsa05) &
       write(*,*) 'cintsatot =', ictot, ifac, ifaq, cintsa, cintsa05, cintsa125

           a4aaeros(ictot,ifac,ifaq)=aaeros
           a4aaerol(ictot,ifac,ifaq)=aaerol
           a4vaeros(ictot,ifac,ifaq)=vaeros
           a4vaerol(ictot,ifac,ifaq)=vaerol

        end do  ! lin

    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(a4cintbg(ictot,ifac,ifaq)<=0.0_r8) then
      write(*,*) 'a4cintbg =', ictot, ifac, ifaq, a4cintbg(ictot,ifac,ifaq)
      write(*,*) 'Error in initialization of a4cintbg'
      stop
     endif
    enddo
    enddo
    enddo

!        write(*,*)'mode 4 ok'


!       Mode 0, BC(ax)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 11
       
          read(9+ifil,996) kcomp, cintbg, cintbg05, cintbg125, &
           aaeros, aaerol, vaeros, vaerol

!         no ictot-, ifac-, ifbc- or ifaq-dependency for this mode,
!         since BC(ax) is purely externally mixed 

           a0cintbg=cintbg
           a0cintbg05=cintbg05
           a0cintbg125=cintbg125

           a0aaeros=aaeros
           a0aaerol=aaerol
           a0vaeros=vaeros
           a0vaerol=vaerol

!        write(*,*)'mode 0 ok'

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc


  993 format(I2,3e10.3,f5.2,19e10.3)
  994 format(I2,e10.3,19e10.3)
!4  995 format(I2,2e10.3,19e10.3)
  995 format(I2,3e10.3,19e10.3)
  996 format(I2,7e11.4)

      do ifil=10,20
        close (ifil)
      end do 

      return
end subroutine initdryp

