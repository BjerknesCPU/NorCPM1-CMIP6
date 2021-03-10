subroutine initlogn

!  Created for CAM3 by Trude Storelvmo, Fall 2007.
!  This subroutine reads the tabulated parameters for "best lognormal fits"
!  of the aerosol size distribution wrt CCN activation as calculated by Alf Kirkevaag.
!

   use shr_kind_mod, only: r8 => shr_kind_r8
   use aerosoldef
   use opttab, only: cat,fac,fbc,faq,cate
   use const

   implicit none

   integer kcomp, ictot, ifac, ifbc, ifaq, irk, istdv
   integer ic, ifil, lin

   open(20,file='inputdata/atm/cam/camoslo/aertab/logntilp1.out' &   ! SO4(n/Ait)
          ,form='formatted',status='old')
   open(21,file='inputdata/atm/cam/camoslo/aertab/logntilp2.out' &   ! BC(n/Ait)
          ,form='formatted',status='old')
   open(22,file='inputdata/atm/cam/camoslo/aertab/logntilp3.out' &   ! OC(n/Ait)
          ,form='formatted',status='old')
   open(13,file='inputdata/atm/cam/camoslo/aertab/logntilp4.out' &   ! BC&OC(n/Ait)
          ,form='formatted',status='old')
   open(14,file='inputdata/atm/cam/camoslo/aertab/logntilp5.out' &   ! SO4(Ait75)
          ,form='formatted',status='old')
   open(15,file='inputdata/atm/cam/camoslo/aertab/logntilp6.out' &   ! MINACC
          ,form='formatted',status='old')
   open(16,file='inputdata/atm/cam/camoslo/aertab/logntilp7.out' &   ! MINCOA
          ,form='formatted',status='old')
   open(17,file='inputdata/atm/cam/camoslo/aertab/logntilp8.out' &   ! SEASF
          ,form='formatted',status='old')
   open(18,file='inputdata/atm/cam/camoslo/aertab/logntilp9.out' &   ! SEASACC
          ,form='formatted',status='old')
   open(19,file='inputdata/atm/cam/camoslo/aertab/logntilp10.out' &  ! SEASCOA
          ,form='formatted',status='old')


!       Modes 5 to 10 (SO4(ait75) and mineral and seasalt-modes + cond./coag./aq.)
!  ************************************************************************

         do ifil = 5,10
          do lin = 1,1296   ! 6**4 entries
           read(9+ifil,995) calog(ifil,lin) &
             ,fraclog5to10(ifil,lin), fabclog5to10(ifil,lin), fraqlog5to10(ifil,lin) &
             ,rk5to10(ifil,lin), stdv5to10(ifil,lin), kcomp

	do ic=1,6
	 if(calog(ifil,lin).eq.cat(kcomp,ic)) then
	  ictot=ic
	  goto 21
	 endif
	end do
   21 continue

	do ic=1,6
	 if(fraclog5to10(ifil,lin).eq.fac(ic)) then
	  ifac=ic
	  goto 31
	 endif
	end do
        PRINT*,'ifac not found',fraclog5to10(ifil,lin)
   31 continue

	do ic=1,6
	 if(fabclog5to10(ifil,lin).eq.fbc(ic)) then
	  ifbc=ic
	  goto 41
	 endif
	end do
        PRINT*,'ifbc not found',fabclog5to10(ifil,lin)
   41 continue

	do ic=1,6
	 if(fraqlog5to10(ifil,lin).eq.faq(ic)) then
	  ifaq=ic
	  goto 51
	 endif
	end do
        PRINT*,'ifaq not found',fraqlog5to10(ifil,lin)
   51 continue

        rrr(kcomp,ictot,ifac,ifbc,ifaq) = rk5to10(ifil,lin)
	sss(kcomp,ictot,ifac,ifbc,ifaq) = stdv5to10(ifil,lin)

          end do   ! lin
         end do    ! ifil

    do kcomp=5,10
    do ifac=1,6
    do ifbc=1,6
    do ictot=1,6
     if(rrr(kcomp,ictot,ifac,ifbc,ifaq)<=0.0_r8) then
      write(*,*) 'rrr =',kcomp,ictot,ifac,ifbc,ifaq,rrr(kcomp,ictot,ifac,ifbc,ifaq)
      write(*,*) 'Error in initialization of rrr'
      stop
     endif
     if(sss(kcomp,ictot,ifac,ifbc,ifaq)<=0.0_r8) then
      write(*,*) 'sss =',ictot,ifac,ifbc,ifaq,sss(kcomp,ictot,ifac,ifbc,ifaq)
      write(*,*) 'Error in initialization of sss'
      stop
     endif

    enddo
    enddo
    enddo
    enddo

!       Modes 1 to 3 (SO4/BC/OC + condesate from H2SO4)
!  ************************************************************************

         do ifil = 1,3
          do lin = 1,16   ! 16 entries   
           read(19+ifil,993) calog1to3(ifil,lin), rk1to3(ifil,lin), &
             stdv1to3(ifil,lin), kcomp 

	do ic=1,16
	 if(calog1to3(ifil,lin).eq.cate(kcomp,ic)) then
	  ictot=ic
	  goto 71
	 endif
	end do
   71 continue

        sss1to3(kcomp,ictot) = stdv1to3(ifil,lin)
        rrr1to3(kcomp,ictot) = rk1to3(ifil,lin)

          end do   ! lin
         end do    ! ifil

    do kcomp=1,3
    do ictot=1,16
     if(sss1to3(kcomp,ictot)<=0.0_r8) then
      write(*,*) 'sss1to3 =',  ictot, sss1to3(kcomp,ictot)
      write(*,*) 'Error in initialization of sss1to3'
      stop
     endif
     if(rrr1to3(kcomp,ictot)<=0.0_r8) then
      write(*,*) 'rrr1to3 =', ictot, rrr1to3(kcomp,ictot)
      write(*,*) 'Error in initialization of rrr1to3'
      stop
     endif
    enddo
    enddo

!       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!  ************************************************************************

         ifil = 4
          do lin = 1,576   ! 16 entries   
           read(9+ifil,994) calog4(lin) &
             ,fraclog4(lin), fraqlog4(lin), rk4(lin), stdv4(lin), kcomp


	do ic=1,16
	 if(calog4(lin).gt.0.9999_r8*cate(kcomp,ic).and. &
	    calog4(lin).lt.1.0001_r8*cate(kcomp,ic)) then
	  ictot=ic
	  goto 91
	 endif
	end do
   91 continue

	do ic=1,6
	 if(fraclog4(lin).eq.fac(ic)) then
	  ifac=ic
	  goto 101
	 endif
	end do
  101 continue

	do ic=1,6
	 if(fraqlog4(lin).eq.faq(ic)) then
          ifaq=ic
          goto 111
	 endif
	end do
  111 continue

        rrr4(ictot,ifac,ifaq) = rk4(lin)
	sss4(ictot,ifac,ifaq) = stdv4(lin)

          end do   ! lin

    do ifac=1,6
    do ifaq=1,6
    do ictot=1,16
     if(rrr4(ictot,ifac,ifaq)<=0.0_r8) then
      write(*,*) 'rrr4 =',ictot,ifac,ifaq,rrr4(ictot,ifac,ifaq)
      write(*,*) 'Error in initialization of rrr4'
      stop

     endif
     if(sss4(ictot,ifac,ifaq)<=0.0_r8) then
      write(*,*) 'sss4 =',ictot,ifac,ifaq,sss4(ictot,ifac,ifaq)
      write(*,*) 'Error in initialization of sss4'
      stop
     endif



    enddo
    enddo
    enddo

        do ifil=13,22
          close (ifil)
        end do 


  993 format(3(x,e12.5),x,I3)
  994 format(5(x,e12.5),x,I3)
  995 format(6(x,e12.5),x,I3)

	return
	end








