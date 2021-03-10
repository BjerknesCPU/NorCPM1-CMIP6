subroutine initccn

!  (Named initabc in CCM-Oslo and CAM2)
!  Tabulating Alf Kirkevaags 'ccnk'-files, to save computing time.
!  Preliminary version created by Jon Egill Kristjansson on 19 November 1998.
!  Re-coded for CAM by Trude Storelvmo and Alf Kirkevaag, 2004 and 2005.
!  Modified for new aerosol schemes by Alf Kirkevaag in January 2006.

   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab
   use const
   use cam_logfile,  only: iulog
   implicit none

   integer kcomp, isup, ictot, ifac, ifbc, ifaq
   integer ic, ifil, lin

   open(40,file='inputdata/atm/cam/camoslo/aertab/ccnk1.out' &   ! SO4(n/Ait)
          ,form='formatted',status='old')
   open(41,file='inputdata/atm/cam/camoslo/aertab/ccnk2.out' &   ! BC(n/Ait)
          ,form='formatted',status='old')
   open(42,file='inputdata/atm/cam/camoslo/aertab/ccnk3.out' &   ! OC(n/Ait)
          ,form='formatted',status='old')
!2oct07   open(13,file='inputdata/atm/cam/camoslo/aertab/ccnk4.out' &   ! BC&OC(n/Ait)
   open(43,file='inputdata/atm/cam/camoslo/aertab/ccnk4.out' &   ! BC&OC(n/Ait)
          ,form='formatted',status='old')
   open(44,file='inputdata/atm/cam/camoslo/aertab/ccnk5.out' &   ! SO4(Ait75)
          ,form='formatted',status='old')
   open(45,file='inputdata/atm/cam/camoslo/aertab/ccnk6.out' &   ! MINACC
          ,form='formatted',status='old')
   open(46,file='inputdata/atm/cam/camoslo/aertab/ccnk7.out' &   ! MINCOA
          ,form='formatted',status='old')
   open(47,file='inputdata/atm/cam/camoslo/aertab/ccnk8.out' &   ! SEASF
          ,form='formatted',status='old')
   open(48,file='inputdata/atm/cam/camoslo/aertab/ccnk9.out' &   ! SEASACC
          ,form='formatted',status='old')
   open(49,file='inputdata/atm/cam/camoslo/aertab/ccnk10.out' &  ! SEASCOA
          ,form='formatted',status='old')


        write(iulog,*)'In subroutine initabc' 

!       Modes 5 to 10 (SO4(ait75) and mineral and seasalt-modes + cond./coag./aq.)
!  ************************************************************************

         do ifil = 5,10

          do lin = 1,11664   ! 9*6**4 entries

           read(39+ifil,993) kcomp, sup(ifil,lin), catot(ifil,lin) &
             ,frac5to10(ifil,lin), fabc5to10(ifil,lin), fraq5to10(ifil,lin) &
             ,fccn5to10(ifil,lin)

        do ic=1,9
	 if(sup(ifil,lin).eq.S(ic)) then
	   isup=ic
	   goto 11
	 endif
       end do
   11 continue

	do ic=1,6
	 if(catot(ifil,lin).eq.cat(kcomp,ic)) then
	  ictot=ic
	  goto 21
	 endif
	end do
   21 continue

	do ic=1,6
	 if(frac5to10(ifil,lin).eq.fac(ic)) then
	  ifac=ic
	  goto 31
	 endif
	end do
   31 continue

	do ic=1,6
	 if(fabc5to10(ifil,lin).eq.fbc(ic)) then
	  ifbc=ic
	  goto 41
	 endif
	end do
   41 continue

	do ic=1,6
	 if(fraq5to10(ifil,lin).eq.faq(ic)) then
	  ifaq=ic
	  goto 51
	 endif
	end do
   51 continue

        fff(kcomp,isup,ictot,ifac,ifbc,ifaq) = fccn5to10(ifil,lin)

!        write(iulog,*)'kcomp, fff =', kcomp, fff(kcomp,isup,ictot,ifac,ifbc,ifaq)  
          end do   ! lin
         end do    ! ifil


    do kcomp=5,10
    do isup=1,9
    do ifac=1,6
    do ifbc=1,6
    do ictot=1,6
!     if(fff(kcomp,isup,ictot,ifac,ifbc,ifaq)<=0.0_r8) then
     if(fff(kcomp,isup,ictot,ifac,ifbc,ifaq).lt.0._r8) then
      write(iulog,*) 'fff =',isup,ictot,ifac,ifbc,ifaq,fff(kcomp,isup,ictot,ifac,ifbc,ifaq)
      write(iulog,*) 'Error in initialization of fff'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 5-10 ok' 

!       Modes 1 to 3 (SO4/BC/OC + condesate from H2SO4)
!  ************************************************************************

         do ifil = 1,3
          do lin = 1,144   ! 9*16 entries   
           read(39+ifil,994) kcomp, sup(ifil,lin), catot1to3(ifil,lin) &
             ,fccn1to3(ifil,lin)

        do ic=1,9
	 if(sup(ifil,lin).eq.S(ic)) then
	   isup=ic
	   goto 61
	 endif
        end do
   61 continue

	do ic=1,16
	 if(catot1to3(ifil,lin).eq.cate(kcomp,ic)) then
	  ictot=ic
	  goto 71
	 endif
	end do
   71 continue

        fff1to3(kcomp,isup,ictot) = fccn1to3(ifil,lin)

!        write(iulog,*)'kcomp, fff1to3 =', kcomp, fff1to3(kcomp,isup,ictot)  

          end do   ! lin
         end do    ! ifil

    do kcomp=1,3
    do isup=1,9
    do ictot=1,16
     if(fff1to3(kcomp,isup,ictot).lt.0._r8) then
      write(iulog,*) 'fff1to3 =', isup, ictot, fff1to3(kcomp,isup,ictot)
      write(iulog,*) 'Error in initialization of fff1to3'
      stop
     endif
    enddo
    enddo
    enddo
!        write(iulog,*)'mode 1-3 ok' 


!       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!  ************************************************************************

         ifil = 4
!4          do lin = 1,864   ! 9*16*6 entries   
          do lin = 1,5184   ! 9*16*6*6 entries   
           read(39+ifil,995) kcomp, sup(ifil,lin), catot4(lin) &
!4             ,frac4(lin), fccn4(lin)
             ,frac4(lin), fraq4(lin), fccn4(lin)

        do ic=1,9
	 if(sup(ifil,lin).eq.S(ic)) then
	   isup=ic
	   goto 81
	 endif
        end do
   81 continue

	do ic=1,16
	 if(catot4(lin).gt.0.9999_r8*cate(kcomp,ic).and. &
	    catot4(lin).lt.1.0001_r8*cate(kcomp,ic)) then
	  ictot=ic
	  goto 91
	 endif
	end do
   91 continue

	do ic=1,6
	 if(frac4(lin).eq.fac(ic)) then
	  ifac=ic
	  goto 101
	 endif
	end do
  101 continue

	do ic=1,6
	 if(fraq4(lin).eq.faq(ic)) then
          ifaq=ic
          goto 111
	 endif
	end do
  111 continue

!4        fff4(isup,ictot,ifac) = fccn4(lin)
        fff4(isup,ictot,ifac,ifaq) = fccn4(lin)

!        write(iulog,*)'s, c, fac, fff4 =', isup,ictot,ifac,fff4(isup,ictot,ifac,ifaq)  

          end do   ! lin

    do isup=1,9
    do ifac=1,6
    do ifaq=1,6
    do ictot=1,16
     if(fff4(isup,ictot,ifac,ifaq).lt.0._r8) then
      write(iulog,*) 'fff4 =',isup,ictot,ifac,ifaq,fff4(isup,ictot,ifac,ifaq)
      write(iulog,*) 'Error in initialization of fff4'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

!        write(iulog,*)'mode 4 ok' 

        do ifil=40,49
          close (ifil)
        end do 

  993 format(i3,x,e13.6,5(x,e12.5))
  994 format(i3,x,e13.6,2(x,e12.5))
!4  995 format(i3,x,e13.6,3(x,e12.5))
  995 format(i3,x,e13.6,4(x,e12.5))
	return
	end








