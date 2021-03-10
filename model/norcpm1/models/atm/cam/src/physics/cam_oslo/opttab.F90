module opttab

! Purpose: To read in look-up tables and calculate optical properties for the aerosols
!    For subroutine interpol:


  use shr_kind_mod, only: r8 => shr_kind_r8
  use cam_logfile,  only: iulog
  implicit none

  private 
  save


  ! Interfaces
  public initopt


!     Array bounds in the tabulated optical parameters  
   integer, public, parameter :: nbands=12    ! number of aerosol spectral bands
   integer, public, parameter :: nbmodes=10   ! number of aerosol background modes (wrt internal mixing)
   integer, public, parameter :: nbmp1=11     ! number of first non-background mode
   integer, public, parameter :: nmodes=14    ! total number of aerosol modes

   real(r8), public, dimension(6) :: fac = (/ 0.0_r8, 0.1_r8,  0.3_r8, 0.5_r8, 0.7_r8, 0.999_r8    /)
   real(r8), public, dimension(6) :: fbc = (/ 0.0_r8, 0.01_r8, 0.1_r8, 0.3_r8, 0.7_r8, 0.999_r8    /)
   real(r8), public, dimension(6) :: faq = (/ 0.0_r8, 0.25_r8, 0.5_r8, 0.75_r8,0.85_r8,1.0_r8      /)
   real(r8), public, dimension(10) :: rh = (/ 0.0_r8, 0.37_r8, 0.47_r8,0.65_r8,0.75_r8,            &
!test                                      0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8,0.98_r8             /)
                                      0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8,0.995_r8             /)
   real(r8), public, dimension(9) :: S   = (/ 1.0005_r8, 1.001_r8,  1.0015_r8, 1.002_r8,           &
                                      1.0025_r8, 1.003_r8,  1.005_r8,  1.008_r8,  1.01_r8  /)  


        real(r8), public, dimension(5:10,6) :: cat = reshape ( (/ &
   1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8, & 
   5.e-4_r8 , 0.01_r8  , 0.02_r8  , 1.e-4_r8 , 0.005_r8 , 0.02_r8  , &
   2.e-3_r8 , 0.05_r8  , 0.1_r8   , 6.e-4_r8 , 0.025_r8 , 0.1_r8   , &
   0.01_r8  , 0.2_r8   , 0.5_r8   , 2.5e-3_r8, 0.1_r8   , 0.5_r8   , &
   0.04_r8  , 0.8_r8   , 2.0_r8   , 1.e-2_r8 , 0.4_r8   , 2.0_r8   , &
   0.15_r8  , 4.0_r8   , 8.0_r8   , 3.5e-2_r8, 2.0_r8   , 8.0_r8     &
                                                      /), (/6,6/) )

        real(r8), public, dimension(4,16) :: cate = reshape ( (/ &
   1.e-10_r8, 1.e-10_r8, 1.e-10_r8, 1.e-10_r8*1.904e-3_r8, & 
   1.e-5_r8 , 1.e-5_r8 , 1.e-4_r8 , 0.01_r8*1.904e-3_r8  , &
   2.e-5_r8 , 2.e-5_r8 , 2.e-4_r8 , 0.05_r8*1.904e-3_r8  , &
   4.e-5_r8 , 4.e-5_r8 , 4.e-4_r8 , 0.1_r8*1.904e-3_r8   , &
   8.e-5_r8 , 8.e-5_r8 , 8.e-4_r8 , 0.2_r8*1.904e-3_r8   , &
   1.5e-4_r8, 1.5e-4_r8, 1.5e-3_r8, 0.4_r8*1.904e-3_r8   , &
   3.e-4_r8 , 3.e-4_r8 , 3.e-3_r8 , 0.7_r8*1.904e-3_r8   , &
   6.e-4_r8 , 6.e-4_r8 , 6.e-3_r8 , 1.0_r8*1.904e-3_r8   , &
   1.2e-3_r8, 1.2e-3_r8, 1.2e-2_r8, 1.5_r8*1.904e-3_r8   , &
   2.5e-3_r8, 2.5e-3_r8, 2.5e-2_r8, 2.5_r8*1.904e-3_r8   , &
   5.e-3_r8 , 5.e-3_r8 , 5.e-2_r8 , 5.0_r8*1.904e-3_r8   , &
   1.e-2_r8 , 1.e-2_r8 , 0.1_r8   , 10.0_r8*1.904e-3_r8  , &
   2.e-2_r8 , 2.e-2_r8 , 0.2_r8   , 25.0_r8*1.904e-3_r8  , &
   4.e-2_r8 , 4.e-2_r8 , 0.4_r8   , 50.0_r8*1.904e-3_r8  , &
   8.e-2_r8 , 8.e-2_r8 , 0.8_r8   , 100.0_r8*1.904e-3_r8 , &
   0.15_r8  , 0.15_r8  , 1.5_r8   , 500.0_r8*1.904e-3_r8 /), (/4,16/) )
   

  real(r8), public :: om1to3(12,10,16,3)
  real(r8), public :: g1to3(12,10,16,3)
  real(r8), public :: be1to3(12,10,16,3)
  real(r8), public :: ke1to3(12,10,16,3)

  real(r8), public :: om4(12,10,16,6,6)
  real(r8), public :: g4(12,10,16,6,6)
  real(r8), public :: be4(12,10,16,6,6)
  real(r8), public :: ke4(12,10,16,6,6)

  real(r8), public,dimension(12) :: om0
  real(r8), public :: g0(12)
  real(r8), public :: be0(12)
  real(r8), public :: ke0(12)

  real(r8), public :: om5to10(12,10,6,6,6,6,5:10)
  real(r8), public :: g5to10(12,10,6,6,6,6,5:10)
  real(r8), public :: be5to10(12,10,6,6,6,6,5:10)
  real(r8), public :: ke5to10(12,10,6,6,6,6,5:10)

  real(r8), public :: e,eps
  parameter (e=2.718281828_r8, eps=1.0e-30_r8)
  


 contains

subroutine initopt

!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1. Files also 
!   involved by this modification: interpol1.F and opttab.h.
!   Modified for new aerosol schemes by Alf Kirkevaag in January 
!   2006.
!---------------------------------------------------------------

!   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer


!   implicit none

!     Tabulating the 'kcomp'-files to save computing time.


      integer kcomp, iwl, irelh, ictot, ifac, ifbc, ifaq
      integer ic, ifil, lin
      real(r8) catot, relh, frac, fabc, fraq
      real(r8) ssa, ass, ext, spext
      real(r8) rh2(10)
      real(r8) :: eps2 = 1.e-2_r8
      real(r8) :: eps4 = 1.e-4_r8
      real(r8) :: eps6 = 1.e-6_r8
      real(r8) :: eps7 = 1.e-7_r8

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      open(40,file='inputdata/atm/cam/camoslo/aertab/kcomp1.out' &
             ,form='formatted',status='old')
      open(41,file='inputdata/atm/cam/camoslo/aertab/kcomp2.out' &
             ,form='formatted',status='old')
      open(42,file='inputdata/atm/cam/camoslo/aertab/kcomp3.out' &
             ,form='formatted',status='old')
      open(43,file='inputdata/atm/cam/camoslo/aertab/kcomp4.out' &
             ,form='formatted',status='old')
      open(44,file='inputdata/atm/cam/camoslo/aertab/kcomp5.out' &
             ,form='formatted',status='old')
      open(45,file='inputdata/atm/cam/camoslo/aertab/kcomp6.out' &
             ,form='formatted',status='old')
      open(46,file='inputdata/atm/cam/camoslo/aertab/kcomp7.out' &
             ,form='formatted',status='old')
      open(47,file='inputdata/atm/cam/camoslo/aertab/kcomp8.out' &
             ,form='formatted',status='old')
      open(48,file='inputdata/atm/cam/camoslo/aertab/kcomp9.out' &
             ,form='formatted',status='old')
      open(49,file='inputdata/atm/cam/camoslo/aertab/kcomp10.out'& 
             ,form='formatted',status='old')
      open(50,file='inputdata/atm/cam/camoslo/aertab/kcomp0.out'& 
             ,form='formatted',status='old')


!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

      do ifil = 5,10
        do lin = 1,155520   ! 12*10*6*6*6*6 

          read(39+ifil,993) kcomp, iwl, relh, catot, frac, fabc, fraq, &
                          ssa, ass, ext, spext


!pga feil i antall siffer i tabellene:
          if(abs(relh-0.99_r8)<eps6) then          
            relh=0.995_r8
          endif    

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

          om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ssa    
          g5to10 (iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ass
          be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ext    ! unit km^-1
          ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, g  =', kcomp, g5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, be =', kcomp, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
!      write(iulog,*) 'kcomp, ke =', kcomp, ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp) 
        end do
      end do


    do kcomp=5,10
    do iwl=1,12
    do irelh=1,10
    do ictot=1,6
    do ifac=1,6
    do ifaq=1,6
     if(be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
      write(iulog,*) 'be5to10 =', iwl, irelh, ictot, ifac, ifbc, ifaq, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
      write(iulog,*) 'Error in initialization of be5to10'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 5-10 ok' 


!       Modes 1 to 3 (SO4/BC/OC + condesate from H2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
      do ifil = 1,3
        do lin = 1,1920   ! 12*10*16 

          read(39+ifil,994) kcomp, iwl, relh, catot, &
                          ssa, ass, ext, spext
!pga feil i antall siffer i tabellene:
          if(abs(relh-0.99_r8)<eps6) then          
            relh=0.995_r8
          endif    

       	  do ic=1,10
	   if(abs(relh-rh(ic))<eps4) then
	    irelh=ic
	    goto 61
	   endif
	  end do
   61     continue

 	  do ic=1,16
!skumelt	   if(abs(catot-cate(kcomp,ic))<eps7) then
	   if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
	    ictot=ic
	    goto 71
	   endif
	  end do
   71     continue

          om1to3(iwl,irelh,ictot,kcomp)=ssa    
          g1to3 (iwl,irelh,ictot,kcomp)=ass
          be1to3(iwl,irelh,ictot,kcomp)=ext    ! unit km^-1
          ke1to3(iwl,irelh,ictot,kcomp)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om1to3(iwl,irelh,ictot,kcomp)
!      write(iulog,*) 'kcomp, g  =', kcomp, g1to3(iwl,irelh,ictot,kcomp)
!      write(iulog,*) 'kcomp, be =', kcomp, be1to3(iwl,irelh,ictot,kcomp)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke1to3(iwl,irelh,ictot,kcomp)
!      if(ifil==1) write(iulog,*) 'iwl,irelh,ictot,kcomp,ke =', iwl,irelh,ictot,kcomp,ke1to3(iwl,irelh,ictot,kcomp)

        end do
      end do
    do kcomp=1,3
    do iwl=1,12
    do irelh=1,10
    do ictot=1,16
     if(be1to3(iwl,irelh,ictot,kcomp)<=0.0_r8) then
      write(iulog,*) 'be1to3 =', iwl, irelh, ictot, be1to3(iwl,irelh,ictot,kcomp)
      write(iulog,*) 'Error in initialization of be1to3'
      stop
     endif
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 1-3 ok' 

!       Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 4
        do lin = 1,69120   ! 12*10*16*6*6 

          read(39+ifil,995) kcomp, iwl, relh, catot, frac, fraq, &
                          ssa, ass, ext, spext
!pga feil i antall siffer i tabellene:
          if(abs(relh-0.99_r8)<eps6) then          
            relh=0.995_r8
          endif    

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

          om4(iwl,irelh,ictot,ifac,ifaq)=ssa    
          g4 (iwl,irelh,ictot,ifac,ifaq)=ass
          be4(iwl,irelh,ictot,ifac,ifaq)=ext    ! unit km^-1
          ke4(iwl,irelh,ictot,ifac,ifaq)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om4(iwl,irelh,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, g  =', kcomp, g4(iwl,irelh,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, be =', kcomp, be4(iwl,irelh,ictot,ifac,ifaq)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke4(iwl,irelh,ictot,ifac,ifaq)
        end do

    do iwl=1,12
    do irelh=1,10
    do ictot=1,16
    do ifac=1,6
    do ifaq=1,6
     if(be4(iwl,irelh,ictot,ifac,ifaq)<=0.0_r8) then
      write(iulog,*) 'be4 =', iwl, irelh, ictot, ifac, ifaq, be4(iwl,irelh,ictot,ifac,ifaq)
      write(iulog,*) 'Error in initialization of be4'
      stop
     endif
    enddo
    enddo
    enddo
    enddo
    enddo

        write(iulog,*)'mode 4 ok' 

!       Mode 0, BC(ax)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

        ifil = 11
        do lin = 1,12   ! 12 

          read(39+ifil,996) kcomp, iwl, relh, &
                          ssa, ass, ext, spext
          om0(iwl)=ssa    
          g0 (iwl)=ass
          be0(iwl)=ext    ! unit km^-1
          ke0(iwl)=spext  ! unit m^2/g

!      write(iulog,*) 'kcomp, om =', kcomp, om0(iwl)
!      write(iulog,*) 'kcomp, g  =', kcomp, g0(iwl)
!      write(iulog,*) 'kcomp, be =', kcomp, be0(iwl)
!      write(iulog,*) 'kcomp, ke =', kcomp, ke0(iwl)

        end do

    do iwl=1,12
     if(be0(iwl)<=0.0_r8) then
      write(iulog,*) 'be0 =', iwl, be0(iwl)
      write(iulog,*) 'Error in initialization of be0'
      stop
     endif
    enddo

        write(iulog,*)'mode 0 ok' 

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

  993 format(2I3,f7.2,3(x,e10.3),f7.2,4(x,e12.5))
  994 format(2I3,f7.2,x,e10.3,4(x,e12.5))
  995 format(2I3,f7.2,3(x,e10.3),4(x,e12.5))
  996 format(2I3,f7.2,4(x,e12.5))

      do ifil=40,50
        close (ifil)
      end do 
      return
end subroutine initopt


end module opttab

