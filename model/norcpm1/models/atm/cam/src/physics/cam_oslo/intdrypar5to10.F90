subroutine intdrypar5to10 (lchnk, ncol, Nnatk, Camk, xfac, xfbc, xfaq, & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,   & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,   &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
           cknorm,cknlt05,ckngt125)

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
!   use inpgraer
   use opttab,   only: cate, cat, fac, faq, fbc, rh, nbmodes, nmodes, nbmp1 

   implicit none

#include <aerodry.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
!
! Input-Output arguments
!
   real(r8), intent(inout) :: xfac(pcols,pver,nbmodes) ! modal (OC+BC)/(SO4+BC+OC)
   real(r8), intent(inout) :: xfbc(pcols,pver,nbmodes) ! modal BC/(OC+BC)
   real(r8), intent(inout) :: xfaq(pcols,pver,nbmodes) ! modal SO4(aq)/SO4
   real(r8), intent(inout) :: &
     cknorm(pcols,pver,0:nmodes), cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
!
! Output arguments: Modal mass concentrations (cint), area (aaero) and volume (vaero)
! (for AeroCom determination of particle effective radii) of each constituent. cint*05 
! and cint*125 are  for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.  
!
   real(r8), intent(out) :: &
     cintbg(pcols,pver,0:nbmodes), cintbg05(pcols,pver,0:nbmodes), cintbg125(pcols,pver,0:nbmodes), & 
     cintbc(pcols,pver,0:nbmodes), cintbc05(pcols,pver,0:nbmodes), cintbc125(pcols,pver,0:nbmodes), & 
     cintoc(pcols,pver,0:nbmodes), cintoc05(pcols,pver,0:nbmodes), cintoc125(pcols,pver,0:nbmodes), &
     cintsc(pcols,pver,0:nbmodes), cintsc05(pcols,pver,0:nbmodes), cintsc125(pcols,pver,0:nbmodes), &
     cintsa(pcols,pver,0:nbmodes), cintsa05(pcols,pver,0:nbmodes), cintsa125(pcols,pver,0:nbmodes), &
     aaeros(pcols,pver,0:nbmodes), aaerol(pcols,pver,0:nbmodes),                                    &
     vaeros(pcols,pver,0:nbmodes), vaerol(pcols,pver,0:nbmodes)
!
!---------------------------Local variables-----------------------------
!
      real(r8) a, b, e, eps, catot, fabc, fraq, xct(pcols,pver)

      real(r8) & 
        arr111, arr112, arr121, arr122, arr211, arr212, arr221, arr222, &
        arr11, arr12, arr21, arr22, arre1, arre2  
      real(r8) &
        arr1111, arr1112, arr1121, arr1122, arr1211, arr1212, arr1221, &
        arr1222, arr2111, arr2112, arr2121, arr2122, arr2211, arr2212, & 
        arr2221, arr2222
      real(r8) &
        arr11111, arr11112, arr11121, arr11122, arr11211, arr11212, & 
        arr11221, arr11222, arr12111, arr12112, arr12121, arr12122, &
        arr12211, arr12212, arr12221, arr12222, arr21111, arr21112, &
        arr21121, arr21122, arr21211, arr21212, arr21221, arr21222, &
        arr22111, arr22112, arr22121, arr22122, arr22211, arr22212, &
        arr22221, arr22222

      integer i, ierr, ictot, ifac, ifbc, ifaq, kcomp, k, icol
      integer ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
        ifbc1(pcols,pver), ifbc2(pcols,pver), ifaq1(pcols,pver),    &
        ifaq2(pcols,pver)
!      Temporary storage of often used array elements
      integer t_ict1, t_ict2, t_ifa1, t_ifa2
      integer t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fbc1, t_fbc2, t_xfbc
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xct,  t_cat1, t_cat2
      real(r8) esssf1, esssf2, esssf3, esssf4, esssf5, esssf6, esssf7, &
        esssf8, esssf9, esssf10, esssf3x6x9x10

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before kcomp-loop'

!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

        do kcomp=5,10

!      initialize output fields
      do k=1,pver
         do icol=1,ncol
        cintbg(icol,k,kcomp)=0.0_r8
        cintbg05(icol,k,kcomp)=0.0_r8
        cintbg125(icol,k,kcomp)=0.0_r8
        cintbc(icol,k,kcomp)=0.0_r8
        cintbc05(icol,k,kcomp)=0.0_r8
        cintbc125(icol,k,kcomp)=0.0_r8
        cintoc(icol,k,kcomp)=0.0_r8
        cintoc05(icol,k,kcomp)=0.0_r8
        cintoc125(icol,k,kcomp)=0.0_r8
        cintsc(icol,k,kcomp)=0.0_r8
        cintsc05(icol,k,kcomp)=0.0_r8
        cintsc125(icol,k,kcomp)=0.0_r8
        cintsa(icol,k,kcomp)=0.0_r8
        cintsa05(icol,k,kcomp)=0.0_r8
        cintsa125(icol,k,kcomp)=0.0_r8
        aaeros(icol,k,kcomp)=0.0_r8
        aaerol(icol,k,kcomp)=0.0_r8
        vaeros(icol,k,kcomp)=0.0_r8
        vaerol(icol,k,kcomp)=0.0_r8
         end do
       end do

!      write(*,*) 'Before x-loop'
      do k=1,pver
         do icol=1,ncol

          if(Nnatk(icol,k,kcomp)>0.0_r8) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cat(kcomp,1)),cat(kcomp,6))
          xfac(icol,k,kcomp) = min(max(xfac(icol,k,kcomp),fac(1)),fac(6))
          xfbc(icol,k,kcomp) = min(max(xfbc(icol,k,kcomp),fbc(1)),fbc(6))
          xfaq(icol,k,kcomp) = min(max(xfaq(icol,k,kcomp),faq(1)),faq(6))

!      write(*,*) 'Before cat-loop', kcomp
      do ictot=1,5
            if(xct(icol,k)>=cat(kcomp,ictot).and. &
            xct(icol,k)<=cat(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k,kcomp)>=fac(ifac).and. &
             xfac(icol,k,kcomp)<=fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

!      write(*,*) 'Before fbc-loop', kcomp
      do ifbc=1,5
            if(xfbc(icol,k,kcomp)>=fbc(ifbc).and. &
             xfbc(icol,k,kcomp)<=fbc(ifbc+1)) then
             ifbc1(icol,k)=ifbc
             ifbc2(icol,k)=ifbc+1
            endif
      end do ! ifbc

!      write(*,*) 'Before faq-loop', kcomp
      do ifaq=1,5
            if(xfaq(icol,k,kcomp)>=faq(ifaq).and. &
            xfaq(icol,k,kcomp)<=faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp)>0.0_r8) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifb1 = ifbc1(icol,k)
      t_ifb2 = ifbc2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)
      t_cat1 = cat(kcomp,t_ict1)
      t_cat2 = cat(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_fbc1 = fbc(t_ifb1)
      t_fbc2 = fbc(t_ifb2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k,kcomp)
      t_xfbc = xfbc(icol,k,kcomp)
      t_xfaq = xfaq(icol,k,kcomp)

      esssf1 = (t_faq2-t_xfaq)
      esssf2 = (t_xfaq-t_faq1)
      esssf3 = 1.0_r8/(t_faq2-t_faq1)
      esssf4 = (t_fbc2-t_xfbc)
      esssf5 = (t_xfbc-t_fbc1)
      esssf6 = 1.0_r8/(t_fbc2-t_fbc1)
      esssf7 = (t_fac2-t_xfac)
      esssf8 = (t_xfac-t_fac1)
      esssf9 = 1.0_r8/(t_fac2-t_fac1)
      esssf10= 1.0_r8/(t_cat2-t_cat1)
      esssf3x6x9x10 = esssf3*esssf6*esssf9*esssf10

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do i=1,19  ! variable number

       if(i==1) then
      arr11111=a5to10cintbg(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbg(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbg(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbg(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbg(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbg(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbg(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbg(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbg(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbg(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbg(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbg(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbg(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbg(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbg(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbg(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==2) then
      arr11111=a5to10cintbg05(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbg05(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbg05(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbg05(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbg05(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbg05(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbg05(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbg05(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbg05(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbg05(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbg05(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbg05(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbg05(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbg05(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbg05(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbg05(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==3) then
      arr11111=a5to10cintbg125(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbg125(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbg125(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbg125(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbg125(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbg125(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbg125(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbg125(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbg125(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbg125(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbg125(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbg125(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbg125(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbg125(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbg125(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbg125(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==4) then
      arr11111=a5to10cintbc(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbc(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbc(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbc(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbc(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbc(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbc(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbc(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbc(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbc(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbc(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbc(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbc(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbc(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbc(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbc(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==5) then
      arr11111=a5to10cintbc05(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbc05(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbc05(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbc05(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbc05(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbc05(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbc05(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbc05(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbc05(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbc05(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbc05(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbc05(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbc05(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbc05(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbc05(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbc05(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==6) then
      arr11111=a5to10cintbc125(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintbc125(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintbc125(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintbc125(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintbc125(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintbc125(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintbc125(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintbc125(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintbc125(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintbc125(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintbc125(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintbc125(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintbc125(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintbc125(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintbc125(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintbc125(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==7) then
      arr11111=a5to10cintoc(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintoc(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintoc(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintoc(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintoc(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintoc(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintoc(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintoc(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintoc(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintoc(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintoc(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintoc(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintoc(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintoc(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintoc(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintoc(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==8) then
      arr11111=a5to10cintoc05(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintoc05(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintoc05(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintoc05(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintoc05(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintoc05(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintoc05(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintoc05(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintoc05(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintoc05(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintoc05(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintoc05(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintoc05(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintoc05(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintoc05(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintoc05(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==9) then
      arr11111=a5to10cintoc125(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintoc125(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintoc125(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintoc125(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintoc125(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintoc125(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintoc125(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintoc125(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintoc125(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintoc125(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintoc125(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintoc125(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintoc125(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintoc125(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintoc125(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintoc125(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==10) then
      arr11111=a5to10cintsc(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsc(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsc(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsc(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsc(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsc(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsc(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsc(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsc(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsc(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsc(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsc(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsc(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsc(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsc(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsc(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==11) then
      arr11111=a5to10cintsc05(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsc05(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsc05(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsc05(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsc05(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsc05(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsc05(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsc05(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsc05(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsc05(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsc05(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsc05(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsc05(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsc05(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsc05(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsc05(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==12) then
      arr11111=a5to10cintsc125(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsc125(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsc125(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsc125(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsc125(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsc125(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsc125(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsc125(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsc125(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsc125(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsc125(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsc125(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsc125(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsc125(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsc125(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsc125(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==13) then
      arr11111=a5to10cintsa(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsa(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsa(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsa(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsa(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsa(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsa(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsa(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsa(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsa(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsa(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsa(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsa(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsa(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsa(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsa(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==14) then
      arr11111=a5to10cintsa05(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsa05(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsa05(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsa05(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsa05(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsa05(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsa05(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsa05(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsa05(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsa05(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsa05(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsa05(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsa05(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsa05(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsa05(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsa05(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==15) then
      arr11111=a5to10cintsa125(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10cintsa125(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10cintsa125(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10cintsa125(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10cintsa125(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10cintsa125(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10cintsa125(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10cintsa125(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10cintsa125(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10cintsa125(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10cintsa125(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10cintsa125(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10cintsa125(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10cintsa125(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10cintsa125(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10cintsa125(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==16) then
      arr11111=a5to10aaeros(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10aaeros(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10aaeros(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10aaeros(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10aaeros(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10aaeros(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10aaeros(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10aaeros(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10aaeros(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10aaeros(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10aaeros(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10aaeros(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10aaeros(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10aaeros(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10aaeros(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10aaeros(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==17) then
      arr11111=a5to10aaerol(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10aaerol(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10aaerol(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10aaerol(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10aaerol(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10aaerol(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10aaerol(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10aaerol(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10aaerol(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10aaerol(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10aaerol(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10aaerol(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10aaerol(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10aaerol(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10aaerol(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10aaerol(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==18) then
      arr11111=a5to10vaeros(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10vaeros(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10vaeros(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10vaeros(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10vaeros(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10vaeros(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10vaeros(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10vaeros(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10vaeros(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10vaeros(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10vaeros(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10vaeros(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10vaeros(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10vaeros(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10vaeros(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10vaeros(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       elseif(i==19) then
      arr11111=a5to10vaerol(t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=a5to10vaerol(t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=a5to10vaerol(t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=a5to10vaerol(t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=a5to10vaerol(t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=a5to10vaerol(t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=a5to10vaerol(t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=a5to10vaerol(t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=a5to10vaerol(t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=a5to10vaerol(t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=a5to10vaerol(t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=a5to10vaerol(t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=a5to10vaerol(t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=a5to10vaerol(t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=a5to10vaerol(t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=a5to10vaerol(t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
       endif

      arr1111=esssf1*arr11111+esssf2*arr11112
      arr1112=esssf1*arr11121+esssf2*arr11122
      arr1121=esssf1*arr11211+esssf2*arr11212
      arr1122=esssf1*arr11221+esssf2*arr11222
      arr1211=esssf1*arr12111+esssf2*arr12112
      arr1212=esssf1*arr12121+esssf2*arr12122
      arr1221=esssf1*arr12211+esssf2*arr12212
      arr1222=esssf1*arr12221+esssf2*arr12222

      arr111=esssf4*arr1111+esssf5*arr1112
      arr112=esssf4*arr1121+esssf5*arr1122
      arr121=esssf4*arr1211+esssf5*arr1212
      arr122=esssf4*arr1221+esssf5*arr1222

      arr11 =esssf7*arr111+esssf8*arr112
      arr12 =esssf7*arr121+esssf8*arr122

!      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf10      

      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf3x6x9x10

!      write(*,*) 'Before array'

       if(i==1) then
         cintbg(icol,k,kcomp)=arre1
       elseif(i==2) then
         cintbg05(icol,k,kcomp)=arre1
       elseif(i==3) then
         cintbg125(icol,k,kcomp)=arre1
       elseif(i==4) then
        cintbc(icol,k,kcomp)=arre1
       elseif(i==5) then
        cintbc05(icol,k,kcomp)=arre1
       elseif(i==6) then
        cintbc125(icol,k,kcomp)=arre1
       elseif(i==7) then
        cintoc(icol,k,kcomp)=arre1
       elseif(i==8) then
        cintoc05(icol,k,kcomp)=arre1
       elseif(i==9) then
        cintoc125(icol,k,kcomp)=arre1
       elseif(i==10) then
        cintsc(icol,k,kcomp)=arre1
       elseif(i==11) then
        cintsc05(icol,k,kcomp)=arre1
       elseif(i==12) then
        cintsc125(icol,k,kcomp)=arre1
       elseif(i==13) then
        cintsa(icol,k,kcomp)=arre1
       elseif(i==14) then
        cintsa05(icol,k,kcomp)=arre1
       elseif(i==15) then
        cintsa125(icol,k,kcomp)=arre1
       elseif(i==16) then
        aaeros(icol,k,kcomp)=arre1
       elseif(i==17) then
        aaerol(icol,k,kcomp)=arre1
       elseif(i==18) then
        vaeros(icol,k,kcomp)=arre1
       elseif(i==19) then
        vaerol(icol,k,kcomp)=arre1
       endif

         end do ! i=1,19 

           endif
         
          cknorm(icol,k,kcomp)  = a5to10cintbg(1,1,1,1,kcomp)
          cknlt05(icol,k,kcomp) = a5to10cintbg05(1,1,1,1,kcomp)
          ckngt125(icol,k,kcomp)= a5to10cintbg125(1,1,1,1,kcomp)

       end do ! icol
      end do ! k

        end do  ! kcomp

      return
end subroutine intdrypar5to10




