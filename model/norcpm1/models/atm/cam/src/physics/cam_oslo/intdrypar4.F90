subroutine intdrypar4 (lchnk, ncol, Nnatk, Camk, xfacin, xfaqin,     & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol)

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab,   only: cate, cat, fac, faq, fbc, rh, nbmodes, nmodes, nbmp1

   implicit none

#include <aerodry.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed (cond+aq) SO4 conc.
   real(r8), intent(in) :: xfacin(pcols,pver)       ! BC/(BC+OC) for he background mode
   real(r8), intent(in) :: xfaqin(pcols,pver)       ! SO4(aq)/SO4
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
      real(r8) a, b, e, eps
      real(r8) xct(pcols,pver), xfac(pcols,pver), xfaq(pcols,pver)

      real(r8) & 
        arr111, arr112, arr121, arr122, arr211, arr212, arr221, arr222, &
        arr11, arr12, arr21, arr22, arre1, arre2  
      real(r8) &
        arr1111, arr1112, arr1121, arr1122, arr1211, arr1212, arr1221, &
        arr1222, arr2111, arr2112, arr2121, arr2122, arr2211, arr2212, & 
        arr2221, arr2222

      integer i, ierr, ictot, ifac, ifaq, kcomp, k, icol
      integer ict1(pcols,pver), ict2(pcols,pver), ifac1(pcols,pver), &
       ifac2(pcols,pver), ifaq1(pcols,pver), ifaq2(pcols,pver)
!      Temporary storage of often used array elements
      integer t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifa1, t_ifa2
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xct,  t_cat1, t_cat2
      real(r8) esssf1, esssf2, esssf3, esssf7, esssf8, esssf9, esssf10

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before kcomp-loop'

!       Mode 4, BC&OC(Ait):

        kcomp=4

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
                 /(Nnatk(icol,k,kcomp)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k),fac(1)),fac(6))
          xfaq(icol,k) = min(max(xfaqin(icol,k),faq(1)),faq(6))

!   write(*,*) 'xct  =', icol, k, xct(icol,k)
!   write(*,*) 'xfac =', icol, k, xfac(icol,k)
!   write(*,*) 'xfaq =', icol, k, xfaq(icol,k)

!      write(*,*) 'Before cat-loop', kcomp
      do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k)<=cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

!ok   write(*,*) 'xct =', xct(icol,k), ict1(icol,k), ict2(icol,k)

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k)<=fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

!ok   write(*,*) 'xfac =', xfac(icol,k), ifac1(icol,k), ifac2(icol,k)     

!      write(*,*) 'Before faq-loop', kcomp
      do ifaq=1,5
            if(xfaq(icol,k)>=faq(ifaq).and. &
            xfaq(icol,k)<=faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq

!ok   write(*,*) 'xfaq =', xfaq(icol,k), ifaq1(icol,k), ifaq2(icol,k)
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)
      t_xct  = xct(icol,k)
      t_xfac = xfac(icol,k)
      t_xfaq = xfaq(icol,k)

      esssf1 = (t_faq2-t_xfaq)
      esssf2 = (t_xfaq-t_faq1)
      esssf3 = 1.0_r8/(t_faq2-t_faq1)

      esssf7 = (t_fac2-t_xfac)
      esssf8 = (t_xfac-t_fac1)
      esssf9 = 1.0_r8/(t_fac2-t_fac1)
      esssf10= 1.0_r8/(t_cat2-t_cat1)

!       write(*,*) 'esssf1,2,7 =', esssf1, esssf2, esssf7
!       write(*,*) 'esssf8,9,10 =', esssf8, esssf9, esssf10

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do i=1,19  ! variable number

       if(i==1) then
!4      arr111=a4cintbg(t_ict1,t_ifc1)
!4      arr112=a4cintbg(t_ict1,t_ifc2)
!4      arr121=a4cintbg(t_ict2,t_ifc1)
!4      arr122=a4cintbg(t_ict2,t_ifc2)
      arr1111=a4cintbg(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbg(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbg(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbg(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbg(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbg(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbg(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbg(t_ict2,t_ifc2,t_ifa2)
       elseif(i==2) then
!4      arr111=a4cintbg05(t_ict1,t_ifc1)
!4      arr112=a4cintbg05(t_ict1,t_ifc2)
!4      arr121=a4cintbg05(t_ict2,t_ifc1)
!4      arr122=a4cintbg05(t_ict2,t_ifc2)
      arr1111=a4cintbg05(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbg05(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbg05(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbg05(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbg05(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbg05(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbg05(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbg05(t_ict2,t_ifc2,t_ifa2)
       elseif(i==3) then
!4      arr111=a4cintbg125(t_ict1,t_ifc1)
!4      arr112=a4cintbg125(t_ict1,t_ifc2)
!4      arr121=a4cintbg125(t_ict2,t_ifc1)
!4      arr112=a4cintbg125(t_ict2,t_ifc2)
      arr1111=a4cintbg125(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbg125(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbg125(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbg125(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbg125(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbg125(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbg125(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbg125(t_ict2,t_ifc2,t_ifa2)
       elseif(i==4) then
!4      arr111=a4cintbc(t_ict1,t_ifc1)
!4      arr112=a4cintbc(t_ict1,t_ifc2)
!4      arr121=a4cintbc(t_ict2,t_ifc1)
!4      arr112=a4cintbc(t_ict2,t_ifc2)
      arr1111=a4cintbc(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbc(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbc(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbc(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbc(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbc(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbc(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbc(t_ict2,t_ifc2,t_ifa2)
       elseif(i==5) then
!4      arr111=a4cintbc05(t_ict1,t_ifc1)
!4      arr112=a4cintbc05(t_ict1,t_ifc2)
!4      arr121=a4cintbc05(t_ict2,t_ifc1)
!4      arr112=a4cintbc05(t_ict2,t_ifc2)
      arr1111=a4cintbc05(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbc05(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbc05(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbc05(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbc05(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbc05(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbc05(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbc05(t_ict2,t_ifc2,t_ifa2)
       elseif(i==6) then
!4      arr111=a4cintbc125(t_ict1,t_ifc1)
!4      arr112=a4cintbc125(t_ict1,t_ifc2)
!4      arr121=a4cintbc125(t_ict2,t_ifc1)
!4      arr112=a4cintbc125(t_ict2,t_ifc2)
      arr1111=a4cintbc125(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintbc125(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintbc125(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintbc125(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintbc125(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintbc125(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintbc125(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintbc125(t_ict2,t_ifc2,t_ifa2)
       elseif(i>=7.and.i<=9) then  ! no added OC
      arr111=eps
      arr112=eps
      arr121=eps
      arr112=eps
       elseif(i==10) then
!4      arr111=a4cintsc(t_ict1,t_ifc1)
!4      arr112=a4cintsc(t_ict1,t_ifc2)
!4      arr121=a4cintsc(t_ict2,t_ifc1)
!4      arr112=a4cintsc(t_ict2,t_ifc2)
      arr1111=a4cintsc(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsc(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsc(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsc(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsc(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsc(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsc(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsc(t_ict2,t_ifc2,t_ifa2)
       elseif(i==11) then
!4      arr111=a4cintsc05(t_ict1,t_ifc1)
!4      arr112=a4cintsc05(t_ict1,t_ifc2)
!4      arr121=a4cintsc05(t_ict2,t_ifc1)
!4      arr112=a4cintsc05(t_ict2,t_ifc2)
      arr1111=a4cintsc05(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsc05(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsc05(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsc05(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsc05(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsc05(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsc05(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsc05(t_ict2,t_ifc2,t_ifa2)
       elseif(i==12) then
!4      arr111=a4cintsc125(t_ict1,t_ifc1)
!4      arr112=a4cintsc125(t_ict1,t_ifc2)
!4      arr121=a4cintsc125(t_ict2,t_ifc1)
!4      arr112=a4cintsc125(t_ict2,t_ifc2)
      arr1111=a4cintsc125(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsc125(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsc125(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsc125(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsc125(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsc125(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsc125(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsc125(t_ict2,t_ifc2,t_ifa2)
!4       elseif(i>=13.and.i<=15) then  ! no SO4(aq)
!4      arr111=eps
!4      arr112=eps
!4      arr121=eps
!4      arr112=eps
       elseif(i==13) then
      arr1111=a4cintsa(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsa(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsa(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsa(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsa(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsa(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsa(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsa(t_ict2,t_ifc2,t_ifa2)
       elseif(i==14) then
      arr1111=a4cintsa05(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsa05(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsa05(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsa05(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsa05(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsa05(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsa05(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsa05(t_ict2,t_ifc2,t_ifa2)
       elseif(i==15) then
      arr1111=a4cintsa125(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4cintsa125(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4cintsa125(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4cintsa125(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4cintsa125(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4cintsa125(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4cintsa125(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4cintsa125(t_ict2,t_ifc2,t_ifa2)
       elseif(i==16) then
!4      arr111=a4aaeros(t_ict1,t_ifc1)
!4      arr112=a4aaeros(t_ict1,t_ifc2)
!4      arr121=a4aaeros(t_ict2,t_ifc1)
!4      arr112=a4aaeros(t_ict2,t_ifc2)
      arr1111=a4aaeros(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4aaeros(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4aaeros(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4aaeros(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4aaeros(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4aaeros(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4aaeros(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4aaeros(t_ict2,t_ifc2,t_ifa2)
       elseif(i==17) then
!4      arr111=a4aaerol(t_ict1,t_ifc1)
!4      arr112=a4aaerol(t_ict1,t_ifc2)
!4      arr121=a4aaerol(t_ict2,t_ifc1)
!4      arr112=a4aaerol(t_ict2,t_ifc2)
      arr1111=a4aaerol(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4aaerol(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4aaerol(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4aaerol(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4aaerol(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4aaerol(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4aaerol(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4aaerol(t_ict2,t_ifc2,t_ifa2)
       elseif(i==18) then
!4      arr111=a4vaeros(t_ict1,t_ifc1)
!4      arr112=a4vaeros(t_ict1,t_ifc2)
!4      arr121=a4vaeros(t_ict2,t_ifc1)
!4      arr112=a4vaeros(t_ict2,t_ifc2)
      arr1111=a4vaeros(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4vaeros(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4vaeros(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4vaeros(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4vaeros(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4vaeros(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4vaeros(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4vaeros(t_ict2,t_ifc2,t_ifa2)
       elseif(i==19) then
!4      arr111=a4vaerol(t_ict1,t_ifc1)
!4      arr112=a4vaerol(t_ict1,t_ifc2)
!4      arr121=a4vaerol(t_ict2,t_ifc1)
!4      arr112=a4vaerol(t_ict2,t_ifc2)
      arr1111=a4vaerol(t_ict1,t_ifc1,t_ifa1)
      arr1112=a4vaerol(t_ict1,t_ifc1,t_ifa2)
      arr1121=a4vaerol(t_ict1,t_ifc2,t_ifa1)
      arr1122=a4vaerol(t_ict1,t_ifc2,t_ifa2)
      arr1211=a4vaerol(t_ict2,t_ifc1,t_ifa1)
      arr1212=a4vaerol(t_ict2,t_ifc1,t_ifa2)
      arr1221=a4vaerol(t_ict2,t_ifc2,t_ifa1)
      arr1222=a4vaerol(t_ict2,t_ifc2,t_ifa2)
       endif

      arr111=esssf1*arr1111+esssf2*arr1112
      arr112=esssf1*arr1121+esssf2*arr1122
      arr121=esssf1*arr1211+esssf2*arr1212
      arr122=esssf1*arr1221+esssf2*arr1222

      arr11 =esssf7*arr111+esssf8*arr112
      arr12 =esssf7*arr121+esssf8*arr122

      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf3*esssf9*esssf10

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
         
       end do ! icol
      end do ! k


      return
end subroutine intdrypar4




