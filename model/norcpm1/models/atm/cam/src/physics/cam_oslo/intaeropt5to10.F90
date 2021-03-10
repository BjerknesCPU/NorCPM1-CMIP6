subroutine intaeropt5to10 (lchnk, ncol, rhum, Nnatk, Camk, xfacin, xfbcin, xfaqin, &
!           bext440, babs440, bext500, babs500, bext550,  babs550,      &
!           bext670, babs670, bext870, babs870,                         &
!           bebg440, babg440, bebg500, babg500, bebg550,                &
!           bebg670, babg670, bebg870, babg870,                         &
!           bebc440, babc440, bebc500, babc500, bebc550,                &
!           bebc670, babc670, bebc870, babc870,                         &
!           beoc440, baoc440, beoc500, baoc500, beoc550,                &
!           beoc670, baoc670, beoc870, baoc870,                         &
!           besu440, basu440, besu500, basu500, besu550,                &
!           besu670, basu670, besu870, basu870,                         &
           bext440, bext500, bext550, bext670, bext870,                &
           bebg440, bebg500, bebg550, bebg670, bebg870,                &
           bebc440, bebc500, bebc550, bebc670, bebc870,                &
           beoc440, beoc500, beoc550, beoc670, beoc870,                &
           besu440, besu500, besu550, besu670, besu870,                &
           babs440, babs500, babs550, babs670, babs870,                &
           bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1,             &
           beoc550lt1, beoc550gt1, besu550lt1, besu550gt1,             &
           backsc550, babg550, babc550, baoc550, basu550) 


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab,   only: cate, cat, fac, faq, fbc, rh, nbmodes, nmodes 

   implicit none

#include <aerocopt.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: rhum(pcols,pver)         ! level relative humidity (fraction)
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: xfacin(pcols,pver,nbmodes) ! modal (OC+BC)/(SO4+BC+OC)
   real(r8), intent(in) :: xfbcin(pcols,pver,nbmodes) ! modal BC/(OC+BC)
   real(r8), intent(in) :: xfaqin(pcols,pver,nbmodes) ! modal SO4(aq)/SO4
!
! Output arguments: Modal total and absorption extiction coefficients (for AeroCom)
! for 550nm (1) and 865nm (2), and for r<1um (lt1) and r>1um (gt1).
!
   real(r8), intent(out) :: &
     bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes), &
     bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes), &
     bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes), &
     bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes), &
     bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes), &
     bebg440(pcols,pver,0:nbmodes), & ! babg440(pcols,pver,0:nbmodes), &
     bebg500(pcols,pver,0:nbmodes), & ! babg500(pcols,pver,0:nbmodes), &
     bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes), &
     bebg670(pcols,pver,0:nbmodes), & ! babg670(pcols,pver,0:nbmodes), &
     bebg870(pcols,pver,0:nbmodes), & ! babg870(pcols,pver,0:nbmodes), &
     bebc440(pcols,pver,0:nbmodes), & ! babc440(pcols,pver,0:nbmodes), &
     bebc500(pcols,pver,0:nbmodes), & ! babc500(pcols,pver,0:nbmodes), &
     bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes), &
     bebc670(pcols,pver,0:nbmodes), & ! babc670(pcols,pver,0:nbmodes), &
     bebc870(pcols,pver,0:nbmodes), & ! babc870(pcols,pver,0:nbmodes), &
     beoc440(pcols,pver,0:nbmodes), & ! baoc440(pcols,pver,0:nbmodes), &
     beoc500(pcols,pver,0:nbmodes), & ! baoc500(pcols,pver,0:nbmodes), &
     beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes), &
     beoc670(pcols,pver,0:nbmodes), & ! baoc670(pcols,pver,0:nbmodes), &
     beoc870(pcols,pver,0:nbmodes), & ! baoc870(pcols,pver,0:nbmodes), &
     besu440(pcols,pver,0:nbmodes), & ! basu440(pcols,pver,0:nbmodes), &
     besu500(pcols,pver,0:nbmodes), & ! basu500(pcols,pver,0:nbmodes), &
     besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes), &
     besu670(pcols,pver,0:nbmodes), & ! basu670(pcols,pver,0:nbmodes), &
     besu870(pcols,pver,0:nbmodes), & ! basu870(pcols,pver,0:nbmodes), &
     bebg550lt1(pcols,pver,0:nbmodes), bebg550gt1(pcols,pver,0:nbmodes), &
     bebc550lt1(pcols,pver,0:nbmodes), bebc550gt1(pcols,pver,0:nbmodes), &
     beoc550lt1(pcols,pver,0:nbmodes), beoc550gt1(pcols,pver,0:nbmodes), &
     besu550lt1(pcols,pver,0:nbmodes), besu550gt1(pcols,pver,0:nbmodes), &
     backsc550(pcols,pver,0:nbmodes)
!
!---------------------------Local variables-----------------------------
!
      real(r8) a, b, e, eps
      real(r8) xct(pcols,pver), xrh(pcols,pver), xfac(pcols,pver,nbmodes), &
        xfbc(pcols,pver,nbmodes), xfaq(pcols,pver,nbmodes)
      real(r8) &
       arr111, arr112, arr121, arr122, arr211, arr212, arr221, arr222, &
!       arr11, arr12, arr21, arr22, arre1, arre2, arr1(50)  
       arr11, arr12, arr21, arr22, arre1, arre2, arr1(38)  
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

      integer i, iv, ierr, irelh, ictot, ifac, ifbc, ifaq, kcomp, k, icol
      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver),&
       ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
       ifbc1(pcols,pver), ifbc2(pcols,pver), ifaq1(pcols,pver),    &
       ifaq2(pcols,pver)

!      Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
      integer t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fbc1, t_fbc2, t_xfbc
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xrh, t_xct, t_rh1, t_rh2
      real(r8)    t_cat1, t_cat2
      real(r8) esssf1, esssf2, esssf3, esssf4, esssf5, esssf6, esssf7, &
       esssf8, esssf9, esssf10, esssf3x6x9x10

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before xrh-loop'
      do k=1,pver
        do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
        end do
      end do

!      write(*,*) 'Before rh-loop', kcomp
      do irelh=1,9
        do k=1,pver
          do icol=1,ncol
           if(xrh(icol,k).ge.rh(irelh).and. &
             xrh(icol,k).le.rh(irelh+1)) then
             irh1(icol,k)=irelh
             irh2(icol,k)=irelh+1
           endif
          end do
        end do
      end do

!      write(*,*) 'Before kcomp-loop'

!       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

        do kcomp=5,10        

!        write(*,*) 'kcomp = ', kcomp 

!     initialize all output fields
      do k=1,pver
        do icol=1,ncol
         bext440(icol,k,kcomp)=0.0_r8 
         babs440(icol,k,kcomp)=0.0_r8 
         bext500(icol,k,kcomp)=0.0_r8 
         babs500(icol,k,kcomp)=0.0_r8 
         bext550(icol,k,kcomp)=0.0_r8 
         babs550(icol,k,kcomp)=0.0_r8 
         bext670(icol,k,kcomp)=0.0_r8 
         babs670(icol,k,kcomp)=0.0_r8 
         bext870(icol,k,kcomp)=0.0_r8 
         babs870(icol,k,kcomp)=0.0_r8 
         bebg440(icol,k,kcomp)=0.0_r8 
!         babg440(icol,k,kcomp)=0.0_r8 
         bebg500(icol,k,kcomp)=0.0_r8 
!         babg500(icol,k,kcomp)=0.0_r8 
         bebg550(icol,k,kcomp)=0.0_r8 
         babg550(icol,k,kcomp)=0.0_r8 
         bebg670(icol,k,kcomp)=0.0_r8 
!         babg670(icol,k,kcomp)=0.0_r8 
         bebg870(icol,k,kcomp)=0.0_r8 
!         babg870(icol,k,kcomp)=0.0_r8 
         bebc440(icol,k,kcomp)=0.0_r8 
!         babc440(icol,k,kcomp)=0.0_r8 
         bebc500(icol,k,kcomp)=0.0_r8 
!         babc500(icol,k,kcomp)=0.0_r8 
         bebc550(icol,k,kcomp)=0.0_r8 
         babc550(icol,k,kcomp)=0.0_r8 
         bebc670(icol,k,kcomp)=0.0_r8 
!         babc670(icol,k,kcomp)=0.0_r8 
         bebc870(icol,k,kcomp)=0.0_r8 
!         babc870(icol,k,kcomp)=0.0_r8 
         beoc440(icol,k,kcomp)=0.0_r8 
!         baoc440(icol,k,kcomp)=0.0_r8 
         beoc500(icol,k,kcomp)=0.0_r8 
!         baoc500(icol,k,kcomp)=0.0_r8 
         beoc550(icol,k,kcomp)=0.0_r8 
         baoc550(icol,k,kcomp)=0.0_r8 
         beoc670(icol,k,kcomp)=0.0_r8 
!         baoc670(icol,k,kcomp)=0.0_r8 
         beoc870(icol,k,kcomp)=0.0_r8 
!         baoc870(icol,k,kcomp)=0.0_r8 
         besu440(icol,k,kcomp)=0.0_r8 
!         basu440(icol,k,kcomp)=0.0_r8 
         besu500(icol,k,kcomp)=0.0_r8 
!         basu500(icol,k,kcomp)=0.0_r8 
         besu550(icol,k,kcomp)=0.0_r8 
         basu550(icol,k,kcomp)=0.0_r8 
         besu670(icol,k,kcomp)=0.0_r8 
!         basu670(icol,k,kcomp)=0.0_r8 
         besu870(icol,k,kcomp)=0.0_r8 
!         basu870(icol,k,kcomp)=0.0_r8 
         bebg550lt1(icol,k,kcomp)=0.0_r8 
         bebg550gt1(icol,k,kcomp)=0.0_r8 
         bebc550lt1(icol,k,kcomp)=0.0_r8 
         bebc550gt1(icol,k,kcomp)=0.0_r8 
         beoc550lt1(icol,k,kcomp)=0.0_r8 
         beoc550gt1(icol,k,kcomp)=0.0_r8 
         besu550lt1(icol,k,kcomp)=0.0_r8 
         besu550gt1(icol,k,kcomp)=0.0_r8 
         backsc550(icol,k,kcomp)=0.0_r8 
        end do
      end do

!      write(*,*) 'Before x-loop'
      do k=1,pver
         do icol=1,ncol
          if(Nnatk(icol,k,kcomp).gt.0) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cat(kcomp,1)),cat(kcomp,6))
          xfac(icol,k,kcomp) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfbc(icol,k,kcomp) = min(max(xfbcin(icol,k,kcomp),fbc(1)),fbc(6))
          xfaq(icol,k,kcomp) = min(max(xfaqin(icol,k,kcomp),faq(1)),faq(6))

      do ictot=1,5
            if(xct(icol,k).ge.cat(kcomp,ictot).and. &
            xct(icol,k).le.cat(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

      do ifac=1,5
            if(xfac(icol,k,kcomp).ge.fac(ifac).and. &
             xfac(icol,k,kcomp).le.fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

      do ifbc=1,5
            if(xfbc(icol,k,kcomp).ge.fbc(ifbc).and. &
             xfbc(icol,k,kcomp).le.fbc(ifbc+1)) then
             ifbc1(icol,k)=ifbc
             ifbc2(icol,k)=ifbc+1
            endif
      end do ! ifbc

      do ifaq=1,5
            if(xfaq(icol,k,kcomp).ge.faq(ifaq).and. &
            xfaq(icol,k,kcomp).le.faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp).gt.0) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifb1 = ifbc1(icol,k)
      t_ifb2 = ifbc2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cat(kcomp,t_ict1)
      t_cat2 = cat(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_fbc1 = fbc(t_ifb1)
      t_fbc2 = fbc(t_ifb2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)

      t_xrh  = xrh(icol,k)
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

!         do iv=1,54  ! variable number
         do iv=1,38  ! variable number

      arr11111=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr11112=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr11121=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr11122=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr11211=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr11212=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr11221=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr11222=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr12111=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr12112=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr12121=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr12122=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr12211=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr12212=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr12221=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr12222=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr21111=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr21112=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr21121=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr21122=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr21211=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr21212=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr21221=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr21222=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      arr22111=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      arr22112=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      arr22121=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      arr22122=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      arr22211=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      arr22212=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      arr22221=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      arr22222=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

      arr1111=esssf1*arr11111+esssf2*arr11112
      arr1112=esssf1*arr11121+esssf2*arr11122
      arr1121=esssf1*arr11211+esssf2*arr11212
      arr1122=esssf1*arr11221+esssf2*arr11222
      arr1211=esssf1*arr12111+esssf2*arr12112
      arr1212=esssf1*arr12121+esssf2*arr12122
      arr1221=esssf1*arr12211+esssf2*arr12212
      arr1222=esssf1*arr12221+esssf2*arr12222
      arr2111=esssf1*arr21111+esssf2*arr21112
      arr2112=esssf1*arr21121+esssf2*arr21122
      arr2121=esssf1*arr21211+esssf2*arr21212
      arr2122=esssf1*arr21221+esssf2*arr21222
      arr2211=esssf1*arr22111+esssf2*arr22112
      arr2212=esssf1*arr22121+esssf2*arr22122
      arr2221=esssf1*arr22211+esssf2*arr22212
      arr2222=esssf1*arr22221+esssf2*arr22222

      arr111=esssf4*arr1111+esssf5*arr1112
      arr112=esssf4*arr1121+esssf5*arr1122
      arr121=esssf4*arr1211+esssf5*arr1212
      arr122=esssf4*arr1221+esssf5*arr1222
      arr211=esssf4*arr2111+esssf5*arr2112
      arr212=esssf4*arr2121+esssf5*arr2122
      arr221=esssf4*arr2211+esssf5*arr2212
      arr222=esssf4*arr2221+esssf5*arr2222

      arr11 =esssf7*arr111+esssf8*arr112
      arr12 =esssf7*arr121+esssf8*arr122
      arr21 =esssf7*arr211+esssf8*arr212
      arr22 =esssf7*arr221+esssf8*arr222

      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf3x6x9x10
      arre2  =((t_cat2-t_xct)*arr21+(t_xct-t_cat1)*arr22)*esssf3x6x9x10   

!      write(*,*) 'Before array'
!old      arr1=((t_rh2-t_xrh)*arre1+(t_xrh-t_rh1)*arre2) &
!old                         /(t_rh2-t_rh1)    
      arr1(iv)=((t_rh2-t_xrh)*arre1+(t_xrh-t_rh1)*arre2) &
                         /(t_rh2-t_rh1)    
!      write(*,*) arr1

!         end do ! iv=1,54 
         end do ! iv=1,38 
 
         bext440(icol,k,kcomp)=arr1(1)
         bext500(icol,k,kcomp)=arr1(2)
         bext670(icol,k,kcomp)=arr1(3)
         bext870(icol,k,kcomp)=arr1(4)
         bebg440(icol,k,kcomp)=arr1(5)
         bebg500(icol,k,kcomp)=arr1(6)
         bebg670(icol,k,kcomp)=arr1(7)
         bebg870(icol,k,kcomp)=arr1(8)
         bebc440(icol,k,kcomp)=arr1(9)
         bebc500(icol,k,kcomp)=arr1(10)
         bebc670(icol,k,kcomp)=arr1(11)
         bebc870(icol,k,kcomp)=arr1(12)
         beoc440(icol,k,kcomp)=arr1(13)
         beoc500(icol,k,kcomp)=arr1(14)
         beoc670(icol,k,kcomp)=arr1(15)
         beoc870(icol,k,kcomp)=arr1(16)
         besu440(icol,k,kcomp)=arr1(17)
         besu500(icol,k,kcomp)=arr1(18)
         besu670(icol,k,kcomp)=arr1(19)
         besu870(icol,k,kcomp)=arr1(20)
         babs440(icol,k,kcomp)=arr1(21)
         babs500(icol,k,kcomp)=arr1(22)
         babs550(icol,k,kcomp)=arr1(23)
         babs670(icol,k,kcomp)=arr1(24)
         babs870(icol,k,kcomp)=arr1(25)
         bebg550lt1(icol,k,kcomp)=arr1(26)
         bebg550gt1(icol,k,kcomp)=arr1(27)
         bebc550lt1(icol,k,kcomp)=arr1(28)
         bebc550gt1(icol,k,kcomp)=arr1(29)
         beoc550lt1(icol,k,kcomp)=arr1(30)
         beoc550gt1(icol,k,kcomp)=arr1(31)
         besu550lt1(icol,k,kcomp)=arr1(32)
         besu550gt1(icol,k,kcomp)=arr1(33)
         backsc550(icol,k,kcomp)=arr1(34)
         babg550(icol,k,kcomp)=arr1(35)
         babc550(icol,k,kcomp)=arr1(36)
         baoc550(icol,k,kcomp)=arr1(37)
         basu550(icol,k,kcomp)=arr1(38)
         bebg550(icol,k,kcomp)=arr1(26)+arr1(27)
         bebc550(icol,k,kcomp)=arr1(28)+arr1(29)
         beoc550(icol,k,kcomp)=arr1(30)+arr1(31)
         besu550(icol,k,kcomp)=arr1(32)+arr1(33)
         bext550(icol,k,kcomp)=bebg550(icol,k,kcomp)+bebc550(icol,k,kcomp) &
                              +beoc550(icol,k,kcomp)+besu550(icol,k,kcomp)

           endif
         
!     if(beocgt1(icol,k,kcomp)>beoc1(icol,k,kcomp)) then
!       write(*,*) '5to10,kcomp,beocgt1,beoc1=', kcomp, beocgt1(icol,k,kcomp), beoc1(icol,k,kcomp)  
!     endif

       end do ! icol
      end do ! k

        end do  ! kcomp

      return
end subroutine intaeropt5to10




