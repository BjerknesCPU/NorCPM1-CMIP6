subroutine intaeropt4 (lchnk, ncol, rhum, mplus10, Nnatk, Camk, xfacin, xfaqin, &
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
   integer, intent(in) :: mplus10                   ! mode number (0) or number + 10 (1)
   real(r8), intent(in) :: rhum(pcols,pver)         ! level relative humidity (fraction)
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed (cond+aq) SO4 conc.
   real(r8), intent(in) :: xfacin(pcols,pver)       ! BC/(BC+OC) for the background mode
   real(r8), intent(in) :: xfaqin(pcols,pver)       ! modal SO4(aq)/SO4
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
      real(r8) xct(pcols,pver), xrh(pcols,pver), xfac(pcols,pver), xfaq(pcols,pver)

      real(r8) &
       arr111, arr112, arr121, arr122, arr211, arr212, arr221, arr222, &
!       arr11, arr12, arr21, arr22, arre1, arre2, arr1(54)  
       arr11, arr12, arr21, arr22, arre1, arre2, arr1(38)  
      real(r8) &
       arr1111, arr1112, arr1121, arr1122, arr1211, arr1212, arr1221, &
       arr1222, arr2111, arr2112, arr2121, arr2122, arr2211, arr2212, &
       arr2221, arr2222

      integer i, iv, ierr, irelh, ictot, ifac, ifaq, kcomp, k, icol, kc10
      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver),&
       ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver), &
       ifaq1(pcols,pver), ifaq2(pcols,pver)

!      Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2,  t_ifa1, t_ifa2
      real(r8)    t_faq1, t_faq2, t_xfaq
      real(r8)    t_fac1, t_fac2, t_xfac
      real(r8)    t_xrh, t_xct, t_rh1, t_rh2
      real(r8)    t_cat1, t_cat2
      real(r8) esssf1, esssf2, esssf3, esssf7, esssf8, esssf9, esssf10

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

!       BC&OC(Ait) mode: 

        do kcomp=4,4        

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

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

!      write(*,*) 'Before x-loop'
      do k=1,pver
         do icol=1,ncol

!x          if(Nnatk(icol,k,kcomp).gt.0) then
          if(Nnatk(icol,k,kc10).gt.0) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
!x                 /(Nnatk(icol,k,kcomp)+eps),cate(kcomp,1)),cate(kcomp,16))
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k),fac(1)),fac(6))
          xfaq(icol,k) = min(max(xfaqin(icol,k),faq(1)),faq(6))

      do ictot=1,15
            if(xct(icol,k).ge.cate(kcomp,ictot).and. &
            xct(icol,k).le.cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

      do ifac=1,5
            if(xfac(icol,k).ge.fac(ifac).and. &
             xfac(icol,k).le.fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

      do ifaq=1,5
            if(xfaq(icol,k).ge.faq(ifaq).and. &
            xfaq(icol,k).le.faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
!x           if(Nnatk(icol,k,kcomp).gt.0) then
           if(Nnatk(icol,k,kc10).gt.0) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_ifa1 = ifaq1(icol,k)
      t_ifa2 = ifaq2(icol,k)

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)

      t_xrh  = xrh(icol,k)
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

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

!         do iv=1,54  ! variable number
         do iv=1,38  ! variable number

      arr1111=bep4(iv,t_irh1,t_ict1,t_ifc1,t_ifa1)
      arr1112=bep4(iv,t_irh1,t_ict1,t_ifc1,t_ifa2)
      arr1121=bep4(iv,t_irh1,t_ict1,t_ifc2,t_ifa1)
      arr1122=bep4(iv,t_irh1,t_ict1,t_ifc2,t_ifa2)
      arr1211=bep4(iv,t_irh1,t_ict2,t_ifc1,t_ifa1)
      arr1212=bep4(iv,t_irh1,t_ict2,t_ifc1,t_ifa2)
      arr1221=bep4(iv,t_irh1,t_ict2,t_ifc2,t_ifa1)
      arr1222=bep4(iv,t_irh1,t_ict2,t_ifc2,t_ifa2)
      arr2111=bep4(iv,t_irh2,t_ict1,t_ifc1,t_ifa1)
      arr2112=bep4(iv,t_irh2,t_ict1,t_ifc1,t_ifa2)
      arr2121=bep4(iv,t_irh2,t_ict1,t_ifc2,t_ifa1)
      arr2122=bep4(iv,t_irh2,t_ict1,t_ifc2,t_ifa2)
      arr2211=bep4(iv,t_irh2,t_ict2,t_ifc1,t_ifa1)
      arr2212=bep4(iv,t_irh2,t_ict2,t_ifc1,t_ifa2)
      arr2221=bep4(iv,t_irh2,t_ict2,t_ifc2,t_ifa1)
      arr2222=bep4(iv,t_irh2,t_ict2,t_ifc2,t_ifa2)

      arr111=esssf1*arr1111+esssf2*arr1112
      arr112=esssf1*arr1121+esssf2*arr1122
      arr121=esssf1*arr1211+esssf2*arr1212
      arr122=esssf1*arr1221+esssf2*arr1222
      arr211=esssf1*arr2111+esssf2*arr2112
      arr212=esssf1*arr2121+esssf2*arr2122
      arr221=esssf1*arr2211+esssf2*arr2212
      arr222=esssf1*arr2221+esssf2*arr2222

      arr11 =esssf7*arr111+esssf8*arr112
      arr12 =esssf7*arr121+esssf8*arr122
      arr21 =esssf7*arr211+esssf8*arr212
      arr22 =esssf7*arr221+esssf8*arr222

      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf3*esssf9*esssf10
      arre2  =((t_cat2-t_xct)*arr21+(t_xct-t_cat1)*arr22)*esssf3*esssf9*esssf10 

!      write(*,*) 'Before array'
!old      arr1=((t_rh2-t_xrh)*arre1+(t_xrh-t_rh1)*arre2) &
!old                         /(t_rh2-t_rh1)    
      arr1(iv)=((t_rh2-t_xrh)*arre1+(t_xrh-t_rh1)*arre2) &
                         /(t_rh2-t_rh1)    
!      write(*,*) arr1

!         end do ! iv=1,50 
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
         
       end do ! icol
      end do ! k

        end do  ! kcomp

      return
end subroutine intaeropt4




