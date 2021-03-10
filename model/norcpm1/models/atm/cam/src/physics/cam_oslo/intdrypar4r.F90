!#include <params.h>

subroutine intdrypar4r (lchnk, ncol, Nnatk, xfacin, xfacnin, & 
           aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)

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
   real(r8), intent(in) :: xfacin(pcols,pver)       ! Aitken BC/(BC+OC)
   real(r8), intent(in) :: xfacnin(pcols,pver)      ! n-mode BC/(BC+OC)
!
! Input-Output arguments
!
   real(r8), intent(inout) :: &
     cknorm(pcols,pver,0:nmodes), cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
   real(r8), intent(inout) :: &
            aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes), &
            vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes)
!
!
!---------------------------Local variables-----------------------------
!
      real(r8) a, b, e, eps
      real(r8) xfac(pcols,pver)

      real(r8) & 
        arr111, arr112, arr121, arr122, arr211, arr212, arr221, arr222, &
        arr11, arr12, arr21, arr22, arre1, arre2  

      integer i, ierr, ictot, ifac, kcomp, k, icol
!feil!      integer ict1(pcols,pver),  xfac(pcols,pver), &
      integer ict1(pcols,pver), ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver)
!      Temporary storage of often used array elements
      integer t_ifc1, t_ifc2
      real(r8) t_fac1, t_fac2, t_xfac
      real(r8) esssf7, esssf8, esssf9, esssf10

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!       Mode 4, BC&OC(Ait):

        kcomp=4

      do k=1,pver
         do icol=1,ncol

          if(Nnatk(icol,k,kcomp)>0.0_r8) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          xfac(icol,k) = min(max(xfacin(icol,k),fac(1)),fac(6))

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k)<=fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp)>0.0_r8) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_xfac = xfac(icol,k)

      esssf7 = (t_fac2-t_xfac)
      esssf8 = (t_xfac-t_fac1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do i=1,3  ! variable number

       if(i==1) then
!4      arr111=a4cintbg(1,t_ifc1)
!4      arr112=a4cintbg(1,t_ifc2)
      arr111=a4cintbg(1,t_ifc1,1)
      arr112=a4cintbg(1,t_ifc2,1)
       elseif(i==2) then
!4      arr111=a4cintbg05(1,t_ifc1)
!4      arr112=a4cintbg05(1,t_ifc2)
      arr111=a4cintbg05(1,t_ifc1,1)
      arr112=a4cintbg05(1,t_ifc2,1)
       elseif(i==3) then
!4      arr111=a4cintbg125(1,t_ifc1)
!4      arr112=a4cintbg125(1,t_ifc2)
      arr111=a4cintbg125(1,t_ifc1,1)
      arr112=a4cintbg125(1,t_ifc2,1)
       endif

      arre1 =esssf7*arr111+esssf8*arr112

!      write(*,*) 'Before array'

       if(i==1) then
         cknorm(icol,k,kcomp)=arre1
       elseif(i==2) then
         cknlt05(icol,k,kcomp)=arre1
       elseif(i==3) then
         ckngt125(icol,k,kcomp)=arre1
       endif

         end do ! i=1,3 

           endif
         
       end do ! icol
      end do ! k


!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc


!       Mode 14, BC&OC(n):

        kcomp=14

!      initialize output fields
      do k=1,pver
         do icol=1,ncol
        aaerosn(icol,k,kcomp)=0.0_r8
        aaeroln(icol,k,kcomp)=0.0_r8
        vaerosn(icol,k,kcomp)=0.0_r8
        vaeroln(icol,k,kcomp)=0.0_r8
         end do
       end do

      do k=1,pver
         do icol=1,ncol

          if(Nnatk(icol,k,kcomp)>0.0_r8) then
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
          xfac(icol,k) = min(max(xfacnin(icol,k),fac(1)),fac(6))

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k)<=fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac
!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
           endif

          end do ! icol
        end do ! k


        do k=1,pver 
          do icol=1,ncol
         
           if(Nnatk(icol,k,kcomp)>0.0_r8) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing
      t_ifc1 = ifac1(icol,k)
      t_ifc2 = ifac2(icol,k)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_xfac = xfac(icol,k)

      esssf7 = (t_fac2-t_xfac)
      esssf8 = (t_xfac-t_fac1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do i=1,7  ! variable number

       if(i==1) then
!4      arr111=a4cintbg(1,t_ifc1)
!4      arr112=a4cintbg(1,t_ifc2)
      arr111=a4cintbg(1,t_ifc1,1)
      arr112=a4cintbg(1,t_ifc2,1)
       elseif(i==2) then
!4      arr111=a4cintbg05(1,t_ifc1)
!4      arr112=a4cintbg05(1,t_ifc2)
      arr111=a4cintbg05(1,t_ifc1,1)
      arr112=a4cintbg05(1,t_ifc2,1)
       elseif(i==3) then
!4      arr111=a4cintbg125(1,t_ifc1)
!4      arr112=a4cintbg125(1,t_ifc2)
      arr111=a4cintbg125(1,t_ifc1,1)
      arr112=a4cintbg125(1,t_ifc2,1)
       elseif(i==4) then
!4      arr111=a4aaeros(1,t_ifc1)
!4      arr112=a4aaeros(1,t_ifc2)
      arr111=a4aaeros(1,t_ifc1,1)
      arr112=a4aaeros(1,t_ifc2,1)
       elseif(i==5) then
!4      arr111=a4aaerol(1,t_ifc1)
!4      arr112=a4aaerol(1,t_ifc2)
      arr111=a4aaerol(1,t_ifc1,1)
      arr112=a4aaerol(1,t_ifc2,1)
       elseif(i==6) then
!4      arr111=a4vaeros(1,t_ifc1)
!4      arr112=a4vaeros(1,t_ifc2)
      arr111=a4vaeros(1,t_ifc1,1)
      arr112=a4vaeros(1,t_ifc2,1)
       elseif(i==7) then
!4      arr111=a4vaerol(1,t_ifc1)
!4      arr112=a4vaerol(1,t_ifc2)
      arr111=a4vaerol(1,t_ifc1,1)
      arr112=a4vaerol(1,t_ifc2,1)
       endif

      arre1 =esssf7*arr111+esssf8*arr112

!      write(*,*) 'Before array'

       if(i==1) then
         cknorm(icol,k,kcomp)=arre1
       elseif(i==2) then
         cknlt05(icol,k,kcomp)=arre1
       elseif(i==3) then
         ckngt125(icol,k,kcomp)=arre1
       elseif(i==4) then
        aaerosn(icol,k,kcomp)=arre1
       elseif(i==5) then
        aaeroln(icol,k,kcomp)=arre1
       elseif(i==6) then
        vaerosn(icol,k,kcomp)=arre1
       elseif(i==7) then
        vaeroln(icol,k,kcomp)=arre1
       endif

         end do ! i=1,7 

           endif
         
       end do ! icol
      end do ! k



      return
end subroutine intdrypar4r

