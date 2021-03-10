subroutine intdrypar1to3 (lchnk, ncol, Nnatk, Camk,                  & 
           cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, & 
           cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
           cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol, &
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
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
!
! Input-Output arguments
!
   real(r8), intent(inout) :: &
     aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes), &
     vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes), &
     cknorm(pcols,pver,0:nmodes), cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
!
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
      real(r8) a, b, e, eps, catot, xct(pcols,pver)

      real(r8) arr11, arr12, arr21, arr22, arre1, arre2  

      integer i, ierr, ictot, kcomp, k, icol
      integer ict1(pcols,pver), ict2(pcols,pver)
!      Temporary storage of often used array elements
      integer t_ict1, t_ict2
      real(r8) t_xct,  t_cat1, t_cat2
      real(r8) esssf10

      parameter (e=2.718281828_r8, eps=1.0e-60_r8)


!      write(*,*) 'Before kcomp-loop'

!       Modes 1-3,  SO4(Ait), BC(Ait) and OC(Ait):

        do kcomp=1,3

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

!      write(*,*) 'Before cat-loop', kcomp
      do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k)<=cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot
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
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_xct  = xct(icol,k)

      esssf10= 1.0_r8/(t_cat2-t_cat1)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

         do i=1,19  ! variable number

       if(i==1) then
      arr11=a1to3cintbg(t_ict1,kcomp)
      arr12=a1to3cintbg(t_ict2,kcomp)
       elseif(i==2) then
      arr11=a1to3cintbg05(t_ict1,kcomp)
      arr12=a1to3cintbg05(t_ict2,kcomp)
       elseif(i==3) then
      arr11=a1to3cintbg125(t_ict1,kcomp)
      arr12=a1to3cintbg125(t_ict2,kcomp)
       elseif(i>=4.and.i<=9) then
      arr11=eps
      arr12=eps
       elseif(i==10) then
      arr11=a1to3cintsc(t_ict1,kcomp)
      arr12=a1to3cintsc(t_ict2,kcomp)
       elseif(i==11) then
      arr11=a1to3cintsc05(t_ict1,kcomp)
      arr12=a1to3cintsc05(t_ict2,kcomp)
       elseif(i==12) then
      arr11=a1to3cintsc125(t_ict1,kcomp)
      arr12=a1to3cintsc125(t_ict2,kcomp)
       elseif(i>=13.and.i<=15) then
      arr11=eps
      arr12=eps
       elseif(i==16) then
      arr11=a1to3aaeros(t_ict1,kcomp)
      arr12=a1to3aaeros(t_ict2,kcomp)
       elseif(i==17) then
      arr11=a1to3aaerol(t_ict1,kcomp)
      arr12=a1to3aaerol(t_ict2,kcomp)
       elseif(i==18) then
      arr11=a1to3vaeros(t_ict1,kcomp)
      arr12=a1to3vaeros(t_ict2,kcomp)
       elseif(i==19) then
      arr11=a1to3vaerol(t_ict1,kcomp)
      arr12=a1to3vaerol(t_ict2,kcomp)
       endif

      arre1  =((t_cat2-t_xct)*arr11+(t_xct-t_cat1)*arr12)*esssf10

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
         
          cknorm(icol,k,kcomp)  = a1to3cintbg(1,kcomp)
          cknlt05(icol,k,kcomp) = a1to3cintbg05(1,kcomp)
          ckngt125(icol,k,kcomp)= a1to3cintbg125(1,kcomp)

       end do ! icol
      end do ! k

        end do  ! kcomp

!      Dry parameters for externally mixed modes modes 11-13,  
!      SO4(n), BC(n) and OC(n):

        do kcomp=11,13

        do k=1,pver 
          do icol=1,ncol
           cknorm(icol,k,kcomp)  = a1to3cintbg(1,kcomp-10)
           cknlt05(icol,k,kcomp) = a1to3cintbg05(1,kcomp-10)
           ckngt125(icol,k,kcomp)= a1to3cintbg125(1,kcomp-10)
           aaerosn(icol,k,kcomp) = a1to3aaeros(1,kcomp-10)
           aaeroln(icol,k,kcomp) = a1to3aaerol(1,kcomp-10)
           vaerosn(icol,k,kcomp) = a1to3vaeros(1,kcomp-10)
           vaeroln(icol,k,kcomp) = a1to3vaerol(1,kcomp-10)
         end do ! icol
        end do ! k

        end do  ! kcomp


      return
end subroutine intdrypar1to3
