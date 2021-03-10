module optinterpol

! Purpose: To read in look-up tables and calculate optical properties for the aerosols
!    For subroutine interpol:


  use shr_kind_mod, only: r8 => shr_kind_r8
  use opttab
  implicit none

  private 
  save

  public interpol0
  public interpol1to3
  public interpol4
  public interpol5to10  


 contains

subroutine interpol0 (lchnk, ncol, daylight, Nnatk, omega, gass, bex, ske)

!   Rewritten for CAM by Alf Kirkevaag in Feb. 2004, and 
!   modified for new aerosol schemes in January 2006.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   logical, intent(in) :: daylight(pcols)           ! only daylight calculations if .true.
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, kcomp, k, icol


      kcomp=0

        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kcomp,i)=0.0_r8
              gass(icol,k,kcomp,i)=0.0_r8
              bex(icol,k,kcomp,i)=0.0_r8
              ske(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol
!           if(Nnatk(icol,k,kcomp)>0.0_r8) then
           if(daylight(icol).and.Nnatk(icol,k,kcomp)>0.0_r8) then

      do i=1,nbands   ! i = wavelength index

        omega(icol,k,kcomp,i)=om0(i)
        gass(icol,k,kcomp,i)=g0(i) 
        bex(icol,k,kcomp,i)=be0(i)
        ske(icol,k,kcomp,i)=ke0(i)

      end do          ! i

           endif
         
          end do ! icol
        end do ! k

      return
end subroutine interpol0


subroutine interpol1to3 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, Camk, &
                         omega, gass, bex, ske, cxstot)

!     Optimized for speed by Arild Burud/NoSerC, June 2002
!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1.
!   Files also involved by this modification: initopt.F
!   and opttab.h. 
!---------------------------------------------------------------
!   Rewritten for CAM by Alf Kirkevaag in Feb. 2004, and 
!   modified (interpol1->interpol1to3) for new aerosol schemes 
!   in November 2005.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   integer, intent(in) :: mplus10                   ! mode number (0) or number + 10 (1)
   logical, intent(in) :: daylight(pcols)           ! only daylight calculations if .true.

   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: rhum(pcols,pver)         ! level relative humidity (fraction)
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)            ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, kcomp, k, icol, kc10
      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver), &
        ict2(pcols,pver)
      real(r8) a, b, relh, catot, camdiff
      real(r8) xct(pcols,pver), cxs(pcols,pver), xrh(pcols,pver)

!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2
      real(r8) t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
      real(r8) esssf10
      real(r8) om11, om12, om21, om22, ome1, ome2  
      real(r8) g11, g12, g21, g22, ge1, ge2 
      real(r8) be11, be12, be21, be22, bex1, bex2
      real(r8) ke11, ke12, ke21, ke22, ske1, ske2



!      write(*,*) 'Before xrh-loop'
      do k=1,pver
        do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
        end do 
      end do

!      write(*,*) 'Before rh-loop'
      do irelh=1,9
       do k=1,pver
        do icol=1,ncol
           if(xrh(icol,k) >= rh(irelh).and. &
             xrh(icol,k)<=rh(irelh+1)) then
             irh1(icol,k)=irelh
             irh2(icol,k)=irelh+1
           endif
         end do
       end do
      end do

      do k=1,pver
        do icol=1,ncol
          xct(icol,k)    = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
        do kcomp=1,3

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

      do k=1,pver
        do icol=1,ncol

!            if(Nnatk(icol,k,kc10)>0.0_r8) then
!g            if(daylight(icol).and.Nnatk(icol,k,kc10)>0.0_r8) then
            if(daylight(icol)) then

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kc10)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kcomp, cxs =', k, kc10, cxs(icol,k)
!          write(*,*) 'k, kcomp, xct =', k, kc10, xct(icol,k)
!        endif

!      write(*,*) 'Before cate-loop', kc10
      do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k) <= cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

            endif ! daylight

          end do ! icol
        end do ! k


!      write(*,*) 'Before init-loop', kc10    ! Bruk For All istedet?
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kc10,i)=0.0_r8
              gass(icol,k,kc10,i)=0.0_r8
              bex(icol,k,kc10,i)=0.0_r8
              ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!           if(Nnatk(icol,k,kc10)>0.0_r8) then
!g           if(daylight(icol).and.Nnatk(icol,k,kc10)>0.0_r8) then
           if(daylight(icol)) then

!      Collect all the vector elements into temporary storage
!      to avoid cache conflicts and excessive cross-referencing

      t_irh1 = irh1(icol,k)
      t_irh2 = irh2(icol,k)
      t_ict1 = ict1(icol,k)
      t_ict2 = ict2(icol,k)

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2

      t_xrh  = xrh(icol,k)
      t_xct  = xct(icol,k)

      esssf10= 1.0_r8/(t_cat2-t_cat1)


      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

      om11=om1to3(i,t_irh1,t_ict1,kcomp)
      om12=om1to3(i,t_irh1,t_ict2,kcomp)
      om21=om1to3(i,t_irh2,t_ict1,kcomp)
      om22=om1to3(i,t_irh2,t_ict2,kcomp)

      ome1  =((t_cat2-t_xct)*om11+(t_xct-t_cat1)*om12)*esssf10
      ome2  =((t_cat2-t_xct)*om21+(t_xct-t_cat1)*om22)*esssf10

!      write(*,*) 'Before omega'
      omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) omega(icol,k,kc10,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

      g11=g1to3(i,t_irh1,t_ict1,kcomp)
      g12=g1to3(i,t_irh1,t_ict2,kcomp)
      g21=g1to3(i,t_irh2,t_ict1,kcomp)
      g22=g1to3(i,t_irh2,t_ict2,kcomp)

      ge1  =((t_cat2-t_xct)*g11+(t_xct-t_cat1)*g12)*esssf10
      ge2  =((t_cat2-t_xct)*g21+(t_xct-t_cat1)*g22)*esssf10

!      write(*,*) 'Before gass'
      gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

      be11=be1to3(i,t_irh1,t_ict1,kcomp)
      be12=be1to3(i,t_irh1,t_ict2,kcomp)
      be21=be1to3(i,t_irh2,t_ict1,kcomp)
      be22=be1to3(i,t_irh2,t_ict2,kcomp)

      bex1  =((t_cat2-t_xct)*be11+(t_xct-t_cat1)*be12)*esssf10      
      bex2  =((t_cat2-t_xct)*be21+(t_xct-t_cat1)*be22)*esssf10   
!osdec05
      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)
!osdec05

!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
       bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i


      do i=8,8            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

      ke11=ke1to3(i,t_irh1,t_ict1,kcomp)
      ke12=ke1to3(i,t_irh1,t_ict2,kcomp)
      ke21=ke1to3(i,t_irh2,t_ict1,kcomp)
      ke22=ke1to3(i,t_irh2,t_ict2,kcomp)

      ske1  =((t_cat2-t_xct)*ke11+(t_xct-t_cat1)*ke12)*esssf10
      ske2  =((t_cat2-t_xct)*ke21+(t_xct-t_cat1)*ke22)*esssf10
!osdec05
      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)
!osdec05

!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

           endif
         
!     write(*,*) 'kc10, bex =', icol,k,kc10,bex(icol,k,kc10,11) 

          end do ! icol
        end do ! k
        end do  ! kcomp

      return
end subroutine interpol1to3




subroutine interpol4 (lchnk, ncol, daylight, rhum, mplus10, Nnatk, Camk, &
                      xfacin, xfaqin, omega, gass, bex, ske, cxstot)

!     Optimized for speed by Arild Burud/NoSerC, June 2002
!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1.
!   Files also involved by this modification: initopt.F
!   and opttab.h. 
!---------------------------------------------------------------
!   Rewritten for CAM by Alf Kirkevaag in Feb. 2004, and 
!   modified (interpol1->interpol4) for new aerosol schemes 
!   in November 2005.


   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                     ! chunk identifier
   integer, intent(in) :: ncol                      ! number of atmospheric columns
   integer, intent(in) :: mplus10                   ! mode number (0) or number + 10 (1)
   logical, intent(in) :: daylight(pcols)           ! only daylight calculations if .true.

   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes) ! modal internally mixed (cond+aq) SO4 conc.
   real(r8), intent(in) :: rhum(pcols,pver)         ! level relative humidity (fraction)
   real(r8), intent(in) :: xfacin(pcols,pver)       ! BC/(BC+OC) for the background mode 
   real(r8), intent(in) :: xfaqin(pcols,pver)       ! = Cso4a2/(Cso4a1+Cso4a2)
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)    ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, ifaq, kcomp, k, icol, kc10
      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver), &
        ifaq1(pcols,pver), ifaq2(pcols,pver)
      real(r8) a, b, relh, catot, fraq, camdiff
      real(r8) xrh(pcols,pver), xct(pcols,pver), xfac(pcols,pver), &
        xfaq(pcols,pver), cxs(pcols,pver)
 
!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifa1, t_ifa2
      real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, &
        t_cat1, t_cat2, t_faq1, t_faq2, t_xfaq
      real(r8) esssf1, esssf2, esssf3, esssf7, esssf8, esssf9, esssf10
      real(r8) om111, om112, om121, om122, om211, om212, om221, om222, &
        om11, om12, om21, om22, ome1, ome2  
      real(r8) g111, g112, g121, g122, g211, g212, g221, g222, &
        g11, g12, g21, g22, ge1, ge2 
      real(r8) be111, be112, be121, be122, be211, be212, be221, be222, &
        be11, be12, be21, be22, bex1, bex2
      real(r8) ke111, ke112, ke121, ke122, ke211, ke212, &
        ke221, ke222, ke11, ke12, ke21, ke22, ske1, ske2
      real(r8) &
        om1111, om1112, om1121, om1122, om1211, om1212, om1221,om1222,  &
        om2111, om2112, om2121, om2122, om2211, om2212, om2221, om2222, & 
        g1111, g1112, g1121, g1122, g1211, g1212, g1221, g1222,         &
        g2111, g2112, g2121, g2122, g2211, g2212, g2221, g2222,         &
        be1111, be1112, be1121, be1122, be1211, be1212, be1221, be1222, & 
        be2111, be2112, be2121, be2122, be2211, be2212, be2221, be2222, &
        ke1111, ke1112, ke1121, ke1122, ke1211, ke1212, ke1221, ke1222, &
        ke2111, ke2112, ke2121, ke2122, ke2211, ke2212, ke2221, ke2222 



!      write(*,*) 'Before xrh-loop'
      do k=1,pver
        do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
        end do 
      end do

!      write(*,*) 'Before rh-loop'
      do irelh=1,9
       do k=1,pver
        do icol=1,ncol
           if(xrh(icol,k) >= rh(irelh).and. &
             xrh(icol,k)<=rh(irelh+1)) then
             irh1(icol,k)=irelh
             irh2(icol,k)=irelh+1
           endif
         end do
       end do
      end do

      do k=1,pver
        do icol=1,ncol
          xct(icol,k)    = 0.0_r8
          xfac(icol,k)   = 0.0_r8
          xfaq(icol,k)   = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
        do kcomp=4,4

           if(mplus10==0) then
             kc10=kcomp
           else
             kc10=kcomp+10
           endif 

      do k=1,pver
        do icol=1,ncol

!            if(Nnatk(icol,k,kc10)>0.0_r8) then
!g            if(daylight(icol).and.Nnatk(icol,k,kc10)>0.0_r8) then
            if(daylight(icol)) then

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kc10)+eps),cate(kcomp,1)),cate(kcomp,16))
          xfac(icol,k) = min(max(xfacin(icol,k),fac(1)),fac(6))
          xfaq(icol,k) = min(max(xfaqin(icol,k),faq(1)),faq(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kc10)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kc10, cxs =', k, kc10, cxs(icol,k)
!        endif

!      write(*,*) 'Before cate-loop', kc10
      do ictot=1,15
            if(xct(icol,k)>=cate(kcomp,ictot).and. &
            xct(icol,k) <= cate(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k)>=fac(ifac).and. &
             xfac(icol,k) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

!      write(*,*) 'Before faq-loop', kcomp
      do ifaq=1,5
            if(xfaq(icol,k) >= faq(ifaq).and. &
            xfaq(icol,k) <= faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq

            endif ! daylight

          end do ! icol
        end do ! k


!      write(*,*) 'Before init-loop', kc10
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kc10,i)=0.0_r8
              gass(icol,k,kc10,i)=0.0_r8
              bex(icol,k,kc10,i)=0.0_r8
              ske(icol,k,kc10,i)=0.0_r8
            end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!           if(Nnatk(icol,k,kc10)>0.0_r8) then
!g           if(daylight(icol).and.Nnatk(icol,k,kc10)>0.0_r8) then
           if(daylight(icol)) then

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

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

      t_rh1  = rh(t_irh1)
      t_rh2  = rh(t_irh2)
      t_cat1 = cate(kcomp,t_ict1)
      t_cat2 = cate(kcomp,t_ict2)
      t_fac1 = fac(t_ifc1)
      t_fac2 = fac(t_ifc2)
      t_faq1 = faq(t_ifa1)
      t_faq2 = faq(t_ifa2)

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2=',t_fac1,t_fac2

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


      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

      om1111=om4(i,t_irh1,t_ict1,t_ifc1,t_ifa1)
      om1112=om4(i,t_irh1,t_ict1,t_ifc1,t_ifa2)
      om1121=om4(i,t_irh1,t_ict1,t_ifc2,t_ifa1)
      om1122=om4(i,t_irh1,t_ict1,t_ifc2,t_ifa2)
      om1211=om4(i,t_irh1,t_ict2,t_ifc1,t_ifa1)
      om1212=om4(i,t_irh1,t_ict2,t_ifc1,t_ifa2)
      om1221=om4(i,t_irh1,t_ict2,t_ifc2,t_ifa1)
      om1222=om4(i,t_irh1,t_ict2,t_ifc2,t_ifa2)
      om2111=om4(i,t_irh2,t_ict1,t_ifc1,t_ifa1)
      om2112=om4(i,t_irh2,t_ict1,t_ifc1,t_ifa2)
      om2121=om4(i,t_irh2,t_ict1,t_ifc2,t_ifa1)
      om2122=om4(i,t_irh2,t_ict1,t_ifc2,t_ifa2)
      om2211=om4(i,t_irh2,t_ict2,t_ifc1,t_ifa1)
      om2212=om4(i,t_irh2,t_ict2,t_ifc1,t_ifa2)
      om2221=om4(i,t_irh2,t_ict2,t_ifc2,t_ifa1)
      om2222=om4(i,t_irh2,t_ict2,t_ifc2,t_ifa2)

      om111=esssf1*om1111+esssf2*om1112
      om112=esssf1*om1121+esssf2*om1122
      om121=esssf1*om1211+esssf2*om1212
      om122=esssf1*om1221+esssf2*om1222
      om211=esssf1*om2111+esssf2*om2112
      om212=esssf1*om2121+esssf2*om2122
      om221=esssf1*om2211+esssf2*om2212
      om222=esssf1*om2221+esssf2*om2222

      om11 =esssf7*om111+esssf8*om112
      om12 =esssf7*om121+esssf8*om122
      om21 =esssf7*om211+esssf8*om212
      om22 =esssf7*om221+esssf8*om222

      ome1  =((t_cat2-t_xct)*om11+(t_xct-t_cat1)*om12)*esssf3*esssf9*esssf10
      ome2  =((t_cat2-t_xct)*om21+(t_xct-t_cat1)*om22)*esssf3*esssf9*esssf10

!      write(*,*) 'Before omega'
      omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

      g1111=g4(i,t_irh1,t_ict1,t_ifc1,t_ifa1)
      g1112=g4(i,t_irh1,t_ict1,t_ifc1,t_ifa2)
      g1121=g4(i,t_irh1,t_ict1,t_ifc2,t_ifa1)
      g1122=g4(i,t_irh1,t_ict1,t_ifc2,t_ifa2)
      g1211=g4(i,t_irh1,t_ict2,t_ifc1,t_ifa1)
      g1212=g4(i,t_irh1,t_ict2,t_ifc1,t_ifa2)
      g1221=g4(i,t_irh1,t_ict2,t_ifc2,t_ifa1)
      g1222=g4(i,t_irh1,t_ict2,t_ifc2,t_ifa2)
      g2111=g4(i,t_irh2,t_ict1,t_ifc1,t_ifa1)
      g2112=g4(i,t_irh2,t_ict1,t_ifc1,t_ifa2)
      g2121=g4(i,t_irh2,t_ict1,t_ifc2,t_ifa1)
      g2122=g4(i,t_irh2,t_ict1,t_ifc2,t_ifa2)
      g2211=g4(i,t_irh2,t_ict2,t_ifc1,t_ifa1)
      g2212=g4(i,t_irh2,t_ict2,t_ifc1,t_ifa2)
      g2221=g4(i,t_irh2,t_ict2,t_ifc2,t_ifa1)
      g2222=g4(i,t_irh2,t_ict2,t_ifc2,t_ifa2)

      g111=esssf1*g1111+esssf2*g1112
      g112=esssf1*g1121+esssf2*g1122
      g121=esssf1*g1211+esssf2*g1212
      g122=esssf1*g1221+esssf2*g1222
      g211=esssf1*g2111+esssf2*g2112
      g212=esssf1*g2121+esssf2*g2122
      g221=esssf1*g2211+esssf2*g2212
      g222=esssf1*g2221+esssf2*g2222

      g11 =esssf7*g111+esssf8*g112
      g12 =esssf7*g121+esssf8*g122
      g21 =esssf7*g211+esssf8*g212
      g22 =esssf7*g221+esssf8*g222

      ge1  =((t_cat2-t_xct)*g11+(t_xct-t_cat1)*g12)*esssf3*esssf9*esssf10
      ge2  =((t_cat2-t_xct)*g21+(t_xct-t_cat1)*g22)*esssf3*esssf9*esssf10

!      write(*,*) 'Before gass'
      gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

      be1111=be4(i,t_irh1,t_ict1,t_ifc1,t_ifa1)
      be1112=be4(i,t_irh1,t_ict1,t_ifc1,t_ifa2)
      be1121=be4(i,t_irh1,t_ict1,t_ifc2,t_ifa1)
      be1122=be4(i,t_irh1,t_ict1,t_ifc2,t_ifa2)
      be1211=be4(i,t_irh1,t_ict2,t_ifc1,t_ifa1)
      be1212=be4(i,t_irh1,t_ict2,t_ifc1,t_ifa2)
      be1221=be4(i,t_irh1,t_ict2,t_ifc2,t_ifa1)
      be1222=be4(i,t_irh1,t_ict2,t_ifc2,t_ifa2)
      be2111=be4(i,t_irh2,t_ict1,t_ifc1,t_ifa1)
      be2112=be4(i,t_irh2,t_ict1,t_ifc1,t_ifa2)
      be2121=be4(i,t_irh2,t_ict1,t_ifc2,t_ifa1)
      be2122=be4(i,t_irh2,t_ict1,t_ifc2,t_ifa2)
      be2211=be4(i,t_irh2,t_ict2,t_ifc1,t_ifa1)
      be2212=be4(i,t_irh2,t_ict2,t_ifc1,t_ifa2)
      be2221=be4(i,t_irh2,t_ict2,t_ifc2,t_ifa1)
      be2222=be4(i,t_irh2,t_ict2,t_ifc2,t_ifa2)

      be111=esssf1*be1111+esssf2*be1112
      be112=esssf1*be1121+esssf2*be1122
      be121=esssf1*be1211+esssf2*be1212
      be122=esssf1*be1221+esssf2*be1222
      be211=esssf1*be2111+esssf2*be2112
      be212=esssf1*be2121+esssf2*be2122
      be221=esssf1*be2211+esssf2*be2212
      be222=esssf1*be2221+esssf2*be2222

      be11 =esssf7*be111+esssf8*be112
      be12 =esssf7*be121+esssf8*be122
      be21 =esssf7*be211+esssf8*be212
      be22 =esssf7*be221+esssf8*be222

      bex1  =((t_cat2-t_xct)*be11+(t_xct-t_cat1)*be12)*esssf3*esssf9*esssf10      
      bex2  =((t_cat2-t_xct)*be21+(t_xct-t_cat1)*be22)*esssf3*esssf9*esssf10      

!osdec05
      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)
!osdec05
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
       bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

!      if(bex(icol,k,kc10,8)<1.e-20_r8) then
!        write(*,995) 'bex(8)=', kc10, t_xrh, t_xct, t_xfac, t_xfaq, bex(icol,k,kc10,8)
!        write(*,996) 'esssf =', esssf1, esssf2, esssf7, esssf8, esssf9, esssf10    
!        write(*,996) 'bees  =', be111, be112, be121, be122, be211, be212, be221, be222
!      endif


      do i=8,8            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

      ke1111=ke4(i,t_irh1,t_ict1,t_ifc1,t_ifa1)
      ke1112=ke4(i,t_irh1,t_ict1,t_ifc1,t_ifa2)
      ke1121=ke4(i,t_irh1,t_ict1,t_ifc2,t_ifa1)
      ke1122=ke4(i,t_irh1,t_ict1,t_ifc2,t_ifa2)
      ke1211=ke4(i,t_irh1,t_ict2,t_ifc1,t_ifa1)
      ke1212=ke4(i,t_irh1,t_ict2,t_ifc1,t_ifa2)
      ke1221=ke4(i,t_irh1,t_ict2,t_ifc2,t_ifa1)
      ke1222=ke4(i,t_irh1,t_ict2,t_ifc2,t_ifa2)
      ke2111=ke4(i,t_irh2,t_ict1,t_ifc1,t_ifa1)
      ke2112=ke4(i,t_irh2,t_ict1,t_ifc1,t_ifa2)
      ke2121=ke4(i,t_irh2,t_ict1,t_ifc2,t_ifa1)
      ke2122=ke4(i,t_irh2,t_ict1,t_ifc2,t_ifa2)
      ke2211=ke4(i,t_irh2,t_ict2,t_ifc1,t_ifa1)
      ke2212=ke4(i,t_irh2,t_ict2,t_ifc1,t_ifa2)
      ke2221=ke4(i,t_irh2,t_ict2,t_ifc2,t_ifa1)
      ke2222=ke4(i,t_irh2,t_ict2,t_ifc2,t_ifa2)

      ke111=esssf1*ke1111+esssf2*ke1112
      ke112=esssf1*ke1121+esssf2*ke1122
      ke121=esssf1*ke1211+esssf2*ke1212
      ke122=esssf1*ke1221+esssf2*ke1222
      ke211=esssf1*ke2111+esssf2*ke2112
      ke212=esssf1*ke2121+esssf2*ke2122
      ke221=esssf1*ke2211+esssf2*ke2212
      ke222=esssf1*ke2221+esssf2*ke2222

      ke11 =esssf7*ke111+esssf8*ke112
      ke12 =esssf7*ke121+esssf8*ke122
      ke21 =esssf7*ke211+esssf8*ke212
      ke22 =esssf7*ke221+esssf8*ke222

      ske1  =((t_cat2-t_xct)*ke11+(t_xct-t_cat1)*ke12)*esssf3*esssf9*esssf10
      ske2  =((t_cat2-t_xct)*ke21+(t_xct-t_cat1)*ke22)*esssf3*esssf9*esssf10
!osdec05
      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)
!osdec05

!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
      endif

      end do ! i

           endif
         
          end do ! icol
        end do ! k

        end do  ! kcomp

  995 format(A10,I3,4(x,e10.3),x,e12.5)
  996 format(A10,8(x,e10.3))

      return
end subroutine interpol4




subroutine interpol5to10 (lchnk, ncol, daylight, rhum, Nnatk, Camk, &
             xfacin, xfbcin, xfaqin, omega, gass, bex, ske, cxstot)

!     Optimized for speed by Arild Burud/NoSerC, June 2002
!---------------------------------------------------------------
!   Modified by Egil Storen/NoSerC July 2002.
!   The sequence of the indices in arrays om1, g1, be1 and ke1
!   (common block /tab1/) has been rearranged to avoid cache
!   problems while running subroutine interpol1.
!   Files also involved by this modification: initopt.F
!   and opttab.h. 
!---------------------------------------------------------------
!   Rewritten for CAM by Alf Kirkevaag in Feb. 2004, and 
!   modified (interpol1->interpol6to10) for new aerosol schemes 
!   in October 2005.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8

   implicit none


!
! Input arguments
!
   integer, intent(in) :: lchnk                       ! chunk identifier
   integer, intent(in) :: ncol                        ! number of atmospheric columns
   logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.

   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration  
   real(r8), intent(in) :: Camk(pcols,pver,nbmodes)   ! modal internally mixed SO4+BC+OC conc.
   real(r8), intent(in) :: rhum(pcols,pver)           ! level relative humidity (fraction)
   real(r8), intent(in) :: xfacin(pcols,pver,nbmodes) ! modal (OC+BC)/(SO4+BC+OC)
   real(r8), intent(in) :: xfbcin(pcols,pver,nbmodes) ! modal BC/(OC+BC)
   real(r8), intent(in) :: xfaqin(pcols,pver,nbmodes) ! modal SO4(aq)/SO4
!
!
! Input-Output arguments
!
   real(r8), intent(inout) :: cxstot(pcols,pver)    ! excess internally mixed mass ("table overshoot") 
!
!
! Output arguments
!
   real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
   real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
   real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
   real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
!
!---------------------------Local variables-----------------------------
!
      integer i, ierr, irelh, ictot, ifac, ifbc, ifaq, kcomp, k, icol
      integer irh1(pcols,pver), irh2(pcols,pver), ict1(pcols,pver), &
        ict2(pcols,pver), ifac1(pcols,pver), ifac2(pcols,pver),     &
        ifbc1(pcols,pver), ifbc2(pcols,pver), ifaq1(pcols,pver),    &
        ifaq2(pcols,pver)
      real(r8) a, b, relh, catot, fabc, fraq, camdiff
      real(r8) xrh(pcols,pver), xct(pcols,pver), xfac(pcols,pver,nbmodes), &
        xfbc(pcols,pver,nbmodes), xfaq(pcols,pver,nbmodes), cxs(pcols,pver)

!     Temporary storage of often used array elements
      integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2,         &
       t_ifb1, t_ifb2, t_ifc1, t_ifc2
      real(r8) t_faq1, t_faq2, t_xfaq, t_fbc1, t_fbc2, t_xfbc, t_fac1, &
       t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
      real(r8) esssf1, esssf2, esssf3, esssf4, esssf5, esssf6, esssf7, &
        esssf8, esssf9, esssf10, esssf3x6x9x10
      real(r8) om111, om112, om121, om122, om211, om212, om221, om222, &
        om11, om12, om21, om22, ome1, ome2  
      real(r8) g111, g112, g121, g122, g211, g212, g221, g222, &
        g11, g12, g21, g22, ge1, ge2 
      real(r8) be111, be112, be121, be122, be211, be212, be221, be222, &
        be11, be12, be21, be22, bex1, bex2
      real(r8) ke111, ke112, ke121, ke122, ke211, ke212, &
        ke221, ke222, ke11, ke12, ke21, ke22, ske1, ske2
      real(r8) &
        om1111, om1112, om1121, om1122, om1211, om1212, om1221,om1222,  &
        om2111, om2112, om2121, om2122, om2211, om2212, om2221, om2222, & 
        g1111, g1112, g1121, g1122, g1211, g1212, g1221, g1222,         &
        g2111, g2112, g2121, g2122, g2211, g2212, g2221, g2222,         &
        be1111, be1112, be1121, be1122, be1211, be1212, be1221, be1222, & 
        be2111, be2112, be2121, be2122, be2211, be2212, be2221, be2222, &
        ke1111, ke1112, ke1121, ke1122, ke1211, ke1212, ke1221, ke1222, &
        ke2111, ke2112, ke2121, ke2122, ke2211, ke2212, ke2221, ke2222 
      real(r8) &
        om11111, om11112, om11121, om11122, om11211, om11212, om11221,  &
        om11222, om12111, om12112, om12121, om12122, om12211, om12212,  &
        om12221, om12222, om21111, om21112, om21121, om21122, om21211,  &
        om21212, om21221, om21222, om22111, om22112, om22121, om22122,  &
        om22211, om22212, om22221, om22222,                             &
        g11111, g11112, g11121, g11122, g11211, g11212, g11221, g11222, &
        g12111, g12112, g12121, g12122, g12211, g12212, g12221, g12222, &
        g21111, g21112,g21121, g21122, g21211, g21212, g21221, g21222,  &
        g22111, g22112, g22121, g22122, g22211, g22212, g22221, g22222, &
        be11111, be11112, be11121, be11122, be11211, be11212, be11221,  &
        be11222, be12111, be12112, be12121, be12122, be12211, be12212,  &
        be12221, be12222, be21111, be21112, be21121, be21122, be21211,  &
        be21212, be21221, be21222, be22111, be22112, be22121, be22122,  &
        be22211, be22212, be22221, be22222,                             &
        ke11111, ke11112, ke11121, ke11122, ke11211, ke11212, ke11221,  &
        ke11222, ke12111, ke12112, ke12121, ke12122, ke12211, ke12212,  &
        ke12221, ke12222, ke21111, ke21112, ke21121, ke21122, ke21211,  &
        ke21212, ke21221, ke21222, ke22111, ke22112, ke22121, ke22122,  &
        ke22211, ke22212, ke22221, ke22222



!      write(*,*) 'Before xrh-loop'
      do k=1,pver
        do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
        end do 
      end do

!      write(*,*) 'Before rh-loop'
      do irelh=1,9
       do k=1,pver
        do icol=1,ncol
           if(xrh(icol,k) >= rh(irelh).and. &
             xrh(icol,k)<=rh(irelh+1)) then
             irh1(icol,k)=irelh
             irh2(icol,k)=irelh+1
           endif
         end do
       end do
      end do

      do k=1,pver
        do icol=1,ncol
          xct(icol,k)  = 0.0_r8
        end do
      end do

!      write(*,*) 'Before kcomp-loop'
!        do kcomp=1,nbmodes        
!        do kcomp=6,10
        do kcomp=5,10

      do k=1,pver
        do icol=1,ncol

            xfac(icol,k,kcomp) = 0.0_r8
            xfbc(icol,k,kcomp) = 0.0_r8
            xfaq(icol,k,kcomp) = 0.0_r8

!            if(Nnatk(icol,k,kcomp)>0.0_r8) then
!g            if(daylight(icol).and.Nnatk(icol,k,kcomp)>0.0_r8) then
            if(daylight(icol)) then

          xct(icol,k)  = min(max(Camk(icol,k,kcomp) &
                 /(Nnatk(icol,k,kcomp)+eps),cat(kcomp,1)),cat(kcomp,6))
          xfac(icol,k,kcomp) = min(max(xfacin(icol,k,kcomp),fac(1)),fac(6))
          xfbc(icol,k,kcomp) = min(max(xfbcin(icol,k,kcomp),fbc(1)),fbc(6))
          xfaq(icol,k,kcomp) = min(max(xfaqin(icol,k,kcomp),faq(1)),faq(6))
          camdiff=Camk(icol,k,kcomp)-xct(icol,k) &
                         *(Nnatk(icol,k,kcomp)+eps)
          cxs(icol,k)=max(0.0_r8,camdiff)
          cxstot(icol,k)=cxstot(icol,k)+cxs(icol,k)

!        if(icol.eq.1) then
!          write(*,*) 'dir: k, kcomp, cxs =', k, kcomp, cxs(icol,k)
!        endif

!      write(*,*) 'Before cat-loop', kcomp
      do ictot=1,5
            if(xct(icol,k)>=cat(kcomp,ictot).and. &
            xct(icol,k) <= cat(kcomp,ictot+1)) then
             ict1(icol,k)=ictot
             ict2(icol,k)=ictot+1
            endif
      end do ! ictot

!      write(*,*) 'Before fac-loop', kcomp
      do ifac=1,5
            if(xfac(icol,k,kcomp)>=fac(ifac).and. &
             xfac(icol,k,kcomp) <= fac(ifac+1)) then
             ifac1(icol,k)=ifac
             ifac2(icol,k)=ifac+1
            endif
      end do ! ifac

!      write(*,*) 'Before fbc-loop', kcomp
      do ifbc=1,5
            if(xfbc(icol,k,kcomp) >= fbc(ifbc).and. &
             xfbc(icol,k,kcomp) <= fbc(ifbc+1)) then
             ifbc1(icol,k)=ifbc
             ifbc2(icol,k)=ifbc+1
            endif
      end do ! ifbc

!      write(*,*) 'Before faq-loop', kcomp
      do ifaq=1,5
            if(xfaq(icol,k,kcomp) >= faq(ifaq).and. &
            xfaq(icol,k,kcomp) <= faq(ifaq+1)) then
             ifaq1(icol,k)=ifaq
             ifaq2(icol,k)=ifaq+1
            endif
      end do ! ifaq

            endif ! daylight

          end do ! icol
        end do ! k


!      write(*,*) 'Before init-loop', kcomp    ! Bruk For All istedet?
        do i=1,nbands
          do icol=1,ncol
            do k=1,pver
              omega(icol,k,kcomp,i)=0.0_r8
              gass(icol,k,kcomp,i)=0.0_r8
              bex(icol,k,kcomp,i)=0.0_r8
              ske(icol,k,kcomp,i)=0.0_r8
            end do
          end do
        end do
         
        do k=1,pver
          do icol=1,ncol

!           if(Nnatk(icol,k,kcomp)>0.0_r8) then
!g           if(daylight(icol).and.Nnatk(icol,k,kcomp)>0.0_r8) then
           if(daylight(icol)) then

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

!      write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
!      write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
!      write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
!      write(*,*) 't_ifb1,t_ifb2=',t_ifb1,t_ifb2
!      write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

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

!      write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
!      write(*,*) 't_fac1,t_fac2,t_fbc1,t_fbc2=',t_fac1,t_fac2,t_fbc1,t_fbc2
!      write(*,*) 't_faq1,t_faq2=',t_faq1,t_faq2

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


      do i=1,nbands            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  single scattering albedo:

      om11111=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      om11112=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      om11121=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      om11122=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      om11211=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      om11212=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      om11221=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      om11222=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      om12111=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      om12112=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      om12121=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      om12122=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      om12211=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      om12212=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      om12221=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      om12222=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      om21111=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      om21112=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      om21121=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      om21122=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      om21211=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      om21212=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      om21221=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      om21222=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      om22111=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      om22112=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      om22121=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      om22122=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      om22211=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      om22212=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      om22221=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      om22222=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

      om1111=esssf1*om11111+esssf2*om11112
      om1112=esssf1*om11121+esssf2*om11122
      om1121=esssf1*om11211+esssf2*om11212
      om1122=esssf1*om11221+esssf2*om11222
      om1211=esssf1*om12111+esssf2*om12112
      om1212=esssf1*om12121+esssf2*om12122
      om1221=esssf1*om12211+esssf2*om12212
      om1222=esssf1*om12221+esssf2*om12222
      om2111=esssf1*om21111+esssf2*om21112
      om2112=esssf1*om21121+esssf2*om21122
      om2121=esssf1*om21211+esssf2*om21212
      om2122=esssf1*om21221+esssf2*om21222
      om2211=esssf1*om22111+esssf2*om22112
      om2212=esssf1*om22121+esssf2*om22122
      om2221=esssf1*om22211+esssf2*om22212
      om2222=esssf1*om22221+esssf2*om22222

      om111=esssf4*om1111+esssf5*om1112
      om112=esssf4*om1121+esssf5*om1122
      om121=esssf4*om1211+esssf5*om1212
      om122=esssf4*om1221+esssf5*om1222
      om211=esssf4*om2111+esssf5*om2112
      om212=esssf4*om2121+esssf5*om2122
      om221=esssf4*om2211+esssf5*om2212
      om222=esssf4*om2221+esssf5*om2222

      om11 =esssf7*om111+esssf8*om112
      om12 =esssf7*om121+esssf8*om122
      om21 =esssf7*om211+esssf8*om212
      om22 =esssf7*om221+esssf8*om222

      ome1  =((t_cat2-t_xct)*om11+(t_xct-t_cat1)*om12)*esssf3x6x9x10      
      ome2  =((t_cat2-t_xct)*om21+(t_xct-t_cat1)*om22)*esssf3x6x9x10   

!      write(*,*) 'Before omega'
      omega(icol,k,kcomp,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                          /(t_rh2-t_rh1)    
!      write(*,*) omega(icol,k,kcomp,i)

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  asymmetry factor   

      g11111=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      g11112=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      g11121=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      g11122=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      g11211=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      g11212=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      g11221=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      g11222=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      g12111=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      g12112=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      g12121=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      g12122=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      g12211=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      g12212=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      g12221=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      g12222=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      g21111=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      g21112=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      g21121=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      g21122=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      g21211=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      g21212=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      g21221=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      g21222=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      g22111=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      g22112=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      g22121=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      g22122=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      g22211=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      g22212=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      g22221=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      g22222=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

      g1111=esssf1*g11111+esssf2*g11112
      g1112=esssf1*g11121+esssf2*g11122
      g1121=esssf1*g11211+esssf2*g11212
      g1122=esssf1*g11221+esssf2*g11222
      g1211=esssf1*g12111+esssf2*g12112
      g1212=esssf1*g12121+esssf2*g12122
      g1221=esssf1*g12211+esssf2*g12212
      g1222=esssf1*g12221+esssf2*g12222
      g2111=esssf1*g21111+esssf2*g21112
      g2112=esssf1*g21121+esssf2*g21122
      g2121=esssf1*g21211+esssf2*g21212
      g2122=esssf1*g21221+esssf2*g21222
      g2211=esssf1*g22111+esssf2*g22112
      g2212=esssf1*g22121+esssf2*g22122
      g2221=esssf1*g22211+esssf2*g22212
      g2222=esssf1*g22221+esssf2*g22222

      g111=esssf4*g1111+esssf5*g1112
      g112=esssf4*g1121+esssf5*g1122
      g121=esssf4*g1211+esssf5*g1212
      g122=esssf4*g1221+esssf5*g1222
      g211=esssf4*g2111+esssf5*g2112
      g212=esssf4*g2121+esssf5*g2122
      g221=esssf4*g2211+esssf5*g2212
      g222=esssf4*g2221+esssf5*g2222

      g11 =esssf7*g111+esssf8*g112
      g12 =esssf7*g121+esssf8*g122
      g21 =esssf7*g211+esssf8*g212
      g22 =esssf7*g221+esssf8*g222

      ge1  =((t_cat2-t_xct)*g11+(t_xct-t_cat1)*g12)*esssf3x6x9x10      
      ge2  =((t_cat2-t_xct)*g21+(t_xct-t_cat1)*g22)*esssf3x6x9x10      

!      write(*,*) 'Before gass'
      gass(icol,k,kcomp,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                   /(t_rh2-t_rh1)    

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol extinction   

      be11111=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      be11112=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      be11121=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      be11122=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      be11211=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      be11212=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      be11221=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      be11222=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      be12111=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      be12112=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      be12121=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      be12122=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      be12211=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      be12212=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      be12221=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      be12222=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      be21111=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      be21112=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      be21121=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      be21122=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      be21211=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      be21212=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      be21221=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      be21222=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      be22111=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      be22112=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      be22121=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      be22122=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      be22211=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      be22212=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      be22221=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      be22222=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

      be1111=esssf1*be11111+esssf2*be11112
      be1112=esssf1*be11121+esssf2*be11122
      be1121=esssf1*be11211+esssf2*be11212
      be1122=esssf1*be11221+esssf2*be11222
      be1211=esssf1*be12111+esssf2*be12112
      be1212=esssf1*be12121+esssf2*be12122
      be1221=esssf1*be12211+esssf2*be12212
      be1222=esssf1*be12221+esssf2*be12222
      be2111=esssf1*be21111+esssf2*be21112
      be2112=esssf1*be21121+esssf2*be21122
      be2121=esssf1*be21211+esssf2*be21212
      be2122=esssf1*be21221+esssf2*be21222
      be2211=esssf1*be22111+esssf2*be22112
      be2212=esssf1*be22121+esssf2*be22122
      be2221=esssf1*be22211+esssf2*be22212
      be2222=esssf1*be22221+esssf2*be22222

      be111=esssf4*be1111+esssf5*be1112
      be112=esssf4*be1121+esssf5*be1122
      be121=esssf4*be1211+esssf5*be1212
      be122=esssf4*be1221+esssf5*be1222
      be211=esssf4*be2111+esssf5*be2112
      be212=esssf4*be2121+esssf5*be2122
      be221=esssf4*be2211+esssf5*be2212
      be222=esssf4*be2221+esssf5*be2222

      be11 =esssf7*be111+esssf8*be112
      be12 =esssf7*be121+esssf8*be122
      be21 =esssf7*be211+esssf8*be212
      be22 =esssf7*be221+esssf8*be222

      bex1  =((t_cat2-t_xct)*be11+(t_xct-t_cat1)*be12)*esssf3x6x9x10      
      bex2  =((t_cat2-t_xct)*be21+(t_xct-t_cat1)*be22)*esssf3x6x9x10      
!osdec05
      bex1=max(bex1,1.e-30_r8)
      bex2=max(bex2,1.e-30_r8)
!osdec05
!      write(*,*) 'Before bex'
      if(t_xrh <= 0.37_r8) then
       bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
        b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
        bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
      endif

      end do ! i


      do i=8,8            ! i = wavelength index

!ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
!  aerosol specific extinction   

      ke11111=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      ke11112=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      ke11121=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      ke11122=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      ke11211=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      ke11212=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      ke11221=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      ke11222=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      ke12111=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      ke12112=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      ke12121=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      ke12122=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      ke12211=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      ke12212=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      ke12221=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      ke12222=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
      ke21111=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
      ke21112=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
      ke21121=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
      ke21122=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
      ke21211=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
      ke21212=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
      ke21221=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
      ke21222=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
      ke22111=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
      ke22112=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
      ke22121=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
      ke22122=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
      ke22211=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
      ke22212=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
      ke22221=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
      ke22222=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

      ke1111=esssf1*ke11111+esssf2*ke11112
      ke1112=esssf1*ke11121+esssf2*ke11122
      ke1121=esssf1*ke11211+esssf2*ke11212
      ke1122=esssf1*ke11221+esssf2*ke11222
      ke1211=esssf1*ke12111+esssf2*ke12112
      ke1212=esssf1*ke12121+esssf2*ke12122
      ke1221=esssf1*ke12211+esssf2*ke12212
      ke1222=esssf1*ke12221+esssf2*ke12222
      ke2111=esssf1*ke21111+esssf2*ke21112
      ke2112=esssf1*ke21121+esssf2*ke21122
      ke2121=esssf1*ke21211+esssf2*ke21212
      ke2122=esssf1*ke21221+esssf2*ke21222
      ke2211=esssf1*ke22111+esssf2*ke22112
      ke2212=esssf1*ke22121+esssf2*ke22122
      ke2221=esssf1*ke22211+esssf2*ke22212
      ke2222=esssf1*ke22221+esssf2*ke22222

      ke111=esssf4*ke1111+esssf5*ke1112
      ke112=esssf4*ke1121+esssf5*ke1122
      ke121=esssf4*ke1211+esssf5*ke1212
      ke122=esssf4*ke1221+esssf5*ke1222
      ke211=esssf4*ke2111+esssf5*ke2112
      ke212=esssf4*ke2121+esssf5*ke2122
      ke221=esssf4*ke2211+esssf5*ke2212
      ke222=esssf4*ke2221+esssf5*ke2222

      ke11 =esssf7*ke111+esssf8*ke112
      ke12 =esssf7*ke121+esssf8*ke122
      ke21 =esssf7*ke211+esssf8*ke212
      ke22 =esssf7*ke221+esssf8*ke222

      ske1  =((t_cat2-t_xct)*ke11+(t_xct-t_cat1)*ke12)*esssf3x6x9x10      
      ske2  =((t_cat2-t_xct)*ke21+(t_xct-t_cat1)*ke22)*esssf3x6x9x10      
!osdec05
      ske1=max(ske1,1.e-30_r8)
      ske2=max(ske2,1.e-30_r8)
!osdec05

!      write(*,*) 'Before ske'
      if(t_xrh <= 0.37_r8) then
        ske(icol,k,kcomp,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
            /(t_rh2-t_rh1)    
      else
        a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
       b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
        ske(icol,k,kcomp,i)=e**(a*t_xrh+b)
      endif

      end do ! i

           endif
         
          end do ! icol
        end do ! k

        end do  ! kcomp

      return
end subroutine interpol5to10




end module optinterpol
