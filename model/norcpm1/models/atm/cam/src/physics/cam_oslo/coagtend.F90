      subroutine coagtend(lchnk,   ncol,    ncnst,   icefrac, &
                   landfrac,     ocnfrac,   t,       q,        cldfrc, &
                   pdel ,dqdt,    dotend,      nrmodes,     ndrops         ) 

! Calculate the coagulation of small aerosols with larger particles and 
! cloud droplets. Only particles smaller that dry radius of 
! 40 nm is assumed to have an efficient coagulation with other particles.

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use constituents, only: pcnst
use const
use cam_history,  only: outfld
use mass,         only: cmidry
use aerosoldef

implicit none
 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: ncnst                ! number of tracers 
!   real(r8), intent(in) :: clat(pcols)          ! Latitude (radians)
   real(r8), intent(in) :: icefrac(pcols)       ! sea ice fraction (fraction)
   real(r8), intent(in) :: landfrac(pcols)      ! land fraction (fraction)
   real(r8), intent(in) :: ocnfrac(pcols)       ! ocean fraction (fraction)
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature (K)
   real(r8), intent(in) :: q(pcols,pver,pcnst) ! TMR including moisture
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction
   real(r8), intent(in) :: pdel(pcols,pver)  ! Delta p
   real(r8), intent(inout) :: dqdt(pcols,pver,pcnst)  ! TMR tendency array
   logical,  intent(inout) :: dotend(pcnst)   
   real(r8), intent(in) :: nrmodes(pcols,pver,pcnst) ! number concentration in each mode
   real(r8), intent(in) :: ndrops(pcols,pver) ! droplet number


! local
   integer :: i,k,m,nsiz 
   real(r8) :: cmi2d(pcols)
   real(r8) :: fice(pcols,pver)     ! Fraction of cwat that is ice. Copied from Iversen and Seland 2003
   real(r8) :: colcoagbc(pcols,pver),colcoagpom(pcols,pver),colcoags4(pcols,pver)
   real(r8) :: clcoagbc(pcols,pver)
   real(r8) :: clcoagpom(pcols,pver)
   real(r8) :: clcoagax(pcols,pver)
   real(r8) :: clcoags4(pcols,pver)

!   dotend(:)=.false.
   dqdt(:,:,2:pcnst)=0._r8
   

    do k=1,pver
      do i=1,ncol
! If warmer than -10 degrees C then water phase
!
         if (t(i,k) > 273.16_r8) fice(i,k) = 0.0_r8
!
! If colder than 0 degrees C but warmer than -25 C mixed phase
!
         if (t(i,k) <= 273.16_r8 .and. t(i,k) >= 248.16_r8) then
            fice(i,k) =(273.16_r8-t(i,k)) / 25.0_r8
         end if
!
! If colder than -25 degrees C then ice phase
!
         if (t(i,k) < 248.16_r8) fice(i,k) = 1.0_r8

      end do
   end do	
   colcoagbc(:,:)=0._r8
   colcoagpom(:,:)=0._r8
   colcoags4(:,:)=0._r8
   clcoagbc(:,:)=0._r8
   clcoagpom(:,:)=0._r8
   clcoagax(:,:)=0._r8
   clcoags4(:,:)=0._r8
!  Explicit coagulation
   do k=1,pver
       do i=1,ncol  
        do nsiz=0,imax	    

         colcoagbc(i,k)=colcoagbc(i,k)+1.e-12_r8*normnk(6,nsiz)*Kp12(6,nsiz)*nrmodes(i,k,l_dst_a2)
         colcoagbc(i,k)=colcoagbc(i,k)+1.e-12_r8*normnk(7,nsiz)*Kp12(7,nsiz)*nrmodes(i,k,l_dst_a3)
         colcoagbc(i,k)=colcoagbc(i,k)+1.e-12_r8*normnk(8,nsiz)*Kp12(8,nsiz)*nrmodes(i,k,l_ss_a1)
         colcoagbc(i,k)=colcoagbc(i,k)+1.e-12_r8*normnk(9,nsiz)*Kp12(9,nsiz)*nrmodes(i,k,l_ss_a2)
         colcoagbc(i,k)=colcoagbc(i,k)+1.e-12_r8*normnk(10,nsiz)*Kp12(10,nsiz)*nrmodes(i,k,l_ss_a3)

         colcoagpom(i,k)=colcoagpom(i,k)+1.e-12_r8*normnk(6,nsiz)*Kp12oc(6,nsiz)*nrmodes(i,k,l_dst_a2)
         colcoagpom(i,k)=colcoagpom(i,k)+1.e-12_r8*normnk(7,nsiz)*Kp12oc(7,nsiz)*nrmodes(i,k,l_dst_a3)
         colcoagpom(i,k)=colcoagpom(i,k)+1.e-12_r8*normnk(8,nsiz)*Kp12oc(8,nsiz)*nrmodes(i,k,l_ss_a1)
         colcoagpom(i,k)=colcoagpom(i,k)+1.e-12_r8*normnk(9,nsiz)*Kp12oc(9,nsiz)*nrmodes(i,k,l_ss_a2)
         colcoagpom(i,k)=colcoagpom(i,k)+1.e-12_r8*normnk(10,nsiz)*Kp12oc(10,nsiz)*nrmodes(i,k,l_ss_a3)

         colcoags4(i,k)=colcoags4(i,k)+1.e-12_r8*normnk(6,nsiz)*Kp12s4(6,nsiz)*nrmodes(i,k,l_dst_a2)
         colcoags4(i,k)=colcoags4(i,k)+1.e-12_r8*normnk(7,nsiz)*Kp12s4(7,nsiz)*nrmodes(i,k,l_dst_a3)
         colcoags4(i,k)=colcoags4(i,k)+1.e-12_r8*normnk(8,nsiz)*Kp12s4(8,nsiz)*nrmodes(i,k,l_ss_a1)
         colcoags4(i,k)=colcoags4(i,k)+1.e-12_r8*normnk(9,nsiz)*Kp12s4(9,nsiz)*nrmodes(i,k,l_ss_a2)
         colcoags4(i,k)=colcoags4(i,k)+1.e-12_r8*normnk(10,nsiz)*Kp12s4(10,nsiz)*nrmodes(i,k,l_ss_a3)
!ak2feb-
	end do

!	clcoagbc(i,k)=1.e-12_r8*Kp12(9,42)*ndrops(i,k)
!	clcoagpom(i,k)=1.e-12_r8*Kp12oc(9,42)*ndrops(i,k)*cldfrc(i,k)
        clcoags4(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)* &
           1.e-12_r8*Kp12s4(3,42)*ndrops(i,k)
        clcoags4(i,k)=min(clcoags4(i,k),2.8e-4_r8)
        colcoags4(i,k)=min(colcoags4(i,k),2.8e-4_r8)

        clcoagbc(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)* &
           1.e-12_r8*Kp12(3,42)*ndrops(i,k)

        clcoagbc(i,k)=min(clcoagbc(i,k),2.8e-4_r8)

        clcoagpom(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)* &
           1.e-12_r8*Kp12oc(3,42)*ndrops(i,k)

        clcoagpom(i,k)=min(clcoagpom(i,k),2.8e-4_r8)

        clcoagax(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)* &
           1.e-12_r8*Kp12ax(3,42)*ndrops(i,k)

        clcoagax(i,k)=min(clcoagax(i,k),2.8e-4_r8)

     end do
   end do




          dotend(l_so4_n)=.true.
          dotend(l_so4_na)=.true.
          dotend(l_so4_a1)=.true.
          dotend(l_so4_a2)=.true.
          dotend(l_so4_ac)=.true.
          dotend(l_bc_n)=.true.
          dotend(l_bc_ni)=.true.
          dotend(l_bc_ax)=.true.
          dotend(l_bc_a)=.true.
          dotend(l_bc_ai)=.true.
          dotend(l_bc_ac)=.true.
!          dotend(l_om_n)=.true.
          dotend(l_om_ni)=.true.
!          dotend(l_om_a)=.true.
          dotend(l_om_ai)=.true.
          dotend(l_om_ac)=.true.

          do k=1,pver
             do i=1,ncol
!	        clcoagbc(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)*clcoag(1)
!	        clcoagpom(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)*clcoag(2)
!	        clcoagax(i,k)=(1._r8-fice(i,k))*cldfrc(i,k)*clcoag(3)
	        
                dqdt(i,k,l_so4_n) = -colcoags4(i,k)*q(i,k,l_so4_n)- &
             clcoags4(i,k)*q(i,k,l_so4_n)

                dqdt(i,k,l_so4_na) = -colcoagpom(i,k)*q(i,k,l_so4_na)- &
             clcoagpom(i,k)*q(i,k,l_so4_na)     

                dqdt(i,k,l_so4_a2) = &
                clcoags4(i,k)*q(i,k,l_so4_n)+ &
             clcoagpom(i,k)*q(i,k,l_so4_a1)+ &
             clcoagpom(i,k)*q(i,k,l_so4_na)          

                dqdt(i,k,l_so4_a1) = -colcoagpom(i,k)*q(i,k,l_so4_a1)- &
             clcoagpom(i,k)*q(i,k,l_so4_a1)

                dqdt(i,k,l_so4_ac) = colcoagpom(i,k)*q(i,k,l_so4_a1)+ &
                  colcoags4(i,k)*q(i,k,l_so4_n)+ &
                  colcoagpom(i,k)*q(i,k,l_so4_na)
	        
!+ &
!             (1-fice(i,k))*clcoag(3)*cldfrc(i,k)*q(i,k,l_so4_a1)

                dqdt(i,k,l_bc_n) = -colcoagbc(i,k)*q(i,k,l_bc_n)- &
             clcoagbc(i,k)*q(i,k,l_bc_n)

                dqdt(i,k,l_bc_ni) = -colcoagpom(i,k)*q(i,k,l_bc_ni)- &
             clcoagpom(i,k)*q(i,k,l_bc_ni)

                dqdt(i,k,l_bc_ax) = - &
             clcoagax(i,k)*q(i,k,l_bc_ax)

	        dqdt(i,k,l_bc_a) = -colcoagpom(i,k)*q(i,k,l_bc_a)- &
             clcoagpom(i,k)*q(i,k,l_bc_a)

              dqdt(i,k,l_bc_ai) = -colcoagpom(i,k)*q(i,k,l_bc_ai)- &
             clcoagpom(i,k)*q(i,k,l_bc_ai)
                 
                dqdt(i,k,l_bc_ac) = colcoagbc(i,k)*q(i,k,l_bc_n)+ &
                clcoagbc(i,k)*q(i,k,l_bc_n) + &
                clcoagax(i,k)*q(i,k,l_bc_ax)+ &
                                  colcoagpom(i,k)*q(i,k,l_bc_ni)+ &
                clcoagpom(i,k)*q(i,k,l_bc_ni)+&
                                  colcoagpom(i,k)*q(i,k,l_bc_a)+ &
                clcoagpom(i,k)*q(i,k,l_bc_a)+&
                                  colcoagpom(i,k)*q(i,k,l_bc_ai)+ &
                clcoagpom(i,k)*q(i,k,l_bc_ai)

!                dqdt(i,k,l_om_n) = -colcoagpom(i,k)*q(i,k,l_om_n)- &
!             clcoagpom(i,k)*q(i,k,l_om_n)

                dqdt(i,k,l_om_ni) = -colcoagpom(i,k)*q(i,k,l_om_ni)- &
             clcoagpom(i,k)*q(i,k,l_om_ni)

!	        dqdt(i,k,l_om_a) = -colcoagpom(i,k)*q(i,k,l_om_a)- &
!             clcoagpom(i,k)*q(i,k,l_om_a)

	        dqdt(i,k,l_om_ai) = -colcoagpom(i,k)*q(i,k,l_om_ai)- &
             clcoagpom(i,k)*q(i,k,l_om_ai)

                dqdt(i,k,l_om_ac) = colcoagpom(i,k)*q(i,k,l_om_ni)+ &
             clcoagpom(i,k)*q(i,k,l_om_ni)+&
                                  colcoagpom(i,k)*q(i,k,l_om_ai)+ &
             clcoagpom(i,k)*q(i,k,l_om_ai)


             end do
          end do     
!   end if       
!   end do

!   call outfld('COAGBC      ',colcoagbc           ,pcols   ,lchnk   )
!   call outfld('COAGPOM     ',colcoagpom          ,pcols   ,lchnk   )
!   call outfld('COAGS4     ',colcoags4          ,pcols   ,lchnk   )
!   call outfld('CLCOAGBC    ',clcoagbc            ,pcols   ,lchnk   )
!   call outfld('CLCOAGOC   ',clcoagpom           ,pcols   ,lchnk   )
!   call outfld('CLCOAGAX    ',clcoagax            ,pcols   ,lchnk   )
!   call outfld('CLCOAGS4    ',clcoags4            ,pcols   ,lchnk   )


!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_ac),cmi2d)
!      call outfld('PCSO4AC ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_a2),cmi2d)
!      call outfld('PCSO4A2 ',cmi2d   ,pcols   ,lchnk) 

!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_n),cmi2d)
!      call outfld('LCSO4N ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_na),cmi2d)
!      call outfld('LCSO4NA ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_a1),cmi2d)
!      call outfld('LCSO4A1 ',cmi2d   ,pcols   ,lchnk) 

!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_n),cmi2d)
!      call outfld('LCBCN ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_om_n),cmi2d)
!      call outfld('LCOMN ',cmi2d   ,pcols   ,lchnk) 

!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_a),cmi2d)
!      call outfld('LCBCA ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_ax),cmi2d)
!      call outfld('LCBCAX ',cmi2d   ,pcols   ,lchnk) 

!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_ni),cmi2d)
!      call outfld('LCBCNI ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_om_ni),cmi2d)
!      call outfld('LCOMNI ',cmi2d   ,pcols   ,lchnk) 

!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_ai),cmi2d)
!      call outfld('LCBCAI ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_om_ai),cmi2d)
!      call outfld('LCOMAI ',cmi2d   ,pcols   ,lchnk) 


   return
   end
