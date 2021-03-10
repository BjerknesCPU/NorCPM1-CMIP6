      subroutine calccol(lchnk,   ncol, q,pdel,rhoair,wetdepflx,cam_out) 

! Calculate columns of total (not mode)
! concentrations and production/loss processes

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use constituents, only: pcnst,wetflxnam
  use cam_history,  only: outfld
  use mass,         only: cmidry
  use aerosoldef
  use camsrfexch_types, only: cam_out_t
  use aero_to_srf, only: set_srf_wetdep
  use physconst, only: rair
  implicit none

 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: q(pcols,pver,pcnst) ! TMR including moisture
   real(r8), intent(in) :: pdel(pcols,pver)  ! Delta p
   real(r8), intent(in) :: rhoair(pcols,pver)  ! Air density
!   real(r8), intent(in) :: drydepflx(pcols,ppcnst) ! Dry dep to the ground
   real(r8), intent(in) :: wetdepflx(pcols,pcnst) ! Wet dep to the ground
   type(cam_out_t), intent(inout) :: cam_out
! local
   integer :: i,k,m,mm 
   real(r8) :: cmi2d(pcols)
   real(r8) :: tdms(pcols,pver)
!   real(r8) :: tmsa(pcols,pver),wetmsa(pcols)
   real(r8) :: tso2(pcols,pver),wetso2(pcols)
   real(r8) :: tso4(pcols,pver),wetso4(pcols)
   real(r8) :: tbc(pcols,pver)
   real(r8) :: tpom(pcols,pver),wetpom(pcols),wetbc(pcols),wetdust(pcols)
   real(r8) :: tdust(pcols,pver)
   real(r8) :: tsalt(pcols,pver),drysalt(pcols),wetsalt(pcols)
   real(r8) :: codms(pcols,pver),coso2(pcols,pver),coso4(pcols,pver)
   real(r8) :: cobc(pcols,pver),copom(pcols,pver)
   real(r8) :: coss(pcols,pver),codust(pcols,pver)

   cmi2d(:)=0._r8
   do k=1,pver
      do i=1,ncol
    	tso4(i,k) = q(i,k,l_so4_n)+q(i,k,l_so4_na)+q(i,k,l_so4_a1)+ & 
       q(i,k,l_so4_a2)+q(i,k,l_so4_ac)+q(i,k,l_so4_pr)
        coso4(i,k)=tso4(i,k)*rhoair(i,k)
	tbc(i,k) = q(i,k,l_bc_n)+q(i,k,l_bc_ax)+   &
   q(i,k,l_bc_ni)+q(i,k,l_bc_a)+q(i,k,l_bc_ai)+q(i,k,l_bc_ac)
        cobc(i,k)=tbc(i,k)*rhoair(i,k)
        tpom(i,k) = q(i,k,l_om_ni)+  &
   q(i,k,l_om_ai)+q(i,k,l_om_ac)
        copom(i,k)=tpom(i,k)*rhoair(i,k)
        tdms(i,k) = q(i,k,l_dms)
        codms(i,k)=tdms(i,k)*(29._r8/46._r8)
	tso2(i,k) = q(i,k,l_so2)
        coso2(i,k)=tso2(i,k)*(29._r8/32._r8)
!        tmsa(i,k) = q(i,k,l_msa)
	tdust(i,k)= q(i,k,l_dst_a2)+q(i,k,l_dst_a3)
        codust(i,k)=tdust(i,k)*rhoair(i,k)
	tsalt(i,k)= q(i,k,l_ss_a1)+q(i,k,l_ss_a2)+q(i,k,l_ss_a3)
        coss(i,k)=tsalt(i,k)*rhoair(i,k)
      end do
   end do

   do i=1,ncol  

! Calculation of dry deposition moved to genaero_intr.F90

!     drymsa(i) = drydepflx(i,l_msa)
!     dryso2(i) = drydepflx(i,l_so2)
!     dryso4(i) = drydepflx(i,l_so4_n)+drydepflx(i,l_so4_na)   &
!+ drydepflx(i,l_so4_a1)+drydepflx(i,l_so4_a2)+drydepflx(i,l_so4_ac)+drydepflx(i,l_so4_pr)
!     drybc(i) = drydepflx(i,l_bc_n)+drydepflx(i,l_bc_ax) &
!+ drydepflx(i,l_bc_ni)+drydepflx(i,l_bc_a)+drydepflx(i,l_bc_ai)+drydepflx(i,l_bc_ac)
!     drypom(i) = drydepflx(i,l_om_n)+drydepflx(i,l_om_ni)&
!+drydepflx(i,l_om_ai)+drydepflx(i,l_om_ac)
!     drydust(i)= drydepflx(i,l_dst_a2)+drydepflx(i,l_dst_a3)
!     drysalt(i)= drydepflx(i,l_ss_a1)+drydepflx(i,l_ss_a2)+drydepflx(i,l_ss_a3)


!     wetmsa(i) = wetdepflx(i,l_msa)
     wetso2(i) = wetdepflx(i,l_so2)
     wetso4(i) = wetdepflx(i,l_so4_n)+wetdepflx(i,l_so4_na)+wetdepflx(i,l_so4_a1) &
+wetdepflx(i,l_so4_a2)+wetdepflx(i,l_so4_ac)+wetdepflx(i,l_so4_pr)
     wetbc(i) = wetdepflx(i,l_bc_n)+wetdepflx(i,l_bc_ax) &
+wetdepflx(i,l_bc_ni)+wetdepflx(i,l_bc_a)+wetdepflx(i,l_bc_ai)+wetdepflx(i,l_bc_ac)
     wetpom(i) = wetdepflx(i,l_om_ni) &
+wetdepflx(i,l_om_ai)+wetdepflx(i,l_om_ac)
     wetdust(i)= wetdepflx(i,l_dst_a2)+wetdepflx(i,l_dst_a3)
     wetsalt(i)= wetdepflx(i,l_ss_a1)+wetdepflx(i,l_ss_a2)+wetdepflx(i,l_ss_a3)

   end do

    call set_srf_wetdep(wetdepflx, cam_out)

!   Column burden for each component
!     do m=1,pcnst
!        call cmidry( lchnk,ncol,pdel, q(:,:,1), q(:,:,m), cmi2d)     
!        call outfld(colname(m),cmi2d   ,pcols   ,lchnk   )
!     end do
!     Deposition for each component
!      do m=1,ppcnst
!	call outfld(dryflxnam(m),drydepflx(:,m),pcols,lchnk)
!	call outfld(wetflxnam(m),wetdepflx(:,m),pcols,lchnk)
!      end do
!#endif
#ifdef SHORTRUN
       do m=1,ncui
         mm=ixac-1+m
	call outfld(wetflxnam(mm),wetdepflx(:,mm),pcols,lchnk)
       end do	
#endif
!  Concentrations
      call outfld('DMSCO ',codms   ,pcols   ,lchnk   )
!      call outfld('MSA ',tmsa   ,pcols   ,lchnk   )
      call outfld('SO2CO ',coso2   ,pcols   ,lchnk   )
      call outfld('SO4 ',coso4   ,pcols   ,lchnk   )
      call outfld('BC ',cobc   ,pcols   ,lchnk   )
      call outfld('POM ',copom   ,pcols   ,lchnk   )
      call outfld('DUST ',codust   ,pcols   ,lchnk   )
      call outfld('SS ',coss   ,pcols   ,lchnk   )
! Deposition fluxes

      call outfld('WET_SO2 ',wetso2   ,pcols   ,lchnk   )
      call outfld('WET_SO4 ',wetso4   ,pcols   ,lchnk   )
      call outfld('WET_BC ',wetbc   ,pcols   ,lchnk   )
      call outfld('WET_POM ',wetpom   ,pcols   ,lchnk   )
      call outfld('WET_DUST ',wetdust   ,pcols   ,lchnk   )
      call outfld('WET_SS ',wetsalt   ,pcols   ,lchnk   )

! Column integrals


      call cmidry( lchnk,ncol,pdel, q(:,:,1), tdms, cmi2d)     
      call outfld('C_DMS ',cmi2d   ,pcols   ,lchnk   )
!      call cmidry( lchnk,ncol,pdel, q(:,:,1), tmsa, cmi2d)     
!      call outfld('C_MSA ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tso2, cmi2d)     
      call outfld('C_SO2 ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tso4, cmi2d)     
      call outfld('C_SO4 ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tbc, cmi2d)     
      call outfld('C_BC ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tpom, cmi2d)     
      call outfld('C_POM ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tdust, cmi2d)     
      call outfld('C_DUST ',cmi2d   ,pcols   ,lchnk   )
      call cmidry( lchnk,ncol,pdel, q(:,:,1), tsalt, cmi2d)     
      call outfld('C_SS ',cmi2d   ,pcols   ,lchnk   )

   return
   end
