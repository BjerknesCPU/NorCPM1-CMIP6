      subroutine condtend(lchnk,   ncol,    ncnst,   t,       q, &
       cldfrc, pdel ,dqdt,    dotend,      nrmodes, dt,loch2so4) 

! Calculate the sulphate nucleation rate, and condensation rate of 
!aerosols used for parameterising the transfer of externally mixed 
!aitken mode particles into an external mixture.
! Note the parameterisation for conversion of externally mixed particles 
!used the h2so4 lifetime onto the particles, and not a given 
! increase in particle radius. Will be improved in future versions of the model

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
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature (K)
   real(r8), intent(in) :: pdel(pcols,pver)  ! Delta p
   real(r8), intent(in) :: q(pcols,pver,pcnst) ! TMR including moisture
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction
   real(r8), intent(inout) :: dqdt(pcols,pver,pcnst)  ! TMR tendency array
   logical,  intent(inout) :: dotend(pcnst)   
   real(r8), intent(in) :: nrmodes(pcols,pver,pcnst) ! number concentration in each mode
   real(r8), intent(in) :: dt                   ! Time step
   real(r8), intent(in) :: loch2so4(pcols,pver) ! h2so4 produced from gas-phase oxidation of so2 or implicit from MSA


! local
   integer :: i,k,m,nsiz 
   real(r8) :: cmi2d(pcols)
   real(r8) :: fice(pcols,pver)     ! Fraction of cwat that is ice. Copied from Iversen and Seland 2003
   real(r8) :: condh2so4(pcols,pver)
   real(r8) :: conds4(pcols,pver)
   real(r8) :: condbc(pcols,pver)
   real(r8) :: condom(pcols,pver)
   real(r8) :: condbcom(pcols,pver)
   real(r8) :: condbck(pcols,pver)
   real(r8) :: condaits4(pcols,pver),condaitbc(pcols,pver)
   real(r8) :: condaitbcp(pcols,pver)
   real(r8) :: condbcax(pcols,pver)
   real(r8) :: conds4pr(pcols,pver)
   real(r8) :: nrmodesait(pcols,pver)
   real(r8) :: nuclso4(pcols,pver)
!   real(r8) :: Drmp(0:100)
   real(r8) ::cmid2d(pcols)


   real(r8) ::  h2so4new(pcols,pver)
   real(r8) ::  dh2so4(pcols,pver)
   real(r8) :: nbcpom(pcols,pver)
#ifdef CONDFIX
   real(r8) :: n_so4_monolayers_age
   real(r8) :: dr_so4_monolayers_age 
   real(r8) :: area_core,volume_shell,volume_monolayer,frac_transfer
#endif 
!   dotend(:)=.false.
   dqdt(:,:,2:pcnst)=0._r8
   


#ifdef CONDFIX  
! thickness of the so4 monolayers (m)
! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3 as in MAM
   ! Assumed number of monolayers
    n_so4_monolayers_age = 3.0_r8
    dr_so4_monolayers_age = n_so4_monolayers_age * 4.76e-10_r8
#endif
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
    condh2so4(:,:)=0._r8
    conds4(:,:)=0._r8
    condbc(:,:)=0._r8
    condom(:,:)=0._r8
    condbcom(:,:)=0._r8
    nbcpom(:,:)=0._r8
    condbck(:,:)=0._r8
    condaits4(:,:)=0._r8
    condaitbc(:,:)=0._r8
    condaitbcp(:,:)=0._r8
    condbcax(:,:)=0._r8
    conds4pr(:,:)=0._r8

!  Explicit coagulation

!   do nsiz=0,imax
!	Drmp(nsiz)=rp(nsiz)*Dmp(nsiz)
!   end do	
    do k=1,pver
      do i=1,ncol

        do nsiz=0,imax	    
         conds4(i,k)= conds4(i,k)+normnk(1,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_so4_n)
         condaits4(i,k)= condaits4(i,k)+normnk(4,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_so4_na)
!
         condbc(i,k)= condbc(i,k)+normnk(2,nsiz)*rp(nsiz)*Dmp_bcn(nsiz)*nrmodes(i,k,l_bc_n)
         condaitbc(i,k)= condaitbc(i,k)+normnk(4,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_bc_a)
!
!         condom(i,k)= condom(i,k)+normnk(3,nsiz)*rp(nsiz)*Dmp_omn(nsiz)*nrmodes(i,k,l_om_n)
!         condaitpom(i,k)= condaitpom(i,k)+normnk(4,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_om_a)
!
         condbcom(i,k)= condbcom(i,k)+normnk(4,nsiz)*rp(nsiz)*Dmp_ni(nsiz)* &
      (nrmodes(i,k,l_bc_ni)+nrmodes(i,k,l_om_ni))
         condaitbcp(i,k)= condaitbcp(i,k)+normnk(4,nsiz)*rp(nsiz)*Dmp(nsiz)* &
       (nrmodes(i,k,l_bc_ai)+nrmodes(i,k,l_om_ai))
!
!ak         condbcax(i,k)= condbcax(i,k)+normnk(5,nsiz)*rp(nsiz)*Dmp_bcn(nsiz)*nrmodes(i,k,l_bc_ax)
         condbcax(i,k)= condbcax(i,k)+normnk(0,nsiz)*rp(nsiz)*Dmp_bcn(nsiz)*nrmodes(i,k,l_bc_ax)
!
!ak         conds4pr(i,k) = conds4pr(i,k)+normnk(11,nsiz)*rp(nsiz)*Dmp(nsiz)* & 
         conds4pr(i,k) = conds4pr(i,k)+normnk(5,nsiz)*rp(nsiz)*Dmp(nsiz)* & 
        nrmodes(i,k,l_so4_pr)
!
         condbck(i,k)= condbck(i,k)+normnk(6,nsiz)*rp(nsiz)*Dmp_dst(nsiz)*nrmodes(i,k,l_dst_a2)
         condbck(i,k)= condbck(i,k)+normnk(7,nsiz)*rp(nsiz)*Dmp_dst(nsiz)*nrmodes(i,k,l_dst_a3)
         condbck(i,k)= condbck(i,k)+normnk(8,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_ss_a1)
         condbck(i,k)= condbck(i,k)+normnk(9,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_ss_a2)
         condbck(i,k)= condbck(i,k)+normnk(10,nsiz)*rp(nsiz)*Dmp(nsiz)*nrmodes(i,k,l_ss_a3)	
	end do

	   condaits4(i,k)=1.e-12_r8*condaits4(i,k)
	   condaitbc(i,k)=1.e-12_r8*condaitbc(i,k)
	   condaitbcp(i,k)=1.e-12_r8*condaitbcp(i,k)
	   condh2so4(i,k)=1.e-12_r8*condh2so4(i,k)
	   conds4(i,k)=1.e-12_r8*conds4(i,k)
	   condbc(i,k)=1.e-12_r8*condbc(i,k)
	   condom(i,k)=1.e-12_r8*condom(i,k)
	   condbcom(i,k)=1.e-12_r8*condbcom(i,k)
	   condbck(i,k)=1.e-12_r8*condbck(i,k)
           condbcax(i,k)=1.e-12_r8*condbcax(i,k)
           conds4pr(i,k)=1.e-12_r8*conds4pr(i,k)

!
	   condh2so4(i,k)=conds4(i,k)+condbc(i,k)+condom(i,k)+condbcom(i,k)+ &
           condaits4(i,k)+condaitbc(i,k)+condaitbcp(i,k)+ &
           condbcax(i,k)+conds4pr(i,k)+condbck(i,k)    

           condbc(i,k)=min(condbc(i,k),2.8e-4_r8)
	   condom(i,k)=min(condom(i,k),2.8e-4_r8)
	   condbcom(i,k)=min(condbcom(i,k),2.8e-4_r8)
           nbcpom(i,k)=nrmodes(i,k,l_bc_ni)+nrmodes(i,k,l_om_ni)
     end do
   end do
!   call outfld('CONDH2S4     ',condh2so4        ,pcols   ,lchnk   )
!   call outfld('CONDS4       ',conds4           ,pcols   ,lchnk   )
!   call outfld('CONDBC       ',condbc           ,pcols   ,lchnk   )
!   call outfld('CONDPOM      ',condom           ,pcols   ,lchnk   )
!   call outfld('CONDBCP      ',condbcom         ,pcols   ,lchnk   )
!   call outfld('CONDBCK      ',condbck         ,pcols   ,lchnk   )
!   call outfld('CONDAS4       ',condaits4           ,pcols   ,lchnk   )
!   call outfld('CONDABC       ',condaitbc           ,pcols   ,lchnk   )
!   call outfld('CONDABCP       ',condaitbcp           ,pcols   ,lchnk   )
!   call outfld('CONDAPOM      ',condaitpom           ,pcols   ,lchnk   )
!   call outfld('CONDBCAX      ',condbcax           ,pcols   ,lchnk   )
!   call outfld('CONDS4PR      ',conds4pr           ,pcols   ,lchnk   )

!     do k=1,pver
!       do i=1 ,ncol
!           if (colcoag(i,k).gt.0) &
!          write(6,*) i,k,colcoag(i,k)
!       end do
!     end do
! !  do m=1,ppcnst
!     if(species_class(m)==spec_class_aerosol) then
!
!            dotend(m)=.true.      

          dotend(l_so4_n)=.true.
          dotend(l_so4_a1)=.true.
!          dotend(l_h2so4)=.true.
          dotend(l_so4_na)=.true.

          dotend(l_bc_n)=.true.
          dotend(l_bc_ni)=.true.
          dotend(l_bc_a)=.true.
          dotend(l_bc_ai)=.true.
#ifdef CONDFIX
          dotend(l_bc_ax)=.true.
#endif 
!          dotend(l_om_n)=.true.
          dotend(l_om_ni)=.true.
!          dotend(l_om_a)=.true.
          dotend(l_om_ai)=.true.

          do k=1,pver
             do i=1,ncol
#ifdef CONDFIX
!                h2so4new(i,k)=q(i,k,l_h2so4)*exp(-condh2so4(i,k)*dt)
!                dh2so4(i,k)=(q(i,k,l_h2so4)-h2so4new(i,k))/dt
! OS 061015 Replaced Euler forward with Euler backward                
! Old forward h2so4new(i,k)=loch2so4(i,k)-condh2so4(i,k)*dt*loch2so4(i,k)
                   
                h2so4new(i,k)=loch2so4(i,k)/(1._r8+condh2so4(i,k)*dt)
#else 
!                h2so4new(i,k)=q(i,k,l_h2so4)*exp(-condh2so4(i,k)*dt)
!                dh2so4(i,k)=(q(i,k,l_h2so4)-h2so4new(i,k))/dt
                h2so4new(i,k)=loch2so4(i,k)-condh2so4(i,k)*dt*loch2so4(i,k)
!                h2so4new(i,k)=q(i,k,l_h2so4)-condh2so4(i,k)*dt*q(i,k,l_h2so4)
                h2so4new(i,k)=max(h2so4new(i,k),0._r8)
#endif 
                dh2so4(i,k)=(loch2so4(i,k)-h2so4new(i,k))/dt
                dqdt(i,k,l_so4_a1) = dh2so4(i,k)


#ifdef CONDFIX 
!              mode number 0, bc_ax

               area_core= nrmodes(i,k,l_bc_ax)*normnsurf(0) 
               volume_shell=dt*condbcax(i,k)*h2so4new(i,k)/rhopart(l_so4_a1)
               volume_monolayer=area_core*dr_so4_monolayers_age
               frac_transfer=min((volume_shell/max(volume_monolayer,1.e-30_r8)),0.999_r8)
               dqdt(i,k,l_bc_ax)=-frac_transfer*q(i,k,l_bc_ax)/dt

!              mode 1 so4_n
             
               area_core= nrmodes(i,k,l_so4_n)*normnsurf(1) 
               volume_shell=dt*conds4(i,k)*h2so4new(i,k)/rhopart(l_so4_a1)
               volume_monolayer=area_core*dr_so4_monolayers_age
               frac_transfer=min((volume_shell/max(volume_monolayer,1.e-30_r8)),0.999_r8)
               dqdt(i,k,l_so4_n)=(h2so4new(i,k)-frac_transfer*q(i,k,l_so4_n))/dt
               dqdt(i,k,l_so4_na)=frac_transfer*q(i,k,l_so4_n)/dt

!              mode 2 bc_n

               area_core= nrmodes(i,k,l_bc_n)*normnsurf(2) 
               volume_shell=dt*condbc(i,k)*h2so4new(i,k)/rhopart(l_so4_a1)
               volume_monolayer=area_core*dr_so4_monolayers_age
               frac_transfer=min((volume_shell/max(volume_monolayer,1.e-30_r8)),0.999_r8)
               dqdt(i,k,l_bc_n)=-frac_transfer*q(i,k,l_bc_n)/dt


!               dqdt(i,k,l_bc_a)=-dqdt(i,k,l_bc_ax)-dqdt(i,k,l_bc_n)
! OS 161015 Started by converting bc_ax into bc_a since both are fossil fuel.
! bc_a on the other hand has bc_n as basis thus too small. Moved into bc_ai

               dqdt(i,k,l_bc_a)=-dqdt(i,k,l_bc_n)


! Mode 3   bc_ni+om_ni

               area_core= (nrmodes(i,k,l_bc_ni)+nrmodes(i,k,l_om_ni))*normnsurf(3) 
               volume_shell=dt*condbcom(i,k)*h2so4new(i,k)/rhopart(l_so4_a1)
               volume_monolayer=area_core*dr_so4_monolayers_age
               frac_transfer=min((volume_shell/max(volume_monolayer,1.e-30_r8)),0.999_r8)
               dqdt(i,k,l_bc_ni)=-frac_transfer*q(i,k,l_bc_ni)/dt
               dqdt(i,k,l_bc_ai)=-dqdt(i,k,l_bc_ni)-dqdt(i,k,l_bc_ax)

               dqdt(i,k,l_om_ni)=-frac_transfer*q(i,k,l_om_ni)/dt
               dqdt(i,k,l_om_ai)=-dqdt(i,k,l_om_ni)

#else 

                dqdt(i,k,l_so4_n) =  -conds4(i,k)*q(i,k,l_so4_n)+ &
                 h2so4new(i,k)/dt
                dqdt(i,k,l_so4_na) = conds4(i,k)*q(i,k,l_so4_n)
        
!               dqdt(i,k,l_h2so4)=-dh2so4(i,k)-h2so4new(i,k)/dt

                dqdt(i,k,l_bc_n)  =  -condbc(i,k)*q(i,k,l_bc_n)
                dqdt(i,k,l_bc_a)  =   condbc(i,k)*q(i,k,l_bc_n)
!                dqdt(i,k,l_om_n)  =  -condom(i,k)*q(i,k,l_om_n)
!                dqdt(i,k,l_om_a)  =   condom(i,k)*q(i,k,l_om_n)
                dqdt(i,k,l_bc_ni) =  -condbcom(i,k)*q(i,k,l_bc_ni)
                dqdt(i,k,l_bc_ai) =   condbcom(i,k)*q(i,k,l_bc_ni)
                dqdt(i,k,l_om_ni) =  -condbcom(i,k)*q(i,k,l_om_ni)
                dqdt(i,k,l_om_ai) =   condbcom(i,k)*q(i,k,l_om_ni)
!+ &
!                condom(i,k)*q(i,k,l_om_n)
        
#endif 

                nuclso4(i,k)=h2so4new(i,k)/dt
             end do
          end do     

      call outfld('NUSO4N ',nuclso4   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), nuclso4,cmi2d)
!      call outfld('PSO4N ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_a1),cmi2d)
!      call outfld('PSO4A1 ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_so4_na),cmi2d)
!      call outfld('PSO4NA ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_a),cmi2d)
!      call outfld('PBCA ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), -dqdt(:,:,l_om_n),cmi2d)
!      call outfld('POMA ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_bc_ai),cmi2d)
!      call outfld('PBCAI ',cmi2d   ,pcols   ,lchnk) 
!      call cmidry(lchnk,ncol,pdel,q(:,:,1), dqdt(:,:,l_om_ai),cmi2d)
!      call outfld('POMAI ',cmi2d   ,pcols   ,lchnk) 

   return
   end
