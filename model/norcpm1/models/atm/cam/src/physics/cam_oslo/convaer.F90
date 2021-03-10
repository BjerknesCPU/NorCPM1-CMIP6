	subroutine convaer(lchnk,ncol,k,rho,qm,Cnso4,Cas75,Cnbc,Cnoc,&
!#define DIAGNCDNC
#ifdef DIAGNCDNC
                           Cabce,Caintmix, f_c,f_bc,f_aq,Nnatk,fnbc,faitbc)
#else
                           Cabce,Caintmix, f_c,f_bc,f_aq,Nnatk,fnbc,faitbc, &
!TS++
		Caitso4, Caitbc, Caitocbc, rhobcocait, C_dst2, C_dst3, C_ss1, C_ss2, C_ss3, &
		rhobcocn, Cnbcioc, Cnocibc)
!TS--
#endif

!
!	Converting aerosol mixing ratios into mass concentrations appropriate  
!       for the CCN-calculations. Original version (with monthly input data)
!       created by Jon Egill Kristjansson on 17 November 1998. Rewritten for 
!       CAM2 with interactive life cycle scheme for (so4, BC and OC) by 
!       Trude Storelvmo, 2004/2005, and for CAM3 (Cam-Oslo, including prognostic 
!       sea-salt and dust aerosols, and number concentrations) by Alf Kirkevaag, 
!       August 2005. Modified for new aerosol schemes by Alf Kirkevaag in January 
!       2006.


        use shr_kind_mod, only: r8 => shr_kind_r8
	use pmgrid
	use ppgrid
        use constituents, only: pcnst, cnst_get_ind
        use const
        use opttab
        use aerosoldef
        use physconst, only: pi


	implicit none

!#include <chemspecies.h>
!#include <depvar.h>

        integer, intent(in) :: lchnk              ! chunk identifier
        integer, intent(in) :: ncol               ! number of atmospheric columns
        integer, intent(in) :: k                  ! model level

	real(r8), intent(in) :: rho(pcols)
        real(r8), intent(in) :: qm(pcols,pver,pcnst) ! Common aerosol (only!) tracers for indirect and direct calculations

        real(r8), intent(out) :: Cnso4(pcols)	  ! SO4(n) (ug/m3)	
        real(r8), intent(out) :: Cas75(pcols)     ! SO4(Ait75) (ug/m3)
        real(r8), intent(out) :: Cnbc(pcols)	  ! BC(n) (ug/m3)	
        real(r8), intent(out) :: Cnoc(pcols)      ! OC(n) (ug/m3)
        real(r8), intent(out) :: Cabce(pcols)	  ! BC(ax) (ug/m3)
        real(r8), intent(out) :: fnbc(pcols)      ! = Cbc/(Cbc+Coc) for BC&OC(n)
        real(r8), intent(out) :: faitbc(pcols)    ! = Cbc/(Cbc+Coc) for BC&OC(Ait)
        real(r8) :: Cnbcioc(pcols)          	  ! BC in BC&OC(n) (ug/m3)	
        real(r8) :: Cnocibc(pcols)                ! OC in BC&OC(n) (ug/m3)
        real(r8) :: Caitso4(pcols)	  ! SO4(Ait) (ug/m3)	
        real(r8) :: Caitbc(pcols)        	  ! BC(Ait) (ug/m3)	
        real(r8) :: Caitoc(pcols)                 ! OC(Ait) (ug/m3)
        real(r8) :: Caitbcioc(pcols)              ! BC in BC&OC(Ait) (ug/m3)	
        real(r8) :: Caitocibc(pcols)              ! OC in BC&OC(Ait) (ug/m3)

        real(r8), intent(out) :: Caintmix(pcols)  ! Total added so4, BC and OC (ug/m3)  
        real(r8), intent(out) :: f_c(pcols)       ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
        real(r8), intent(out) :: f_bc(pcols)      ! = Cbc/(Cbc+Coc)
        real(r8), intent(out) :: f_aq(pcols)      ! = Cso4a2/(Cso4a1+Cso4a2+Cso4ac)
!        real(r8), intent(out) :: Nnatk(pcols,nmodes)
        real(r8), intent(out) :: Nnatk(pcols,0:nmodes)

#ifndef DIAGNCDNC
        real(8)  Caitocbc(pcols)
#endif
        real(8)  efact_bcocn(pcols), rhobcocn(pcols)
        real(8)  efact_bcocait(pcols), rhobcocait(pcols)
        real(r8) C_dst1(pcols), C_dst2(pcols), C_dst3(pcols)
        real(r8) C_ss1(pcols) , C_ss2(pcols) , C_ss3(pcols)

	integer i, kcomp

	real(r8) rhofac(pcols)  ! Unit conversion factor (for kg/kg --> ug/m3)
	real(r8) totant(pcols)  ! Int. mixed (cond./coag./aq.) SO4+BC+OC concentration


!       Converting from mass mixing ratios to mass concentrations
	do i=1, ncol

	 rhofac(i) = 1.e9_r8*rho(i)

	 totant(i) = qm(i,k,l_bc_ac)+qm(i,k,l_om_ac) &
	            +3._r8*(qm(i,k,l_so4_a1)+qm(i,k,l_so4_a2)+qm(i,k,l_so4_ac))
	 Caintmix(i) = rhofac(i)*totant(i)

         f_c(i)   = min((qm(i,k,l_bc_ac)+qm(i,k,l_om_ac))/(totant(i)+eps),0.999_r8)
	 f_bc(i)  = min(qm(i,k,l_bc_ac)/(qm(i,k,l_bc_ac)+qm(i,k,l_om_ac)+eps),0.999_r8)
         f_aq(i)  = qm(i,k,l_so4_a2) &
                   /(qm(i,k,l_so4_a1)+qm(i,k,l_so4_a2)+qm(i,k,l_so4_ac)+eps)

	 Cnso4(i)  = rhofac(i)*qm(i,k,l_so4_n)*3._r8   ! SO4(n) mass
	 Cas75(i)  = rhofac(i)*qm(i,k,l_so4_pr)*3._r8  ! SO4(Ait75) mass
	 Cnbc(i)   = rhofac(i)*qm(i,k,l_bc_n)       ! BC(n) mass
	 !Cnoc(i)   = rhofac(i)*qm(i,k,l_om_n)       ! OC(n) mass !mode removed
	 Cnoc(i)   = eps                            ! OC(n) mass
	 Cabce(i)  = rhofac(i)*qm(i,k,l_bc_ax)      ! BC(ax) mass
	 Cnbcioc(i)= rhofac(i)*qm(i,k,l_bc_ni)      ! BC mass for BC&OC(n) mode
	 Cnocibc(i)= rhofac(i)*qm(i,k,l_om_ni)      ! OC mass for BC&OC(n) mode 
         rhobcocn(i)= (Cnbcioc(i)+Cnocibc(i)) &     ! mass density for the BC&OC(n) mode 
                    /(Cnbcioc(i)/rhopart(l_bc_ni)+Cnocibc(i)/rhopart(l_om_ni)+eps)
         fnbc(i)   = Cnbcioc(i)/(Cnbcioc(i)+Cnocibc(i)+eps)         

	 Caitso4(i)= rhofac(i)*qm(i,k,l_so4_na)*3._r8 ! SO4(Ait) mass (without cond. SO4) 
   	 Caitbc(i) = rhofac(i)*qm(i,k,l_bc_a)      ! BC(Ait) mass
   	 !Caitoc(i) = rhofac(i)*qm(i,k,l_om_a)      ! OC(Ait) mass !corinna: mode no longer present
   	 Caitbcioc(i) = rhofac(i)*qm(i,k,l_bc_ai)  ! BC mass for BC&OC(Ait) mode
   	 Caitocibc(i) = rhofac(i)*qm(i,k,l_om_ai)  ! OC mass for BC&OC(Ait) mode
#ifndef DIAGNCDNC
	 Caitocbc(i) = Caitocibc(i)+Caitbcioc(i)   ! TS: Total mass for BC&OC(Ait) mode 
#endif
         rhobcocait(i)= (Caitbcioc(i)+Caitocibc(i)) &  ! mass density for the BC&OC(Ait) mode 
                    /(Caitbcioc(i)/rhopart(l_bc_ni)+Caitocibc(i)/rhopart(l_om_ni)+eps)
         faitbc(i) = Caitbcioc(i)/(Caitbcioc(i)+Caitocibc(i)+eps)

!        the dst_1 mode does not exist in the aerocomB emissions:
         C_dst1(i) = 0.0_r8 
         C_dst2(i) = qm(i,k,l_dst_a2)*rhofac(i)
         C_dst3(i) = qm(i,k,l_dst_a3)*rhofac(i)
         C_ss1(i)  = qm(i,k,l_ss_a1)*rhofac(i)
         C_ss2(i)  = qm(i,k,l_ss_a2)*rhofac(i)
         C_ss3(i)  = qm(i,k,l_ss_a3)*rhofac(i)

	end do

        do i=1,ncol
!       Number concentration for the externally mixed fractal BC a-mode
          Nnatk(i,0) = Cabce(i)*efact_bcax
!       Number concentrations for the internally mixed SO4, BC and OC Aitken modes
          Nnatk(i,1) = Caitso4(i)*efact_so4na
          Nnatk(i,2) = Caitbc(i)*efact_bca
          !Nnatk(i,3) = Caitoc(i)*efact_oma !corinna: mode no longer present
          Nnatk(i,3) =  0._r8 
          efact_bcocait(i) = 1.e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_ai)))**2) &
                         *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ai))**3*rhobcocait(i)+eps) 
          Nnatk(i,4) = (Caitbcioc(i)+Caitocibc(i))*efact_bcocait(i)
!       Number concentration for the internally mixed SO4(Ait75) mode
!          Nnatk(i,5) = 0.0_r8                                                           ! SO4(Ait75)=0 naa !!!!!!!!!
          Nnatk(i,5) = Cas75(i)*efact_so4pr
!       Number concentrations for the internally mixed sea-salt and mineral modes
          Nnatk(i,6) = C_dst2(i)*efact_dst2
          Nnatk(i,7) = C_dst3(i)*efact_dst3
          Nnatk(i,8) = C_ss1(i)*efact_ss1
          Nnatk(i,9) = C_ss2(i)*efact_ss2
          Nnatk(i,10)= C_ss3(i)*efact_ss3
!         and for the externally mixed SO4, BC and OC n-modes
          Nnatk(i,11)= Cnso4(i)*efact_so4n
          Nnatk(i,12)= Cnbc(i)*efact_bcn
          Nnatk(i,13)= Cnoc(i)*efact_omn
          efact_bcocn(i) = 1.e-15_r8/(e**(4.5_r8*(log(sgpart(l_om_ni)))**2) &
                         *(4.0_r8/3.0_r8)*pi*(effsize(l_om_ni))**3*rhobcocn(i)+eps) 
          Nnatk(i,14)= (Cnbcioc(i)+Cnocibc(i))*efact_bcocn(i)
        end do

!        do kcomp=0,14 
!         do i=1,ncol
!          if(Nnatk(i,kcomp).ne.0._r8) &
!           write(*,*) 'lon, lev, Nnatk', i, k, Nnatk(i,kcomp)
!         end do
!        end do

!	write(*,*) 'End of convaer'

	return
	end






