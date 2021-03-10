! Calculates DMS and SO2 gas chemistry rates
! To a large eextent build on a simular rooutine in MATCH
! Ø Seland 2010


!-----------------------------------------------------------------------
subroutine gaschem(lchnk,   ncol,    ncnst,   dt,       &
                   q,       t,       pmid,    landfrac, &
                   pdel, condprod,   s2prod,  msaprod,  &
                   dqdt,    dotend,  loh,h2so4_gasprod                      )

!-----------------------------------------------------------------------
!
! computes TMR (tracer mixing ratio) tendencies for gas phase chemistry
!    by "integrating" the chemistry ODE's then 
!
! this version does sulfur chemistry only using offline oxidants
!    DMS oxidation to MSA and SO2
!    SO2 oxidation to H2SO4
!
! the gas-phase reaction coding follows that used in MIRAGE-1
!
! the real work is done in subr chemtend
!
!-----------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use physconst,    only: rair,mwdry,avogad
use mass,         only: cmidry
use aerosoldef

implicit none

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   integer,  intent(in) :: ncnst                ! number of tracers 
   real(r8), intent(in) :: dt                   ! Time step
   real(r8), intent(in) :: q(pcols,pver,ncnst) ! TMR including moisture
                                                ! (TMR = tracer mixing ratio)
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
   real(r8), intent(in) :: landfrac(pcols)      ! Land area fraction
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure between layers 
   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! TMR tendency array
   real(r8), intent(inout) :: condprod(pcols) ! condensation production (kg/m2/s)
   real(r8), intent(inout) :: s2prod(pcols) ! dms-so2
   real(r8), intent(inout) :: msaprod(pcols) ! dms-so2
   logical,  intent(out) :: dotend(ncnst)   ! identifies the species that
                                             ! tendencies are computed for


   real(r8), intent(in) :: loh(pcols,pver) ! OH concentration
   real(r8), intent(inout) :: h2so4_gasprod(pcols,pver) ! h2so4 produced in one timestep

!   real(r8), intent(out) :: h2so4_gasprod(pcols,pver,2)
                                             ! h2so4 production tendency

! local 
   integer :: i, k, m
   real(r8) :: cavogadro
  real(r8) :: caircell(pcols,pver),aircon(pcols,pver),rhoair(pcols,pver)
  real(r8) :: rk(pcols,pver,3)
  real(r8) :: so4prod(pcols,pver)
  real(r8) :: dms2msa(pcols,pver)
  real(r8) :: kD12(pcols,pver)

! set tendency flags
   dotend(:) = .false.
   dotend(l_dms) = .true.
   dotend(l_so2) = .true.
!   dotend(l_msa) = .true.
!   dotend(l_h2so4) = .true.
   dotend(l_om_ni)=.true.
   dotend(l_so4_n) = .true.
   dotend(l_so4_a1) = .true.

!   dotend(l_h2so4) = .true.
!   if (l_h2o2 > 0) dotend(l_h2o2) = .true.

! zero out tendencies (this is probably not necessary!!)
   dqdt(:,:,2:ncnst) = 0.0_r8


! Convert constants to cgs prefix c
       cavogadro = avogad/1000._r8
        
! Nomenclature:
!    t         = ambient temperature (K)
!    rhoair    = air density (kilograms/m^3 space)
!    aircon    = total air molar concentration (moles air/cm^3 space)
!    caircell  = total air molecular concentration (molecules/cm^3 space)
!    avogadro  = molecules/kmol
!    cavogadro = molecules/mol

   do k=1,pver
      do i=1,ncol   
         rhoair(i,k)=pmid(i,k)/(rair*t(i,k))
! Going from kgm to gcm
         aircon(i,k)=0.001_r8*rhoair(i,k)/mwdry
 	caircell(i,k) = aircon(i,k) * cavogadro
       end do
   end do

   call rxnrate( lchnk, ncol, t ,caircell,rk)         

!   call condrate(lchnk, ncol, kD12)

   call chemtend( lchnk, ncol, ncnst, landfrac, q, dqdt, dt ,rk,loh,so4prod,dms2msa)

!   h2so4_gasprod(1:ncol,:,1) = dqdt(1:ncol,:,l_h2so4)
!   h2so4_gasprod(1:ncol,:,2) =    q(1:ncol,:,l_h2so4)
   do k=1,pver
     do i=1,ncol
       h2so4_gasprod(i,k)=so4prod(i,k)*dt
     end do
   end do     

   call cmidry( lchnk,ncol,pdel, q(:,:,1), so4prod(:,:), condprod(:))

   call cmidry( lchnk,ncol,pdel, q(:,:,1), -dqdt(:,:,l_dms)-dms2msa(:,:), s2prod(:))

   call cmidry( lchnk,ncol,pdel, q(:,:,1), dms2msa(:,:), msaprod(:))

   return
   end subroutine gaschem

!-----------------------------------------------------------------------
subroutine rxnrate( lchnk, ncol, te,caircell,rk )
!-----------------------------------------------------------------------
!
! computes gas phase reaction rate constants
!
!-----------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid

implicit none
!include 'chemrates.h'

   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric column
   real(r8), intent(in) :: te(pcols,pver)  ! Temperature in Kelvin
   real(r8), intent(in) :: caircell(pcols,pver) 
   real(r8), intent(out) :: rk(pcols,pver,3)  
! local
        
        real(r8) :: troe
        external troe
 
        integer  :: i,k
        real(r8) :: rk0
        real(r8) :: rki
        real(r8) :: fc
        real(r8) :: rk0num
        real(r8) :: rk0den

   do k=1,pver
      do i=1,ncol   
!  sulfur reactions
          rk0 = 3.0e-31_r8*caircell(i,k)*(te(i,k)/300._r8)**(-3.3_r8)
          rki = 1.5e-12_r8
          fc  = 0.6_r8
        rk(i,k,1) = troe(rk0,rki,fc)

        rk(i,k,2) = 1.2e-11_r8*exp(-260.0_r8/te(i,k))
          rk0num = te(i,k)*exp(-234.0_r8/te(i,k)) + 8.46e-10*exp(+7230.0_r8/te(i,k)) &
                    + 2.68e-10_r8*exp(+7810.0_r8/te(i,k))
          rk0den = 1.04e+11_r8*te(i,k) + 88.1_r8*exp(+7460.0_r8/te(i,k))
          rk0 = rk0num/rk0den
          rk(i,k,3) = max( 0.0_r8, rk0 - rk(i,k,2) )
       end do
     end do
     return
     end subroutine rxnrate


!---------------------------------------------------------------------
!      FUNCTION TROE - calculates three-body rate coefficients,
!                      troe, given values of k0, ki, and Fc.
!
!      History:  Created Fall, 1992 by Rick D. Saylor 
!
!---------------------------------------------------------------------
       function troe(rk0,rki,fc)
use shr_kind_mod, only: r8 => shr_kind_r8
implicit none

       real(r8) troe
! rki, rk0, fc - parameters in Troe expression for 3-body reaction
       real (r8) rk0,rki,fc 
       real(r8) rn,rt,arg,f

       rn = 0.75_r8 - 1.27_r8*log10(fc)
       rt = rk0/rki
       arg = (log10(rt)/rn)*(log10(rt)/rn) + 1._r8
       f = fc**(1.0_r8/arg)
       troe = rk0*rki*f/(rk0+rki)


       return
       end function troe


!-----------------------------------------------------------------------
subroutine chemtend( lchnk, ncol, ncnst, landfrac,q, dqdt, dtc,rk ,loh,so4prod,msaprod)
!-----------------------------------------------------------------------
!
! this routine actually computes TMR (tracer mixing ratio) tendencies 
!    for gas phase chemistry
!
! this version does sulfur chemistry only using offline oxidants
!    DMS oxidation to MSA and SO2
!
!-----------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use time_manager, only: get_nstep
use aerosoldef

implicit none
!include 'chemrates.h'  
!include 'chemspecies.h' 

   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric column
   integer, intent(in) :: ncnst                 ! number of tracers
   real(r8), intent(in) :: landfrac(pcols)      ! Land area fraction  
   real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array 
   real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array
   real(r8), intent(in) :: dtc                  ! Time step for chemistry
   real(r8), intent(in) :: rk(pcols,pver,3)   
   real(r8), intent(in) :: loh(pcols,pver) ! OH concentration
   real(r8), intent(out) :: so4prod(pcols,pver)
   real(r8), intent(out) :: msaprod(pcols,pver)

! local
   integer  :: i, k, nstep
   real(r8) :: ohp, alpha, indtc, invair
   real(r8) :: pdms, ddms, pso2, dso2, ph2so4, dh2so4,pmsa, dmsa,pso4
   real(r8) :: dmsnew, dmsold, so2new, so2old
!   real(r8) :: h2so4old, h2so4new
!   real(r8) :: msaold, msanew
   real(r8) :: fdiurn_oh(pcols)
   real(r8) :: fno3(pcols,pver)
   real(r8) :: fso4nuc
   real(r8), parameter :: tnudgeh2o2_inv = 1.0_r8/(10.0_r8*86400.0_r8)
   real(r8), parameter :: dirchem = 2.8e-4_r8
   real(r8),parameter :: gaslife = 2.8e-4_r8
   

   call oxi_diurnal_var( lchnk, ncol, dtc, landfrac, fdiurn_oh ,fno3)

   nstep = get_nstep()

! alpha = MSA yield from rxn 39
   alpha=0.4_r8

   indtc=1._r8/dtc

   fso4nuc=0.05_r8
   do k=1,pver
     do i=1,ncol
        ohp = loh(i,k)*fdiurn_oh(i)
! DMS
        dmsold=q(i,k,l_dms)
        dmsnew = dmsold/(1._r8+dtc*(rk(i,k,2)+rk(i,k,3))*ohp)
        dmsnew = dmsnew - dtc*fno3(i,k)*dirchem*dmsold
        dqdt(i,k,l_dms) = (dmsnew-dmsold)*indtc
! SO2
        so2old = q(i,k,l_so2)
        pso2 = (rk(i,k,2)+alpha*rk(i,k,3))*dmsnew*ohp
        dso2 = rk(i,k,1)*ohp
        so2new = (so2old + pso2*dtc)/(1._r8 + dso2*dtc)
        so2new = so2new + dtc*fno3(i,k)*dirchem*dmsold
        dqdt(i,k,l_so2)=(so2new-so2old)*indtc


! H2SO4
!        h2so4old = q(i,k,l_h2so4)
        ph2so4 = rk(i,k,1)*so2new*ohp
        so4prod(i,k)=ph2so4
!        dh2so4 = 0._r8
!gaslife        
!        h2so4new = (h2so4old + ph2so4*dtc)/(1._r8 + dh2so4*dtc)
!        dqdt(i,k,l_h2so4) = (h2so4new-h2so4old)*indtc
! SO4

!        pso4 = gaslife*h2so4new
!        dqdt(i,k,l_so4_n)=fso4nuc*pso4
!        dqdt(i,k,l_so4_a1)=(1._r8-fso4nuc)*pso4
         
!        h2so4old=q(i,k,l_h2so4)        
!        ph2so4 = rk(i,k,1)*so2new*ohp
!        h2so4new = h2so4old + ph2so4*dtc
!        dqdt(i,k,l_h2so4) = (h2so4new-h2so4old)*indtc
! MSA
!        msaold= q(i,k,l_msa)
        pmsa = (1._r8 - alpha)*rk(i,k,3)*dmsnew*ohp
!        msanew = msaold + pmsa*dtc
!        dqdt(i,k,l_h2so4)=dqdt(i,k,l_h2so4)+pmsa
        msaprod(i,k)=pmsa

! MSA=CH3SO3H  Atmoic weight 96
        dqdt(i,k,l_om_ni)= 3._r8*pmsa
!        dqdt(i,k,l_msa) = (msanew-msaold)*indtc

! H2O2 - only if it is prognosed (in which case l_h2o2 > 0)
!    when nstep==0, set tendency to bring the MR to lh2o2(i,k)
!       (this should be moved to the initialization section)
!    when nstep>0, tendency = gas-phase production term (offline)
!       plus weak nudging toward lh2o2(i,k)
!
!        if (l_h2o2 > 0) then
!           if (nstep == 0) then
!              dqdt(i,k,l_h2o2) = (lh2o2(i,k) - q(i,k,l_h2o2))*indtc
!           else
!              dqdt(i,k,l_h2o2) = max( lgprh2o2(i,k)*fdiurn_oh(i), 0.0_r8 )   &
!                               + (lh2o2(i,k) - q(i,k,l_h2o2))*tnudgeh2o2_inv
!           end if
!        end if

     end do    ! "i = 1, ncol"
   end do      ! "k = 1, pver"

return
end subroutine chemtend



!-----------------------------------------------------------------------
subroutine oxi_diurnal_var( lchnk, ncol, dtc, landfrac, fdiurn_oh ,fno3)
!-----------------------------------------------------------------------
!
! this routine computes diurnal modulation factors to be applied to
!    the offline OH MRs and the offline H2O2 production rates
!
! the current (very simple) algorithm sets the factors to zero during night
!    and a constant value during day, such that the daily average is 1.0
!    (i.e., a step function)
!
! the sunrisesetxx routine (part of the MIRAGE-1 gas chemistry) is used
!    to obtain sunrise and set times
!
!-----------------------------------------------------------------------
use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use physconst, only: pi
use phys_grid, only: get_rlat_all_p, get_rlon_all_p
use time_manager, only: get_curr_date

implicit none

   integer, intent(in) :: lchnk             ! chunk identifier
   integer, intent(in) :: ncol              ! number of atmospheric column
   real(r8), intent(in) :: dtc              ! Time step for chemistry

   real(r8), intent(in) :: landfrac(pcols)      ! Land area fraction 

   real(r8), intent(out) :: fdiurn_oh(pcols)   ! diurnal modulation factor 

   real(r8), intent(out) :: fno3(pcols,pver) ! Rapid reaction with no3 at night
                                               ! for offline oh MR

! local
   integer :: i                            ! column index
   integer :: iriseset                     ! sunrise/set flag
   integer :: iday, imon, iyr, jyr         ! date stuff
   integer :: j                            ! working var
   integer :: ncsec                        ! time stuff

   real(r8) :: deglat, deglon              ! lat and long (degrees)
   real(r8) :: solardec                    ! solar declination (degrees)
   real(r8) :: sum                         ! working vars
   real(r8) :: trise, tset                 ! sunrise and set times (h then d)
   real(r8) :: tlight                      ! amount of daylight (d)
   real(r8) :: trisej, tsetj               ! working vars
   real(r8) :: t1, t2, ta, tb              ! working vars
   real(r8) :: rlats(pcols), rlons(pcols)  ! latitude & longitude (radians)



  fno3(:,:)=0._r8

! adjust century so that year is between 1950 and 2049, as the "solar"
! routines may balk at other years
   call get_curr_date( iyr, imon, iday, ncsec )
   jyr = mod( iyr, 100 ) + 1900
   if (jyr < 1950) jyr = jyr + 100
   if (jyr > 2049) jyr = jyr - 100

! get lat and lon
   call get_rlat_all_p( lchnk, ncol, rlats )
   call get_rlon_all_p( lchnk, ncol, rlons )

   do i = 1, ncol
      deglat = rlats(i)*180.0_r8/pi
      deglat = max( -89.9999_r8, min( +89.9999_r8, deglat ) )
      deglon = rlons(i)*180.0_r8/pi

! get sunrise and sunset times in UTC hours
      call sunrisesetxx( deglon, deglat, jyr, imon, iday,   &
                iriseset, trise, tset, solardec )



! convert rise/set times to days
! compute tlight = amount of daylight
! handle case of all day or night
      if (iriseset > 0) then
         trise = trise/24.0_r8
         tset  = tset/24.0_r8
         tlight = tset - trise
         if (tlight < 0.0_r8) then
            tset = tset + 1.0_r8
            tlight = tlight + 1.0_r8
         end if
      else
         trise = 0.0_r8
         if (abs(deglat+solardec) .ge. 90.0_r8) then
            tset = 1.0_r8
         else
            tset = 0.0_r8
         end if
         tlight = tset - trise
      end if

! if all day or all night (or very close to it), set fdiurn = 1.0_r8
      if ((tlight .ge. 0.99_r8) .or. (tlight .le. 0.01_r8)) then
         fdiurn_oh(i) = 1.0_r8
! otherwise determine overlap between current timestep and daylight times
! to account for all overlap possibilities, need to try this 
! with rise/set times shifted by +/- 1 day 
      else
         t1 = ncsec/86400.0_r8
         t2 = t1 + dtc/86400.0_r8
         sum = 0.0_r8
         do j = -1, 1
            trisej = trise + j
            tsetj  = trisej + tlight
            ta = max( t1, trisej )
            tb = min( t2, tsetj )
            sum = sum + max( tb-ta, 0.0_r8 )/tlight
         end do
         fdiurn_oh(i) = sum/(t2-t1)
      end if
   
!  NO3+DMS direct reaction over NH land and at nighttime 
 
   if ((fdiurn_oh(i).lt.0.5_r8.or.tlight.le.0.01_r8).and.deglat.gt.10._r8) then
   
   fno3(i,pver) = 1._r8*landfrac(i)
   fno3(i,pver-1) = 0.75_r8*landfrac(i)
   fno3(i,pver-2) = 0.5_r8*landfrac(i)
   fno3(i,pver-3) = 0.25_r8*landfrac(i)


   end if

   end do   ! "i = 1, ncol"

return
end subroutine oxi_diurnal_var





