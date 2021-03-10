! CAM-Oslo Ø Seland Calculates aqueous phase chemistry and scavenging of SO2
      subroutine aqgas(lchnk,   ncol,    dt,       &
                   qstate,       t,       pmid,    cldfrc,       &
                   pdel,    conicw,           &
                   cmfdqr,  precs,   conds,   evaps,             &
                   wetdepflx,                  &
                   aqprod,  fracis,               &
                   lo3,lh2o2,totcond,ixcldliq) 



use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
!use physconst, only: avogadro,rair,mwdry
use constituents, only: pcnst
use mass,         only: cmidry
use cam_history,  only: outfld
use aerosoldef

implicit none
 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
!   integer,  intent(in) :: ncnst                ! number of tracers 
   real(r8), intent(in) :: dt                   ! Time step
   real(r8), intent(inout) :: qstate(pcols,pver,pcnst) ! TMR including moisture
   real(r8), intent(in) :: t(pcols,pver)        ! Temperature in Kelvin
   real(r8), intent(in) :: pmid(pcols,pver)     ! Air pressure in Pa
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction

!   real(r8), intent(in) :: icwmr1 (pcols,pver)  ! in cloud water mixing ration for zhang scheme
!   real(r8), intent(in) :: icwmr2 (pcols,pver)  ! in cloud water mixing ration for hack  scheme
   real(r8), intent(in) :: pdel(pcols,pver)     ! pressure between layers
   real(r8), intent(in) :: conicw(pcols,pver)   ! convective cloud water
   real(r8), intent(in) :: cmfdqr(pcols,pver)   ! rate of production of convective precip
   real(r8), intent(in) :: precs(pcols,pver)    ! rate of production of stratiform precip
   real(r8), intent(in) :: conds(pcols,pver)    ! rate of production of condensate
   real(r8), intent(in) :: evaps(pcols,pver)    ! rate of evaporation of precip

   real(r8), intent(inout) :: wetdepflx(pcols,pcnst) ! constituent wetdep fluxes (kg/m2/s)
   real(r8), intent(inout) :: aqprod(pcols) ! aqueous production (kg/m2/s)
!   real(r8), intent(out) :: dqdt(pcols,pver,pcnst)  ! TMR tendency array
!   logical,  intent(out) :: dotend(pcnst)   
   real(r8), intent(inout) :: fracis(pcols,pver,pcnst)
   real(r8), intent(in) :: lo3(pcols,pver) ! O3 concentration
   real(r8), intent(in) :: lh2o2(pcols,pver) ! H2O2 concentration
   real(r8), intent(in)  :: totcond(pcols,pver) ! Total condensate
   integer , intent(in) :: ixcldliq ! Number of constituents in aerosol parameterisation
! local

   real(r8) :: khet(pcols,pver)
   real(r8) :: s2tend(pcols,pver),s4tend(pcols,pver)
   real(r8) :: scavt(pcols,pver),evsa2(pcols,pver),scavso2(pcols,pver)
   integer :: i,k 
!   real(r8) :: fice(pcols,pver)     ! Fraction of cwat that is ice. Copied from Iversen and Seland 2003
   real(r8) :: frach2(pcols,pver)
!   real(r8) :: totcond(pcols,pver)
   real(r8) :: dh2o2(pcols,pver)

   s2tend(:,:)=0._r8
   s4tend(:,:)=0._r8
!   totcond(:,:)=0._r8
   dh2o2(:,:)=0._r8


   call cldy_camuio(lchnk,   ncol,    pcnst,          &
                 qstate(:,:,l_so2), t,       pmid,    cldfrc,     &
                 totcond, conicw,  khet,&
                 lo3,         qstate(:,:,l_qh2o2) ,frach2) 



   call outfld('KHET'   , khet    , pcols, lchnk)
   call outfld('CH2O2'   , lh2o2    , pcols, lchnk)



!   call outfld('MAXPR'  , maxprec    , pcols, lchnk)
!   call outfld('FPREC'  , fprec    , pcols, lchnk)
   call  wetdeps2(lchnk,   ncol,    dt,       &
                 qstate(:,:,l_so2), pdel,    cldfrc,  conicw ,&
                 qstate(:,:,ixcldliq),   cmfdqr,  precs,   conds,    &
                 evaps,  khet,     &
                 fracis(:,:,l_so2), scavt,evsa2)     
    do k=1,pver
      do i=1,ncol

        s2tend(i,k)=-khet(i,k)*cldfrc(i,k)*qstate(i,k,l_so2)
        s4tend(i,k)=max((khet(i,k)*cldfrc(i,k)*qstate(i,k,l_so2)+scavt(i,k)),0._r8) & 
           +evsa2(i,k)
        scavso2(i,k)=s4tend(i,k)+s2tend(i,k)
        dh2o2(i,k)=frach2(i,k)*s2tend(i,k)
        dh2o2(i,k)=max(dh2o2(i,k),-qstate(i,k,l_qh2o2)/dt)
      end do
    end do   

       call cmidry( lchnk,ncol,pdel, qstate(:,:,1), scavso2, wetdepflx(:,l_so2))    
       call cmidry( lchnk,ncol,pdel, qstate(:,:,1), s4tend, aqprod(:))

   do k=1,pver
     do i=1,ncol
       qstate(i,k,l_so2)=qstate(i,k,l_so2)+s2tend(i,k)*dt
       qstate(i,k,l_so4_a2)=qstate(i,k,l_so4_a2)+s4tend(i,k)*dt
       qstate(i,k,l_qh2o2)=qstate(i,k,l_qh2o2)+dh2o2(i,k)*dt
     end do
   end do

   return
   end
