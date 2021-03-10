     subroutine wetdeps2(lchnk,   ncol,    dt,                         &
                 tracer, pdel,    cldfrc,   conicw,                    &
                 cwat,   cmfdqr,  precs,   conds,    evaps,  khet,     &
                 fracis, scavt,   evsa2) 
! Wet scavenging of SO2, in-cloud only

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use physconst,     only: gravit
implicit none



! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: dt                   ! Time step
   real(r8), intent(in) :: tracer(pcols,pver)   !  
   real(r8), intent(in) :: pdel(pcols,pver)     ! Dp between layers
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction
   real(r8), intent(in) :: conicw(pcols,pver)   ! convective cloud water
   real(r8), intent(in) :: cwat(pcols,pver)     ! Cloud water
   real(r8), intent(in) :: cmfdqr(pcols,pver)   ! rate of production of convective precip
   real(r8), intent(in) :: precs(pcols,pver)    ! rate of production of stratiform precip
   real(r8), intent(in) :: conds(pcols,pver)    ! rate of production of condensate
   real(r8), intent(in) :: evaps(pcols,pver)    ! rate of evaporation of precip

!   real(r8), intent(in) :: icwmr1 (pcols,pver)  ! in cloud water mixing ration for zhang scheme
!   real(r8), intent(in) :: icwmr2 (pcols,pver)  ! in cloud water mixing ration for hack  scheme
   real(r8), intent(in) :: khet(pcols,pver)  ! Oxidation rate
   real(r8), intent(out) :: fracis(pcols,pver)! fraction of species not scavenged               
   real(r8), intent(out) :: scavt(pcols,pver)    ! scavenging tend 
    real(r8), intent(out) :: evsa2(pcols,pver)

 


!c local variables

      integer :: i                 ! chunk index
      integer :: k                 ! z index

      real(r8) :: fracev               ! fraction of precip from above that is evaporating
      real(r8) :: fracp                ! fraction of cloud water converted to precip
      real(r8) :: omsm                 ! 1 - (a small number)
      real(r8) :: pdog                 ! work variable (pdel/gravit)
      real(r8) :: precabc(pcols)       ! conv precip from above (work array)
      real(r8) :: precabs(pcols)       ! strat precip from above (work array)
      real(r8) :: precbl               ! precip falling out of level (work array)
      real(r8) :: rat                  ! ratio of amount available to amount removed
      real(r8) :: scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) :: scavabc(pcols)       ! scavenged tracer flux from above (work array)
      real(r8) :: srcc                 ! tend for convective rain
      real(r8) :: srcs                 ! tend for stratiform rain
      real(r8) :: srct                 ! work variable
!      real(r8) :: fins                 ! fraction of rem. rate by strat rain
!      real(r8) :: finc                 ! fraction of rem. rate by conv. rain
      real(r8) :: srcs1                ! work variable
      real(r8) :: srcs2                ! work variable
      real(r8) :: cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) :: cldmabc(pcols)       ! maximum cloud at or above this level
      real(r8) :: odds                 ! limit on removal rate (proportional to prec)
      real(r8) :: winclo(pcols,pver)
      real(r8) :: wbeclo(pcols,pver)
! ------------------------------------------------------------------------
      omsm = 1._r8-1.e-10_r8          ! used to prevent roundoff errors below zero



!c this section of code is for highly soluble aerosols,
! OS and also relative highly soluble gases


!c keep track of the amount of stratiform precip falling from above
!c also the amount of tracer falling from above
      scavt(:,:)=0._r8
      evsa2(:,:)=0._r8
      do i = 1,ncol
         precabs(i) = 0._r8
         precabc(i) = 0._r8
         scavab(i) = 0._r8
         scavabc(i) = 0._r8
         cldmabs(i) = 0._r8
         cldmabc(i) = 0._r8
      end do


      do k = 1,pver
         do i = 1,ncol

!            winclo(i,k) = 1._r8
	    winclo(i,k) = 2.e6_r8*khet(i,k)/1000._r8
            winclo(i,k) = min(winclo(i,k),0.5_r8)
	    wbeclo(i,k) = 0._r8
!            weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice
!c	    if (winclo(i,k).gt.0._r8) 
!c     &	write(6,*) 'wetdeps2 ',i,k,weight,winclo(i,k),khet(i,k)

            pdog = pdel(i,k)/gravit

!c calculate the fraction of strat precip from above 
!c                 which evaporates within this layer
!            fracev = evaps(i,k)*pdel(i,k)/gravit  &
!                 /max(1.e-12_r8,precabs(i))

             fracev = evaps(i,k)*pdel(i,k)/gravit &
                      /max(1.e-12_r8,precabs(i))
             ! trap to ensure reasonable ratio bounds
             fracev = max(0._r8,min(1._r8,fracev))

!c now do the convective scavenging

!c set odds proportional to fraction of the grid box that is swept by the 
!c precipitation =precabc/rhoh20*(area of sphere projected on plane
!c                                /volume of sphere)*dt
!c assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
!c unless the fraction of the area that is cloud is less than odds, in which
!c case use the cloud fraction (assumes precabs is in kg/m2/s)
!c is really: precabs*3/4/1000./1e-3*dt
!c here I use .1 from Balkanski


!c use a local rate of convective rain production for incloud scav
!cos            odds=max(min(1._r8,
!cos     $           cmfdqr(i,k)*pdel(i,k)/gravit*.1*dt),0._r8)


!Test on what happens cmfdqr is grid-values instead of incloud values

              odds = (cmfdqr(i,k)/max(cldfrc(i,k),1.e-4_r8))*dt/max(1.e-8_r8,conicw(i,k))
              odds = max(min(1._r8,odds),0._r8)       
              srcs1 =  winclo(i,k)*cldfrc(i,k)*odds*tracer(i,k)/dt
! end test

!            odds=max(min(1._r8,
!     $           cmfdqr(i,k)*pdel(i,k)/gravit*winclo(i,k)*dt),0._r8)
!            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1._r8-weight)
!     $           /dt 
!!
!
!
!
!c scavenge below layers of precip
            cldmabc(i) = max(cldfrc(i,k),cldmabc(i))

!cos            odds=max(
!cos     $             min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8)*.1*dt),
!cos     $             0._r8)
            odds=max(        &
                  min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8)  &
                  *dt),0._r8)
            srcs2 = wbeclo(i,k)*cldmabc(i)*odds*tracer(i,k)    &
                /dt
            srcc = srcs1 + srcs2
!            write(6,*) i,k,srcs1,srcc,
!            finc = srcs1/(srcc + 1.e-36_r8)
!c            if (lat.eq.34.and.i.eq.101) then
!c               write (6,*) ' conv in wetdepa ' 
!c               write (6,67)
!c     $              k, odds, srcs1, srcs2, precabc(i),
!c     $              precabc(i)/max(cldmabc(i),1.e-5_r8)
!c            endif

!c now do the stratiform scavenging


            cldmabs(i) = max(cldfrc(i,k),cldmabs(i))


!c for layers with significant stratiform rain
!c fracp is the fraction of cloud water converted to precip
            fracp =             &
                max(0._r8,        &
!cos     $               min(precs(i,k)*dt
                    min(precs(i,k)*dt  &
!c     $                   /max(cwat(i,k)+conds(i,k)*dt,1.e-12_r8),
                        /max(cwat(i,k),1.e-12_r8),      &
                        1._r8)                          & 
                   )                


!c assume the corresponding amnt of tracer is removed
!c++mcb -- remove cldc; change cldt to cldv 
            srcs1 = winclo(i,k)*0.3_r8*cldfrc(i,k)*fracp*tracer(i,k)/dt  
                   ! scavenge only the liquid phase
!c--mcb

!c now calc the amount removed by rain from above

!c same calc as above but for strat rain
!cos            odds=max(
!cos     $             min(1._r8,precabs(i)/max(cldmabs(i),1.e-5_r8)*.1*dt),
!cos     $             0._r8)
            odds=max(                                        &
                  min(1._r8,precabs(i)/max(cldmabs(i),1.e-5_r8)   &
                  *dt),0._r8)

!c scavenging due to rain from above

            srcs2 =wbeclo(i,k)*(0.3_r8*cldmabs(i)*odds)          &
                *tracer(i,k)/dt

!c total strat scavenged 
            srcs = srcs1 + srcs2
! last here

!c fraction taken by incloud processes
!            fins=srcs1/(srcs + 1.e-100_r8)

!c make sure we dont take out more than is there
!c ratio of amount available to amount removed

            rat = cldfrc(i,k)*tracer(i,k)/  &  
       max(dt*(srcc+srcs),1.e-36_r8)    

!c          if (rat.lt.1.and.rat.gt.0.and.tracer(i,k).gt.1.e-13_r8) 
!c     &	write(6,*) 'rat ',i,lat,k,rat,tracer(i,k),srcc,srcs
!c	  if (rat2.lt.1.and.rat2.gt.0.and.tracer(i,k).gt.1.e-13_r8) 
!c     &  write(6,*) 'rat2 ',i,lat,k,rat2,cldv(i,k),tracer(i,k)  
            if (rat.lt.1._r8) then
               srcs = srcs*rat
               srcc = srcc*rat
            endif
            srct = (srcc+srcs)*omsm
!c	if (lat.eq.45.and.i.eq.14) then
!c	write(6,*) 'srct ',srct,srcc,srcs,srcs1,srcs2
!c	end if
!c fraction that is not removed within the cloud
!c (assumed to be interstitial, and subject to convective transport)


            fracis(i,k) =             &
                1._r8 -                 &
                max(0._r8,              &
                    min(1._r8,          & 
                        dt*srct/max(cldfrc(i,k)*tracer(i,k),1.e-36_r8)  &
                       )             &
                   )

!c
!c tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            scavt(i,k) = -srct     
!                + fracev*scavab(i)*gravit/pdel(i,k)
            evsa2(i,k)=fracev*scavab(i)*gravit/pdel(i,k)
!            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

!c now keep track of what we have scavenged from above by stratiform rain
            scavab(i) = scavab(i)*(1-fracev) + srcs*pdel(i,k)/gravit
            precabs(i) = precabs(i)     &
                + (precs(i,k) - evaps(i,k))*pdel(i,k)/gravit
!            scavabc(i) = scavabc(i) + srcc*pdel(i,k)/gravit
            precabc(i) = precabc(i)     &
                + (cmfdqr(i,k))*pdel(i,k)/gravit
!            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)/gravit
         end do
      end do


      return
      end
