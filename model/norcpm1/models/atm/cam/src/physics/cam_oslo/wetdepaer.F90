     subroutine wetdepaer (lchnk,ncol,deltat, pdel, cldfrc, &
       cmfdqr, precs, conds,                                &
       evaps, fice, cwat, tracer,conicw , wincon,           &
       scavt, fracis,winclo,wbeclo,wbeconclo)   

! wet scavenging of particles in CAM-Oslo  Ø Seland 2010


use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid
use physconst,     only: gravit
implicit none



! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: deltat                   ! Time step
   real(r8), intent(in) :: tracer(pcols,pver)   !  
   real(r8), intent(in) :: pdel(pcols,pver)     ! Dp between layers
   real(r8), intent(in) :: cldfrc(pcols,pver)   ! Volume of cloud fraction
   real(r8), intent(in) :: fice(pcols,pver)     ! Fraction of cwat that is ice
   real(r8), intent(in) :: cwat(pcols,pver)     ! Cloud water
   real(r8), intent(in) :: cmfdqr(pcols,pver)   ! rate of production of convective precip
   real(r8), intent(in) :: precs(pcols,pver)    ! rate of production of stratiform precip
   real(r8), intent(in) :: conds(pcols,pver)    ! rate of production of condensate
   real(r8), intent(in) :: evaps(pcols,pver)    ! rate of evaporation of precip

   real(r8), intent(in) :: conicw(pcols,pver)   ! convective cloud water

   real(r8), intent(in) :: winclo
   real(r8), intent(in) :: wbeclo
   real(r8), intent(in) :: wincon
   real(r8), intent(in) :: wbeconclo 

   real(r8), intent(inout) :: fracis(pcols,pver)! fraction of species not scavenged               
   real(r8), intent(out) :: scavt(pcols,pver)    ! scavenging tend 



! local variables

      integer :: i                 ! x index
      integer :: k                 ! z index

      real(r8) :: fracev               ! fraction of precip from above that is evaporating
      real(r8) :: fracp                ! fraction of cloud water converted to precip
      real(r8) :: gafrac               ! fraction of tracer in gas phasea
      real(r8) :: hconst               ! henry's law solubility constant when equation is expressed
                                ! in terms of mixing ratios
      real(r8) :: omsm                 ! 1 - (a small number)
      real(r8) :: pdog                 ! work variable (pdel/grav)
      real(r8) :: precabc(pcols)       ! conv precip from above (work array)
      real(r8) :: precabs(pcols)       ! strat precip from above (work array)
      real(r8) :: precbl               ! precip falling out of level (work array)
      real(r8) :: rat                  ! ratio of amount available to amount removed
      real(r8) :: scavab(pcols)        ! scavenged tracer flux from above (work array)
      real(r8) :: scavabc(pcols)       ! scavenged tracer flux from above (work array)
      real(r8) :: srcc                 ! tend for convective rain
      real(r8) :: srcs                 ! tend for stratiform rain
      real(r8) :: srct                 ! work variable
      real(r8) :: tracab(pcols)        ! column integrated tracer amount
!      real(r8) :: vfall                ! fall speed of precip
      real(r8) :: fins                 ! fraction of rem. rate by strat rain
      real(r8) :: finc                 ! fraction of rem. rate by conv. rain
      real(r8) :: srcs1                ! work variable
      real(r8) :: srcs2                ! work variable
      real(r8) :: cldmabs(pcols)       ! maximum cloud at or above this level
      real(r8) :: cldmabc(pcols)       ! maximum cloud at or above this level
      real(r8) :: odds                 ! limit on removal rate (proportional to prec)
      real(r8) :: rdel
      real(r8) :: rgravit

! ------------------------------------------------------------------------
      omsm = 1._r8-1.e-10_r8          ! used to prevent roundoff errors below zero

! assume 4 m/s fall speed currently (should be improved)
!      vfall = 4._r8


! this section of code is for highly soluble aerosols,
! the assumption is that within the cloud that
! all the tracer is in the cloud water

! for convective clouds, we dont know the cloud water amount. 
!                        and assume that all the cloud water falls out on each time step
!                        therefore all the tracer goes with it.

! for stratiform clouds, the fraction of cloud water converted to precip defines
!                        the amount of tracer which is pulled out.


! keep track of the amount of stratiform precip falling from above
! also the amount of tracer falling from above
      do i = 1,ncol
         precabs(i) = 0._r8
         precabc(i) = 0._r8
         scavab(i) = 0._r8
         scavabc(i) = 0._r8
         tracab(i) = 0._r8
         cldmabs(i) = 0._r8
         cldmabc(i) = 0._r8
      end do
      rgravit=1._r8/gravit
      rdel=1._r8/deltat

      do k = 1,pver
         do i = 1,ncol


!            tc     = t(i,k) - 273.16_r8
!            weight = max(0._r8,min(-tc*0.04_r8,1.0_r8)) ! fraction of condensate that is ice

            pdog = pdel(i,k)*rgravit


! calculate the fraction of strat precip from above 
!                 which evaporates within this layer
            fracev = evaps(i,k)*pdel(i,k)/gravit &
                     /max(1.e-12_r8,precabs(i))
            ! trap to ensure reasonable ratio bounds
            fracev = max(0._r8,min(1._r8,fracev))


! now do the convective scavenging

! set odds proportional to fraction of the grid box that is swept by the 
! precipitation =precabc/rhoh20*(area of sphere projected on plane
!                                /volume of sphere)*deltat
! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
! unless the fraction of the area that is cloud is less than odds, in which
! case use the cloud fraction (assumes precabs is in kg/m2/s)
! is really: precabs*3/4/1000./1e-3*deltat
! here I use .1 from Balkanski

! Picked up from the wetdep subroutine in standard CAM and used instead
! of the values in CCM-OSLO
            ! use a local rate of convective rain production for incloud scav
            !odds=max(min(1._r8, &
            !     cmfdqr(i,k)*pdel(i,k)/grav*0.1_r8*deltat),0._r8)
            !++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            !            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1._r8-weight) &
            ! srcs1 =  cldv(i,k)*odds*tracer(i,k)*(1._r8-weight) &
            !srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1._r8-weight) &
            !         /deltat




! use a local rate of convective rain production for incloud scav
!cos            odds=max(min(1._r8,
!cos     $           cmfdqr(i,k)*pdel(i,k)/gravit*.1*deltat),0._r8)


!c++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
!            srcs1 =  cldc(i,k)*odds*tracer(i,k)*(1._r8-weight)
!            srcs1 =  cldt(i,k)*odds*tracer(i,k)*(1._r8-weight)

!            srcs1 =  1._r8*odds*tracer(i,k)*(1._r8-fice(i,k))   &
!                *rdel 
!Test on what happens cmfdqr is grid-values instead of incloud values
              odds = (cmfdqr(i,k)/max(cldfrc(i,k),1.e-4_r8))*deltat/max(1.e-8_r8,conicw(i,k))
              odds = max(min(1._r8,odds),0._r8)       
              srcs1 =  wincon*cldfrc(i,k)*odds*tracer(i,k)*(1._r8-fice(i,k)) &
                  *rdel 

!     $           cmfdqr(i,k)*pdel(i,k)/gravit*winclo*deltat),0._r8)
!c--mcb

! scavenge below layers of precip
!c++mcb -- change cldc to cldt; change cldt to cldv (9/17/96)
            cldmabc(i) = max(cldfrc(i,k),cldmabc(i))
            odds=max(        &
                  min(1._r8,precabc(i)/max(cldmabc(i),1.e-5_r8)  &
                  *deltat),0._r8)                      
            srcs2 = wbeconclo*cldmabc(i)*odds*tracer(i,k)             &
                *rdel
            srcc = srcs1 + srcs2
            finc = srcs1/(srcc + 1.e-36_r8)

! now do the stratiform scavenging

!c++mcb -- change cldt to cldv (9/17/96)
            cldmabs(i) = max(cldfrc(i,k),cldmabs(i))

! for layers with significant stratiform rain
! fracp is the fraction of cloud water converted to precip
           fracp =            &
                max(0._r8,       &
!cos     $               min(precs(i,k)*deltat
                    min(precs(i,k)*deltat  &
!     $                   /max(cwat(i,k)+conds(i,k)*deltat,1.e-12_r8),
                        /max(cwat(i,k),1.e-12_r8),   &
                        1._r8)   &
                   )

! assume the corresponding amnt of tracer is removed
!c++mcb -- remove cldc; change cldt to cldv 
!            srcs1 = (cldt(i,k)-cldc(i,k))*fracp*tracer(i,k)/deltat
!            srcs1 = cldt(i,k)*fracp*tracer(i,k)/deltat
            srcs1 = winclo*0.3_r8*cldfrc(i,k)*fracp*tracer(i,k)*rdel   &
                *(1._r8-fice(i,k))   ! scavenge only the liquid phase
!c--mcb

! now calc the amount removed by rain from above

! same calc as above but for strat rain
            odds=max(    &
                  min(1._r8,precabs(i)/max(cldmabs(i),1.e-5_r8)  &
                  *deltat),0._r8)

! scavenging due to rain from above

!               srcs2 =((cldt(i,k)-cldc(i,k))*odds)
            srcs2 =wbeclo*(0.3_r8*cldmabs(i)*odds)*tracer(i,k)*rdel


! total strat scavenged 
            srcs = srcs1 + srcs2

! fraction taken by incloud processes
            fins=srcs1/(srcs + 1.e-36_r8)

! make sure we dont take out more than is there
! ratio of amount available to amount removed
!            rat = tracer(i,k)/max(deltat*(srcc+srcs),1.e-100_r8)
            rat = cldfrc(i,k)*tracer(i,k)/     &
        max(deltat*(srcc+srcs),1.e-36_r8)    
!          if (rat.lt.1._r8and.rat.gt.0._r8and.tracer(i,k).gt.1.e-13_r8) 
!     &	write(6,*) 'rat ',i,lat,k,rat,tracer(i,k),srcc,srcs
!	  if (rat2.lt.1._r8and.rat2.gt.0._r8and.tracer(i,k).gt.1.e-13_r8) 
!     &  write(6,*) 'rat2 ',i,lat,k,rat,rat2,cldt(i,k),tracer(i,k)  
            if (rat.lt.1._r8) then
               srcs = srcs*rat
               srcc = srcc*rat
            end if

            srct = (srcc+srcs)*omsm

! fraction that is not removed within the cloud
! (assumed to be interstitial, and subject to convective transport)
! if not grid-val
            fracis(i,k) =         &
                1._r8 -              &
                max(0._r8,           &
                    min(1._r8,       &
                        deltat*srct/max(cldfrc(i,k)*tracer(i,k),1.e-36_r8) &
                       )          &
                   )      



! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
            scavt(i,k) = -srct &
                + fracev*scavab(i)*gravit/pdel(i,k)

!            iscavt(i,k) = -(srcc*finc + srcs*fins)*omsm

! now keep track of what we have scavenged from above by stratiform rain
            scavab(i) = scavab(i)*(1-fracev) + srcs*pdel(i,k)*rgravit
            precabs(i) = precabs(i) &
                + (precs(i,k) - evaps(i,k))*pdel(i,k)*rgravit
!            scavabc(i) = scavabc(i) + srcc*pdel(i,k)*rgravit
            precabc(i) = precabc(i)  &
                + (cmfdqr(i,k))*pdel(i,k)*rgravit
!            tracab(i) = tracab(i) + tracer(i,k)*pdel(i,k)*rgravit
!	if (lat.eq.45._r8and.i.eq.14) then
!            write(6,*) 'wetdepa ',i,k,scavab(i),precabs(i),scavabc(i)
!	    write(6,*) precabc(i),tracab(i),scavt(i,k)
!        end if

         end do
      end do


      return
      end
