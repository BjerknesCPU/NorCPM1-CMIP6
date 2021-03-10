      subroutine getflux_clim
c
c --- ------------------------------------------------------------------
c --- Get climatological fluxes
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_forc.h"
#include "common_clim.h"
c
c --- Set interpolation parameters
      xx=x
      ll1=l1
      ll2=l2
      ll3=l3
      ll4=l4
      ll5=l5
c
c --- ------------------------------------------------------------------
c --- Compute the air-sea fluxes
c --- ------------------------------------------------------------------
c
      call asflux
c
      return
      end
