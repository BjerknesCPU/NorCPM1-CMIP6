      subroutine getflux_syn
c
c --- ------------------------------------------------------------------
c --- Get synoptic fluxes
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
#include "param_atm.h"
#include "common_syn.h"
c
c --- Set interpolation parameter
      xx=real(mod(nstep-1,nstep_in_day)+1)/nstep_in_day
c
c --- ------------------------------------------------------------------
c --- The first time step of a new day, read new forcing fields
c --- ------------------------------------------------------------------
c
      if (mod(nstep,nstep_in_day).eq.1) call rdatm_syn
c
c --- ------------------------------------------------------------------
c --- Compute the air-sea fluxes
c --- ------------------------------------------------------------------
c
      call asflux
c
      return
      end
