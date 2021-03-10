      block data atmdat
c
      use mod_xc
c
      implicit none
c
#include "param_atm.h"
#include "common_atm.h"
c
      data
c
c ---  atm_mval     - value of a point that should not have a physical
c ---                 value
c ---  atm_fval     - value of a point that do not have a physical
c ---                 value but should have one
c ---  atm_ice_csmt - constant determining how much the atm. fields are
c ---                 smoothed over ice covered regions
c ---  atm_rnf_csmt - constant determining how much the atm. runoff is
c ---                 smoothed at the coastal discharge point
c ---  atm_crnf     - runoff adjustment factor
c ---  atm_cswa     - short-wave radiation adjustment factor
#if   defined(NCEP)
     . atm_mval,atm_fval /-9999.,-99999./
     .,atm_ice_csmt,atm_rnf_csmt /2.e14,1.e13/
     .,atm_crnf,atm_cswa /0.82073,0.88340/
#elif defined(ERA)
     . atm_mval,atm_fval /-9999.,-99999./
     .,atm_ice_csmt,atm_rnf_csmt /0.0,1.e13/
     .,atm_crnf,atm_cswa /0.7234,0.9721/
#endif
c
      end
