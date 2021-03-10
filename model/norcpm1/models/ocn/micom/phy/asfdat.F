      block data asfdat
c
      use mod_xc
c
      implicit none
c
#include "common_asf.h"
c
      data
     . zu       /10./     ! m                  meas. height of wind
     .,zt       /10./     ! m                  meas. height of temperature
     .,zq       /10./     ! m                  meas. height of specific humidity
     .,emiss    /.97/     ! --                 emissivity of water
     .,cpair    /1004.7/  ! J / (K kg)         specific heat of dry air
     .,stefanb  /5.67e-8/ ! W / (m^2 K^4)      stefan-boltzman constant
#ifdef CLIM
     .,tciter   /5/       !                    iterations in the computation
                          !                    of transfer coefficients
#endif
#ifdef SYN
     .,tciter   /1/       !                    iterations in the computation
                          !                    of transfer coefficients
#endif
c
      end
