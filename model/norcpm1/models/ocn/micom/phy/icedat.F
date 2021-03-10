      block data icedat
c
      use mod_xc
c
      implicit none
c
#include "common_ice.h"
c
      data
     . albi_f   /.70/     ! --                 max albedo over ice
     .,albi_m   /.60/     ! --                 max albedo over melting ice
     .,albs_f   /.85/     ! --                 albedo over snow
     .,albs_m   /.75/     ! --                 albedo over melting snow
     .,rhoice   /906./    ! kg / m^3           density of ice
     .,rhosnw   /330./    ! kg / m^3           density of snow
     .,rkice    /2.04/    ! w / (m k)          ice conductivity
     .,fusi     /3.02e8/  ! j / m^3            heat of fusion of ice
     .,fuss     /1.10e8/  ! j / m^3            heat of fusion of snow
     .,fice_max /.995/    ! --                 maximum fractional ice cover
     .,tice_m   /273.05/  ! k                  melting point of ice
     .,tsnw_m   /273.15/  ! k                  melting point of snow
     .,hice_nhmn/.50/     ! m                  min. ice thickness northern hemi.
     .,hice_shmn/.30/     ! m                  min. ice thickness southern hemi.
     .,gamma    /2.e-7/   ! 1 / s              snow aging timescale
     .,sice     /6./      ! per mil            salinity of seaice
     .,sref     /34.65/   ! per mil            global ref. surface salinity
     .,rksnw    /.31/     ! w / (m k)          snow conductivity
     .,cwi      /0.006/   ! --                 ice-ocean heat transfer coeff.
     .,cuc      /1.e3/    ! w / (m^2 k)        const. for heat flux
                          !                    assossiated with
                          !                    under-cooled water,
                          !                    resulting in a temp.
                          !                    adjustment of a 20 m
                          !                    mixed layer towards
                          !                    freezing point with an
                          !                    e-folding timescale of
                          !                    approx. 1 day.
     .,cdiff    /1.e3/    ! m^2 / s            horizontal diffusivity
     .,cdfac    /1.  /    !                    cdiff is multiplied by
                          !                    cdfac in channels of one
                          !                    grid point width
c
      end
