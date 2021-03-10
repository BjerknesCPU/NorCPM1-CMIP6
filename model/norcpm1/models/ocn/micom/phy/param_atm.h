c
c --- ------------------------------------------------------------------
c --- Parameters related to atmospheric forcing:
c ---   atm_idm  - zonal dimension of atmospheric grid
c ---   atm_jdm  - meridional dimension of atmospheric grid
c ---   atm_abdm - max. number of runoff discharge basins for each
c ---              atmospheric grid point
c ---   atm_nwgt - number of neighbours uses in the interpolation
c --- ------------------------------------------------------------------
c
      integer atm_idm,atm_jdm,atm_abdm,atm_nwgt
#if   defined(NCEP)
      parameter (atm_idm=192,atm_jdm=94,atm_abdm=10,atm_nwgt=12)
#elif defined(ERA)
      parameter (atm_idm=320,atm_jdm=160,atm_abdm=9,atm_nwgt=12)
#endif
