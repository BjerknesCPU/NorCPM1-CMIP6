c
c --- ------------------------------------------------------------------
c --- Parameters related to nesting:
c ---   idm_o,jdm_o     - grid dimensions of course outer grid
c ---   i1o,i2o,j1o,j2o - subdomain index interval in coarse outer grid
c ---   nbz             - number of grid cells in boundary zone
c ---   cnst            - constant defining boundary zone weights
c --- ------------------------------------------------------------------
c
      integer idm_o,jdm_o,i1o,i2o,j1o,j2o,nbz
      real cnst
      parameter (idm_o=196,jdm_o=360,
     .           i1o=19,i2o=159,j1o=34,j2o=254,
     .           nbz=5,cnst=3.)
