      real function swtfrz(s)
c
c --- computes sea water freezing point temperature [deg. celsius],
c
      implicit none
c
c --- input variables:
c ---   s      - salinity [psu]
c
      real s
c
cmsi#ifdef CCSMCOUPLED
      swtfrz=-1.8
cmsi#else
cmsi      swtfrz=-0.0547*s
cmsi#endif
c
      return
      end
