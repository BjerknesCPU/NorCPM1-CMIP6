      subroutine albw_eval(dangle,lat,cc,albw_d,albw)
c
c --- ------------------------------------------------------------------
c --- compute 24 hrs mean albedo at the marine surface layer
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: lat,cc,albw
      real dangle,albw_d
c
      real pi2,deg,eepsil,fraci,absh2o,s0,decli,sundv,sin2,cos2,stot,
     .     sads,hangle,cosz,srad,sdir,sdif,altdeg,cfac,ssurf,albdir
      integer ifrac,i,j,l,npart
c
c --- set various quantities
c
      pi2=8.*atan(1.)          !        2 times pi
      deg=360./pi2             !        convert from radians to degrees
      eepsil=1.e-9             !        small number
c
      ifrac=24                 !        split each 12 hrs day into ifrac parts
      fraci=1./ifrac           !        1 over ifrac
c
      absh2o=0.09              ! ---    absorption of water and ozone
      s0=1365.                 ! w/m^2  solar constant
c
c --- compute astronomic quantities
c
      decli=.006918+.070257*sin(dangle)   -.399912*cos(dangle)
     $             +.000907*sin(2.*dangle)-.006758*cos(2.*dangle)
     $             +.001480*sin(3.*dangle)-.002697*cos(3.*dangle)
c
      sundv=1.00011+.001280*sin(dangle)   +.034221*cos(dangle)
     $             +.000077*sin(2.*dangle)+.000719*cos(2.*dangle)
c
c$OMP PARALLEL DO PRIVATE(sin2,cos2,stot,sads,hangle,cosz,srad,sdir,
c$OMP+ sdif,altdeg,cfac,ssurf,albdir)
      do j=1,jj
        do l=1,isp(j)
        do i=max(1,ifp(j,l)),min(ii,ilp(j,l))
c
c --- --- compute astronomic quantities
c
          sin2=sin(lat(i,j)/deg)*sin(decli)
          cos2=cos(lat(i,j)/deg)*cos(decli)
c
c --- --- split each day into ifrac parts, and compute the solar
c --- --- radiance for each part. by assuming symmetry of the irradiance
c --- --- about noon, it is sufficient to compute the irradiance for the
c --- --- first 12 hrs of the (24 hrs) day (mean for the first 12 hrs
c --- --- equals then the mean for the last 12 hrs)
c
          stot=0.
          sads=0.
c
          do npart=1,ifrac
c
c --- ----- hour angle, in radians
            hangle=pi2*(npart-.5)*fraci*.5
c
c --- ----- cosine of the zenith angle
            cosz=min(1.,max(0.,sin2+cos2*cos(hangle)))
c
c --- ----- extraterrestrial radiation
            srad =s0*sundv*cosz
c
c --- ----- direct radiation component
            sdir=srad*0.7**min(100.,1./(cosz+eepsil))
c
c --- ----- diffusive radiation component
            sdif=((1.-absh2o)*srad-sdir)*.5
c
c --- ----- solar noon altitude in degrees
            altdeg=max(0.,asin(sin2+cos2))*deg
c
c --- ----- cloudiness correction
            cfac=(1.-0.62*cc(i,j)+0.0019*altdeg)
c
            ssurf=(sdir+sdif)*cfac+eepsil
            stot=stot+ssurf
c
c --- ----- albedo for direct light
            albdir=0.03*exp(0.742*acos(cosz)**2.866)
c
c --- ----- radiation weighted sum of direct albedo
            sads=sads+albdir*ssurf
c
          enddo
c
c --- --- daily mean albedo over water
          albw(i,j)=(1.-cc(i,j))*sads/stot+cc(i,j)*albw_d
c
        enddo
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
