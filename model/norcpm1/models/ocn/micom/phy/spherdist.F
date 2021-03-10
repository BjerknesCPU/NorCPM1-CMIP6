      real function spherdist(r,lon1,lat1,lon2,lat2)
c
c --- Computes the distance between geo. pos. lon1,lat1 and lon2,lat2 on
c --- a sphere with radius r
c
      real r,lon1,lat1,lon2,lat2
c
      real rad
      parameter(rad=1.74532925199432958e-02)
c
      real lambda,phi,x1,y1,z1,x2,y2,z2
c
      phi=lon1*rad
      lambda=lat1*rad
      x1=cos(lambda)*cos(phi) 
      y1=cos(lambda)*sin(phi) 
      z1=sin(lambda)
c
      phi=lon2*rad
      lambda=lat2*rad
      x2=cos(lambda)*cos(phi) 
      y2=cos(lambda)*sin(phi) 
      z2=sin(lambda)
c
      spherdist=r*acos(min(1.,max(-1.,x1*x2+y1*y2+z1*z2)))
c
      return
      end
