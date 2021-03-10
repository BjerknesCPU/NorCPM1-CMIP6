      subroutine uvrotr2g(missing_value,fill_value,plat,plon,
     .                    umask,vmask,u,v)
c
c --- ------------------------------------------------------------------
c --- Rotate zonal/meridional vector components at C-grid scalar points
c --- to model grid components at C-grid velocity points
c --- ------------------------------------------------------------------
c
      use mod_xc
c
      implicit none
c
      real missing_value,fill_value
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: plat,plon,u,v
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: umask,vmask
c
      real latlim,rad
      parameter (latlim=87.,rad=1.74532925199432958e-02)
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: utmp,vtmp
      real ulat,vlat,dlat,dlon,psi
      integer i,j
c
c --- ------------------------------------------------------------------
c --- rotate vector components
c --- ------------------------------------------------------------------
c
c$OMP PARALLEL DO PRIVATE(ulat,dlat,dlon,psi,vlat)
      do j=1,jj
        do i=1,ii
c
c --- --- u-component
          if (umask(i,j).eq.0) then
            utmp(i,j)=missing_value
          else
            ulat=.5*(plat(i,j)+plat(i-1,j))
            if (abs(ulat).gt.latlim) then
              utmp(i,j)=fill_value
            else
              dlat=plat(i,j)-plat(i-1,j)
              dlon=mod(plon(i  ,j)+360.,360.)
     .            -mod(plon(i-1,j)+360.,360.)
              if (abs(dlon+360.).lt.abs(dlon)) dlon=dlon+360.
              if (abs(dlon-360.).lt.abs(dlon)) dlon=dlon-360.
              psi=atan2(dlat,cos(ulat*rad)*dlon)
              utmp(i,j)=.5*((u(i,j)+u(i-1,j))*cos(psi)
     .                     +(v(i,j)+v(i-1,j))*sin(psi))
            endif
          endif
c
c --- --- v-component
          if (vmask(i,j).eq.0) then
            vtmp(i,j)=missing_value
          else
            vlat=.5*(plat(i,j)+plat(i,j-1))
            if (abs(vlat).gt.latlim) then
              vtmp(i,j)=fill_value
            else
              dlat=plat(i,j)-plat(i,j-1)
              dlon=mod(plon(i,j  )+360.,360.)
     .            -mod(plon(i,j-1)+360.,360.)
              if (abs(dlon+360.).lt.abs(dlon)) dlon=dlon+360.
              if (abs(dlon-360.).lt.abs(dlon)) dlon=dlon-360.
              psi=atan2(dlat,cos(vlat*rad)*dlon)
              vtmp(i,j)=.5*((u(i,j)+u(i,j-1))*cos(psi)
     .                     +(v(i,j)+v(i,j-1))*sin(psi))
            endif
          endif
c
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1,jj
        do i=1,ii
          u(i,j)=utmp(i,j)
          v(i,j)=vtmp(i,j)
        enddo
      enddo
c$OMP END PARALLEL DO
c
c --- ------------------------------------------------------------------
c --- extrapolate values to velocity points near pole singularity
c --- ------------------------------------------------------------------
c
      call fill_global(missing_value,fill_value,halo_uv,u)
      call fill_global(missing_value,fill_value,halo_vv,v)
c
      return
      end
