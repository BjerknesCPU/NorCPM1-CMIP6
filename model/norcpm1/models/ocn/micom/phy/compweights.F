      subroutine compweights(nxs,nys,nxt,nyt,n_weight,
     .                       slat,slon,tlat,tlon,
     .                       weight,i_weight,j_weight)
c
c --- ------------------------------------------------------------------
c --- Compute weights for gaussian interpolation. The weights are
c --- adapted to the coarsest of the source and target grid at each
c --- target grid location. The source grid is assumed to be cyclic in
c --- the i-direction, while the target grid is assumed to be cyclic in
c --- the j-direction.
c
c --- Created 08.11.02 by mats.bentsen@nersc.no
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer nxs,nys,nxt,nyt,n_weight
      integer*2, dimension(n_weight,nxt,nyt) :: i_weight,j_weight
      real, dimension(nxs,nys) :: slat,slon
      real, dimension(nxt,nyt) :: tlat,tlon
      real*4, dimension(n_weight,nxt,nyt) :: weight
c
      real, dimension(0:n_weight) :: r_weight
      real pi,min_d,d,r,l2
      integer nw_2,is,js,it,jt,iso,jso,m,n,i,j,ii,ism1,isp1,jsm1,jsp1,
     .        itm1,itp1,jtm1,jtp1
c
      real spherdist
      external spherdist
c
      pi=4.*atan(1.)
      r_weight(0)=-1.
      nw_2=n_weight/2+1
c
      is=1
      js=1
c
      do jt=1,nyt
        do it=1,nxt
c
          min_d=spherdist(1.,slon(is,js),slat(is,js),
     .                       tlon(it,jt),tlat(it,jt))
c
 100      iso=is
          jso=js
c
          i=mod(iso-2+nxs,nxs)+1
          d=spherdist(1.,slon(i,jso),slat(i,jso),
     .                   tlon(it,jt),tlat(it,jt))
          if (d.lt.min_d) then
            is=i
            js=jso
            min_d=d
          endif
          i=mod(iso,nxs)+1
          d=spherdist(1.,slon(i,jso),slat(i,jso),
     .                   tlon(it,jt),tlat(it,jt))
          if (d.lt.min_d) then
            is=i
            js=jso
            min_d=d
          endif
          j=max(jso-1,1)
          d=spherdist(1.,slon(iso,j),slat(iso,j),
     .                   tlon(it,jt),tlat(it,jt))
          if (d.lt.min_d) then
            is=iso
            js=j
            min_d=d
          endif
          j=min(jso+1,nys)
          d=spherdist(1.,slon(iso,j),slat(iso,j),
     .                   tlon(it,jt),tlat(it,jt))
          if (d.lt.min_d) then
            is=iso
            js=j
            min_d=d
          endif
c
          if (is.ne.iso.or.js.ne.jso) goto 100
c
          do m=1,n_weight
            r_weight(m)=999999.
          enddo
c
          do j=min(nys-nw_2*2,max(1,       js-nw_2)),
     .         min(nys       ,max(nw_2*2+1,js+nw_2))
            do ii=is-nw_2,is+nw_2
              i=mod(ii-1+nxs,nxs)+1
              r=spherdist(1.,slon(i ,j ),slat(i ,j ),
     .                       tlon(it,jt),tlat(it,jt))
              m=n_weight+1
 10           m=m-1
              if (r.lt.r_weight(m)) goto 10
              m=m+1
              if (m.le.n_weight) then
                do n=n_weight-1,m,-1
                  r_weight(n+1)=r_weight(n)
                  i_weight(n+1,it,jt)=i_weight(n,it,jt)
                  j_weight(n+1,it,jt)=j_weight(n,it,jt)
                enddo
                r_weight(m)=r
                i_weight(m,it,jt)=i
                j_weight(m,it,jt)=j
              endif
            enddo
          enddo
c
          itm1=min(nxt-2,max(1,it-1))
          itp1=min(nxt  ,max(3,it+1))
          jtm1=mod(jt-2+nyt,nyt)+1
          jtp1=mod(jt,nyt)+1
          ism1=mod(is-2+nxs,nxs)+1
          isp1=mod(is,nxs)+1
          jsm1=min(nys-2,max(1,js-1))
          jsp1=min(nys  ,max(3,js+1))
c
          l2=.25*max( spherdist(1.,tlon(itm1,jt),tlat(itm1,jt),
     .                             tlon(itp1,jt),tlat(itp1,jt))
     .               *spherdist(1.,tlon(it,jtm1),tlat(it,jtm1),
     .                             tlon(it,jtp1),tlat(it,jtp1))
     .              , spherdist(1.,slon(ism1,js ),slat(ism1,js ),
     .                             slon(isp1,js ),slat(isp1,js ))
     .               *spherdist(1.,slon(is ,jsm1),slat(is ,jsm1),
     .                             slon(is ,jsp1),slat(is ,jsp1)))/pi
c
          do m=1,n_weight
            weight(m,it,jt)=max(exp(-.5*r_weight(m)*r_weight(m)/l2),
     .                          1.e-9)
          enddo
c
        enddo
      enddo
c
#ifdef ARCTIC
      do it=1,nxt
        do m=1,n_weight
          weight(m,it,nyt)=weight(m,nxt-it+1,nyt-1)
          i_weight(m,it,nyt)=i_weight(m,nxt-it+1,nyt-1)
          j_weight(m,it,nyt)=j_weight(m,nxt-it+1,nyt-1)
        enddo
      enddo
#endif
c
      return
      end
