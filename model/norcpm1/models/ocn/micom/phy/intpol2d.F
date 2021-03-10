      subroutine intpol2d(nxs,nys,nxt,nyt,n_weight,
     .                    weight,i_weight,j_weight,
     .                    missing_value,fill_value,
     .                    smask,tmask,sdata,tdata)
c
c --- ------------------------------------------------------------------
c --- Gaussian interpolation of 2D fields using precomputed weights.
c
c --- Created 08.11.02 by mats.bentsen@nersc.no
c --- ------------------------------------------------------------------
c
      implicit none
c
      integer nxs,nys,nxt,nyt,n_weight
      real*4, dimension(n_weight,nxt,nyt) :: weight
      integer*2, dimension(n_weight,nxt,nyt) :: i_weight,j_weight
      real missing_value,fill_value
      integer, dimension(nxs,nys) :: smask
      integer, dimension(nxt,nyt) :: tmask
      real, dimension(nxs,nys) :: sdata
      real, dimension(nxt,nyt) :: tdata
c
      real w_sum,d_sum,w
      integer it,jt,m,n,is,js
c
c --- interpolate
c$OMP PARALLEL DO PRIVATE(w_sum,d_sum,n,is,js,w)
      do jt=1,nyt
        do it=1,nxt
          if (tmask(it,jt).eq.0) then
            tdata(it,jt)=missing_value
          else
            w_sum=0.
            d_sum=0.
            do n=1,n_weight
              is=i_weight(n,it,jt)
              js=j_weight(n,it,jt)
              if (smask(is,js).eq.1.and.
     .            sdata(is,js).ne.missing_value) then
                w=weight(n,it,jt)
                w_sum=w_sum+w
                d_sum=d_sum+sdata(is,js)*w
              endif
            enddo
            if (w_sum.eq.0.) then
              tdata(it,jt)=fill_value
            else
              tdata(it,jt)=d_sum/w_sum
            endif
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
