      subroutine intpol_nest(nxs,nys,nxt,nyt,i1s,i2s,j1s,j2s,
     .                       missing_value,fill_value,
     .                       smask,tmask,sdata,tdata)
c
c --- Interpolate variables from coarse source grid to a fine target
c --- grid, assuming that each source grid cell is divided into exactly
c --- nfac x nfac grid cells of the target grid with nfac being an
c --- integer.
c
      implicit none
c
c --- Input arguments
c ---   nxs,nys         - grid dimensions of the source grid
c ---   nxt,nyt         - grid dimensions of target grid
c ---   i1s,i2s,j1s,j2s - subdomain index interval relative to source
c ---                     grid
c ---   missing_value   - value of missing data
c ---   fill_value      - value of data to be filled
c ---   smask           - source data mask
c ---   tmask           - target data mask
c ---   sdata           - source data
c ---   tdata           - target data
c
      integer nxs,nys,nxt,nyt,i1s,i2s,j1s,j2s,
     .        smask(nxs,nys),tmask(nxt,nyt)
      real missing_value,fill_value,sdata(nxs,nys),tdata(nxt,nyt)
c
      integer nfac,iit,jjt,iis,jjs,it,jt,is,js
      real w(3,3,(nxt-3)/(i2s-i1s+1),(nxt-3)/(i2s-i1s+1)),wi(3),wj(3),
     .     q,d_sum,w_sum
c
      nfac=(nxt-3)/(i2s-i1s+1)
c
c --- compute weights
      q=.5*(nfac+1)
      do jjt=1,nfac
        wj(1)=max(0.,q-jjt)
        wj(3)=max(0.,jjt-q)
        wj(2)=nfac-wj(1)-wj(3)
        do iit=1,nfac
          wi(1)=max(0.,q-iit)
          wi(3)=max(0.,iit-q)
          wi(2)=nfac-wi(1)-wi(3)
          do jjs=1,3
            do iis=1,3
              w(iis,jjs,iit,jjt)=wi(iis)*wj(jjs)
            enddo
          enddo
        enddo
      enddo
c
c --- interpolate
c$OMP PARALLEL DO PRIVATE(iit,jjt,d_sum,w_sum,is,js)
      do jt=1,nyt
        jjt=modulo(jt-2,nfac)+1
        do it=1,nxt
          if (tmask(it,jt).eq.0) then
            tdata(it,jt)=missing_value
          else
            iit=modulo(it-2,nfac)+1
            d_sum=0.
            w_sum=0.
            do jjs=1,3
              js=mod((jt+nfac-2)/nfac+j1s+jjs-3+nys-1,nys)+1
              do iis=1,3
                is=(it+nfac-2)/nfac+i1s+iis-3
                if (smask(is,js).eq.1.and.
     .              sdata(is,js).ne.missing_value) then
                  d_sum=d_sum+sdata(is,js)*w(iis,jjs,iit,jjt)
                  w_sum=w_sum+w(iis,jjs,iit,jjt)
                endif
              enddo
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
c --- fill missing values
c
      call fill_global(nxt,nyt,missing_value,fill_value,tdata)
c
      end
