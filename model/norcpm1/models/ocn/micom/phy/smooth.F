      subroutine smooth(mbdy,w,field,mask,
     .                  scuy,scvx,scuxi,scvyi,scp2i)
c
c --- smooth the field where the mask=1 in a scale independent and
c --- conservative manner
c
      use mod_xc
c
      implicit none
c
      integer mbdy
      real w
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) ::
     .  field,scuy,scvx,scuxi,scvyi,scp2i
      integer, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: mask
c
      real, dimension(1-nbdy:idm+nbdy,1-nbdy:jdm+nbdy) :: uflx,vflx
      integer i,j
c
c$OMP PARALLEL DO
      do j=1-mbdy+1,jj+mbdy
        do i=1-mbdy+1,ii+mbdy
          if (mask(i-1,j).eq.1.and.mask(i,j).eq.1) then
            uflx(i,j)=scuy(i,j)*scuxi(i,j)*(field(i-1,j)-field(i,j))
          else
            uflx(i,j)=0.
          endif
          if (mask(i,j-1).eq.1.and.mask(i,j).eq.1) then
            vflx(i,j)=scvx(i,j)*scvyi(i,j)*(field(i,j-1)-field(i,j))
          else
            vflx(i,j)=0.
          endif
        enddo
      enddo
c$OMP END PARALLEL DO
c
c$OMP PARALLEL DO
      do j=1-mbdy+1,jj+mbdy-1
        do i=1-mbdy+1,ii+mbdy-1
          if (mask(i,j).eq.1)
     .      field(i,j)=field(i,j)-w*scp2i(i,j)*(uflx(i+1,j)-uflx(i,j)
     .                                         +vflx(i,j+1)-vflx(i,j))
        enddo
      enddo
c$OMP END PARALLEL DO
c
      return
      end
