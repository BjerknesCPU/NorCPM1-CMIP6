!ut subroutine koagsub(kcomp)

subroutine koagsub(kcomp,rhob,rk)

!      Modified by Arild Burud/NoSerC - April 2002
!      - dmpxr introduced, is used in parmix()
!      - a**n.0 changed to a**n for speedup
!      - a**0.5 changed to sqrt(a) 
        
      use shr_kind_mod, only: r8 => shr_kind_r8
      use const
      use physconst,    only: pi
      implicit none

      real(r8) diff1(0:101), diff2(0:101), g12(0:101), &
           g1(0:101), g2(0:101), c12(0:101), c1(0:101), c2(0:101), & 
           mfv1(0:101), mfv2(0:101)
      integer i,kcomp
!inn+
!      real(r8), intent(in) :: rk(14)
!      real(r8), intent(in) :: rhob(14)
!      real(r8), intent(in) :: rk(10)
!      real(r8), intent(in) :: rhob(10)
      real(r8), intent(in) :: rk(0:10)
      real(r8), intent(in) :: rhob(0:10)
!inn-
      real (r8) moc,mbc,ms4,noc,nbc,ns4


!     coagulation coefficient for SO4 (Brownian, Fuchs form)
      do i=0,imax

        c1(i)=4.786e4_r8/sqrt(rhob(kcomp)*rp(i)**3)
        c2(i)=4.786e4_r8/sqrt(rhob(1)*rk(1)**3)
        c12(i)=sqrt(c1(i)**2+c2(i)**2)
        diff1(i)=(11.64_r8/rp(i))*(5.0_r8+0.253_r8/rp(i)+0.024_r8/rp(i)**2     &
         +0.00457_r8/rp(i)**3)/(5.0_r8-0.0633_r8/rp(i)+0.0446_r8/rp(i)**2)
        diff2(i)=(11.64_r8/rk(1))*(5.0_r8+0.253_r8/rk(1)+0.024_r8/rk(1)**2  &
         +0.00457_r8/rk(1)**3)/(5.0_r8-0.0633_r8/rk(1)+0.0446_r8/rk(1)**2)
        mfv1(i)=8.0_r8*diff1(i)/(pi*c1(i))
        mfv2(i)=8.0_r8*diff2(i)/(pi*c2(i))
        g1(i)=((2.0_r8*rp(i)+mfv1(i))**3     &
         -(4.0_r8*rp(i)**2+mfv1(i)**2)**1.5_r8) &
         /(6.0_r8*rp(i)*mfv1(i))-2.0_r8*rp(i)
        g2(i)=((2.0_r8*rk(1)+mfv2(i))**3        &
         -(4.0_r8*rk(1)**2+mfv2(i)**2)**1.5_r8)    &
         /(6.0_r8*rk(1)*mfv2(i))-2.0_r8*rk(1)
        g12(i)=sqrt(g1(i)**2+g2(i)**2)
        Kp12s4(kcomp,i)=4*pi*(rp(i)+rk(1))*(diff1(i)+diff2(i))  &
         /((rp(i)+rk(1))/(rp(i)+rk(1)+g12(i))                  &
          +(4.0_r8/c12(i))*(diff1(i)+diff2(i))/(rk(1)+rp(i)))
      enddo


!     coagulation coefficient for OC (Brownian, Fuchs form)
      do i=0,imax
        c1(i)=4.786e4_r8/sqrt(rhob(kcomp)*rp(i)**3)
        c2(i)=4.786e4_r8/sqrt(rhob(3)*rk(3)**3)
        c12(i)=sqrt(c1(i)**2+c2(i)**2)
        diff1(i)=(11.64_r8/rp(i))*(5.0_r8+0.253_r8/rp(i)+0.024_r8/rp(i)**2     &
         +0.00457_r8/rp(i)**3)/(5.0_r8-0.0633_r8/rp(i)+0.0446_r8/rp(i)**2)
        diff2(i)=(11.64_r8/rk(3))*(5.0_r8+0.253_r8/rk(3)+0.024_r8/rk(3)**2  &
         +0.00457_r8/rk(3)**3)/(5.0_r8-0.0633_r8/rk(3)+0.0446_r8/rk(3)**2)
        mfv1(i)=8.0_r8*diff1(i)/(pi*c1(i))
        mfv2(i)=8.0_r8*diff2(i)/(pi*c2(i))
        g1(i)=((2.0_r8*rp(i)+mfv1(i))**3     &
         -(4.0_r8*rp(i)**2+mfv1(i)**2)**1.5_r8) &
         /(6.0_r8*rp(i)*mfv1(i))-2.0_r8*rp(i)
        g2(i)=((2.0_r8*rk(3)+mfv2(i))**3        &
         -(4.0_r8*rk(3)**2+mfv2(i)**2)**1.5_r8)    &
         /(6.0_r8*rk(3)*mfv2(i))-2.0_r8*rk(3)
        g12(i)=sqrt(g1(i)**2+g2(i)**2)
        Kp12oc(kcomp,i)=4*pi*(rp(i)+rk(3))*(diff1(i)+diff2(i))  &
         /((rp(i)+rk(3))/(rp(i)+rk(3)+g12(i))                  &
          +(4.0_r8/c12(i))*(diff1(i)+diff2(i))/(rk(3)+rp(i)))
      enddo


!     coagulation coefficient for BC_AX (Brownian, Fuchs form)
      do i=0,imax
        c1(i)=4.786e4_r8/sqrt(rhob(kcomp)*rp(i)**3)
!        c2(i)=4.786e4_r8/sqrt(rhob(5)*rk(5)**3)
        c2(i)=4.786e4_r8/sqrt(rhob(0)*rk(0)**3)
        c12(i)=sqrt(c1(i)**2+c2(i)**2)
        diff1(i)=(11.64_r8/rp(i))*(5.0_r8+0.253_r8/rp(i)+0.024_r8/rp(i)**2     &
         +0.00457_r8/rp(i)**3)/(5.0_r8-0.0633_r8/rp(i)+0.0446_r8/rp(i)**2)
!        diff2(i)=(11.64_r8/rk(5))*(5.0_r8+0.253_r8/rk(5)+0.024_r8/rk(5)**2  &
!         +0.00457_r8/rk(5)**3)/(5.0_r8-0.0633_r8/rk(5)+0.0446_r8/rk(5)**2)
        diff2(i)=(11.64_r8/rk(0))*(5.0_r8+0.253_r8/rk(0)+0.024_r8/rk(0)**2  &
         +0.00457_r8/rk(0)**3)/(5.0_r8-0.0633_r8/rk(0)+0.0446_r8/rk(0)**2)
        mfv1(i)=8.0_r8*diff1(i)/(pi*c1(i))
        mfv2(i)=8.0_r8*diff2(i)/(pi*c2(i))
        g1(i)=((2.0_r8*rp(i)+mfv1(i))**3     &
         -(4.0_r8*rp(i)**2+mfv1(i)**2)**1.5_r8) &
         /(6.0_r8*rp(i)*mfv1(i))-2.0_r8*rp(i)
!        g2(i)=((2.0_r8*rk(5)+mfv2(i))**3        &
!         -(4.0_r8*rk(5)**2+mfv2(i)**2)**1.5_r8)    &
!         /(6.0_r8*rk(5)*mfv2(i))-2.0_r8*rk(5)
        g2(i)=((2.0_r8*rk(0)+mfv2(i))**3        &
         -(4.0_r8*rk(0)**2+mfv2(i)**2)**1.5_r8)    &
         /(6.0_r8*rk(0)*mfv2(i))-2.0_r8*rk(0)
        g12(i)=sqrt(g1(i)**2+g2(i)**2)
!        Kp12ax(kcomp,i)=4*pi*(rp(i)+rk(5))*(diff1(i)+diff2(i))  &
!         /((rp(i)+rk(5))/(rp(i)+rk(5)+g12(i))                  &
!          +(4.0_r8/c12(i))*(diff1(i)+diff2(i))/(rk(5)+rp(i)))
        Kp12ax(kcomp,i)=4*pi*(rp(i)+rk(0))*(diff1(i)+diff2(i))  &
         /((rp(i)+rk(0))/(rp(i)+rk(0)+g12(i))                  &
          +(4.0_r8/c12(i))*(diff1(i)+diff2(i))/(rk(0)+rp(i)))
      enddo


!     coagulation coefficient for BC (Brownian, Fuchs form),
!     and diffusion coefficient for sulphate as H2SO4
      do i=0,imax

        c1(i)=4.786e4_r8/sqrt(rhob(kcomp)*rp(i)**3)
        c2(i)=4.786e4_r8/sqrt(rhob(2)*rk(2)**3)
        c12(i)=sqrt(c1(i)**2+c2(i)**2)
        diff1(i)=(11.64_r8/rp(i))*(5.0_r8+0.253_r8/rp(i)+0.024_r8/rp(i)**2     &
         +0.00457_r8/rp(i)**3)/(5.0_r8-0.0633_r8/rp(i)+0.0446_r8/rp(i)**2)
        diff2(i)=(11.64_r8/rk(2))*(5.0_r8+0.253_r8/rk(2)+0.024_r8/rk(2)**2  &
         +0.00457_r8/rk(2)**3)/(5.0_r8-0.0633_r8/rk(2)+0.0446_r8/rk(2)**2)
        mfv1(i)=8.0_r8*diff1(i)/(pi*c1(i))
        mfv2(i)=8.0_r8*diff2(i)/(pi*c2(i))
        g1(i)=((2.0_r8*rp(i)+mfv1(i))**3       &
         -(4.0_r8*rp(i)**2+mfv1(i)**2)**1.5_r8)   &
         /(6.0_r8*rp(i)*mfv1(i))-2.0_r8*rp(i)
        g2(i)=((2.0_r8*rk(2)+mfv2(i))**3          &
         -(4.0_r8*rk(2)**2+mfv2(i)**2)**1.5_r8)      &
         /(6.0_r8*rk(2)*mfv2(i))-2.0_r8*rk(2)
        g12(i)=sqrt(g1(i)**2+g2(i)**2)
        Kp12(kcomp,i)=4*pi*(rp(i)+rk(2))*(diff1(i)+diff2(i))  &
         /((rp(i)+rk(2))/(rp(i)+rk(2)+g12(i))                &
          +(4.0_r8/c12(i))*(diff1(i)+diff2(i))/(rk(2)+rp(i)))
        Dmp(i)=diff/(rp(i)/(rp(i)+mfv)+4.0_r8*diff/(th*rp(i)))

        Dmp_dst(i)=diff/(rp(i)/(rp(i)+mfv)+4.0_r8*diff/(0.3_r8*th*rp(i)))
        Dmp_bcn(i)=diff/(rp(i)/(rp(i)+mfv)+4.0_r8*diff/(0.3_r8*th*rp(i)))
        Dmp_omn(i)=diff/(rp(i)/(rp(i)+mfv)+4.0_r8*diff/(0.7_r8*th*rp(i)))
        Dmp_ni(i) =diff/(rp(i)/(rp(i)+mfv)+4.0_r8*diff/(0.5_r8*th*rp(i)))

!        if(kcomp==1.or.kcomp>=8) then   ! SO4 or Sea-salt
        if(kcomp==1.or.kcomp==5.or.kcomp>=8) then  ! SO4 or Sea-salt
          Dmpa(kcomp,i)=Dmp(i)
        elseif(kcomp==2) then           ! BC 
          Dmpa(kcomp,i)=Dmp_bcn(i)
        elseif(kcomp==3) then           ! OC
          Dmpa(kcomp,i)=Dmp_omn(i)
        elseif(kcomp==4) then           ! OC&BC
          Dmpa(kcomp,i)=Dmp_ni(i)
        else                            ! Dust
          Dmpa(kcomp,i)=Dmp_dst(i)
        endif

!       if(kcomp.eq.3) &
!       write(6,*) i,rp(i),1.e-12_r8*rp(i)*Dmp(i),1.e-12_r8*rp(i)*Dmp_ni(i)
      enddo

      return
end subroutine koagsub
