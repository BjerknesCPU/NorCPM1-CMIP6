subroutine modalapp(lchnk,ncol,Nnatk,Ca,f_c,f_bc,f_aq,Cam,fcm,fbcm,faqm)

!     Calculation of the apportionment of internally mixed SO4, BC and OC 
!     mass between the various background mineral and sea-salt modes. Separated 
!     from pmxsub into a independent subroutine by Alf Kirkevåg on September 
!     12'th, 2005.
!     Erroneous ic (irc=i instead of irc=ir) corrected September 14'th 2005 
!     Modified for new aerosol schemes by Alf Kirkevaag in January 2006: Now
!     also Aitken-modes are subject to condensation of H2SO4, and both n and
!     Aitken modes may coagulate onto the mineral/sea-salt background aerosol.

   use ppgrid
   use shr_kind_mod, only: r8 => shr_kind_r8
   use opttab
   use const

   implicit none


!
! Input arguments
!
   integer, intent(in)  :: lchnk                     ! chunk identifier
   integer, intent(in)  :: ncol                      ! number of atmospheric columns
   real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes)! aerosol mode number concentration  
   real(r8), intent(in) :: Ca(pcols,pver)            ! internally mixed mass, tot=SO4+OC+BC   
   real(r8), intent(in) :: f_c(pcols,pver)           ! mass fraction (OC+BC)/tot
   real(r8), intent(in) :: f_bc(pcols,pver)          ! mass fraction BC/(OC+BC)
   real(r8), intent(in) :: f_aq(pcols,pver)          ! mass fraction SO4(aq)/SO4
!
! Output arguments
!
   real(r8), intent(out) :: Cam(pcols,pver,nbmodes)  ! modal internal mass, tot=SO4+BC+OC 
   real(r8), intent(out) :: fcm(pcols,pver,nbmodes)  ! modal mass fraction (OC+BC)/tot (not valid for kcomp=4)
   real(r8), intent(out) :: fbcm(pcols,pver,nbmodes) ! modal mass fraction BC/(OC+BC)
   real(r8), intent(out) :: faqm(pcols,pver,nbmodes) ! modal mass fraction SO4(aq)/SO4
!
! Local variables
!
   real(r8) NrD(pcols,pver,nbmodes), NK12(pcols,pver,nbmodes), NK12oc(pcols,pver,nbmodes), & 
            Nag(pcols,pver,nbmodes), NrDtot(pcols,pver), NK12tot(pcols,pver),  &
            NK12octot(pcols,pver), Nagtot(pcols,pver)
   real(r8) fcondk(nbmodes), fcoagk(nbmodes), fcoagock(nbmodes), faqk(nbmodes) 
   real(r8) cabck(nbmodes), caock(nbmodes), caqsk(nbmodes), ccondsk(nbmodes)
   real(r8) es000(nbmodes,imax),es001(nbmodes,imax),es002(nbmodes,imax),es003(nbmodes,imax)

   integer i, k, kcomp, ir, icol


      do i=1,10
        do ir=1,imax
	   es003(i,ir)=nk(i,ir)*dlogr
           es001(i,ir)=es003(i,ir)*Dmpa(i,ir)*rp(ir)
!	   es002(i,ir)=es003(i,ir)*Kp12(i,ir)     ! Kp12    not used. Instead we use a common:
	   es000(i,ir)=es003(i,ir)*Kp12oc(i,ir)   ! Kp12oc  with r based on 40nm (BC/OC/SO4)
	enddo
        do k=1,pver
         do icol=1,ncol               !                              Ait.   a and c
!g            NrD(icol,k,i)=0.0_r8    ! H2SO4 condensation: on modes 1-4 and 5-10
            NrD(icol,k,i)=0.0_r8      ! H2SO4 cond./coag.:  on modes 1-4  /  5-10
            NK12(icol,k,i)=0.0_r8     ! BC coagulation:     on modes         5-10
            NK12oc(icol,k,i)=0.0_r8   ! OC coagulation:     on modes         5-10
            Nag(icol,k,i)=0.0_r8      ! (NH4)2SO4 aq. chem: on modes   4 and 5-10
         enddo
        enddo 
      enddo ! i

! Mode 0 is a purely externally mixed BC(ax) mode, and is therefore skipped here.

      do i=1,3
! Modes 1-3: only cond. SO4 is to be internally mixed with these modes
        do k=1,pver
          do icol=1,ncol
!g            if(Nnatk(icol,k,i)>0.0_r8) then      
             do ir=1,imax
               NrD(icol,k,i)=NrD(icol,k,i)+Nnatk(icol,k,i)*es001(i,ir)
             enddo
!g            endif   ! Nnatk>0
          end do
        end do
      enddo ! i

      do i=4,4
! Mode 4: cond. and aq. SO4 is to be internally mixed with this mode
        do k=1,pver
          do icol=1,ncol
!g            if(Nnatk(icol,k,i)>0.0_r8) then      
             do ir=1,imax
               NrD(icol,k,i)=NrD(icol,k,i)+Nnatk(icol,k,i)*es001(i,ir)
             enddo
             Nag(icol,k,i)=Nnatk(icol,k,i)*nk(i,irc)*l10rp(ir)
             do ir=irc+1,imax
               Nag(icol,k,i)=Nag(icol,k,i)+Nnatk(icol,k,i)*es003(i,ir)
             enddo
!g            endif   ! Nnatk>0
          end do
        end do
      enddo ! i

      do i=5,10
! Modes 5-10: assuming that SO4 from condensation is lumped together with SO4 
! from coagulation, opposite of the case for modes 1-4.
        do k=1,pver
          do icol=1,ncol
!g            if(Nnatk(icol,k,i)>0.0_r8) then      
             do ir=1,imax
!g               NrD(icol,k,i)=NrD(icol,k,i)+Nnatk(icol,k,i)*es001(i,ir)    ! --> coag istedet!!!
!f               NK12s(icol,k,i)=NK12s(icol,k,i)+Nnatk(icol,k,i)*es000(i,ir)
               NrD(icol,k,i)=NrD(icol,k,i)+Nnatk(icol,k,i)*es000(i,ir)      ! lump coag. and cond.
               NK12(icol,k,i)=NK12(icol,k,i)+Nnatk(icol,k,i)*es000(i,ir)
               NK12oc(icol,k,i)=NK12oc(icol,k,i)+Nnatk(icol,k,i)*es000(i,ir)
             enddo
             Nag(icol,k,i)=Nnatk(icol,k,i)*nk(i,irc)*l10rp(ir)
             do ir=irc+1,imax
               Nag(icol,k,i)=Nag(icol,k,i)+Nnatk(icol,k,i)*es003(i,ir)
             enddo
!g            endif   ! Nnatk>0
          end do
        end do
      enddo ! i

      do k=1,pver
        do icol=1,ncol
          NrDtot(icol,k)=0.0_r8
          NK12tot(icol,k)=0.0_r8
          NK12octot(icol,k)=0.0_r8
          Nagtot(icol,k)=0.0_r8
        end do
      enddo  

! Then the contribution (0 or not) from each mode is added to the total process
! specific masses:
      do i=1,10
        do k=1,pver
         do icol=1,ncol
!g           if(Nnatk(icol,k,i)>0.0_r8) then      
            NrDtot(icol,k)=NrDtot(icol,k)+NrD(icol,k,i)
            NK12tot(icol,k)=NK12tot(icol,k)+NK12(icol,k,i)
            NK12octot(icol,k)=NK12octot(icol,k)+NK12oc(icol,k,i)
            Nagtot(icol,k)=Nagtot(icol,k)+Nag(icol,k,i)
!g           endif
           Cam(icol,k,i)=0.0_r8
           fcm(icol,k,i)=0.0_r8
           fbcm(icol,k,i)=0.0_r8
           faqm(icol,k,i)=0.0_r8
         end do
        enddo  
      enddo  

! And finally the contribution from each mode relative to the totals are calculated,
! assuming that the apprtionment of mass for the first iteration (in time) is representative
! for the whole apportionment process (which is ok for small and moderate masses added):
      do i=1,10
        do k=1,pver
         do icol=1,ncol
!g           if(Nnatk(icol,k,i)>0.0_r8) then      

            fcondk(i)=NrD(icol,k,i)/(NrDtot(icol,k)+eps)
            fcoagk(i)=NK12(icol,k,i)/(NK12tot(icol,k)+eps)
            fcoagock(i)=NK12oc(icol,k,i)/(NK12octot(icol,k)+eps)
            faqk(i)=Nag(icol,k,i)/(Nagtot(icol,k)+eps)

            cabck(i)=fcoagk(i)*f_c(icol,k)*f_bc(icol,k)*Ca(icol,k)
            caock(i)=fcoagock(i)*f_c(icol,k)*(1.0_r8-f_bc(icol,k))*Ca(icol,k)
            caqsk(i)=faqk(i)*f_aq(icol,k)*(1.0_r8-f_c(icol,k))*Ca(icol,k)
            ccondsk(i)=fcondk(i)*(1.0_r8-f_aq(icol,k))*(1.0_r8-f_c(icol,k))*Ca(icol,k)

            Cam(icol,k,i)=cabck(i)+caock(i)+caqsk(i)+ccondsk(i)
            fcm(icol,k,i)=(cabck(i)+caock(i))/(Cam(icol,k,i)+eps)
            fbcm(icol,k,i)=cabck(i)/(cabck(i)+caock(i)+eps)
            faqm(icol,k,i)=caqsk(i)/(caqsk(i)+ccondsk(i)+eps)

!g           endif
         end do
        end do
      enddo  

      return
end subroutine modalapp
