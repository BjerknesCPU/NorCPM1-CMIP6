subroutine modalapp2d(lchnk,ncol,nlons,ind,Nnatk,Ca,f_c,f_bc,f_aq,Cam,fcm,fbcm,faqm)

!     Calculation of the apportionment of internally mixed SO4, BC and OC 
!     mass between the various background mineral and sea-salt modes. Separated 
!     from pmxsub into a independent subroutine by Alf Kirkevåg on September 
!     12'th, 2005, and converted to 2D for use in parmix on September 15'th.
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
   integer, intent(in)  :: lchnk                ! chunk identifier
   integer, intent(in)  :: ncol                 ! number of atmospheric columns
   integer, intent(in)  :: nlons                ! number of cloudy chunks 
   integer, intent(in)  :: ind(pcols)           ! array of cloudy chunks
   real(r8), intent(in) :: Nnatk(pcols,0:nmodes)  ! aerosol mode number concentration  
   real(r8), intent(in) :: Ca(pcols)            ! internally mixed mass, tot=SO4+OC+BC   
   real(r8), intent(in) :: f_c(pcols)           ! mass fraction (OC+BC)/tot
   real(r8), intent(in) :: f_bc(pcols)          ! mass fraction BC/(OC+BC)
   real(r8), intent(in) :: f_aq(pcols)          ! mass fraction SO4(aq)/SO4
!
! Output arguments
!
   real(r8), intent(out) :: Cam(pcols,nbmodes)  ! modal internal mass, tot=SO4+BC+OC 
   real(r8), intent(out) :: fcm(pcols,nbmodes)  ! modal mass fraction (OC+BC)/tot
   real(r8), intent(out) :: fbcm(pcols,nbmodes) ! modal mass fraction BC/(OC+BC)
   real(r8), intent(out) :: faqm(pcols,nbmodes) ! modal mass fraction SO4(aq)/SO4
!
! Local variables
!
   real(r8) NrD(pcols,nbmodes), NK12(pcols,nbmodes), NK12oc(pcols,nbmodes), & 
            Nag(pcols,nbmodes), NrDtot(pcols), NK12tot(pcols),  &
            NK12octot(pcols), Nagtot(pcols)
   real(r8) fcondk(nbmodes), fcoagk(nbmodes), fcoagock(nbmodes), faqk(nbmodes) 
   real(r8) cabck(nbmodes), caock(nbmodes), caqsk(nbmodes), ccondsk(nbmodes)
   real(r8) es000(nbmodes,imax),es001(nbmodes,imax),es002(nbmodes,imax),es003(nbmodes,imax)

   integer i, kcomp, ir, icol, long


      do i=1,10
        do ir=1,imax
	   es003(i,ir)=nk(i,ir)*dlogr
	   es001(i,ir)=es003(i,ir)*Dmpa(i,ir)*rp(ir)
!	   es002(i,ir)=es003(i,ir)*Kp12(i,ir)     ! Kp12    not used. Instead we use a common:
	   es000(i,ir)=es003(i,ir)*Kp12oc(i,ir)   ! Kp12oc  with r based on 40nm is used (BC/OC/SO4)
	enddo
        do icol=1,ncol             !                              Ait.   a and c
!g           NrD(icol,i)=0.0_r8    ! H2SO4 condensation: on modes 1-4 and 5-10
           NrD(icol,i)=0.0_r8      ! H2SO4 cond./coag.:  on modes 1-4  /  5-10
           NK12(icol,i)=0.0_r8     ! BC coagulation:     on modes         5-10
           NK12oc(icol,i)=0.0_r8   ! OC coagulation:     on modes         5-10  
           Nag(icol,i)=0.0_r8      ! (NH4)2SO4 aq. chem: on modes   4 and 5-10
        end do
     enddo ! i
     do icol=1,ncol
       NrDtot(icol)=0.0_r8  
       NK12tot(icol)=0.0_r8
       NK12octot(icol)=0.0_r8
       Nagtot(icol)=0.0_r8
     enddo  

!    mode 0 is a purely externally mixed BC(ax) mode, and is therefore skipped here.

     do i=1,3
! Modes 1-3: only cond. SO4 is to be internally mixed with these modes
!        do icol=1,ncol
        do long=1,nlons
        icol=ind(long)
!g          if(Nnatk(icol,i)>0.0_r8) then
           do ir=1,imax
             NrD(icol,i)=NrD(icol,i)+Nnatk(icol,i)*es001(i,ir)
           enddo
!g          endif   ! Nnatk>0
        end do
     enddo ! i

     do i=4,4
! Mode 4: cond. and aq. SO4 is to be internally mixed with this mode
!        do icol=1,ncol
        do long=1,nlons
        icol=ind(long)
!g          if(Nnatk(icol,i)>0.0_r8) then
           do ir=1,imax
             NrD(icol,i)=NrD(icol,i)+Nnatk(icol,i)*es001(i,ir)
           enddo
           Nag(icol,i)=Nnatk(icol,i)*nk(i,irc)*l10rp(ir)
           do ir=irc+1,imax
             Nag(icol,i)=Nag(icol,i)+Nnatk(icol,i)*es003(i,ir)
           enddo
!g          endif   ! Nnatk>0
        end do
     enddo ! i

     do i=5,10
! Modes 5-10: assuming that SO4 from condensation is lumped together with SO4 
! from coagulation, opposite of the case for modes 1-4.
!        do icol=1,ncol
        do long=1,nlons
        icol=ind(long)
!g          if(Nnatk(icol,i)>0.0_r8) then
           do ir=1,imax
!g             NrD(icol,i)=NrD(icol,i)+Nnatk(icol,i)*es001(i,ir)    ! --> coag istedet!!!
!f             NK12s(icol,i)=NK12s(icol,i)+Nnatk(icol,i)*es000(i,ir)
             NrD(icol,i)=NrD(icol,i)+Nnatk(icol,i)*es000(i,ir)      ! lump coag. and cond.
             NK12(icol,i)=NK12(icol,i)+Nnatk(icol,i)*es000(i,ir)
             NK12oc(icol,i)=NK12oc(icol,i)+Nnatk(icol,i)*es000(i,ir)
           enddo
           Nag(icol,i)=Nnatk(icol,i)*nk(i,irc)*l10rp(ir)
           do ir=irc+1,imax
             Nag(icol,i)=Nag(icol,i)+Nnatk(icol,i)*es003(i,ir)
           enddo
!g          endif   ! Nnatk>0
        end do
     enddo ! i

! Then the contribution (0 or not) from each mode is added to the total process
! specific masses: 
      do i=1,10
!        do icol=1,ncol
        do long=1,nlons
        icol=ind(long)
!g          if(Nnatk(icol,i)>0.0_r8) then      
            NrDtot(icol)=NrDtot(icol)+NrD(icol,i)
            NK12tot(icol)=NK12tot(icol)+NK12(icol,i)
            NK12octot(icol)=NK12octot(icol)+NK12oc(icol,i)
            Nagtot(icol)=Nagtot(icol)+Nag(icol,i)
!g          endif
          Cam(icol,i)=0.0_r8
          fcm(icol,i)=0.0_r8
          fbcm(icol,i)=0.0_r8
          faqm(icol,i)=0.0_r8
        enddo  
      enddo  

! And finally the contribution from each mode relative to the totals are calculated,
! assuming that the apprtionment of mass for the first iteration (in time) is representative
! for the whole apprtionment process (which is ok for small and moderate masses added): 
      do i=1,10
!        do icol=1,ncol
        do long=1,nlons
        icol=ind(long)
!g          if(Nnatk(icol,i)>0.0_r8) then 

            fcondk(i)=NrD(icol,i)/(NrDtot(icol)+eps)
            fcoagk(i)=NK12(icol,i)/(NK12tot(icol)+eps)
            fcoagock(i)=NK12oc(icol,i)/(NK12octot(icol)+eps)
            faqk(i)=Nag(icol,i)/(Nagtot(icol)+eps)

            cabck(i)=fcoagk(i)*f_c(icol)*f_bc(icol)*Ca(icol)
            caock(i)=fcoagock(i)*f_c(icol)*(1.0_r8-f_bc(icol))*Ca(icol)
            caqsk(i)=faqk(i)*f_aq(icol)*(1.0_r8-f_c(icol))*Ca(icol)
            ccondsk(i)=fcondk(i)*(1.0_r8-f_aq(icol))*(1.0_r8-f_c(icol))*Ca(icol)

            Cam(icol,i)=cabck(i)+caock(i)+caqsk(i)+ccondsk(i)
            fcm(icol,i)=(cabck(i)+caock(i))/(Cam(icol,i)+eps)
            fbcm(icol,i)=cabck(i)/(cabck(i)+caock(i)+eps)
            faqm(icol,i)=caqsk(i)/(caqsk(i)+ccondsk(i)+eps)

!g          endif
        end do
      enddo  

      return
end subroutine modalapp2d
