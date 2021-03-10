module mass
! Module for calculating column integrals
  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst, only: gravit
  use ppgrid

  implicit none
  private

  public cmidry
  public cmi1d
  public cmidpdivg

contains

      subroutine cmidry( lchnk,ncol,pdel, sh, as,cmi )  

 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: pdel(pcols,pver)     ! Delta p
   real(r8), intent(in) :: sh(pcols,pver)       ! Specific humidity
   real(r8), intent(in) :: as(pcols,pver)       ! Advected species   
   real(r8), intent(out) :: cmi(pcols)          ! Column mass integral


! local variables

      integer :: i, k
      real(r8) :: rgravit
!-----------------------------------------------------------------------
      rgravit=1._r8/gravit
      do i = 1, ncol
         cmi(i) = 0.0_r8
      end do
      do k = 1, pver
         do i = 1, ncol
            cmi(i) = cmi(i) + as(i,k)*(1._r8-sh(i,k))*pdel(i,k)
         end do
      end do
      do i = 1, ncol
         cmi(i) = cmi(i) * rgravit
      end do
      return
      end subroutine cmidry


      subroutine cmi1d( lchnk,ncol,pdel, sh, as,cmi )  

 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: pdel(pcols,pver)     ! Delta p
   real(r8), intent(in) :: sh(pcols,pver)       ! Specific humidity
   real(r8), intent(in) :: as(pcols)            ! field   
   real(r8), intent(out) :: cmi(pcols)          ! Column mass integral


! local variables

      integer :: i
      real(r8) :: rgravit
!-----------------------------------------------------------------------
      rgravit=1._r8/gravit
      do i = 1, ncol
         cmi(i) = 0.0_r8
      end do
!      do k = 1, pver
         do i = 1, ncol
            cmi(i) = cmi(i) + as(i)*(1._r8-sh(i,pver))*pdel(i,pver)
         end do
!      end do
      do i = 1, ncol
         cmi(i) = cmi(i) * rgravit
      end do
      return
      end subroutine cmi1d


      subroutine cmidpdivg( lchnk,ncol,pdel, as,cmi )  

 

! arguments
   integer,  intent(in) :: lchnk                ! chunk identifier
   integer,  intent(in) :: ncol                 ! number of atmospheric column
   real(r8), intent(in) :: pdel(pcols,pver)     ! Delta p
!   real(r8), intent(in) :: sh(pcols,pver)       ! Specific humidity
   real(r8), intent(in) :: as(pcols,pver)       ! Advected species   
   real(r8), intent(out) :: cmi(pcols)          ! Column mass integral


! local variables

      integer :: i, k
      real(r8) :: rgravit
!-----------------------------------------------------------------------
      rgravit=1._r8/gravit
      do i = 1, ncol
         cmi(i) = 0.0_r8
      end do
      do k = 1, pver
         do i = 1, ncol
            cmi(i) = cmi(i) + as(i,k)*pdel(i,k)
         end do
      end do
      do i = 1, ncol
         cmi(i) = cmi(i) * rgravit
      end do
      return
      end subroutine cmidpdivg



end module mass
