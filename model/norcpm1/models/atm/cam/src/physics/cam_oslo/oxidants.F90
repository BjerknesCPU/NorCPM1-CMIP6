!OS Used the same structure for reading in oxidants as for reading in 
! emissions. Taken from Match.
module oxidants
!----------------------------------------------------------------------- 
! 
! Purpose: Module to read in emissions and converting the emissions 
! from lat,lon to chunks keeping them in arrys
! Emission module.  .
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!
!
   use spmd_utils,  only: masterproc
   use pmgrid,      only: plon, plat
   use ppgrid,      only: pcols, pver, begchunk, endchunk
   use phys_grid,   only: scatter_field_to_chunk, get_ncols_p

   use physconst,   only: rair,mwdry,avogad
!   use uioinpsubs, only: read_1pathname_camuioinp
   use cam_logfile,  only: iulog
   implicit none
   private
!ny
!   save
!ny


! chunk arrays for storing oxidants
! Unit: molecules/cm3
   real(r8), allocatable, dimension(:,:,:,:) :: &
      co3    ! O3 concentrations (col,pver,chnk,month)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      ch2o2   ! h2o2 concentrations (col,pver,chnk,month)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      coh     ! oh concentrations(col,pver,chnk,month)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      cstoh   ! New oh fields only for troposphere 
!oh concentrations(col,pver,chnk,month)
! full pathnames for oxidant data files
   character*120 ::  fnoxid_o3	! o3
   character*120 ::  fnoxid_h2o2    ! h2o2 
   character*120 ::  fnoxid_oh  ! oh
   character*120 ::  fnoxid_stoh  ! Oh from CTM-1 Needed for values in the stratosphere   
   logical f19tof09 ! True if the data must be converted from f19_25 to f09_12
   public get_oxidant_pathnames ! Gets pathnames for oxidation data files
   public oxini              ! Reads emissions data and converts to chunks
   public getoxid              ! Interface to reading oxidations   

contains


 



   subroutine get_oxidant_pathnames
!-----------------------------------------------------------------------
!
! Purpose: Reads oxidation pathnames from mirage2.inp and the type of file
!          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

   use dycore, only: get_resolution

   implicit none
   character(len=32) :: hgrid
!-----------------------------------------------------------------------
#ifdef READINPUTFILE   
!-----------------------------------------------------------------------
   call read_1pathname_camuioinp( fnoxid_o3,  &
      'OXIDANT PATHNAME - O3' )
   if ( masterproc ) write(iulog,9100)  &
      'O3 concentration dataset: ',  &
      trim(fnoxid_o3)

   call read_1pathname_camuioinp( fnoxid_h2o2,  &
      'OXIDANT PATHNAME - H2O2' )
   if ( masterproc ) write(iulog,9100)  &
      'H2O2 concentration dataset: ',  &
      trim(fnoxid_h2o2)

   call read_1pathname_camuioinp( fnoxid_oh,  &
      'OXIDANT PATHNAME - OH' )
   if ( masterproc ) write(iulog,9100)  &
      'OH concentration dataset: ',  &
      trim(fnoxid_oh)



   if (masterproc) write(iulog,*)

9100	format( a / 4x, a )
#endif

fnoxid_o3='inputdata/atm/cam/camoslo/o3_ctm2.FV19.asciiflat'
fnoxid_h2o2='inputdata/atm/cam/camoslo/h2o2_ctm2.FV19.asciiflat'
fnoxid_oh='inputdata/atm/cam/camoslo/oh_ctm2.FV19.asciiflat'
fnoxid_stoh='inputdata/atm/cam/camoslo/oh.FV19.asciiflat'

  f19tof09=.false.

  hgrid = get_resolution()


  if(trim(hgrid)=='0.9x1.25') f19tof09=.true.



   return
   end subroutine get_oxidant_pathnames




   subroutine oxini

   use error_messages, only: alloc_err, handle_ncerr

#if ( defined SPMD )
   use mpishorthand
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
   integer :: i, j,k, lat, m     ! longitude, level, latitude, time indices
   integer :: n                  ! constituent index
   integer :: istat                 ! error return
!  horizontal grid specifier


! temporary global arrays for reading concentration of oxidants
   real(r8) :: xem_xy(plon,plat)
   real(r8) :: xem_xymon(plon,plat,12)
   real(r8) :: xem_xysea(plon,plat,4)
   real(r8) :: xem_xzy(plon,pver,plat)
   real(r8) :: xox_xzymon(plon,pver,plat,12)


!
! allocate chunk arrays to store oxidants
!
  allocate( co3(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'oxini','co3', &
       pcols*pver*(endchunk-begchunk+1)*12 )

  allocate( ch2o2(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'oxini','ch2o2', &
       pcols*pver*(endchunk-begchunk+1)*12 )

  allocate( coh(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'oxini','coh', &
       pcols*pver*(endchunk-begchunk+1)*12 )

  allocate( cstoh(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'oxini','cstoh', &
       pcols*pver*(endchunk-begchunk+1)*12 )



!
! read oxidants to temporary global arrays, convert units, 
!    and scatter to chunk arrays
! SPMD: Master does all the work.  Sends needed info to slaves
!

!xox_xzymon(:,:,:)=0._r8
   xox_xzymon=0.0_r8 
   if (masterproc) then

!   write(iulog,*) fnoxid_o3 
   if (fnoxid_o3 /= 'none') call read_asciiflat_xzyt( 	&
       fnoxid_o3, xox_xzymon, pver, 12, .true. ,f19tof09)
!   xox_xzymon = xox_xzymon 
   end if
   do m=1,12
     call scatter_field_to_chunk( 1,pver,1,plon, xox_xzymon(:,:,:,m),co3(:,:,:,m) )
   end do 

   if (masterproc) then
   xox_xzymon=0.0_r8  
!   write(iulog,*) fnoxid_h2o2 
   if (fnoxid_h2o2 /= 'none') call read_asciiflat_xzyt( 	&
       fnoxid_h2o2, xox_xzymon, pver, 12, .true. ,f19tof09)
!   xox_xzymon = xox_xzymon 
   end if
   do m=1,12
   call scatter_field_to_chunk( 1,pver,1,plon, xox_xzymon(:,:,:,m),ch2o2(:,:,:,m) )
   end do 

   if (masterproc) then
   xox_xzymon=0.0_r8  
!   write(iulog,*) fnoxid_oh 
   if (fnoxid_oh /= 'none') call read_asciiflat_xzyt( 	&
       fnoxid_oh, xox_xzymon, pver, 12, .true. ,f19tof09)
!   xox_xzymon = xox_xzymon 
   end if
   do m=1,12
   call scatter_field_to_chunk( 1,pver,1,plon, xox_xzymon(:,:,:,m),coh(:,:,:,m) )
   end do 


  if (masterproc) then
   xox_xzymon=0.0_r8  
!   write(iulog,*) fnoxid_stoh 
   if (fnoxid_stoh /= 'none') call read_asciiflat_xzyt( 	&
       fnoxid_stoh, xox_xzymon, pver, 12, .true.,f19tof09 )
!   xox_xzymon = xox_xzymon 
   end if
   do m=1,12
   call scatter_field_to_chunk( 1,pver,1,plon,xox_xzymon(:,:,:,m),cstoh(:,:,:,m) )
   end do 






   return
   end subroutine oxini




   subroutine getoxid(lchnk, ncol,mth, t, pmid,lo3,lh2o2,loh,&
           cldv,qh2o2,dqh2o2)



   implicit none

!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns
  integer , intent(in)   :: mth              ! Month number
  real(r8), intent(in)   :: t(pcols,pver)    ! Temperature
  real(r8), intent(in)   :: pmid(pcols,pver) ! Air pressure in Pa

  real(r8), intent(out)  :: lo3(pcols,pver) ! O3 concentration
  real(r8), intent(out)  :: lh2o2(pcols,pver) ! H2O2 concentration from file
  real(r8), intent(out)  :: loh(pcols,pver) ! OH concentration
  real(r8), intent(in)   :: cldv(pcols,pver) ! Local cloud volume
  real(r8), intent(in)   :: qh2o2(pcols,pver) ! H2O2 from tracer-array
  real(r8), intent(out)  :: dqh2o2(pcols,pver) ! Replenishment rate of H2O2
!
! Local variables.
!
   integer :: i, k, m                 ! longitude, level indices
   real(r8) :: cavogadro
  real(r8) caircell(pcols,pver),aircon(pcols,pver),rhoair(pcols,pver)
  real(r8) :: replh2o2(pcols,pver) !Replacement rate of h2o2
  real(r8) :: cldmax(pcols)

!     real(r8) loh(pcols,pver),lo3(pcols,pver),lh2o2(pcols,pver),lgprh2o2(pcols,pver)
!-----------------------------------------------------------------------
! 


! season:  DecJanFeb=1, MarAprMay=2, JunJulAug=3, SepOctNov=4
!      season = (mth + 3)/3
!      if (season > 4) season = season - 4
!      if (season < 1) season = season + 4

! Convert constants to cgs prefix c
       cavogadro = avogad/1000._r8

   do k=1,pver
      do i=1,ncol   
         rhoair(i,k)=pmid(i,k)/(rair*t(i,k))
! Going from kgm to gcm
         aircon(i,k)=0.001_r8*rhoair(i,k)/mwdry
        caircell(i,k) = aircon(i,k) * cavogadro
       end do
   end do

   do i=1,ncol     
     do k=1,pver
       lo3(i,k)=co3(i,k,lchnk,mth)
       lh2o2(i,k)=ch2o2(i,k,lchnk,mth)

      end do	
      do k=1,12
       loh(i,k)=cstoh(i,k,lchnk,mth)
      end do
        do k=13,26
          loh(i,k)=coh(i,k,lchnk,mth)*caircell(i,k)
        end do
   end do	
     
!   do k=1,pver
!     do i=1,ncol
!        lo3(i,k) = lo3(i,k)/caircell(i,k)
!        lh2o2(i,k) = lh2o2(i,k)/caircell(i,k)
!      end do		
!   end do		


    cldmax(:)=0._r8
   do k=1,pver
      do i=1,ncol
        cldmax(i)=max(cldv(i,k),cldmax(i))
!        replh2o2(i,k)=max(1._r8-cldmax(i),0.083_r8)*2.8e-4_r8
        replh2o2(i,k)=max(0.0833_r8,min(1._r8,(1.1_r8-cldmax(i))*(1.1_r8-cldmax(i))))& 
        *2.8e-4_r8
	dqh2o2(i,k)=replh2o2(i,k)*(lh2o2(i,k)-qh2o2(i,k))
      end do
   end do	
   return

end subroutine getoxid




end module oxidants
