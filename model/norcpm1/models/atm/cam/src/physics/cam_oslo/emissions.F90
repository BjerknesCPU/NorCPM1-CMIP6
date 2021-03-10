
module emissions
!----------------------------------------------------------------------- 
! 
! Purpose: Module to read in emissions and converting the emissions 
! from lat,lon to chunks keeping them in arrys
! Emission module.  .
! At present time all emissions file are hard-coded. 
! Will, should, and must be improved
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
!
!

   use spmd_utils, only: masterproc
   use pmgrid,      only: plon, plat
   use ppgrid,      only: pcols, pver, begchunk, endchunk 

   use phys_grid,   only: scatter_field_to_chunk, get_ncols_p
   use abortutils, only: endrun
   use physconst,   only: gravit
   use wrap_nf
   use netcdf
   use aerosoldef      
   use cam_logfile,  only: iulog
   implicit none
   private


! chunk arrays for storing emissions
!   real(r8), allocatable, dimension(:,:,:) :: &
   real(r8), allocatable, dimension(:,:,:) :: &
      em_so2_ff   ! so2 ff emisson rates (col,lev,chnk)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_so2_bb   ! so2 bb emisson rates (col,lev,chnk)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_so2_volc     ! so2 volcanic emisson rates (col,lev,chnk)
#ifdef AEROCOM20
   real(r8), allocatable, dimension(:,:,:,:) :: &
      em_so2_air        ! so2 air emisson rates (col,lev,chnk,month)
#endif
   real(r8), allocatable, dimension(:,:,:) :: &
      em_om_ff      ! om  ff emisson rates (col,lev,chnk)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      em_om_soa      ! om  ff emisson rates (col,lev,chnk,month)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_om_bb      ! om bb emisson rates (col,lev chnk)
#ifdef AEROCOM20
   real(r8), allocatable, dimension(:,:,:,:) :: &
      em_om_air ! om air emisson rates (col,lev,chnk,month)
#endif
   real(r8), allocatable, dimension(:,:,:) :: &
      em_bc_ff        ! bc ff emisson rates (col,lev,chnk)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_bc_bb        ! bc bb emisson rates (col,lev,chnk)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_bc_air        ! bc air emisson rates (col,lev,chnk)

   real(r8), allocatable, dimension(:,:,:) :: &
      em_dms        ! dms emisson rates (col,chnk,day)

   real(r8), allocatable, dimension(:,:,:) :: &
      em_dust1        ! dust 1 emisson rates (col,chnk,day)
   real(r8), allocatable, dimension(:,:,:) :: &
      em_dust2        ! dust 2 emisson rates (col,chnk,day)


   real(r8), allocatable, dimension(:,:,:) :: &
      em_ss1oc        ! seasalt 1 emisson rates (col,chnk,day)
! full pathnames for emissions data files
!   "ff" means fossil fuel
!   "bb" means biomass burning (usually excludes boreal forest fires)
!   character*120 ::  fnemis_dms_ocean	! dms oceanic
   character*120 ::  fnemis_so2_volc    ! so2 volcanic
   character*120 ::  fnemis_so2_ff  ! so2 ff 
   character*120 ::  fnemis_so2_bb  ! so2 bb 
!   character*120 ::  fnemis_so2_ff_elev ! so2 ff elevated
   character*120 ::  fnemis_om_ff 	! om ff
   character*120 ::  fnemis_om_soa 	! om soa
   character*120 ::  fnemis_om_bb 	! om bb
   character*120 ::  fnemis_bc_ff	! bc ff
   character*120 ::  fnemis_bc_bb	! bc bbb
   character*120 ::  fnemis_bc_air	! bc air
#ifdef AEROCOM20
   character*120 ::  fnemis_om_air	! om air
   character*120 ::  fnemis_so2_air	! so2 air
#endif
   character*120 ::  fnemis_dms(12)	! dms
   character*120 ::  fnemis_dust(12)	! dust
   character*120 ::  fnemis_ss(12)	! seasalt

   character*120 ::  fnemis_ff(3)       ! Fossil fuel filenames SO2,BC,OM
   character*120 ::  fnemis_bb(3)       ! Biomass burning filenames
   character*120 ::  fnemis_air(1)      ! Emissions from aircraft filenames
   character*12   ::  field_ff(3)        ! Field name for fossil fuel 
   character*12   ::  field_bb(3)        ! Field name for fossil fuel 
   character*12   ::  field_air(1)       ! Field name for fossil fuel 
   integer  fflev,bblev,airlev          ! Vertical levels in inputfile

   integer, private :: ncid_dms  ! dms emis. dataset id
   integer, private :: ncid_dust  ! dustemis. dataset id
   integer, private :: ncid_ss  ! sesaltemis. dataset id
   integer, private :: dms_id
   integer, private :: dust1_id,dust2_id
   integer, private :: ss1_id,ss2_id,ss3_id
   integer, private :: ff_id(3),bb_id(3),air_id(1)
      
   integer, private :: ncid_ff(3),ncid_bb(3),ncid_air(1)

   logical f19tof09 ! True if the data must be converted from f19_25 to f09_12
   ! Only needed for SOA, BC,POM and SO2 found at both resolutions.

   integer       :: currem_mnd       ! Month used for daily emis data
   public        :: currem_mnd
   integer       :: currem_year      ! Emission year, default 1850
   public        :: currem_year  
   public get_emis_pathnames ! Gets pathnames for emissions data files
   public emini              ! Reads emissions data and converts to chunks
   public emint              ! Decides if a new month read in is needed
   public emvaryear          ! Monthly mean emissions data for datasets that vary with different years
   public emday              ! Daily emission data
   public emdayf09           ! Daily emission data for f09_12 resolution
   public getem              ! Interface to adding emissions onto tendencies   

contains




   subroutine get_emis_pathnames
!-----------------------------------------------------------------------
!
! Purpose: Defines emission files
! Also reads the current resolution
!          
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   use dycore, only: get_resolution
   implicit none
   integer n
   character(len=32) :: hgrid
!-----------------------------------------------------------------------
   
!-----------------------------------------------------------------------

  hgrid = get_resolution()
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_so2ff_1850-2005.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_so2bb_1850-2005.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_bcff_1850-2005.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_bcbb_1850-2005.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcair_150513_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_pomff_1850-2005.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_pombb_1850-2005.nc'
      field_bb(3)='emiss_pom_bb'

#ifdef RCP85
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_so2ff_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_so2bb_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcff_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcbb_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcair_150513_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_pomff_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_pombb_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'

#endif
  if(trim(hgrid)=='0.9x1.25') then
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_so2ff_f09_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_so2bb_f09_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcff_f09_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcbb_f09_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_bcair_f09_041113_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_pomff_f09_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP85_pombb_f09_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'

  end if

#ifdef RCP45
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_so2ff_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_so2bb_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_bcff_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_bcbb_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_bcair_150513_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_pomff_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP45_pombb_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'
#endif


#ifdef RCP26
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_so2ff_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_so2bb_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcff_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcbb_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcair_150513_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_pomff_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_pombb_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'

  if(trim(hgrid)=='0.9x1.25') then
      fnemis_ff(1)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_so2ff_f09_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_so2bb_f09_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcff_f09_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcbb_f09_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_bcair_f09_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_pomff_f09_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='/work/shared/noresm/inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP26_pombb_f09_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'
  end if
#endif

#ifdef RCP60
      fnemis_ff(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_so2ff_1850-2100.nc'
      field_ff(1)='emiss_so2'
      fnemis_bb(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_so2bb_1850-2100.nc'
      field_bb(1)='emiss_so2_bb'
      fnemis_ff(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_bcff_1850-2100.nc'
      field_ff(2)='emiss_bc'
      fnemis_bb(2)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_bcbb_1850-2100.nc'
      field_bb(2)='emiss_bc_bb'
      fnemis_air(1)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_bcair_150513_1850-2100.nc'
      field_air(1)='emiss_bc'
      fnemis_ff(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_pomff_1850-2100.nc'
      field_ff(3)='emiss_pom'
      fnemis_bb(3)='inputdata/atm/cam/camoslo/emissions/IPCC_emiss_RCP60_pombb_1850-2100.nc'
      field_bb(3)='emiss_pom_bb'
#endif



      fflev=2
      bblev=8
      airlev=16

      fnemis_so2_volc='inputdata/atm/cam/camoslo/emissions/aerocom.soxvolc.FV19.asciiflat'
      fnemis_om_soa='inputdata/atm/cam/camoslo/emissions/aerocom.pomsoa.FV19.asciiflat'


  if(trim(hgrid)=='0.9x1.25') then
      fnemis_so2_ff='inputdata/atm/cam/camoslo/emissions/1850.soxff.FV09.asciiflat'
      fnemis_so2_bb='inputdata/atm/cam/camoslo/emissions/1850.soxbb.FV09.asciiflat'
      fnemis_so2_volc='inputdata/atm/cam/camoslo/emissions/aerocom.soxvolc.FV19.asciiflat'
      fnemis_bc_ff='inputdata/atm/cam/camoslo/emissions/1850.bcff.FV09.asciiflat'
      fnemis_bc_bb='inputdata/atm/cam/camoslo/emissions/1850.bcbb.FV09.asciiflat'
      fnemis_bc_air='inputdata/atm/cam/camoslo/emissions/1850.bcair.FV09.asciiflat'
      fnemis_om_ff='inputdata/atm/cam/camoslo/emissions/1850.pomff.FV09.asciiflat'
      fnemis_om_bb='inputdata/atm/cam/camoslo/emissions/1850.pombb.FV09.asciiflat'
      fnemis_om_soa='inputdata/atm/cam/camoslo/emissions/aerocom.pomsoa.FV19.asciiflat'
  end if

      fnemis_dms(1)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200001.nc'
      fnemis_dms(2)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200002.nc'
      fnemis_dms(3)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200003.nc'
      fnemis_dms(4)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200004.nc'
      fnemis_dms(5)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200005.nc'
      fnemis_dms(6)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200006.nc'
      fnemis_dms(7)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200007.nc'
      fnemis_dms(8)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200008.nc'
      fnemis_dms(9)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200009.nc'
      fnemis_dms(10)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200010.nc'
      fnemis_dms(11)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200011.nc'
      fnemis_dms(12)='inputdata/atm/cam/camoslo/emissions/dmsgcmFV19_200012.nc'


      fnemis_dust(1)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200001.nc'
      fnemis_dust(2)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200002.nc'
      fnemis_dust(3)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200003.nc'
      fnemis_dust(4)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200004.nc'
      fnemis_dust(5)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200005.nc'
      fnemis_dust(6)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200006.nc'
      fnemis_dust(7)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200007.nc'
      fnemis_dust(8)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200008.nc'
      fnemis_dust(9)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200009.nc'
      fnemis_dust(10)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200010.nc'
      fnemis_dust(11)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200011.nc'
      fnemis_dust(12)='inputdata/atm/cam/camoslo/emissions/dustgcmFV19_200012.nc'

      fnemis_ss(1)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200001.nc'
      fnemis_ss(2)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200002.nc'
      fnemis_ss(3)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200003.nc'
      fnemis_ss(4)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200004.nc'
      fnemis_ss(5)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200005.nc'
      fnemis_ss(6)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200006.nc'
      fnemis_ss(7)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200007.nc'
      fnemis_ss(8)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200008.nc'
      fnemis_ss(9)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200009.nc'
      fnemis_ss(10)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200010.nc'
      fnemis_ss(11)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200011.nc'
      fnemis_ss(12)='inputdata/atm/cam/camoslo/emissions/saltgcmFV19_200012.nc'

!write(iulog,*) 'fnemis_ss(12) ',fnemis_ss(12)

  f19tof09=.false.

  hgrid = get_resolution()


  if(trim(hgrid)=='0.9x1.25') f19tof09=.true.



   return
   end subroutine get_emis_pathnames






   subroutine emini

   use error_messages, only: alloc_err, handle_ncerr
   use time_manager, only: get_curr_date
#if ( defined SPMD )
   use mpishorthand
#endif
#ifdef CMIP6
   use ppgrid, only : begchunk, endchunk, pver
#endif
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!#include <comctl.h>
!-----------------------------------------------------------------------
!#include <comlun.h>
!-----------------------------------------------------------------------
!#include <aermodes.h>
   integer :: i, j,k, lat, m     ! longitude, level, latitude, time indices
   integer :: n                  ! constituent index
   integer :: istat                 ! error return

! temporary global arrays for reading emission rates
   real(r8) :: xem_xy(plon,plat)
   real(r8) :: xem_xymon(plon,plat,12)
   real(r8) :: xem_xysea(plon,plat,4)
   real(r8) :: xem_xzymon(plon,pver,plat,12)
   real(r8) :: xem_xzy(plon,pver,plat)
   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: emyear ! Emission year

         call get_curr_date(yr, mon, day, ncsec)

!
! allocate chunk arrays to store emissions
!

  allocate( em_so2_ff(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_so2_ff', &
       pcols*pver*(endchunk-begchunk+1))

  allocate( em_so2_bb(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_so2_bb', &
       pcols*pver*(endchunk-begchunk+1))

  allocate( em_so2_volc(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_so2_volc', &
       pcols*(endchunk-begchunk+1)*pver )

  allocate( em_om_ff(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_om_ff', &
       pcols*pver*(endchunk-begchunk+1))

  allocate( em_om_soa(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'emini','em_om_soa', &
       pcols*pver*(endchunk-begchunk+1)*12)

  allocate( em_om_bb(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_om_bb', &
       pcols*pver*(endchunk-begchunk+1))

  allocate( em_bc_ff(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_bc_ff', &
       pcols*pver*(endchunk-begchunk+1))

  allocate( em_bc_bb(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_bc_bb', &
       pcols*pver*(endchunk-begchunk+1))

 allocate( em_bc_air(pcols,pver,begchunk:endchunk), stat=istat )
  call alloc_err( istat, 'emini','em_bc_air', &
       pcols*pver*(endchunk-begchunk+1))

#ifdef AEROCOM20

 allocate( em_so2_air(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'emini','em_so2_air', &
       pcols*pver*(endchunk-begchunk+1)*12 )
 allocate( em_om_air(pcols,pver,begchunk:endchunk,12), stat=istat )
  call alloc_err( istat, 'emini','em_om_air', &
       pcols*pver*(endchunk-begchunk+1)*12 )
#endif

!
! read emissions to temporary global arrays, convert units, 
!    and scatter to chunk arrays
! SPMD: Master does all the work.  Sends needed info to slaves
!

  xem_xzy = 0.0_r8
  if (masterproc) then
   if (fnemis_so2_volc /= 'none') call read_asciiflat_xzyt( 	&
       fnemis_so2_volc, xem_xzy, pver, 1, .false. ,f19tof09)
!   xem_xzy = xem_xzy
  end if

   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzy, em_so2_volc )

 
   xem_xzymon = 0.0_r8
   if (masterproc) then
!   xem_xzymon = 0.0_r8
   if (fnemis_om_soa /= 'none') call read_asciiflat_xzyt( 	&
       fnemis_om_soa, xem_xzymon,  pver, 12, .false. , f19tof09  )
!   xem_xzymon = xem_xzymon 
   end if
   do m=1,12
!   call scatter_field_to_chunk( 1,pver,12,plon, xem_xzymon, em_om_soa )
     call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,m), em_om_soa(:,:,:,m) )
   end do 

#ifdef AEROCOM20
   if (masterproc) then
   xem_xzymon = 0.0_r8
   if (fnemis_om_air /= 'none') call read_asciiflat_xzyt( 	&
       fnemis_om_air, xem_xzymon,  pver, 12, .false. , .false.  )
   end if
   call scatter_field_to_chunk( 1,pver,12,plon, xem_xzymon, em_om_air )

   if (masterproc) then
   xem_xzymon = 0.0_r8
   if (fnemis_so2_air /= 'none') call read_asciiflat_xzyt( 	&
       fnemis_so2_air, xem_xzymon,  pver, 12, .false. , .false.  )
   end if
   call scatter_field_to_chunk( 1,pver,12,plon, xem_xzymon, em_so2_air )


#endif


! Daily emissions

!Allocate  DMS+ dust + SS

      allocate( em_dms(pcols,begchunk:endchunk,31), stat=istat )
      call alloc_err( istat, 'emini','em_dms', &
        pcols*(endchunk-begchunk+1)*31 )

      allocate( em_dust1(pcols,begchunk:endchunk,31), stat=istat )
      call alloc_err( istat, 'emini','em_dust1', &
        pcols*(endchunk-begchunk+1)*31 )

      allocate( em_dust2(pcols,begchunk:endchunk,31), stat=istat )
      call alloc_err( istat, 'emini','em_dust2', &
        pcols*(endchunk-begchunk+1)*31 )

      allocate( em_ss1oc(pcols,begchunk:endchunk,31), stat=istat )
      call alloc_err( istat, 'emini','em_ss1oc', &
         pcols*(endchunk-begchunk+1)*31 )
      currem_mnd=mon
   if (f19tof09) then 
     call emdayf09(mon)
   else
     call emday(mon)
   end if
   emyear=max(yr,1850)
   emyear=min(emyear,2100)
#if defined(AER1850) || defined(AERYR1850)
   emyear=1850
#endif
#if defined(AER2000) || defined(AERYR2000)
   emyear=2000
#endif
   currem_year=emyear
   call emvaryear(mon,emyear)




   return
   end subroutine emini





   subroutine emvaryear(mon,year)

  use error_messages, only: alloc_err, handle_ncerr
   use ioFilemod, only : getfil
   implicit none
!  Arguments
   integer,intent(in)  :: mon  ! Current model month
   integer,intent(in)  :: year  ! Current model year
!   logical,intent(in)  :: init ! Initiate the emission arrays or not.

! local
!   integer :: i, j,k, lat, m     ! longitude, level, latitude, time indices
!   integer :: n                  ! constituent index
   integer ind_decade1            ! Decadal indexes used for interpolation
   integer ind_decade2
   real(r8) frac1,frac2           !  Linear averaging factors
!   Frac1 = (1840+10*ind_decade2)-year)/10
!   Frac2 = 1-frac1 
   integer :: istat                 ! error ret
   real(r8) :: xem_xyzmon1ff(plon,plat,fflev) ! Decadal value for (year-1840)/10 
   real(r8) :: xem_xyzmon2ff(plon,plat,fflev) ! Decadal value for 1+(year-1840)/10

   real(r8) :: xem_xyzmon1bb(plon,plat,bblev) ! Decadal value for (year-1840)/10 
   real(r8) :: xem_xyzmon2bb(plon,plat,bblev) ! Decadal value for 1+(year-1840)/10

   real(r8) :: xem_xyzmon1air(plon,plat,airlev) ! Decadal value for (year-1840)/10 
   real(r8) :: xem_xyzmon2air(plon,plat,airlev) ! Decadal value for 1+(year-1840)/10


   real(r8) :: xem_xzymon(plon,pver,plat,3) ! Averaged values, note change of indexes
   integer :: i,j,k,m,invk
   integer :: nff,nbb,nair
   character(len=256) locfn      ! local filename

   integer londimid              ! netcdf id for longitude dimension
   integer latdimid              ! netcdf id for latitude dimension
   integer levdimid              ! netcdf id for level dimension
!   integer lonid                 ! netcdf id for longitude variable
!   integer latid                 ! netcdf id for latitude variable
!   integer levid                 ! netcdf id for level variable
   integer timeid                ! netcdf id for time variable
   integer lonsiz     ! size of longitude dimension on oxidants dataset
   integer levsiz     ! size of level dimension on oxidants dataset
   integer latsiz     ! size of latitude dimension on oxidants dataset
   integer timsiz     ! size of time dimension on oxidants dataset
   integer cnt5(5)               ! array of counts for each dimension
   integer strt5(5)              ! array of starting indices
! Calculation of decadal index
   ind_decade1=int(real(year-1840,r8)/10._r8)
   ind_decade2=ind_decade1+1
! Calculation of averaging factors
   frac1 = (1840._r8+10._r8*real(ind_decade2,r8)-real(year,r8))/10._r8
   frac2 = 1._r8-frac1

   if (masterproc) then
     xem_xzymon(:,:,:,:)=0._r8
     do nff=1,3
       xem_xyzmon1ff(:,:,:)=0._r8
       xem_xyzmon2ff(:,:,:)=0._r8     
       call getfil(fnemis_ff(nff), locfn,0)

       call wrap_open(locfn, 0, ncid_ff(nff))
       call wrap_inq_dimid( ncid_ff(nff), 'lon', londimid)
       call wrap_inq_dimid( ncid_ff(nff), 'lat', latdimid)
       call wrap_inq_dimid( ncid_ff(nff), 'lev',levdimid  )
       call wrap_inq_dimid( ncid_ff(nff), 'month',timeid  )

       call wrap_inq_dimlen( ncid_ff(nff), londimid, lonsiz   )
       if (lonsiz /= plon) then
         write(iulog,*)'EMVARYEAR: lonsiz=',lonsiz,' must = ',plon
         call endrun
       end if

       call wrap_inq_dimlen( ncid_ff(nff), latdimid, latsiz   )
       if (latsiz /= plat) then
         write(iulog,*)'EMVARYEAR: latsiz=',latsiz,' must = ',plat
         call endrun
       end if

       call wrap_inq_dimlen( ncid_ff(nff), levdimid, levsiz   )
       if (levsiz > fflev) then
         write(iulog,*)'EMVARYEAR: levsiz=',levsiz
         write(iulog,*)'EMVARYEAR: fflev=',fflev
         call endrun
      end if

       call wrap_inq_dimlen( ncid_ff(nff), timeid, timsiz   )
       if (timsiz > 12) then
         write(iulog,*)'EMVARYEAR: timesiz=',timsiz
         call endrun
      end if


       call wrap_inq_varid( ncid_ff(nff), field_ff(nff) , ff_id(nff) )

       strt5(1) = 1
       strt5(2) = 1
       strt5(3) = 1
       strt5(4) = mon
       strt5(5) = ind_decade1
       cnt5(1)  = lonsiz
       cnt5(2)  = latsiz
       cnt5(3)  = levsiz
       cnt5(4)  = 1
       cnt5(5)  = 1
       call wrap_get_vara_realx(ncid_ff(nff),ff_id(nff),strt5,cnt5,xem_xyzmon1ff)
       strt5(5)=ind_decade2
       call wrap_get_vara_realx(ncid_ff(nff),ff_id(nff),strt5,cnt5,xem_xyzmon2ff)
       write(iulog,*) 'Read new fossil fuel em. file ',fnemis_ff(nff)
       write(iulog,*) 'Fraction from decade 1', 10*ind_decade1+1840,frac1
       write(iulog,*) 'Read new month ',mon
!       do m=1,12
         do j=1,plat
           do k=pver,pver-fflev+1,-1
             invk=pver+1-k
             do i=1,plon
               xem_xzymon(i,k,j,nff)=frac1*xem_xyzmon1ff(i,j,invk)+ &
                 frac2*xem_xyzmon2ff(i,j,invk)
             end do
           end do
         end do
!       end do   
     end do  
     do nff=1,3
        call wrap_close(ncid_ff(nff))
     end do
   end if  
   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,1), em_so2_ff )
   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,2), em_bc_ff )
   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,3), em_om_ff )



   if (masterproc) then
     xem_xzymon(:,:,:,:)=0._r8
     do nbb=1,3
       xem_xyzmon1bb(:,:,:)=0._r8
       xem_xyzmon2bb(:,:,:)=0._r8  
       call getfil(fnemis_bb(nbb), locfn,0)
       call wrap_open(locfn, 0, ncid_bb(nbb))
       call wrap_inq_dimid( ncid_bb(nbb), 'lon', londimid)
       call wrap_inq_dimid( ncid_bb(nbb), 'lat', latdimid)
       call wrap_inq_dimid( ncid_bb(nbb), 'month',timeid  )
       call wrap_inq_dimid( ncid_bb(nbb), 'lev',levdimid  )

       call wrap_inq_dimlen( ncid_bb(nbb), londimid, lonsiz   )
       if (lonsiz /= plon) then
         write(iulog,*)'EMVARYEAR: lonsiz=',lonsiz,' must = ',plon
         call endrun
       end if

       call wrap_inq_dimlen( ncid_bb(nbb), latdimid, latsiz   )
       if (latsiz /= plat) then
         write(iulog,*)'EMVARYEAR: latsiz=',latsiz,' must = ',plat
         call endrun
       end if

       call wrap_inq_dimlen( ncid_bb(nbb), timeid, timsiz   )
       if (timsiz > 12) then
         write(iulog,*)'EMVARYEAR: timesiz=',timsiz
         call endrun
      end if

       call wrap_inq_dimlen( ncid_bb(nbb), levdimid, levsiz   )
       if (levsiz > bblev) then
         write(iulog,*)'EMVARYEAR: levsiz=',levsiz
         write(iulog,*)'EMVARYEAR: bblev=',bblev
         call endrun
      end if
       call wrap_inq_varid( ncid_bb(nbb), field_bb(nbb) , bb_id(nbb) )

       strt5(1) = 1
       strt5(2) = 1
       strt5(3) = 1
       strt5(4) = mon
       strt5(5) = ind_decade1
       cnt5(1)  = lonsiz
       cnt5(2)  = latsiz
       cnt5(3)  = levsiz
       cnt5(4)  = 1
       cnt5(5)  = 1
       call wrap_get_vara_realx (ncid_bb(nbb),bb_id(nbb),strt5,cnt5,xem_xyzmon1bb)

       strt5(5)=ind_decade2
       call wrap_get_vara_realx (ncid_bb(nbb),bb_id(nbb),strt5,cnt5,xem_xyzmon2bb)

       write(iulog,*) 'Read biomass em. file ',fnemis_bb(nbb)
       write(iulog,*) 'Fraction from decade 1', 10*ind_decade1+1840,frac1
!       do m=1,12
         do j=1,plat
           do k=pver,pver-bblev+1,-1
             invk=pver+1-k
             do i=1,plon
               xem_xzymon(i,k,j,nbb)=frac1*xem_xyzmon1bb(i,j,invk)+ &
                 frac2*xem_xyzmon2bb(i,j,invk)
             end do
           end do
         end do
!       end do   
     end do
     do nbb=1,3
       call wrap_close(ncid_bb(nbb))
     end do
   end if  

   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,1), em_so2_bb )
   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,2), em_bc_bb )
   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,3), em_om_bb )


   if (masterproc) then
     xem_xzymon(:,:,:,:)=0._r8
     do nair=1,1
       xem_xyzmon1air(:,:,:)=0._r8
       xem_xyzmon2air(:,:,:)=0._r8     
       call getfil(fnemis_air(nair), locfn,0)

       call wrap_open(locfn, 0, ncid_air(nair))
       call wrap_inq_dimid( ncid_air(nair), 'lon', londimid)
       call wrap_inq_dimid( ncid_air(nair), 'lat', latdimid)
       call wrap_inq_dimid( ncid_air(nair), 'month',timeid  )
       call wrap_inq_dimid( ncid_air(nair), 'lev',levdimid  )

       call wrap_inq_dimlen( ncid_air(nair), londimid, lonsiz   )
       if (lonsiz /= plon) then
         write(iulog,*)'EMVARYEAR: lonsiz=',lonsiz,' must = ',plon
         call endrun
       end if

       call wrap_inq_dimlen( ncid_air(nair), latdimid, latsiz   )
       if (latsiz /= plat) then
         write(iulog,*)'EMVARYEAR: latsiz=',latsiz,' must = ',plat
         call endrun
       end if

       call wrap_inq_dimlen( ncid_air(nair), timeid, timsiz   )
       if (timsiz > 12) then
         write(iulog,*)'EMVARYEAR: timesiz=',timsiz
         call endrun
      end if

       call wrap_inq_dimlen( ncid_air(nair), levdimid, levsiz   )
       if (levsiz > airlev) then
         write(iulog,*)'EMVARYEAR: levsiz=',levsiz
         write(iulog,*)'EMVARYEAR: airlev=',airlev
         call endrun
      end if

       call wrap_inq_varid( ncid_air(nair), field_air(nair) , air_id(nair) )

       strt5(1) = 1
       strt5(2) = 1
       strt5(3) = 1
       strt5(4) = mon
       strt5(5) = ind_decade1
       cnt5(1)  = lonsiz
       cnt5(2)  = latsiz
       cnt5(3)  = levsiz
       cnt5(4)  = 1
       cnt5(5)  = 1
       call wrap_get_vara_realx (ncid_air(nair),air_id(nair),strt5,cnt5,xem_xyzmon1air)

       strt5(5)=ind_decade2
       call wrap_get_vara_realx (ncid_air(nair),air_id(nair),strt5,cnt5,xem_xyzmon2air)

       write(iulog,*) 'Aircraft emis. file ',fnemis_air(nair)
       write(iulog,*) 'Fraction from decade 1', 10*ind_decade1+1840,frac1
!       do m=1,12
         do j=1,plat
           do k=pver,pver-airlev+1,-1
             invk=pver+1-k
             do i=1,plon
               xem_xzymon(i,k,j,nair)=frac1*xem_xyzmon1air(i,j,invk)+ &
                 frac2*xem_xyzmon2air(i,j,invk)
             end do
           end do
         end do
!       end do   
     end do
     do nair=1,1
       call wrap_close(ncid_air(nair))
     end do
   end if  

   call scatter_field_to_chunk( 1,pver,1,plon, xem_xzymon(:,:,:,1), em_bc_air )


   return
   end subroutine emvaryear




   subroutine emday(mon)

   use error_messages, only: alloc_err, handle_ncerr
   use ioFilemod, only : getfil

   implicit none

!  Arguments
   integer,intent(in)  :: mon  ! Current model month
!   logical,intent(in)  :: init ! Initiate the emission arrays or not.

! local
!   integer :: i, j,k, lat, m     ! longitude, level, latitude, time indices
!   integer :: n                  ! constituent index
   integer :: istat                 ! error return
   real(r8) :: xem_xyday1(plon,plat,31)
   real(r8) :: xem_xyday2(plon,plat,31)
   real(r8) :: xem_xyday3(plon,plat,31)
   integer :: n
   character(len=256) locfn      ! local filename

   integer londimid              ! netcdf id for longitude dimension
   integer latdimid              ! netcdf id for latitude dimension
   integer levdimid              ! netcdf id for level dimension
!   integer lonid                 ! netcdf id for longitude variable
!   integer latid                 ! netcdf id for latitude variable
!   integer levid                 ! netcdf id for level variable
   integer timeid                ! netcdf id for time variable
   integer lonsiz     ! size of longitude dimension on oxidants dataset
   integer levsiz     ! size of level dimension on oxidants dataset
   integer latsiz     ! size of latitude dimension on oxidants dataset
   integer timsiz     ! size of time dimension on oxidants dataset
   integer cnt3(3)               ! array of counts for each dimension
   integer strt3(3)              ! array of starting indices


if (masterproc) then
      xem_xyday1(:,:,:)=0.0_r8
      call getfil(fnemis_dms(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_dms)
      call wrap_inq_dimid( ncid_dms, 'lon', londimid)
      call wrap_inq_dimid( ncid_dms, 'lat', latdimid)
      call wrap_inq_dimid( ncid_dms, 'time',timeid  )


      call wrap_inq_dimlen( ncid_dms, londimid, lonsiz   )
      if (lonsiz /= plon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dms, latdimid, latsiz   )
      if (latsiz /= plat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dms, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

      call wrap_inq_varid( ncid_dms, 'dmsflux' , dms_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz


     call wrap_get_vara_realx (ncid_dms,dms_id,strt3,cnt3,xem_xyday1)
	write(iulog,*) 'read new dms emission file ',fnemis_dms(mon)
        write(iulog,*) 'Number of days ',timsiz
        call wrap_close(ncid_dms)
endif
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_dms )




if (masterproc) then
      xem_xyday1(:,:,:)=0.0_r8
      xem_xyday2(:,:,:)=0.0_r8

      call getfil(fnemis_dust(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_dust)
      call wrap_inq_dimid( ncid_dust, 'lon', londimid)
      call wrap_inq_dimid( ncid_dust, 'lat', latdimid)
      call wrap_inq_dimid( ncid_dust, 'time',timeid  )


      call wrap_inq_dimlen( ncid_dust, londimid, lonsiz   )
      if (lonsiz /= plon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dust, latdimid, latsiz   )
      if (latsiz /= plat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dust, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

      call wrap_inq_varid( ncid_dust, 'dustflux01' , dust1_id )
      call wrap_inq_varid( ncid_dust, 'dustflux02' , dust2_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz

     call wrap_get_vara_realx (ncid_dust,dust1_id,strt3,cnt3,xem_xyday1)

     call wrap_get_vara_realx (ncid_dust,dust2_id,strt3,cnt3,xem_xyday2)

	write(iulog,*) 'read new dust emission file ',fnemis_dust(mon)
        write(iulog,*) 'Number of days ',timsiz
        call wrap_close(ncid_dust)
     end if

     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_dust1 )
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday2, em_dust2 )



if (masterproc) then

      xem_xyday1(:,:,:)=0.0_r8
      call getfil(fnemis_ss(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_ss)
      call wrap_inq_dimid( ncid_ss, 'lon', londimid)
      call wrap_inq_dimid( ncid_ss, 'lat', latdimid)
      call wrap_inq_dimid( ncid_ss, 'time',timeid  )

      call wrap_inq_dimlen( ncid_ss, londimid, lonsiz   )
      if (lonsiz /= plon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_ss, latdimid, latsiz   )
      if (latsiz /= plat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_ss, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

	write(iulog,*) 'read new ss emission file ',fnemis_ss(mon)
        write(iulog,*) 'Number of days ',timsiz


      call wrap_inq_varid( ncid_ss, 'saltflux01' , ss1_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz

     call wrap_get_vara_realx (ncid_ss,ss1_id,strt3,cnt3,xem_xyday1)
     call wrap_close(ncid_ss)
end if
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_ss1oc )
   return
   end subroutine emday




   subroutine emdayf09(mon)

   use error_messages, only: alloc_err, handle_ncerr
   use ioFilemod, only : getfil
   implicit none
!  Arguments
   integer,intent(in)  :: mon  ! Current model month
!   logical,intent(in)  :: init ! Initiate the emission arrays or not.

! local
   integer :: i, j,k     ! longitude, latitude,level
!   integer :: n                  ! constituent index
   integer :: istat                 ! error return
   real(r8) :: xem_xyday1(plon,plat,31)
   real(r8) :: xem_xyday2(plon,plat,31)
   real(r8) :: xem_xyday3(plon,plat,31)
   real(r8) :: xem_xydayfv19(144,96,31)
   integer,parameter :: fv19lon=144
   integer,parameter :: fv19lat=96
   integer :: n
   character(len=256) locfn      ! local filename

   integer londimid              ! netcdf id for longitude dimension
   integer latdimid              ! netcdf id for latitude dimension
   integer levdimid              ! netcdf id for level dimension
!   integer lonid                 ! netcdf id for longitude variable
!   integer latid                 ! netcdf id for latitude variable
!   integer levid                 ! netcdf id for level variable
   integer timeid                ! netcdf id for time variable
   integer lonsiz     ! size of longitude dimension on oxidants dataset
   integer levsiz     ! size of level dimension on oxidants dataset
   integer latsiz     ! size of latitude dimension on oxidants dataset
   integer timsiz     ! size of time dimension on oxidants dataset
   integer cnt3(3)               ! array of counts for each dimension
   integer strt3(3)              ! array of starting indices
   
if (masterproc) then
      xem_xyday1(:,:,:)=0.0_r8
      xem_xydayfv19(:,:,:)=0.0_r8
      call getfil(fnemis_dms(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_dms)
      call wrap_inq_dimid( ncid_dms, 'lon', londimid)
      call wrap_inq_dimid( ncid_dms, 'lat', latdimid)
      call wrap_inq_dimid( ncid_dms, 'time',timeid  )


      call wrap_inq_dimlen( ncid_dms, londimid, lonsiz   )
      if (lonsiz /= fv19lon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dms, latdimid, latsiz   )
      if (latsiz /= fv19lat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dms, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

      call wrap_inq_varid( ncid_dms, 'dmsflux' , dms_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz


     call wrap_get_vara_realx (ncid_dms,dms_id,strt3,cnt3,xem_xydayfv19)
	write(iulog,*) 'read new dms emission file ',fnemis_dms(mon)
        write(iulog,*) 'Number of days ',timsiz


     do k=1,cnt3(3)
        do j=1,cnt3(2)*2
            do i=1,cnt3(1)*2
              xem_xyday1(i,j,k)= &
        xem_xydayfv19(nint(real(i,r8)/2._r8),nint(real(j,r8)/2._r8),k)
            end do
         end do
     end do
     call wrap_close(ncid_dms)
endif
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_dms )




if (masterproc) then
      xem_xyday1(:,:,:)=0.0_r8
      xem_xydayfv19(:,:,:)=0.0_r8
      xem_xyday2(:,:,:)=0.0_r8

      call getfil(fnemis_dust(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_dust)
      call wrap_inq_dimid( ncid_dust, 'lon', londimid)
      call wrap_inq_dimid( ncid_dust, 'lat', latdimid)
      call wrap_inq_dimid( ncid_dust, 'time',timeid  )


      call wrap_inq_dimlen( ncid_dust, londimid, lonsiz   )
      if (lonsiz /= fv19lon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dust, latdimid, latsiz   )
      if (latsiz /= fv19lat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_dust, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

      call wrap_inq_varid( ncid_dust, 'dustflux01' , dust1_id )
      call wrap_inq_varid( ncid_dust, 'dustflux02' , dust2_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz

     call wrap_get_vara_realx (ncid_dust,dust1_id,strt3,cnt3,xem_xydayfv19)
     do k=1,cnt3(3)
        do j=1,cnt3(2)*2
            do i=1,cnt3(1)*2
              xem_xyday1(i,j,k)= &
        xem_xydayfv19(nint(real(i,r8)/2._r8),nint(real(j,r8)/2._r8),k)
!             write(6,*) i,j,k,nint(real(i,r8)/2._r8),xem_xyday1(i,j,k)
            end do
         end do
     end do

     xem_xydayfv19(:,:,:)=0.0_r8

     call wrap_get_vara_realx (ncid_dust,dust2_id,strt3,cnt3,xem_xydayfv19)

	write(iulog,*) 'read new dust emission file ',fnemis_dust(mon)
        write(iulog,*) 'Number of days ',timsiz

     do k=1,cnt3(3)
        do j=1,cnt3(2)*2
            do i=1,cnt3(1)*2
              xem_xyday2(i,j,k)= &
        xem_xydayfv19(nint(real(i,r8)/2._r8),nint(real(j,r8)/2._r8),k)
!             write(6,*) i,j,k,nint(real(i,r8)/2._r8),xem_xyday1(i,j,k)
            end do
         end do
     end do


        call wrap_close(ncid_dust)
     end if
 
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_dust1 )
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday2, em_dust2 )



if (masterproc) then

      xem_xyday1(:,:,:)=0.0_r8
      xem_xydayfv19(:,:,:)=0.0_r8
      call getfil(fnemis_ss(mon), locfn,0)

      call wrap_open(locfn, 0, ncid_ss)
      call wrap_inq_dimid( ncid_ss, 'lon', londimid)
      call wrap_inq_dimid( ncid_ss, 'lat', latdimid)
      call wrap_inq_dimid( ncid_ss, 'time',timeid  )




      call wrap_inq_dimlen( ncid_ss, londimid, lonsiz   )
      if (lonsiz /= fv19lon) then
         write(iulog,*)'EMDAY: lonsiz=',lonsiz,' must = ',plon
         call endrun
      end if

      call wrap_inq_dimlen( ncid_ss, latdimid, latsiz   )
      if (latsiz /= fv19lat) then
         write(iulog,*)'EMDAY: latsiz=',latsiz,' must = ',plat
         call endrun
      end if

      call wrap_inq_dimlen( ncid_ss, timeid, timsiz   )
      if (timsiz > 31) then
         write(iulog,*)'EMDAY: timesiz=',timsiz
         call endrun
      end if

	write(iulog,*) 'read new ss emission file ',fnemis_ss(mon)
        write(iulog,*) 'Number of days ',timsiz

!   if (fnemis_bc_bb /= 'none') call read_asciiflat_xzyt( 	&
!       fnemis_bc_bb, xem_xyday,  1, 30, .true. )

      call wrap_inq_varid( ncid_ss, 'saltflux01' , ss1_id )
      call wrap_inq_varid( ncid_ss, 'saltflux02' , ss2_id )
      call wrap_inq_varid( ncid_ss, 'saltflux03' , ss3_id )

      strt3(1) = 1
      strt3(2) = 1
      strt3(3) = 1
      cnt3(1)  = lonsiz
      cnt3(2)  = latsiz
      cnt3(3)  = timsiz

     call wrap_get_vara_realx (ncid_ss,ss1_id,strt3,cnt3,xem_xydayfv19)

     do k=1,cnt3(3)
        do j=1,cnt3(2)*2
            do i=1,cnt3(1)*2
              xem_xyday1(i,j,k)= &
        xem_xydayfv19(nint(real(i,r8)/2._r8),nint(real(j,r8)/2._r8),k)
            end do
         end do
     end do

        call wrap_close(ncid_ss)
end if
     call scatter_field_to_chunk( 1,1,31,plon, xem_xyday1, em_ss1oc )
   return
   end subroutine emdayf09







  subroutine emint
	
! Test if new month or new year with respect to emissions 
! Read in new data if needed

   use time_manager, only: get_curr_date
#if ( defined SPMD )
   use mpishorthand
#endif
 implicit none
 integer  :: yr, mon, day, ncsec ! components of a date
 integer emyear
         call get_curr_date(yr, mon, day, ncsec)

      if (currem_mnd.ne.mon) then
        currem_mnd=mon

        if (f19tof09) then 
         call emdayf09(mon)
        else
         call emday(mon)
        end if
 

      emyear=max(yr,1850)
      emyear=min(emyear,2100)
          currem_year=emyear
#ifdef AERYR1850
          emyear=1850
#endif
#ifdef AERYR2000
          emyear=2000
#endif
          call emvaryear(mon,emyear)
      end if         

	return
	end subroutine emint



#ifdef CMIP6
   subroutine getem(state,lchnk, ncol,&
#else 
   subroutine getem(lchnk, ncol,&
#endif
		    sst, ocnfrac, icefrac, u,v,zm,&
                    mth,day,ncnst,doemis,q ,pdel,dqdt)
   use physconst,     only: mwdry
   use constituents,      only: sflxnam
   use cam_history,        only: outfld
   use mass,         only: cmidry
   use const
#ifdef CMIP6
   use physics_types,only : physics_state
   use mo_extfrc, only : extfrc_set
   use mo_srf_emissions,  only : set_srf_emissions
   use chem_mods, only : extcnt, gas_pcnst
   use physconst, only : rga
#endif

   implicit none

#ifdef CMIP6
   type(physics_state), intent(in)    :: state
#endif


!-----------------------------------------------------------------------
!
! Arguments
!
  integer , intent(in)   :: lchnk            ! chunk identifier
  integer , intent(in)   :: ncol             ! number of atmospheric columns
  integer , intent(in)   :: mth              ! Month number
  integer , intent(in)   :: day              ! day number
  integer, intent(in) :: ncnst           ! number of tracers to transport
  logical, intent(inout) :: doemis(ncnst)  ! flag for new emissions
  real(r8), intent(in) :: q(pcols,pver,ncnst)  ! Tracer array including moistur
  real(r8), intent(in) :: pdel(pcols,pver)       ! thickness between interfaces
  real(r8), intent(out) :: dqdt(pcols,pver,ncnst)  ! Tracer tendency array
  real(r8), intent(in) :: sst(pcols), ocnfrac(pcols), icefrac(pcols)
  real(r8), intent(in) :: u(pcols,pver),v(pcols,pver),zm(pcols,pver)



!
! Local variables.
!
   integer :: i, k, m,mm                 ! longitude, level indices
   integer :: season

   real(r8) :: emsfactor(pver)		
! used to be 1/32.1*mwdry*gravit/pdel Now gravit/pdel
   real(r8) :: emcfactor(pver)          
! Used to be 1/12.*mwdry*gravit/pdel Now gravit/pdel
   real(r8) :: embfactor(pver)     ! gravit/pdel 
   real nudge ! nudging coefficient (/s)
   data nudge/1.e-5_r8/
   save nudge
   real(r8) :: s4fac,s43afac
!,s4acfac
   real(r8) :: emsh,emsm,emsl
   real(r8) :: fomxff,fomxbb
   real(r8) :: fbcaxff
   real(r8) :: fbcnff,fbcnbb
   real(r8) :: fsoan
   real(r8) :: emchff,emcmff,emclff,emchbb,emcmbb,emclbb
   real(r8) :: cmi(pcols)   
   real(r8) :: wspd10(pcols)
! Variables required for prognostic sea salt parameterization
      real(r8) :: W, em_ss1, em_ss2, em_ss3	! W = ocean white cap fraction

   real (r8), parameter :: z0=0.0001_r8  ! m roughness length over oceans. 
! Taken from sea-salt emissions in 
!standard cam (aerosol_intr.F90)

#ifdef CMIP6
   real(r8) :: sflx(ncol,gas_pcnst)
   real(r8) :: extfrc(ncol,pver,max(1,extcnt))
   real(r8), parameter :: m2km  = 1.e-3_r8
   real(r8) :: dqdt_cmip6(pcols,pver,ncnst)  ! Tracer tendency array
   real(r8) :: dqdt_bc_air(pcols,pver) 
   real(r8) :: dqdt_so2_volc(pcols,pver) 
#endif

!-----------------------------------------------------------------------
! 

#ifdef CMIP6
   dqdt_cmip6 = 0._r8
   sflx=0._r8  ! needed as BC_NI not set for surface
   call set_srf_emissions( lchnk, sflx(:,:), ncol )

   extfrc=0._r8
   call extfrc_set( lchnk, m2km * state%zi(:ncol,1:pver+1), extfrc, ncol )
#endif



      s4fac=0.025_r8
      s43afac=0.75_r8
      fomxff=0.75_r8
      fomxbb=0.5_r8
      fbcaxff=0.1_r8
      fbcnff=1._r8-fbcaxff
      fbcnbb=0.5_r8
      fsoan=0.2_r8

! season:  DecJanFeb=1, MarAprMay=2, JunJulAug=3, SepOctNov=4
      season = (mth + 3)/3
      if (season > 4) season = season - 4
      if (season < 1) season = season + 4

!	write(iulog,*) 'i getem'

!          doemis(l_om_n) = .true.
          doemis(l_om_ni) = .true.

          doemis(l_bc_n) = .true.
          doemis(l_bc_ax) = .true.
          doemis(l_bc_ni) = .true.

          doemis(l_so2) = .true.
          doemis(l_so4_pr)=.true.
          doemis(l_dms) = .true.
	  
	  doemis(l_dst_a2)=.true.
	  doemis(l_dst_a3)=.true.
	  doemis(l_ss_a1)=.true.
	  doemis(l_ss_a2)=.true.
	  doemis(l_ss_a3)=.true.

      do i=1,ncol
        wspd10(i)= sqrt(u(i,pver)*u(i,pver)+v(i,pver)*v(i,pver))  
        wspd10(i)=wspd10(i)*log(10._r8/z0)/log(zm(i,pver)/z0)
      end do

      do i=1,ncol     

! emission units are [kg/m^2/s]
! tendency units are [kg/kg_air/s]
! emfactor is 1/(kg_air/m^2)
          do k = 1, pver
!             emsfactor(k) = gravit/pdel(i,k)
!             emcfactor(k) = gravit/pdel(i,k)
!             embfactor(k) = gravit/pdel(i,k)
              emsfactor(k) = gravit/(pdel(i,k)*(1._r8-q(i,k,1)))  ! use dry mass consistent with cmidry diagnostics (Ingo Bethke, 2018)  
              emcfactor(k) = emsfactor(k)
              embfactor(k) = emsfactor(k)
          end do
          do k=1,pver
!            dqdt(i,k,l_om_n) = 1.4_r8*em_om_ff(i,k,lchnk,mth) * emcfactor(k)

            dqdt(i,k,l_om_ni) = (1.4_r8*em_om_ff(i,k,lchnk)+ &
              2.6_r8*em_om_bb(i,k,lchnk)+ &
              1.96_r8*em_om_soa(i,k,lchnk,mth))*emcfactor(k)

#ifdef CMIP6
            dqdt_cmip6(i,k,l_om_ni) = (extfrc(i,k,6) + 1.96_r8*em_om_soa(i,k,lchnk,mth))  *emcfactor(k) 
#endif



#ifdef AEROCOM20
           dqdt(i,k,l_bc_n)=fbcnff*em_bc_ff(i,k,lchnk,mth)*emcfactor(k)
           dqdt(i,k,l_bc_ax)=fbcaxff*em_bc_ff(i,k,lchnk,mth)*emcfactor(k)
#else
            
          dqdt(i,k,l_bc_n) = fbcnff*(em_bc_ff(i,k,lchnk) + &
            em_bc_air(i,k,lchnk)) * emcfactor(k)
          dqdt(i,k,l_bc_ax) = fbcaxff*(em_bc_ff(i,k,lchnk)+  &
            em_bc_air(i,k,lchnk)) *emcfactor(k)
#endif

          dqdt(i,k,l_bc_ni) = em_bc_bb(i,k,lchnk) * emcfactor(k)
#ifdef CMIP6 
!         dqdt_cmip6(i,k,l_bc_n) = (extfrc(i,k,3) + fbcnff*em_bc_air(i,k,lchnk))*emcfactor(k) 
          dqdt_cmip6(i,k,l_bc_n) = extfrc(i,k,3) *emcfactor(k) ! aviation already included
!         dqdt_cmip6(i,k,l_bc_ax) = (extfrc(i,k,4) + fbcaxff*em_bc_air(i,k,lchnk))*emcfactor(k) 
          dqdt_cmip6(i,k,l_bc_ax) = extfrc(i,k,4) *emcfactor(k) ! aviation already included
          dqdt_cmip6(i,k,l_bc_ni) = extfrc(i,k,5)*emcfactor(k) 
          dqdt_bc_air(i,k) = em_bc_air(i,k,lchnk)*emcfactor(k) 
#endif
#ifdef AEROCOM20
          dqdt(i,k,l_bc_ni) = dqdt(i,k,l_bc_ni)+   &
           em_bc_air(i,k,lchnk)*emcfactor(k)
          dqdt(i,k,l_om_ni) = dqdt(i,k,l_om_ni)+   &
           em_om_air(i,k,lchnk)*emcfactor(k)
#endif
	 end do

! Spracklen Oceanic OC source of 5.5 Tg / year
           dqdt(i,pver,l_om_ni)=dqdt(i,pver,l_om_ni)+  &
           1.4_r8*76.5_r8*em_ss1oc(i,lchnk,day)*emcfactor(pver)
#ifdef CMIP6
           dqdt_cmip6(i,pver,l_om_ni)=dqdt_cmip6(i,pver,l_om_ni)+  &
           1.4_r8*76.5_r8*em_ss1oc(i,lchnk,day)*emcfactor(pver)
#endif


	  do k=1,pver
	           
	  dqdt(i,k,l_so2) = &
          (1._r8-s4fac)*(em_so2_volc(i,k,lchnk)+0.5_r8*em_so2_ff(i,k,lchnk)+&
             0.5_r8*em_so2_bb(i,k,lchnk))*emsfactor(k)
#ifdef CMIP6
!         dqdt_cmip6(i,k,l_so2) = (1._r8-s4fac)*(em_so2_volc(i,k,lchnk)+0.5_r8*extfrc(i,k,1))*emcfactor(k)
          dqdt_cmip6(i,k,l_so2) =0.5_r8* extfrc(i,k,1)*emcfactor(k) ! volcanoes already included
          dqdt_so2_volc(i,k) = (1._r8-s4fac)*em_so2_volc(i,k,lchnk)*emcfactor(k) 
#endif

	  dqdt(i,k,l_so4_pr) = &
          s4fac*(em_so2_volc(i,k,lchnk)+0.5_r8*em_so2_ff(i,k,lchnk)+&
             0.5_r8*em_so2_bb(i,k,lchnk))*emcfactor(k)
!          dqdt(i,k,l_so4_ac)=s4fac*s4acfac*em_so2_ff(i,k,lchnk)*emsfactor(k)
#ifdef CMIP6
          dqdt_cmip6(i,k,l_so4_pr) = 0.5_r8*extfrc(i,k,2)*emcfactor(k)
#endif

#ifdef AEROCOM20
          dqdt(i,k,l_so2) = dqdt(i,k,l_so2)+ 0.5_r8*0.995_r8 * &
           em_so2_air(i,k,lchnk)*emsfactor(k)
          dqdt(i,k,l_so4_pr) = dqdt(i,k,l_so4_pr)+0.5_r8*0.005_r8 *   &
           em_so2_air(i,k,lchnk)*emsfactor(k)
#endif


	  end do


! 10 % of the coarse seasalt emissions is emitted as accumulation mode
!in order to increase the numbers of CCN

	  dqdt(i,pver,l_dms) = em_dms(i,lchnk,day)*emsfactor(pver)

	  dqdt(i,pver,l_dst_a2) = em_dust1(i,lchnk,day) * embfactor(pver)
	  dqdt(i,pver,l_dst_a3) = em_dust2(i,lchnk,day) * embfactor(pver)

#ifdef CMIP6
	  dqdt_cmip6(i,pver,l_so2) = dqdt_cmip6(i,pver,l_so2) + 0.5_r8*sflx(i,1)*emsfactor(pver)
	  dqdt_cmip6(i,pver,l_so4_pr) = dqdt_cmip6(i,pver,l_so4_pr) + 0.5_r8*sflx(i,2)*emsfactor(pver)
	  dqdt_cmip6(i,pver,l_bc_n) = dqdt_cmip6(i,pver,l_bc_n) + sflx(i,3)*emsfactor(pver)
	  dqdt_cmip6(i,pver,l_bc_ax) = dqdt_cmip6(i,pver,l_bc_ax) + sflx(i,4)*emsfactor(pver)
	  dqdt_cmip6(i,pver,l_bc_ni) = dqdt_cmip6(i,pver,l_bc_ni) + sflx(i,5)*emsfactor(pver) ! should always be zero
	  dqdt_cmip6(i,pver,l_om_ni) = dqdt_cmip6(i,pver,l_om_ni) + sflx(i,6)*emsfactor(pver)
#endif

!#ifdef PROGSSLT
! Purpose: 
! Prognostic emission of ocean sea salt.
! Derived from Martensson et al. (2003) (JGR-108, 4297,doi:10.1029/
! 2002JD002263.
! 
! valid from 20 nm to ~2500 nm dry diameter (based 
! on lab experiment with artificial sea water)
!	
! currently we recommend that it is combined with 
! the parameterisation by Monahan et al. (1986) for 
! dry diameters > 2-3 um even if it lacks 
! temperature dependence (despite that Bowyer et 
! al. (1990) found a similar dependency in the tank
! from Monahan et al. (1986))
!
! you should not apply Monahan et al. (1986) 
! <0.9 um dry diameter (~0.8 um wet radius, 80% RH)

! Calculations of source strength and size distribution 
! calculate whitecap area
          W = (3.84_r8*10._r8**(-6.0_r8))*(wspd10(i)**3.41_r8)
! account for ocean fraction and sea ice. 
	  W = ocnfrac(i) * (1._r8-icefrac(i)) * W

! Temperature dependence of sea salt emissions based on linear fitting
! of Martensson parameterization using the 3 CAM-Oslo log-normal modes
          em_ss1 = W*(sst(i)*ssa1 + ssb1)
          em_ss2 = W*(sst(i)*ssa2 + ssb2)
          em_ss3 = W*(sst(i)*sst(i)*ssc3+sst(i)*ssa3 + ssb3)
!   if (em_ss1.lt.0) write(6,*) 'negative ss1',W,sst(i),ssa1,ssb1
!   if (em_ss2.lt.0) write(6,*) 'negative ss2',W,sst(i),ssa2,ssb2
!   if (em_ss3.lt.0) write(6,*) 'negative ss3',W,sst(i),ssa3,ssb3,273.15_r8*ssa3
   
!          if (wspd10(i).gt.15._r8) then
!          write(6,*) 'SS source'
!          write(6,*)u(i,pver),v(i,pver),wspd10(i)
!          write(6,*) ocnfrac(i),icefrac(i),sst(i)
!          write(6,*) W,ssa1,ssb1,em_ss1
!          write(6,*) ssa2,ssb2,em_ss2
!          write(6,*) ssc3,ssa3,ssb3,em_ss3
!          end if
!          write(6,*) i,W,sst(i),ssa1,ssb1,icefrac(i),ocnfrac(i)
! update sea salt tendencies
	  dqdt(i,pver,l_ss_a1) = em_ss1 * embfactor(pver)
	  dqdt(i,pver,l_ss_a2) = em_ss2 * embfactor(pver)
	  dqdt(i,pver,l_ss_a3) = em_ss3 * embfactor(pver)





!	  dqdt(i,pver,l_ss_a1) = em_ss1(i,lchnk,day) * embfactor(pver)
!orig	  dqdt(i,pver,l_ss_a2) = em_ss2(i,lchnk,day) * embfactor(pver)
!orig	  dqdt(i,pver,l_ss_a3) = em_ss3(i,lchnk,day) * embfactor(pver)
!ostst3	  dqdt(i,pver,l_ss_a2) = em_ss2(i,lchnk,day) * embfactor(pver) + &
!ostst3          0.1_r8*em_ss3(i,lchnk,day) * embfactor(pver)
!ostst3	  dqdt(i,pver,l_ss_a3) = 0.9_r8*em_ss3(i,lchnk,day) * embfactor(pver)
!caktst4
       end do
#ifdef SHORTRUN
       do m=1,ncui
         mm=ixac-1+m
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,mm), cmi)
        call outfld(sflxnam(mm),cmi,pcols,lchnk)

       end do
#endif
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_dms), cmi)
        call outfld('EMI_DMS',cmi,pcols,lchnk)         

        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_so2), cmi)
        call outfld('EMI_SO2',cmi,pcols,lchnk)         
#ifdef CMIP6
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_so2_volc, cmi)
        call outfld('EMI_SO2_VOLC',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_so2), cmi)
        call outfld('EMI_SO2_CMIP6',cmi,pcols,lchnk)         
#endif 
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_so4_n)+ & 
         dqdt(:,:,l_so4_na)+dqdt(:,:,l_so4_a1)+dqdt(:,:,l_so4_a2)+ &
         dqdt(:,:,l_so4_pr)+dqdt(:,:,l_so4_ac), cmi)
        call outfld('EMI_SO4',cmi,pcols,lchnk)         
#ifdef CMIP6
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_so4_n)+ & 
         dqdt(:,:,l_so4_na)+dqdt(:,:,l_so4_a1)+dqdt(:,:,l_so4_a2)+ &
         dqdt_cmip6(:,:,l_so4_pr)+dqdt(:,:,l_so4_ac), cmi)
        call outfld('EMI_SO4_CMIP6',cmi,pcols,lchnk)         
#endif 

        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_bc_n)+ &
        dqdt(:,:,l_bc_ni)+dqdt(:,:,l_bc_a)+dqdt(:,:,l_bc_ai)+ &
        dqdt(:,:,l_bc_ax)+dqdt(:,:,l_bc_ac), cmi)
        call outfld('EMI_BC',cmi,pcols,lchnk)         
#ifdef CMIP6
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_bc_air, cmi)
        call outfld('EMI_BC_AIR',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_bc_n), cmi)
        call outfld('EMI_BC_N',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_bc_ni), cmi)
        call outfld('EMI_BC_NI',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_bc_ax), cmi)
        call outfld('EMI_BC_AX',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_bc_n)+ &
        dqdt_cmip6(:,:,l_bc_ni)+dqdt(:,:,l_bc_a)+dqdt(:,:,l_bc_ai)+ &
        dqdt_cmip6(:,:,l_bc_ax)+dqdt(:,:,l_bc_ac), cmi)
        call outfld('EMI_BC_CMIP6',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_bc_n), cmi)
        call outfld('EMI_BC_N_CMIP6',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_bc_ni), cmi)
        call outfld('EMI_BC_NI_CMIP6',cmi,pcols,lchnk)         
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_bc_ax), cmi)
        call outfld('EMI_BC_AX_CMIP6',cmi,pcols,lchnk)         
#endif 

        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_om_ni)+ &
        dqdt(:,:,l_om_ai)+dqdt(:,:,l_om_ac), cmi)
        call outfld('EMI_POM',cmi,pcols,lchnk)         
#ifdef CMIP6
        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt_cmip6(:,:,l_om_ni)+ &
        dqdt(:,:,l_om_ai)+dqdt(:,:,l_om_ac), cmi)
        call outfld('EMI_POM_CMIP6',cmi,pcols,lchnk)         
#endif 

        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_dst_a2)+ & 
        dqdt(:,:,l_dst_a3), cmi)
        call outfld('EMI_DUST',cmi,pcols,lchnk)         

        call cmidry( lchnk,ncol,pdel, q(:,:,1), dqdt(:,:,l_ss_a1)+ &
        dqdt(:,:,l_ss_a2)+dqdt(:,:,l_ss_a3), cmi)
        call outfld('EMI_SS',cmi,pcols,lchnk)         
	
#ifdef CMIP6
        ! copy CMIP6 tendencies 
        dqdt(:,:,l_bc_n)=dqdt_cmip6(:,:,l_bc_n)
        dqdt(:,:,l_bc_ax)=dqdt_cmip6(:,:,l_bc_ax)
        dqdt(:,:,l_bc_ni)=dqdt_cmip6(:,:,l_bc_ni)
        dqdt(:,:,l_om_ni)=dqdt_cmip6(:,:,l_om_ni)
        dqdt(:,:,l_so2)=dqdt_cmip6(:,:,l_so2)
        dqdt(:,:,l_so4_pr)=dqdt_cmip6(:,:,l_so4_pr)
#endif 
    return
end subroutine getem




end module emissions
