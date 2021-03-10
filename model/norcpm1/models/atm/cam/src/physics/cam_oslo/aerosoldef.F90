
module aerosoldef

!---------------------------------------------------------------------------------
! Module to set up register aerosols indexes, number of gas and particle 
! species and their scavenging rates. Tables for hunidity growth
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use physconst,    only: mwdry,cpair
  use constituents, only: pcnst,cnst_add, cnst_num_avail
  use cam_logfile,  only: iulog
  use camsrfexch_types, only: hub2atm_setopts

  implicit none
  save
  private          ! Make default type private to the module
 
!
! Public interfaces
!
  public aero_register           ! register consituents
  public aero_init               ! register output
  public aer_phys_init
  public inittabrh


  integer, public, parameter :: ncui=23          ! Total number of constituents
  integer,public, parameter :: ncgas=3           ! number of gas species
  integer,public, parameter :: ncaer=20          ! number of aerosol species
!  integer, public parameter :: ndvel_gas = 1     ! Total number of dry deposited gases
  integer, public:: ixac                     ! index of 1st constituent
  integer, public:: ixae                     ! index of 1st aerosol constituent
  integer, dimension(ncui) :: c_i            ! global constituent index


  integer ::ltot_adv_chem	! number of camuio chem/aerosol advected
! species
  integer :: lbgn_adv_chem	! index of 1st  camuio chem/aerosol adv species
  integer :: lend_adv_chem	! index of last camuio chem/aerosol adv species
  integer :: ltot_nad_chem	! number of camuio chem/aerosol non-advected
				!     species
  integer :: lbgn_nad_chem	! index of 1st  camuio chem/aerosol nad species
  integer :: lend_nad_chem	! index of last camuio chem/aerosol nad species
  integer, public, dimension (pcnst) :: species_class	
! indicates species class (nonmirage,cldphysics, aerosol, gas, othermirage )

  integer, public, parameter ::  spec_class_other = 0
  integer, public, parameter ::  spec_class_cldphysics = 1 
  integer, public, parameter ::  spec_class_aerosol = 2 
  integer, public, parameter ::  spec_class_gas = 3 
  integer, public, parameter ::  n_land_type = 11

! following are species indices for individual camuio species
  integer,public :: 					&
  l_so4_n,    l_so4_na,   l_so4_a1, l_so4_a2, l_so4_ac,             &
  l_bc_n,     l_bc_ax,    l_bc_ni,  l_bc_a,  l_bc_ai,l_bc_ac,     &
  l_om_ni,    l_om_ai  ,l_om_ac,               &
  l_so4_pr,   &
  l_dst_a2,   l_dst_a3,              &
  l_ss_a1,    l_ss_a2,    l_ss_a3

  integer,public  ::					&
  l_so2,      l_dms,      l_qh2o2


! Mapping of proxies for dry deposition
  integer, public,dimension (pcnst) :: nvdgas

  real(r8), public,dimension (pcnst) ::scavin
  real(r8), public,dimension (pcnst) ::scavbe
  real(r8), public,dimension (pcnst) ::scavbecon 
  real(r8), public,dimension (pcnst) ::scavconv
  real(r8), public,dimension (pcnst) ::effsize
  real(r8), public,dimension (pcnst) ::rhopart
  real(r8), public,dimension (pcnst) ::sgpart
  real(r8), public,dimension (pcnst) ::mindepvel
  real(r8), public,dimension (pcnst,3) ::predvel

  real(r8), public,dimension (10)      :: rhtab
  real(r8), public,dimension (10,pcnst):: rdivr0(10,pcnst)

    data rhtab/ 0.0_r8, 0.37_r8, 0.47_r8, 0.65_r8, 0.75_r8, 0.80_r8, 0.85_r8, 0.90_r8, 0.95_r8, 0.98_r8 /

contains





!===============================================================================
  subroutine aero_register
!----------------------------------------------------------------------- 
! 
! Register aerosol modes and indices, should be changed to read in values 
! instead of hard-coding it. 
! 
!-----------------------------------------------------------------------

  use mpishorthand


!---------------------------Local workspace-----------------------------
  integer :: i
  integer :: idone_dum
  integer :: iostatval
  integer :: m,mm                                   ! tracer index
  character*3 :: adv_nad_dum
  character*8 :: name_dum, spec_class_dum
!   character*120 :: longname(ncui)
  character*300 :: txtline
  real(r8), parameter :: zero  = 0._r8
!-----------------------------------------------------------------------
!   register the species
! constituent names, short and long
  character(len=8), dimension(ncui), parameter :: c_names = & 
      (/'DMS', 'SO2', 'H2O2','SO4_N','SO4_NA','SO4_A1','SO4_A2', &
      'SO4_AC','SO4_PR','BC_N','BC_AX','BC_NI','BC_A','BC_AI','BC_AC', &
      'OM_NI', 'OM_AI','OM_AC','DST_A2','DST_A3','SS_A1','SS_A2','SS_A3'/)


  character(len=120), dimension(ncui), parameter :: longname = & 
      (/'DMS ', 'SO2 ', 'H2O2','SO4_N','SO4_NA','SO4_A1','SO4_A2', &
      'SO4_AC','SO4_PR','BC_N','BC_AX','BC_NI','BC_A','BC_AI','BC_AC', &
      'OM_NI', 'OM_AI','OM_AC','DST_A2','DST_A3','SS_A1','SS_A2','SS_A3'/)



  do i=1,ncui
    call cnst_add(c_names(i), mwdry, cpair, zero, c_i(i), &  
       longname(i))
  end do

  ixac=c_i(1)
  ixae=ixac+ncgas
  species_class(:)=0
  do i=ixac,ixac+ncgas-1
    species_class(i)=3
  end do
  do i=ixac+ncgas,pcnst
    species_class(i)=2
  end do 


  do m = 1, ncui
    mm = m+ixac-1 
    if (c_names(m) .eq. 'SO2     ') l_so2      = mm
    if (c_names(m) .eq. 'DMS     ') l_dms      = mm
!	if (c_names(m) .eq. 'H2SO4   ') l_h2so4    = mm
!	if (c_names(m) .eq. 'MSA     ') l_msa      = mm
    if (c_names(m) .eq. 'H2O2   ')  l_qh2o2    = mm

    if (c_names(m) .eq. 'SO4_N   ') l_so4_n     = mm
    if (c_names(m) .eq. 'SO4_NA  ') l_so4_na   = mm
    if (c_names(m) .eq. 'SO4_A1  ') l_so4_a1   = mm
    if (c_names(m) .eq. 'SO4_A2  ') l_so4_a2   = mm
    if (c_names(m) .eq. 'SO4_AC  ') l_so4_ac   = mm
    if (c_names(m) .eq. 'SO4_PR  ') l_so4_pr   = mm 

    if (c_names(m) .eq. 'BC_N   ')  l_bc_n      = mm
    if (c_names(m) .eq. 'BC_AX  ') l_bc_ax     = mm
    if (c_names(m) .eq. 'BC_NI  ')  l_bc_ni    = mm
    if (c_names(m) .eq. 'BC_A   ')  l_bc_a     = mm
    if (c_names(m) .eq. 'BC_AI  ')   l_bc_ai   = mm
    if (c_names(m) .eq. 'BC_AC  ')  l_bc_ac    = mm

!	if (c_names(m) .eq. 'OM_N   ') l_om_n       = mm
    if (c_names(m) .eq. 'OM_NI  ') l_om_ni     = mm
    if (c_names(m) .eq. 'OM_AI  ') l_om_ai     = mm
    if (c_names(m) .eq. 'OM_AC  ') l_om_ac     = mm

!	if (c_names(m) .eq. 'DST_A1 ') l_dst_a1   = mm
    if (c_names(m) .eq. 'DST_A2 ') l_dst_a2   = mm
    if (c_names(m) .eq. 'DST_A3 ') l_dst_a3   = mm

    if (c_names(m) .eq. 'SS_A1  ') l_ss_a1     = mm
    if (c_names(m) .eq. 'SS_A2  ') l_ss_a2     = mm
    if (c_names(m) .eq. 'SS_A3  ') l_ss_a3     = mm


  end do

! Set size and deposition properties

  call aer_phys_init
       
!!   set the "species location pointers"
!   call init_chem_species


!   set the aerosol mode definitions
!	write(iulog,*) '*** halting in chemadd_register_cnst ***'
!	write(iulog,*) '*** before call to init_aer_modes    ***'
!	call endrun
!	call init_aer_modes

    return
  end subroutine aero_register


  subroutine aer_phys_init

  use spmd_utils,   only: masterproc
  integer m

!   write pointers values
  if ( masterproc ) then
    write(iulog,*)
    write(iulog,*) 'init_chem_species -- species location pointers'
    write(iulog,*)
    write(iulog,*) 'l_so2      =', l_so2  
    write(iulog,*) 'l_dms      =', l_dms  
!       write(iulog,*) 'l_h2so4    =', l_h2so4
!       write(iulog,*) 'l_msa      =', l_msa  
    write(iulog,*) 'l_qh2o2    =', l_qh2o2  
!	write(iulog,*) 'l_co       =', l_co   
!	write(iulog,*) 'l_h2o2     =', l_h2o2 
!	write(iulog,*) 'l_ch3o2h   =', l_ch3o2h
!	write(iulog,*) 'l_oh       =', l_oh   
!	write(iulog,*) 'l_o3       =', l_o3   

    write(iulog,*)
    write(iulog,*) 'l_so4_n    =', l_so4_n
    write(iulog,*) 'l_so4_na   =', l_so4_na
    write(iulog,*) 'l_so4_a1   =', l_so4_a1
    write(iulog,*) 'l_so4_a2   =', l_so4_a2
    write(iulog,*) 'l_so4_ac   =', l_so4_ac
    write(iulog,*) 'l_so4_pr     =',l_so4_pr

    write(iulog,*)
    write(iulog,*) 'l_bc_n     =', l_bc_n
    write(iulog,*) 'l_bc_ax    =', l_bc_ax
    write(iulog,*) 'l_bc_ni    =', l_bc_ni
    write(iulog,*) 'l_bc_a     =', l_bc_a 
    write(iulog,*) 'l_bc_ai    =', l_bc_ai 
    write(iulog,*) 'l_bc_ac    =', l_bc_ac

    write(iulog,*)
!       write(iulog,*) 'l_om_n     =', l_om_n
    write(iulog,*) 'l_om_ni     =',l_om_ni
    write(iulog,*) 'l_om_ai     =',l_om_ai
    write(iulog,*) 'l_om_ac     =',l_om_ac

!	write(iulog,*) 'l_dst_a1  =', l_dst_a1
    write(iulog,*) 'l_dst_a2  =', l_dst_a2
    write(iulog,*) 'l_dst_a3  =', l_dst_a3

    write(iulog,*)
    write(iulog,*) 'l_ss_a1    =', l_ss_a1
    write(iulog,*) 'l_ss_a2    =', l_ss_a2
    write(iulog,*) 'l_ss_a3    =', l_ss_a3


  end if

  do m=1,pcnst
    nvdgas(m)=-999888777
  end do
  if (l_dms   .gt. 0) nvdgas(l_dms)   = 15
    if (l_qh2o2   .gt. 0) nvdgas(l_qh2o2)   = 15
!    if (l_msa   .gt. 0) nvdgas(l_msa)   = 11
  if (l_so2   .gt. 0) nvdgas(l_so2)   = 1
!    if (l_h2so4 .gt. 0) nvdgas(l_h2so4) = 5
!        if (l_h2o2  .gt. 0) nvdgas(l_h2o2)  = 6

             

  do m=1,pcnst
    scavin(m) = -999888777._r8
    scavbe(m) = -999888777._r8
    effsize(m) = -999888777._r8
    rhopart(m) = -999888777._r8
    sgpart(m) = -999888777._r8
    mindepvel(m) = -999888777._r8
    predvel(m,:) = -999888777._r8
  end do
  scavin(l_so4_n)=0._r8
  scavin(l_so4_na)=0.25_r8
  scavin(l_so4_a1)=0.25_r8
  scavin(l_so4_a2)=1._r8
  scavin(l_so4_ac)=1._r8
  scavin(l_so4_pr)=0.8_r8

  scavin(l_bc_n)=0._r8
  scavin(l_bc_ax)=0._r8
  scavin(l_bc_ni)=0.1_r8
  scavin(l_bc_a) =0.25_r8
  scavin(l_bc_ai)=0.25_r8
  scavin(l_bc_ac)=1._r8

!    scavin(l_om_n)=0._r8
  scavin(l_om_ni)=0.1_r8
  scavin(l_om_ai)=0.25_r8
  scavin(l_om_ac)=1._r8

  scavin(l_dst_a2)=0.25_r8
  scavin(l_dst_a3)=0.25_r8
  scavin(l_ss_a1)=1._r8
  scavin(l_ss_a2)=1._r8
  scavin(l_ss_a3)=1._r8

  scavbe(l_so4_n)=0.04_r8
  scavbe(l_so4_na)=0.02_r8
  scavbe(l_so4_a1)=0.02_r8
  scavbe(l_so4_a2)=0.01_r8
  scavbe(l_so4_ac)=0.02_r8
  scavbe(l_so4_pr)=0.01_r8

    scavbe(l_bc_n)=0.08_r8
    scavbe(l_bc_ax)=0.01_r8
    scavbe(l_bc_ni)=0.02_r8
    scavbe(l_bc_a) =0.02_r8
    scavbe(l_bc_ai)=0.02_r8
    scavbe(l_bc_ac)=0.02_r8

!    scavbe(l_om_n)=0.02_r8
    scavbe(l_om_ni)=0.02_r8
    scavbe(l_om_ai)=0.02_r8
    scavbe(l_om_ac)=0.02_r8
!             scavbe(l_dst_a1)=0.02_r8
    scavbe(l_dst_a2)=0.02_r8
    scavbe(l_dst_a3)=0.2_r8
    scavbe(l_ss_a1)=0.02_r8
    scavbe(l_ss_a2)=0.02_r8
    scavbe(l_ss_a3)=0.5_r8
	    
    do m=1,pcnst
       scavbecon(m)=scavbe(m)/10._r8
    end do
    scavbecon(l_dst_a3)=scavbe(l_dst_a3)
    scavbecon(l_ss_a3)= scavbe(l_ss_a3)
             
    do m=1,pcnst
       scavconv(m)=scavin(m)
    end do
    scavconv(l_so4_a2)=1._r8

    effsize(l_so4_n)=0.0118e-6_r8
    effsize(l_so4_na)=0.04e-6_r8
    effsize(l_so4_a1)=0.04e-6_r8
    effsize(l_so4_a2)=0.1e-6_r8
    effsize(l_so4_ac)=0.1e-6_r8
    effsize(l_so4_pr)=0.075e-6_r8

    effsize(l_bc_n)=0.0118e-6_r8
    effsize(l_bc_ax)=0.1e-6_r8
    effsize(l_bc_ni)=0.04e-6_r8
    effsize(l_bc_a) =0.04e-6_r8
    effsize(l_bc_ai)=0.04e-6_r8
    effsize(l_bc_ac)=0.1e-6_r8

!    effsize(l_om_n)=0.04e-6_r8
    effsize(l_om_ni)=0.04e-6_r8
    effsize(l_om_ai)=0.04e-6_r8
    effsize(l_om_ac)=0.1e-6_r8
    effsize(l_dst_a2)=0.22e-6_r8

    effsize(l_dst_a3)=0.63e-6_r8
    effsize(l_ss_a1)=0.022e-6_r8
    effsize(l_ss_a2)=0.13e-6_r8
    effsize(l_ss_a3)=0.74e-6_r8

!25jan06             sgpart(l_so4_n)=1.59    _r8
    sgpart(l_so4_n)=1.8_r8
    sgpart(l_so4_na)=1.8_r8
    sgpart(l_so4_a1)=1.8_r8
    sgpart(l_so4_a2)=1.59_r8
    sgpart(l_so4_ac)=1.59_r8
    sgpart(l_so4_pr)=1.59_r8

!
             sgpart(l_bc_n)=1.8_r8
    sgpart(l_bc_ax)=1.6_r8
    sgpart(l_bc_ni)=1.8_r8
    sgpart(l_bc_a) =1.8_r8
    sgpart(l_bc_ai)=1.8_r8
    sgpart(l_bc_ac)=1.59_r8
!
!    sgpart(l_om_n)=1.8_r8
    sgpart(l_om_ni)=1.8_r8
    sgpart(l_om_ai)=1.8_r8
    sgpart(l_om_ac)=1.59_r8

    sgpart(l_dst_a2)=1.59_r8
    sgpart(l_dst_a3)=2.00_r8
    sgpart(l_ss_a1)=1.59_r8
    sgpart(l_ss_a2)=1.59_r8
    sgpart(l_ss_a3)=2.00_r8


    rhopart(l_so4_n)=1841._r8
    rhopart(l_so4_na)=1841._r8
    rhopart(l_so4_a1)=1841._r8
    rhopart(l_so4_a2)=1769._r8
    rhopart(l_so4_ac)=1769._r8
    rhopart(l_so4_pr)=1841._r8
    rhopart(l_bc_n)=2000._r8
    rhopart(l_bc_ax)=507._r8
    rhopart(l_bc_ni)=2000._r8
    rhopart(l_bc_a) =2000._r8
    rhopart(l_bc_ai)=2000._r8
    rhopart(l_bc_ac)=2000._r8
!    rhopart(l_om_n)=1500._r8
    rhopart(l_om_ni)=1500._r8
    rhopart(l_om_ai)=1500._r8
    rhopart(l_om_ac)=1500._r8

    rhopart(l_dst_a2)=2600._r8
    rhopart(l_dst_a3)=2600._r8
    rhopart(l_ss_a1)=2200._r8
    rhopart(l_ss_a2)=2200._r8
    rhopart(l_ss_a3)=2200._r8

    predvel(l_dms,:)=0._r8

!    predvel(l_msa,1)=0.003_r8
!    predvel(l_msa,2)=0.007_r8
!    predvel(l_msa,3)=0.004_r8
    predvel(l_so2,1)=0.002_r8
    predvel(l_so2,2)=0.007_r8
    predvel(l_so2,3)=0.003_r8

!    predvel(l_h2so4,:)=0._r8
    predvel(l_qh2o2,:)=0._r8
    predvel(l_so4_n,:)=0.003_r8
    predvel(l_so4_na,:)=0.001_r8
    predvel(l_so4_a1,:)=0.001_r8
    predvel(l_so4_a2,:)=0.001_r8
    predvel(l_so4_ac,:)=0.001_r8
    predvel(l_so4_pr,:)=0.001_r8
    predvel(l_bc_n,:)=0.003_r8
    predvel(l_bc_ax,:)=0.001_r8
    predvel(l_bc_ni,:)=0.001_r8
    predvel(l_bc_a,:) =0.001_r8
    predvel(l_bc_ai,:)=0.001_r8
    predvel(l_bc_ac,:)=0.001_r8

!    predvel(l_om_n,:)=0.001_r8
    predvel(l_om_ni,:)=0.001_r8
    predvel(l_om_ai,:)=0.001_r8
    predvel(l_om_ac,:)=0.001_r8
    predvel(l_dst_a2,:)=0.0005_r8
    predvel(l_dst_a3,:)=0.01_r8
    predvel(l_ss_a1,:)=0.001_r8
    predvel(l_ss_a2,:)=0.001_r8
    predvel(l_ss_a3,:)=0.04_r8
                         
    call inittabrh


! os Coded with compilator flag (aeroslo)
!! tell camsrfexch_types to allocate fv & ram1 -- needed for dry depositon of particles
!       call hub2atm_setopts(aero_dust_in=.true.)


    return
  end subroutine aer_phys_init


  subroutine inittabrh
  ! Tables for hygroscopic growth
	
    integer :: i





    real(r8) :: rr0ss(10),rr0so4(10),rr0bcoc(10)

    data rr0ss / 1.00_r8, 1.00_r8, 1.02_r8, 1.57_r8, 1.88_r8, 1.97_r8, 2.12_r8, 2.35_r8, 2.88_r8, 3.62_r8 /
    data rr0so4 / 1.00_r8, 1.34_r8, 1.39_r8, 1.52_r8, 1.62_r8, 1.69_r8, 1.78_r8, 1.92_r8, 2.22_r8, 2.79_r8 /    
    data rr0bcoc / 1.00_r8, 1.02_r8, 1.03_r8, 1.12_r8, 1.17_r8, 1.20_r8, 1.25_r8, 1.31_r8, 1.46_r8, 1.71_r8 /

    rdivr0(:,:)=1._r8

    do i=1,10
       rdivr0(i,l_so4_n)=rr0so4(i)
       rdivr0(i,l_so4_na)=rr0so4(i)
       rdivr0(i,l_so4_a1)=rr0so4(i)
       rdivr0(i,l_so4_a2)=rr0so4(i)
       rdivr0(i,l_so4_ac)=rr0so4(i)
       rdivr0(i,l_so4_pr)=rr0so4(i)

       rdivr0(i,l_bc_a)=rr0so4(i)

!      rdivr0(i,l_bc_n)=rr0bcoc(i)
       rdivr0(i,l_bc_ni)=rr0bcoc(i)
       rdivr0(i,l_bc_ai)=rr0bcoc(i)
       rdivr0(i,l_bc_ac)=rr0bcoc(i)

!       rdivr0(i,l_om_n)=rr0bcoc(i)
       rdivr0(i,l_om_ni)=rr0bcoc(i)
       rdivr0(i,l_om_ai)=rr0bcoc(i)
       rdivr0(i,l_om_ac)=rr0bcoc(i)

       rdivr0(i,l_ss_a1)=rr0ss(i)
       rdivr0(i,l_ss_a2)=rr0ss(i)
       rdivr0(i,l_ss_a3)=rr0ss(i)	
    end do
    return
  end subroutine inittabrh


!===============================================================================
  subroutine aero_init
!-----------------------------------------------------------------------
! 
! Purpose: Set flux fields
! 
! 
!
! 
!-----------------------------------------------------------------------
!    use history,    only: addfld, add_default, phys_decomp

!----------------------------------------------------------------------- 
!
! Purpose: declare history variables, initialize data sets
! No initalisation in this version. Restart fields not read in here.
!
!-----------------------------------------------------------------------

    use tracers_suite,   only: init_tr, get_tracer_name
    use cam_history,     only: addfld, add_default, phys_decomp
    use ppgrid,          only: pver
    use constituents,    only: cnst_get_ind, cnst_name, cnst_longname, & 
      sflxnam, wetflxnam,numbnam
   use spmd_utils,   only: masterproc
    use drydep_so2,    only : dvel_so2_inti
    character*90  dvel_lnd
    character*90  clim_soilw, season_wes,depvel_file
!   character(len=16)::   drydep_list(:) 
!   local variables

! Local
    integer m, mm
!   character(len=8) :: name   ! constituent name

!   if ( tracers_flag ) then     
     
    dvel_lnd='inputdata/atm/cam/chem/trop_mozart/dvel/regrid_vegetation.nc'
!     dvel_lnd='inputdata/atm/cam/chem/'
     clim_soilw='inputdata/atm/cam/chem/trop_mozart/dvel/clim_soilw.nc'
     season_wes='inputdata/atm/cam/chem/trop_mozart/dvel/season_wes.nc'


    do m = 1,ncui
       mm = ixac-1+m    
!         name = get_tracer_name(m)
!         call cnst_get_ind(name, mm)
       call addfld (cnst_name(mm), 'kg/kg   ', pver, 'A', cnst_longname(mm), phys_decomp)
#ifdef SHORTRUN
       call addfld (sflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
       call addfld (wetflxnam(mm),   'kg/m2/s ',    1, 'A', trim(cnst_name(mm))//' surface flux', phys_decomp)
       call addfld (numbnam(mm),   'm-3 ',    pver, 'A', trim(cnst_name(mm))//' Aerosol number concentration', phys_decomp)

       call addfld (trim(cnst_name(mm))//'DRY','kg/m2/s ',1, 'A','trim(cnst_name(mm)) dry deposition',phys_decomp)

       call addfld (trim(cnst_name(mm))//'DV','kg/m2/s ',1, 'A','trim(cnst_name(mm)) dry deposition',phys_decomp)

#endif
       call add_default (cnst_name(mm), 1, ' ')
#ifdef SHORTRUN
       call add_default (sflxnam(mm),   1, ' ')
       call add_default (wetflxnam(mm),   1, ' ')
       call add_default (numbnam(mm),   1, ' ')
       call add_default (trim(cnst_name(mm))//'DRY', 1, ' ')

       call add_default (trim(cnst_name(mm))//'DV', 1, ' ')
#endif
    end do
    call addfld ('S4AQ ','kg S/m2 s-1  ',1, 'A','Sulphate aq.phase prod' ,phys_decomp)
    call addfld ('S4GA ','kg S/m2 s-1  ',1, 'A','Sulphate gas-phase prod' ,phys_decomp)    
    call addfld ('S2GA ','kg S/m2 s-1  ',1, 'A','SO2 gas-phase prod' ,phys_decomp)   
    call addfld ('MSAGA ','kg S/m2 s-1  ',1, 'A','MSA gas-phase prod' ,phys_decomp)  

    call add_default ('S4AQ', 1, ' ')
    call add_default ('S4GA', 1, ' ')
    call add_default ('S2GA', 1, ' ')
    call add_default ('MSAGA', 1, ' ')


    call dvel_so2_inti(dvel_lnd, clim_soilw, season_wes)
    return
  end subroutine aero_init

!=============================================================================




  end module aerosoldef


