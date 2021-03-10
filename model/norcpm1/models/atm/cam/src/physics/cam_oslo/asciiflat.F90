!$Id: asciiflat.F90,v 1.1 2002/10/20 22:27:02 rce Exp $
!$Log: asciiflat.F90,v $
!Revision 1.1  2002/10/20 22:27:02  rce
!Initial revision
!
!-----------------------------------------------------------------------
	subroutine read_asciiflat_xzyt( filename, field, &
		kmaxd, mmaxd, k_must_match, rfv19 )
!
!   reads multiple times (months or seasons) from an "ascii-flat"
!   emissions file into a real*8 array
!
        use shr_kind_mod, only: r8 => shr_kind_r8
	use pmgrid, only: plon, plat
	use units
        use abortutils, only: endrun
        use cam_logfile,  only: iulog
	implicit none

! arguments
	character*(*), intent(in) :: filename  ! name of asciiflat file
	integer, intent(in) :: kmaxd  ! vertical dimension for field
	integer, intent(in) :: mmaxd  ! time dimension for field
	logical, intent(in) :: k_must_match
			!  if true, asciiflat file must have kmaxd levels
			!  if false, asciiflat file can have fewer levels
        logical, intent(in) :: rfv19
                        !  if true file has dimension fv19lon,fv19lat 
                        ! if false file must have dimension plon,plat
	real(r8), intent(out) :: field(plon,kmaxd,plat,mmaxd)

! local variables
        integer :: iversion	
!   iversion = indicates file format version -- should be 2
        integer :: itot, jtot, ktot
!   itot, jtot, ktot = number of longitude, latitude, vertical grids
        integer :: ntextlines
!   ntextlines = number of informational text lines
        integer :: ixcordinfo, iycordinfo, izcordinfo
!   if positive, file contains grid information
        integer :: lun
        integer :: i, j, k, m

	real(r8):: xcord1, dxcord, ycord1, dycord, zcord1, dzcord
!   coordinate info (origin and increment) for x, y, z

        character*80 textline

   integer,parameter :: fv19lon=144
   integer,parameter :: fv19lat=96
   real(r8) :: fieldfv19(fv19lon,kmaxd,fv19lat)
! get logical unit and open file
	lun = getunit()
        open( lun, file=filename, form='formatted', status='old', err=8100 )

9000    format(a)
9010    format(7i10)
9020    format(8e14.6)


! outermost loop -- read mmaxd times from the file
m_loop:	do m = 1, mmaxd

        read(lun,9010,err=8200) iversion
        if (iversion .gt. 3) then
            write(iulog,9400) trim(filename), 'iversion', 2, iversion
            call endrun()
        end if

        read(lun,9010) itot, jtot, ktot, ntextlines,  &
		ixcordinfo, iycordinfo, izcordinfo
        if (ixcordinfo .gt. 0) read(lun,9020,err=8300) xcord1, dxcord
        if (iycordinfo .gt. 0) read(lun,9020,err=8300) ycord1, dycord
        if (izcordinfo .gt. 0) read(lun,9020,err=8300) zcord1, dzcord

        if (rfv19) then

! check dimensions
        if (itot .ne. fv19lon) then
            write(iulog,9400) trim(filename), 'itot', plon, itot
            call endrun()

        else if (jtot .ne. fv19lat) then
            write(iulog,9400) trim(filename), 'jtot', plat, jtot
            call endrun()

        else if ((k_must_match .eqv. .true.) .and. (ktot .ne. kmaxd)) then
            write(iulog,9400) trim(filename), 'ktot', kmaxd, ktot
	    call endrun()

        else if ((k_must_match .eqv. .false.) .and. (ktot .gt. kmaxd)) then
            write(iulog,9400) trim(filename), 'ktot', kmaxd, ktot
	    call endrun()

	end if

        else


! check dimensions
        if (itot .ne. plon) then
            write(iulog,9400) trim(filename), 'itot', plon, itot
            call endrun()

        else if (jtot .ne. plat) then
            write(iulog,9400) trim(filename), 'jtot', plat, jtot
            call endrun()

        else if ((k_must_match .eqv. .true.) .and. (ktot .ne. kmaxd)) then
            write(iulog,9400) trim(filename), 'ktot', kmaxd, ktot
	    call endrun()

        else if ((k_must_match .eqv. .false.) .and. (ktot .gt. kmaxd)) then
            write(iulog,9400) trim(filename), 'ktot', kmaxd, ktot
	    call endrun()

	end if

        end if ! rfv19

        do i = 1, ntextlines
            read(lun,9000,err=8300) textline
        end do


! asciiflat file has "surface" level first, "top of atmosphere" last
! reverse this order when reading into field array
        if(rfv19) then

	read(lun,9020) (((fieldfv19(i,k,j), i=1,itot), j=1,jtot),  &
		k=kmaxd,kmaxd+1-ktot,-1)

       do j=1,plat
          do k=kmaxd,kmaxd+1-ktot,-1
             do i=1,plon
                 field(i,k,j,m)=  &
                  fieldfv19(nint(real(i,r8)/2._r8),k,nint(real(j,r8)/2._r8))
             end do
          end do
      end do
      else

 
	read(lun,9020) (((field(i,k,j,m), i=1,itot), j=1,jtot),  &
		k=kmaxd,kmaxd+1-ktot,-1)
      end if
        end do m_loop


! done reading -- close file and free unit
	close( lun )
	call freeunit( lun )
	return


!   open error
8100	write(iulog,9100) trim(filename)
	call endrun()

!   read errors
8200	write(iulog,9200) trim(filename)
	call endrun()

8300	write(iulog,9200) trim(filename)
	call endrun()

9100	format( / 'read_asciiflat_xzyt open error for:' / 4x, a )
9200	format( / 'read_asciiflat_xzyt type-1 read error for:' / 4x, a )
9300	format( / 'read_asciiflat_xzyt type-2 read error for:' / 4x, a )
9400	format( / 'read_asciiflat_xzyt mismatcherror for:' / 4x, a /  &
		4x, a, ' expected, found = ', 2i10 )

	end subroutine read_asciiflat_xzyt

