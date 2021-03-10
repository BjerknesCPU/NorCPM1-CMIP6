module m_get_micom_nrens

contains

  integer function get_micom_nrens(nx, ny)
#if defined (QMPI)
    use qmpi , only : stop_mpi, master
#else
    use qmpi_fake
#endif
    implicit none

    integer, intent(in) :: nx, ny

    integer imem 
    logical ex
    character(len=3) :: cmem

    imem = 1
    ex = .true.
    do while (ex)
       write(cmem, '(i3.3)') imem
       inquire(exist = ex, file = 'forecast' // cmem // '.nc')
       if (ex) then
          imem = imem + 1
       end if
    end do
    get_micom_nrens = imem - 1
  end function get_micom_nrens

end module m_get_micom_nrens
