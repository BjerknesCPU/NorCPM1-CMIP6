module distribute

#if defined(QMPI)
  use qmpi
#else
  use qmpi_fake
#endif

  !
  ! public stuff
  !
  integer, public :: my_number_of_iterations, my_first_iteration, my_last_iteration
  integer, dimension(:), allocatable, public :: number_of_iterations, first_iteration, last_iteration
  integer, dimension(:), allocatable, public :: randommap

contains
  subroutine distribute_iterations(nz)
    implicit none

    integer, intent(in) :: nz

    integer :: i, j
    real(8) :: num_procs_real, mean_iterations

    if (.not. allocated(number_of_iterations)) then
       allocate(number_of_iterations(qmpi_num_proc))
    end if
    if (.not. allocated(first_iteration)) then
       allocate(first_iteration(qmpi_num_proc))
    end if
    if (.not. allocated(last_iteration)) then
       allocate(last_iteration(qmpi_num_proc))
    end if

!    if (master) then
!       print *, 'Distribution of iterations:'
!    end if

    num_procs_real = qmpi_num_proc
    mean_iterations = nz / num_procs_real

    j = -1
    if (int(mean_iterations) .eq. mean_iterations) then
       my_number_of_iterations = nz/qmpi_num_proc
       if (master) then
          number_of_iterations(:) = nz / qmpi_num_proc
!          print *, 'All procs get ', number_of_iterations(1), 'iterations'
       endif
       j = qmpi_num_proc
    else
       do i = 1, qmpi_num_proc
          if (i * floor(mean_iterations) +&
               (qmpi_num_proc-i) * ceiling(mean_iterations) .eq. nz) then
             j = i
             exit
          endif
       end do

       if (qmpi_proc_num + 1 .le. j) then
          my_number_of_iterations = floor(mean_iterations)
       else
          my_number_of_iterations = ceiling(mean_iterations)
       endif

       if (master) then
          number_of_iterations(1:j) = floor(mean_iterations)
          number_of_iterations(j+1:qmpi_num_proc) = ceiling(mean_iterations)
          if ((j * floor(mean_iterations) +&
               (qmpi_num_proc - j) * ceiling(mean_iterations)) .ne. nz) then
             print *, 'ERROR in distribute_iteration()'
             stop
          endif
       endif
    endif

    if (master) then
       first_iteration(1) = 1; 
       last_iteration(1) = number_of_iterations(1)
       do i = 2, qmpi_num_proc
          first_iteration(i) = last_iteration(i - 1) + 1 
          last_iteration(i) = first_iteration(i) + number_of_iterations(i) - 1
       end do
    endif

    if (qmpi_proc_num + 1 .le. j) then
       my_first_iteration = qmpi_proc_num*my_number_of_iterations + 1
    else
       my_first_iteration = j * (my_number_of_iterations - 1) +&
            (qmpi_proc_num - j) * my_number_of_iterations + 1
    endif
    my_last_iteration = my_first_iteration + my_number_of_iterations - 1

    !print *, 'I am', qmpi_proc_num, ', my_first_ind =', my_first_iteration,&
    !     ', my_last_ind =', my_last_iteration
  end subroutine distribute_iterations

  subroutine distribute_iterations_field(nz, names, levels)
    implicit none

    integer         , intent(in) :: nz
    character(len=8), dimension(nz), intent(in):: names
    integer         , dimension(nz), intent(in):: levels 

    integer :: i, j1, j2, j3, index_dp, dplen, num_procs_int
    real(8) :: num_procs_real
    real(8) :: mean_iterations1, mean_iterations2, mean_iterations3

    if (.not. allocated(number_of_iterations)) then
       allocate(number_of_iterations(qmpi_num_proc))
    end if
    if (.not. allocated(first_iteration)) then
       allocate(first_iteration(qmpi_num_proc))
    end if
    if (.not. allocated(last_iteration)) then
       allocate(last_iteration(qmpi_num_proc))
    end if

    ! Find the index of field 'dp' for the level 1
    dplen = 0
    do i=1,nz
       if (trim(names(i)) .eq. 'dp') dplen = dplen + 1
    end do
    if (dplen .eq. 0) then
       ! no dp as varibale, just call standard function
       call distribute_iterations(nz)
       return
    else if (qmpi_num_proc .le. 2) then
       print *, 'ERROR in distribute_iteration_field():'
       print *, 'qmpi_num_proc should be greater than 2'
       print *, '  option 1: use more than 2 cpus'
       print *, '  option 2: do not update dp'
       stop
    end if
    do i=1,nz
       if ((trim(names(i)) .eq. 'dp')) then
          index_dp = i
          exit
       end if
    end do

!    if (master) then
!       print *, 'Distribution of iterations:'
!    end if

    ! All DP take one proc
    num_procs_real = (qmpi_num_proc - 1) * (index_dp - 1)&
         /(nz-dplen)
    num_procs_int = max(int(num_procs_real), 1)
    mean_iterations1 =  (index_dp - 1) / real(num_procs_int)
    mean_iterations2 = dplen
    mean_iterations3 =  (nz - dplen - (index_dp - 1))&
         /real(max(qmpi_num_proc - 1 - num_procs_int, 1))

    j1 = -1
    do i = 1, num_procs_int
       if (i * floor(mean_iterations1) + (num_procs_int - i) &
            * ceiling(mean_iterations1) .eq. (index_dp - 1)) then
          j1 = i
          exit
       endif
    end do

    j2 = num_procs_int + 1

    j3 = -1
    do i = j2+1,qmpi_num_proc 
       if ((i - j2) * floor(mean_iterations3) +&
            (qmpi_num_proc - i) * ceiling(mean_iterations3) .eq. &
            (nz - dplen - (index_dp - 1))) then
          j3 = i
          exit
       endif
    end do
    
    if (qmpi_proc_num + 1 .le. j1) then
       my_number_of_iterations = floor(mean_iterations1)
    else if (qmpi_proc_num + 1 .lt. j2) then
       my_number_of_iterations = ceiling(mean_iterations1)
    else if (qmpi_proc_num + 1 .eq. j2) then
       my_number_of_iterations = nint(mean_iterations2)
    else if (qmpi_proc_num + 1 .le. j3) then
       my_number_of_iterations = floor(mean_iterations3)
    else
       my_number_of_iterations = ceiling(mean_iterations3)
    endif

    if (master) then
       number_of_iterations(1:j1) = floor(mean_iterations1)
       number_of_iterations(j1+1:j2-1) = ceiling(mean_iterations1)
       number_of_iterations(j2) = nint(mean_iterations2)
       number_of_iterations(j2+1:j3) = floor(mean_iterations3)
       number_of_iterations(j3+1:qmpi_num_proc) = floor(mean_iterations3)
       if ((j1 * floor(mean_iterations1) + (j2 - 1 - j1) * &
            ceiling(mean_iterations1)) .ne. index_dp-1) then
          print *, j1, floor(mean_iterations1), ceiling(mean_iterations1), j2
          print *, 'ERROR in distribute_iteration_field for j1'
          stop
       endif
       if (((j3 - j2) * floor(mean_iterations3) + (qmpi_num_proc - j3) *&
            ceiling(mean_iterations3)) .ne. (nz - dplen - (index_dp - 1))) then
          print *, 'ERROR in distribute_iteration_field for j3'
          stop
       endif
       if ((j1 * floor(mean_iterations1) +&
            (j2 - 1 - j1) * ceiling(mean_iterations1) +&
            nint(mean_iterations2) + &
            (j3 - j2) * floor(mean_iterations3) +&
            (qmpi_num_proc - j3) * ceiling(mean_iterations3)) .ne. nz) then
          print *, 'ERROR in distribute_iteration_field for j2'
          stop
       endif
    endif

    if (master) then
       first_iteration(1) = 1; 
       last_iteration(1) = number_of_iterations(1)
       do i = 2, qmpi_num_proc
          first_iteration(i) = last_iteration(i - 1) + 1 
          last_iteration(i) = first_iteration(i) + number_of_iterations(i) - 1
       end do
    endif

    if (qmpi_proc_num + 1 .le. j1) then
       my_first_iteration = qmpi_proc_num*my_number_of_iterations + 1
    else if (qmpi_proc_num + 1 .lt. j2) then
       my_first_iteration = j1 * (my_number_of_iterations - 1) +&
            (qmpi_proc_num - j1) * my_number_of_iterations + 1
    else if (qmpi_proc_num + 1 .eq. j2) then 
       my_first_iteration = index_dp
    else if (qmpi_proc_num + 1 .le. j3) then
       my_first_iteration = index_dp + dplen + (qmpi_proc_num - j2) *&
            my_number_of_iterations
    else
       my_first_iteration = index_dp + dplen + &
            (j3 - j2) * (my_number_of_iterations - 1) +&
            (qmpi_proc_num - j3) * my_number_of_iterations
    endif
    my_last_iteration = my_first_iteration + my_number_of_iterations - 1

    !print *, 'I am', qmpi_proc_num, ', my_first_ind =', my_first_iteration,&
    !     ', my_last_ind =', my_last_iteration
  end subroutine distribute_iterations_field

end module distribute

