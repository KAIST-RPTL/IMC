module communication
use VARIABLES
use TH_HEADER
use mpi
implicit none

! MPI parameters
!integer :: max_proc     = 4 ! Maximum processes for Child
integer :: child_process    ! CHILD PROCESS
integer :: intercomm        ! Intercommunication
integer :: c_icore          ! MPI rank for child
integer :: c_ncore          ! MPI_number of cores of child
integer :: c_info           ! MPI info for child; input directory
character(50) :: comm_child ! Command to run child

integer, dimension(1) :: errcodes
real(8) :: snd, rcv

contains

! =========================================================
! This subroutine initiates CHILD process
! Currently (2023/03/30) working on STARTH...
subroutine INIT_CHILD
character(len=50), dimension(2) :: arg
    !arg(1) = trim(comm_child)
    arg(1) = trim(directory)//'/START'
    arg(2) = " "
    c_info = MPI_INFO_NULL
    call MPI_COMM_SPAWN(comm_child, arg, 1, c_info, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
    call COMM_WITH_CHILD
end subroutine

subroutine COMM_WITH_CHILD
integer :: iz, idx
    call MPI_SEND(nth, 3, MPI_INTEGER, 0, 0, intercomm, ierr)
    call MPI_SEND(power_th, nth(1)*nth(2)*nth(3), MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(t_in, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(p_in, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(mflux, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_RECV(t_comm_cool, (nth(1)+1)*(nth(2)+1)*(nth(3)+1), MPI_REAL8, 0, 0, intercomm, MPI_STATUS_IGNORE, ierr)
    call MPI_RECV(rho_comm_cool, (nth(1)+1)*(nth(2)+1)*(nth(3)+1), MPI_REAL8, 0, 0, intercomm, MPI_STATUS_IGNORE, ierr)
end subroutine

subroutine END_CHILD(comm)
    integer :: comm
    call MPI_COMM_FREE(comm, ierr)
end subroutine

subroutine TH_ASSIGN_GRID
    use constants, only : K_B
    integer :: ix, iy, iz, idx
    t_fuel = 0d0
    
    t_clad = 0d0 ! TODO
    
    do iz = 1, nth(3)+1
        do iy = 1, nth(2)+1
            do ix = 1, nth(1)+1
                idx = ix + (iy-1) * (nth(1)+1)  
                t_bulk(ix, iy, iz) = t_comm_cool(iz, idx)
                rho_bulk(ix,iy,iz) = rho_comm_cool(iz, idx)
            enddo
        enddo
        if(icore==score) then
            print *, 'ZVAL', iz, th_iter
            print *, 'BULKTEMP'
            do iy = 1, nth(2)+1
                write(*,'(<nth(1)+1>F10.3)') (t_bulk(ix, iy, iz), ix = 1, nth(1)+1)
            enddo
            print *, 'BULKRHO'
            do iy = 1, nth(2)+1
                write(*,'(<nth(1)+1>F10.3)') (rho_bulk(ix, iy, iz), ix = 1, nth(1)+1)
            enddo
        endif

        t_bulk = t_bulk * K_B

    enddo
end subroutine

end module communication
