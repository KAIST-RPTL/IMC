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
    if(icore==score) then
        call MPI_COMM_SPAWN(comm_child, arg, 1, c_info, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
        print *, 'WRLD', icore, MPI_COMM_WORLD
        call COMM_WITH_CHILD
        call END_CHILD(intercomm)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    print *, 'PROCESS TERM', icore
end subroutine

subroutine COMM_WITH_CHILD
integer :: iz, idx, i, j
    call MPI_SEND(nth, 3, MPI_INTEGER, 0, 0, intercomm, ierr)
    call MPI_SEND(fuel_th, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(cld_th, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(gap_th, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(power_th, nth(1)*nth(2)*nth(3), MPI_REAL8, 0, 0, intercomm, ierr)
    !call MPI_SEND(t_fuel, nth(1)*nth(2)*nth(3), MPI_REAL8, 0,0,intercomm,ierr)
    call MPI_SEND(t_in, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(p_in, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_SEND(mflux, 1, MPI_REAL8, 0, 0, intercomm, ierr)
    call MPI_RECV(t_comm_cool, (nth(1)+1)*(nth(2)+1)*(nth(3)+1), MPI_REAL8, 0, 0, intercomm, MPI_STATUS_IGNORE, ierr)
    call MPI_RECV(t_comm_fuel, nth(1)*nth(2)*nth(3), MPI_REAL8, 0,0,intercomm, MPI_STATUS_IGNORE, ierr)
    call MPI_RECV(rho_comm_cool, (nth(1)+1)*(nth(2)+1)*(nth(3)+1), MPI_REAL8, 0, 0, intercomm, MPI_STATUS_IGNORE, ierr)
    print *, 'IMC: RECV well!'
end subroutine

subroutine END_CHILD(comm)
    integer :: comm
    call MPI_COMM_DISCONNECT(comm, ierr)
end subroutine

subroutine TH_ASSIGN_GRID
    use constants, only : K_B
    integer :: ix, iy, iz, idx
    real(8) :: th_w
!    t_fuel = 900d0
    
    if(th_iter>0)then
        th_w = 0.6d0
    else
        th_w = 1d0
    endif
    t_clad = 300d0 * K_B ! TODO
    
    do iz = 1, nth(3)
        do iy = 1, nth(2)
            do ix = 1, nth(1)
                idx = ix + (iy-1) * (nth(1))  
                t_fuel(ix, iy, iz) = t_comm_fuel(iz, idx) * th_w * K_B& 
                    + t_fuel(ix, iy, iz) * (1d0-th_w)
            enddo
        enddo
        if(icore==score) then
            print *, 'ZVAL', iz, th_iter
            print *, 'FUELTEMP'
            do iy = 1, nth(2)
                write(*,'(<nth(1)+1>F10.3)') (t_fuel(ix, iy, iz)/K_B, ix = 1, nth(1))
            enddo
        endif
    enddo
    do iz = 1, nth(3)+1
        do iy = 1, nth(2)+1
            do ix = 1, nth(1)+1
                idx = ix + (iy-1) * (nth(1)+1)  
                t_bulk(ix, iy, iz) = t_comm_cool(iz, idx) * th_w * K_B&
                    + t_bulk(ix,iy,iz) * (1d0-th_w)
                rho_bulk(ix,iy,iz) = rho_comm_cool(iz, idx) * th_w &
                    + rho_bulk(ix, iy, iz) * (1d0-th_w)
            enddo
        enddo
        if(icore==score) then
            print *, 'BULKTEMP', iz
            do iy = 1, nth(2)+1
                write(*,'(<nth(1)+1>F10.3)') (t_bulk(ix, iy, iz)/K_B, ix = 1, nth(1)+1)
            enddo
        endif
!            print *, 'BULKRHO'
!            do iy = 1, nth(2)+1
!                write(*,'(<nth(1)+1>F10.3)') (rho_bulk(ix, iy, iz), ix = 1, nth(1)+1)
!            enddo
!        endif
    enddo
    
!    t_fuel = t_fuel * K_B
!    t_clad = t_clad * K_B
!    t_bulk = t_bulk * K_B
end subroutine

end module communication
