module communication
use VARIABLES
use mpi
implicit none
logical :: do_child = .false.   ! Existence of child process

! MPI parameters
!integer :: max_proc     = 4 ! Maximum processes for Child
integer :: intercomm        ! Intercommunication
integer :: c_icore          ! MPI rank for child
integer :: c_ncore          ! MPI_number of cores of child
integer :: c_info           ! MPI info for child; input directory
character(50) :: comm_child ! Command to run child

integer, dimension(1) :: errcodes
integer :: snd, rcv

contains

! =========================================================
! This subroutine initiates CHILD process
! Currently (2023/03/30) working on STARTH...
subroutine INIT_CHILD
character(len=50), dimension(2) :: arg
    write(*,*) 'BEGIN'
    !arg(1) = trim(comm_child)
    arg(1) = trim(directory)//'/START'
    arg(2) = " "
    c_info = MPI_INFO_NULL
    call MPI_COMM_SPAWN(comm_child, arg, 1, c_info, 0, MPI_COMM_SELF, intercomm, errcodes, ierr)
    call SND_RCV_CHILD
    write(*,*) 'DONE', intercomm
    call MPI_COMM_RANK(intercomm, c_icore, ierr)
    call MPI_COMM_SIZE(intercomm, c_ncore, ierr)
    write(*,*) c_icore, c_ncore

end subroutine

subroutine SND_RCV_CHILD
    snd = 12345
    write(*,*) 'SENDING...', snd, intercomm, MPI_COMM_SELF
    call MPI_SEND(snd, 1, MPI_INTEGER, 0, 0, intercomm, ierr)
    !call MPI_COMM_JOIN(intercomm, ierr)
    call MPI_RECV(rcv, 1, MPI_INTEGER, 0, 0, intercomm, MPI_STATUS_IGNORE, ierr)

    write(*,*) 'RECEIVED', snd, rcv
end subroutine


end module communication
