module bank_header 
	use mpi 
    use variables, only: latent
	
    implicit none
    type :: Bank
        real(8) :: wgt           ! weight of bank site
        real(8) :: xyz(1:3)        ! location of bank particle
        real(8) :: uvw(1:3)        ! diretional cosines
        real(8) :: E             ! energy for CE
        integer :: G             ! energy group if in MG mode.
        logical :: delayed
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
!        integer, allocatable :: delayedarr(:)
!        real(8), allocatable :: delayedlam(:), nlifearr(:)
        
        ! In case of latent is defined as 'parameter'
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr(1:latent)
      contains 
    end type Bank
    
    type :: PrecBank
        real(8) :: wgt          ! weight of bank site
        real(8) :: xyz(1:3)     ! location of bank particle
        real(8) :: E            ! Incident energy for CE
        integer :: G            ! Incident energy group if in MG mode.
		integer :: idx			! material or isotope index
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
!        integer, allocatable :: delayedarr(:)
!        real(8), allocatable :: delayedlam(:), nlifearr(:)
        
        ! In case of latent is defined as 'parameter'
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr(1:latent)
    end type PrecBank
	
	
	type :: VRCBank 
        real(8) :: wgt          ! weight of bank site
        real(8) :: xyz(1:3)     ! location of bank particle
        real(8) :: uvw(1:3)     ! location of bank particle
        real(8) :: E            ! Incident energy for CE
        integer :: G            ! Incident energy group if in MG mode.
	end type VRCBank
	
	
	
    ! Source and fission bank
    type(Bank), allocatable, target :: source_bank(:)
    type(Bank), allocatable, target :: fission_bank(:)
    type(Bank), allocatable         :: temp_bank(:)
    type(Bank)                      :: thread_bank(20000)
    !$OMP THREADPRIVATE(thread_bank)
    integer                         :: bank_idx
    !$OMP THREADPRIVATE(bank_idx)
    
	
	type(Bank), allocatable :: split_bank(:)
	type(Bank), allocatable :: split_bank_temp(:)
    type(Bank)              :: split_thread(15000)
    !$OMP THREADPRIVATE(split_thread)
    integer                 :: split_idx
    !$OMP THREADPRIVATE(split_idx)
	
	

	type(Bank), allocatable :: delayed_bank(:)
	type(Bank), allocatable :: prompt_bank(:)
	type(Bank), allocatable :: dynamic_bank(:)
    type(Bank)				:: thread_bank_init(15000)
    !$OMP THREADPRIVATE(thread_bank_init)
    integer                 :: init_idx
    !$OMP THREADPRIVATE(init_idx)
	
	
	
	type(PrecBank), allocatable :: prec_bank(:), prec_bank_temp(:), prec_bank_local(:)
    type(PrecBank)              :: prec_thread(15000)
    !$OMP THREADPRIVATE(prec_thread)
    integer                 :: prec_idx
    !$OMP THREADPRIVATE(prec_idx)
	
	
	type(Bank), allocatable :: vrc_bank(:), vrc_bank_temp(:)
    type(Bank)              :: vrc_thread(10000)
    !$OMP THREADPRIVATE(vrc_thread)
    integer                 :: vrc_idx
    !$OMP THREADPRIVATE(vrc_idx)
	
	
	
	integer :: MPI_bank, MPI_precbank!, MPI_vrcbank
	
	contains 
    subroutine MPI_banktype() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:8) :: blocklength, displacement, oldtype 
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 3
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 1
		blocklength(7) = 8
		blocklength(8) = 8
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex + intex
		displacement(6) = 8*realex + intex + logicex
		displacement(7) = 8*realex + intex + logicex + realex
		displacement(8) = 8*realex + intex + logicex + realex + 8*realex
		
		oldtype(0:3) = MPI_double_precision
		oldtype(4)   = MPI_INTEGER
		oldtype(5)   = MPI_logical
		oldtype(6:8) = MPI_double_precision
		
		call MPI_TYPE_STRUCT (9, blocklength, displacement, oldtype,MPI_bank, ierr)
		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	end subroutine 
	
	
    subroutine MPI_precbanktype() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:7) :: blocklength, displacement, oldtype 
		
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 1
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 8
		blocklength(7) = 8
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr) 
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr) 
		displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
		displacement(3) = 5*realex
		displacement(4) = 5*realex + intex 
		displacement(5) = 5*realex + 2*intex
		displacement(6) = 5*realex + 2*intex + realex
		displacement(7) = 5*realex + 2*intex + 9*realex
		
		oldtype(0:2) = MPI_double_precision
		oldtype(3:4)   = MPI_INTEGER
		oldtype(5:7)   = MPI_double_precision
				
		call MPI_TYPE_STRUCT (8, blocklength, displacement, oldtype,MPI_precbank, ierr)
		call MPI_TYPE_COMMIT (MPI_precbank, ierr)
	end subroutine 	

    subroutine MPI_banktype_ifp() 
		integer :: ierr 
		integer :: realex, intex, logicex
		!integer, dimension(0:8) :: blocklength, displacement, oldtype 
        integer, dimension(0:11) :: blocklength, displacement, oldtype
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 3
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 1
		blocklength(7) = 8
		blocklength(8) = 8
        ! ADJOINT MPI_BANK
        blocklength(9) = latent
        blocklength(10)= latent
        blocklength(11)= latent
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex + intex
		displacement(6) = 8*realex + intex + logicex
		displacement(7) = 8*realex + intex + logicex + realex
		displacement(8) = 8*realex + intex + logicex + realex + 8*realex
        displacement(9) = 8*realex + intex + logicex + realex + 8*realex + 8*realex
        displacement(10)= 8*realex + intex + logicex + realex + 8*realex + 8*realex + latent*intex
        displacement(11)= 8*realex + intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex
		
		oldtype(0:3) = MPI_double_precision
		oldtype(4)   = MPI_INTEGER
		oldtype(5)   = MPI_logical
		oldtype(6:8) = MPI_double_precision
        oldtype(9)   = MPI_INTEGER
        oldtype(10)  = MPI_double_precision
        oldtype(11)  = MPI_double_precision
		
		call MPI_TYPE_STRUCT (12, blocklength, displacement, oldtype,MPI_bank, ierr)
		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	end subroutine 
	
	
    subroutine MPI_precbanktype_ifp() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:10) :: blocklength, displacement, oldtype 
		
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 1
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 8
		blocklength(7) = 8
        ! ADJOINT MPI_PRECBANK
        blocklength(8) = latent
        blocklength(9) = latent
        blocklength(10)= latent
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr) 
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr) 
		displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
		displacement(3) = 5*realex
		displacement(4) = 5*realex + intex 
		displacement(5) = 5*realex + 2*intex
		displacement(6) = 5*realex + 2*intex + realex
		displacement(7) = 5*realex + 2*intex + 9*realex
        displacement(8) = 5*realex + 2*intex + 9*realex + 8 * realex
        displacement(9) = 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex
        displacement(10)= 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex + latent*realex
		
		oldtype(0:2) = MPI_double_precision
		oldtype(3:4)   = MPI_INTEGER
		oldtype(5:7)   = MPI_double_precision
        oldtype(8)     = MPI_INTEGER
        oldtype(9)     = MPI_double_precision
        oldtype(10)    = MPI_double_precision
				
		call MPI_TYPE_STRUCT (11, blocklength, displacement, oldtype,MPI_precbank, ierr)
		call MPI_TYPE_COMMIT (MPI_precbank, ierr)
	end subroutine 	
	
end module
