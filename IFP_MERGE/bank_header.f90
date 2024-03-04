module bank_header 
	use mpi 
    use variables, only: latent
	
    implicit none
	
    ! (TSOH-IFP)	type :: Bank
    ! (TSOH-IFP)	    real(8) :: wgt           ! weight of bank site
    ! (TSOH-IFP)	    real(8) :: xyz(1:3)        ! location of bank particle
    ! (TSOH-IFP)	    real(8) :: uvw(1:3)        ! diretional cosines
    ! (TSOH-IFP)	    real(8) :: E             ! energy for CE
    ! (TSOH-IFP)	    integer :: G             ! energy group if in MG mode.
    ! (TSOH-IFP)	    logical :: delayed
	! (TSOH-IFP)		real(8) :: time
	! (TSOH-IFP)		real(8) :: beta(1:8)
	! (TSOH-IFP)		real(8) :: lambda(1:8)        
    ! (TSOH-IFP)	    ! In case of latent is defined as 'parameter'
    ! (TSOH-IFP)	    integer :: delayedarr(1:latent)
    ! (TSOH-IFP)	    real(8) :: delayedlam(1:latent)
    ! (TSOH-IFP)	    real(8) :: nlifearr  (1:latent)
    ! (TSOH-IFP)	  contains 
    ! (TSOH-IFP)	end type Bank
	
	! TSOH-IFP
	TYPE :: BANK
        real(8) :: wgt           ! weight of bank site
        real(8) :: xyz(1:3)      ! location of bank particle
        real(8) :: uvw(1:3)      ! diretional cosines
        real(8) :: E             ! energy for CE
        integer :: G             ! energy group if in MG mode.
		INTEGER :: G_delayed     ! energy group for precursor in MG mode (when delayed == .TRUE.)
        logical :: delayed
		! ********************************** IFP RELATED ********************************************* CAVEAT) IF APPEND BELOW lambda(1:8) doesn't work. Reckon due to memory leakage (but why?)
		INTEGER :: i_Parent  ! Index for parent  from previous cycle
		INTEGER :: i_Current ! Index for current odrder of the stored source
		INTEGER :: dIMT      ! Stores cycle-wise IMT value (delayed group)
		REAL(8) :: dlam      ! Stores cycle-wise lam value (decay constant)
		REAL(8) :: nlife     ! Stores cycle-wise life time (travel time)
		! ********************************************************************************************
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
	END TYPE BANK
    
    type :: PrecBank
        real(8) :: wgt          ! weight of bank site
        real(8) :: xyz(1:3)     ! location of bank particle
        real(8) :: E            ! Incident energy for CE
        integer :: G            ! Incident energy group if in MG mode.
		integer :: idx			! material or isotope index
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
        
        ! (TSOH-IFP): ! In case of latent is defined as 'parameter' / THIS PART WILL BE MODIFIED AFTERWARDS
        ! (TSOH-IFP): integer :: delayedarr(1:latent)
        ! (TSOH-IFP): real(8) :: delayedlam(1:latent)
        ! (TSOH-IFP): real(8) :: nlifearr  (1:latent)
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
	! TSOH-IFP: STORES PREVIOUS CYCLE-WISE SOURCE BANKS
	TYPE(Bank), ALLOCATABLE :: mat_source_bank    (:,:)
	
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
    ! (TSOH-IFP)	subroutine MPI_banktype() 
	! (TSOH-IFP)		integer :: ierr 
	! (TSOH-IFP)		integer :: realex, intex, logicex
	! (TSOH-IFP)		integer, dimension(0:8) :: blocklength, displacement, oldtype 
	! (TSOH-IFP)		
	! (TSOH-IFP)		blocklength(0) = 1
	! (TSOH-IFP)		blocklength(1) = 3
	! (TSOH-IFP)		blocklength(2) = 3
	! (TSOH-IFP)		blocklength(3) = 1
	! (TSOH-IFP)		blocklength(4) = 1
	! (TSOH-IFP)		blocklength(5) = 1
	! (TSOH-IFP)		blocklength(6) = 1
	! (TSOH-IFP)		blocklength(7) = 8
	! (TSOH-IFP)		blocklength(8) = 8
	! (TSOH-IFP)		
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
	! (TSOH-IFP)		displacement(0) = 0
	! (TSOH-IFP)		displacement(1) =   realex
	! (TSOH-IFP)		displacement(2) = 4*realex
	! (TSOH-IFP)		displacement(3) = 7*realex
	! (TSOH-IFP)		displacement(4) = 8*realex
	! (TSOH-IFP)		displacement(5) = 8*realex + intex
	! (TSOH-IFP)		displacement(6) = 8*realex + intex + logicex
	! (TSOH-IFP)		displacement(7) = 8*realex + intex + logicex + realex
	! (TSOH-IFP)		displacement(8) = 8*realex + intex + logicex + realex + 8*realex
	! (TSOH-IFP)		
	! (TSOH-IFP)		oldtype(0:3) = MPI_double_precision
	! (TSOH-IFP)		oldtype(4)   = MPI_INTEGER
	! (TSOH-IFP)		oldtype(5)   = MPI_logical
	! (TSOH-IFP)		oldtype(6:8) = MPI_double_precision
	! (TSOH-IFP)		
	! (TSOH-IFP)		call MPI_TYPE_STRUCT (9, blocklength, displacement, oldtype,MPI_bank, ierr)
	! (TSOH-IFP)		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	! (TSOH-IFP)	end subroutine 
	
	! TSOH-IFP
    subroutine MPI_banktype() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:14) :: blocklength, displacement, oldtype 
		! IN THE SAME ORDER TO THAT OF TYPE
		blocklength(0) = 1       ! real(8)
		blocklength(1) = 3       ! real(8)
		blocklength(2) = 3       ! real(8)
		blocklength(3) = 1       ! real(8)
		blocklength(4) = 1       ! integer
		blocklength(5) = 1       ! INTEGER
		blocklength(6) = 1       ! logical
		! ************* IFP REALTED ****************
		blocklength(7) = 1  ! INTEGER
		blocklength(8) = 1  ! INTEGER
		blocklength(9) = 1  ! INTEGER
		blocklength(10) = 1 ! REAL(8)
		blocklength(11) = 1 ! REAL(8)
		! ******************************************
		blocklength(12) = 1 ! real(8)
		blocklength(13) = 8 ! real(8)
		blocklength(14) = 8 ! real(8)

		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		
		! DISPLACEMENT (메모리상 어디에 위치하는지)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex +   intex
		displacement(6) = 8*realex + 2*intex
		displacement(7) = 8*realex + 2*intex + logicex
		displacement(8) = 8*realex + 2*intex + logicex +   intex
		displacement(9) = 8*realex + 2*intex + logicex + 2*intex
		displacement(10)= 8*realex + 2*intex + logicex + 3*intex
		displacement(11)= 8*realex + 2*intex + logicex + 3*intex +   realex
		displacement(12)= 8*realex + 2*intex + logicex + 3*intex + 2*realex
		displacement(13)= 8*realex + 2*intex + logicex + 3*intex + 3*realex
		displacement(14)= 8*realex + 2*intex + logicex + 3*intex + 3*realex + 8*realex

		! IN THE SAME ORDER TO THAT OF TYPE
		oldtype(0:3) = MPI_double_precision
		oldtype(4)   = MPI_INTEGER
		oldtype(5)   = MPI_INTEGER
		oldtype(6)   = MPI_logical
		oldtype(7:9) = MPI_INTEGER
		oldtype(10)  = MPI_double_precision
		oldtype(11)  = MPI_double_precision
		oldtype(12)  = MPI_double_precision
		oldtype(13)  = MPI_double_precision
		oldtype(14)  = MPI_double_precision

		call MPI_TYPE_STRUCT (15, blocklength, displacement, oldtype, MPI_bank, ierr)
		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	end subroutine MPI_banktype
	
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
	end subroutine MPI_precbanktype 	

    ! (TSOH-IFP)	subroutine MPI_banktype_ifp() 
	! (TSOH-IFP)		integer :: ierr 
	! (TSOH-IFP)		integer :: realex, intex, logicex
    ! (TSOH-IFP)	    integer, dimension(0:11) :: blocklength, displacement, oldtype
	! (TSOH-IFP)		
	! (TSOH-IFP)		blocklength(0) = 1
	! (TSOH-IFP)		blocklength(1) = 3
	! (TSOH-IFP)		blocklength(2) = 3
	! (TSOH-IFP)		blocklength(3) = 1
	! (TSOH-IFP)		blocklength(4) = 1
	! (TSOH-IFP)		blocklength(5) = 1
	! (TSOH-IFP)		blocklength(6) = 1
	! (TSOH-IFP)		blocklength(7) = 8
	! (TSOH-IFP)		blocklength(8) = 8
    ! (TSOH-IFP)	    ! ADJOINT MPI_BANK
    ! (TSOH-IFP)	    blocklength(9) = latent
    ! (TSOH-IFP)	    blocklength(10)= latent
    ! (TSOH-IFP)	    blocklength(11)= latent
	! (TSOH-IFP)		
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
	! (TSOH-IFP)		displacement(0) = 0
	! (TSOH-IFP)		displacement(1) =   realex
	! (TSOH-IFP)		displacement(2) = 4*realex
	! (TSOH-IFP)		displacement(3) = 7*realex
	! (TSOH-IFP)		displacement(4) = 8*realex
	! (TSOH-IFP)		displacement(5) = 8*realex + intex
	! (TSOH-IFP)		displacement(6) = 8*realex + intex + logicex
	! (TSOH-IFP)		displacement(7) = 8*realex + intex + logicex + realex
	! (TSOH-IFP)		displacement(8) = 8*realex + intex + logicex + realex + 8*realex
    ! (TSOH-IFP)	    displacement(9) = 8*realex + intex + logicex + realex + 8*realex + 8*realex
    ! (TSOH-IFP)	    displacement(10)= 8*realex + intex + logicex + realex + 8*realex + 8*realex + latent*intex
    ! (TSOH-IFP)	    displacement(11)= 8*realex + intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex
	! (TSOH-IFP)		
	! (TSOH-IFP)		oldtype(0:3) = MPI_double_precision
	! (TSOH-IFP)		oldtype(4)   = MPI_INTEGER
	! (TSOH-IFP)		oldtype(5)   = MPI_logical
	! (TSOH-IFP)		oldtype(6:8) = MPI_double_precision
    ! (TSOH-IFP)	    oldtype(9)   = MPI_INTEGER
    ! (TSOH-IFP)	    oldtype(10)  = MPI_double_precision
    ! (TSOH-IFP)	    oldtype(11)  = MPI_double_precision
	! (TSOH-IFP)		
	! (TSOH-IFP)		call MPI_TYPE_STRUCT (12, blocklength, displacement, oldtype,MPI_bank, ierr)
	! (TSOH-IFP)		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	! (TSOH-IFP)	end subroutine 
	
    ! (TSOH-IFP)	subroutine MPI_precbanktype_ifp() 
	! (TSOH-IFP)		integer :: ierr 
	! (TSOH-IFP)		integer :: realex, intex, logicex
	! (TSOH-IFP)		integer, dimension(0:10) :: blocklength, displacement, oldtype 
	! (TSOH-IFP)		
	! (TSOH-IFP)		blocklength(0) = 1
	! (TSOH-IFP)		blocklength(1) = 3
	! (TSOH-IFP)		blocklength(2) = 1
	! (TSOH-IFP)		blocklength(3) = 1
	! (TSOH-IFP)		blocklength(4) = 1
	! (TSOH-IFP)		blocklength(5) = 1
	! (TSOH-IFP)		blocklength(6) = 8
	! (TSOH-IFP)		blocklength(7) = 8
    ! (TSOH-IFP)	    ! ADJOINT MPI_PRECBANK
    ! (TSOH-IFP)	    blocklength(8) = latent
    ! (TSOH-IFP)	    blocklength(9) = latent
    ! (TSOH-IFP)	    blocklength(10)= latent
	! (TSOH-IFP)		
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr) 
	! (TSOH-IFP)		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr) 
	! (TSOH-IFP)		displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
	! (TSOH-IFP)		displacement(3) = 5*realex
	! (TSOH-IFP)		displacement(4) = 5*realex + intex 
	! (TSOH-IFP)		displacement(5) = 5*realex + 2*intex
	! (TSOH-IFP)		displacement(6) = 5*realex + 2*intex + realex
	! (TSOH-IFP)		displacement(7) = 5*realex + 2*intex + 9*realex
    ! (TSOH-IFP)	    displacement(8) = 5*realex + 2*intex + 9*realex + 8 * realex
    ! (TSOH-IFP)	    displacement(9) = 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex
    ! (TSOH-IFP)	    displacement(10)= 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex + latent*realex
	! (TSOH-IFP)		
	! (TSOH-IFP)		oldtype(0:2) = MPI_double_precision
	! (TSOH-IFP)		oldtype(3:4)   = MPI_INTEGER
	! (TSOH-IFP)		oldtype(5:7)   = MPI_double_precision
    ! (TSOH-IFP)	    oldtype(8)     = MPI_INTEGER
    ! (TSOH-IFP)	    oldtype(9)     = MPI_double_precision
    ! (TSOH-IFP)	    oldtype(10)    = MPI_double_precision
	! (TSOH-IFP)				
	! (TSOH-IFP)		call MPI_TYPE_STRUCT (11, blocklength, displacement, oldtype,MPI_precbank, ierr)
	! (TSOH-IFP)		call MPI_TYPE_COMMIT (MPI_precbank, ierr)
	! (TSOH-IFP)	end subroutine 	
	
end module
