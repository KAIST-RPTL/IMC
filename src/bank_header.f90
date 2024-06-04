module bank_header 
	use mpi 
    use variables, only: latent, nainfo_src
	
    implicit none
	! --- FOR ST & DP MC SIMULATION
    type :: Bank
        real(8) :: wgt           ! weight of bank site
        real(8) :: xyz(1:3)      ! location of bank particle
        real(8) :: uvw(1:3)      ! diretional cosines
        real(8) :: E             ! energy for CE
        integer :: G             ! energy group if in MG mode.
		INTEGER :: G_delayed     ! energy group for precursor in MG mode (when delayed == .TRUE.) / INCLUDED NEWLY (필요함)
        logical :: delayed
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
        ! In case of latent is defined as 'parameter' / WILL STICK TO THIS APPROACH
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr  (1:latent)
		! DUMMY VARIABLE FOR COMMUNICATION (이거 빼면 MPI 송수신 시 nlifearr 파트에 에러생김)
		REAL(8) :: DUMMY           
    end type Bank
	
	! --- FOR ADJOINT SPECTRUM TALLY (ONLY FOR STEADY-STATE SIMULATION)
	TYPE :: Bank_LONG
        real(8) :: wgt           ! weight of bank site
        real(8) :: xyz(1:3)      ! location of bank particle
        real(8) :: uvw(1:3)      ! diretional cosines
        real(8) :: E             ! energy for CE
        integer :: G             ! energy group if in MG mode.
		INTEGER :: G_delayed     ! energy group for precursor in MG mode (when delayed == .TRUE.) / INCLUDED NEWLY (필요함)
        logical :: delayed
		real(8) :: time
		real(8) :: beta(1:8)
		real(8) :: lambda(1:8)
        ! In case of latent is defined as 'parameter' / WILL STICK TO THIS APPROACH
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr  (1:latent)
		! DUMMY VARIABLE FOR COMMUNICATION (이거 빼면 MPI 송수신 시 nlifearr 파트에 에러생김)
		REAL(8) :: DUMMY  
		! IFP DISTRIBUTION TALLY
		INTEGER :: i_Parent  ! Index for parent  from previous cycle
		INTEGER :: i_Current ! Index for current odrder of the stored source
		INTEGER :: arr_code(1:nainfo_src) ! --- (X,Y,Z,G) for neutron absorption (ADJOINT related) / ENERGY BIN INDEX FOR TALLYING ENERGY SPECTRUM
		REAL(8) :: arr_PwTL(1:nainfo_src) ! --- WEIGHT    for neutron absorption (ADJOINT related)
	END TYPE Bank_LONG
    
	! --- ONLY FOR TDMC / NOT USED FOR ST & DP CALCULATION
    type :: PrecBank
        real(8) :: wgt          ! weight of bank site
        real(8) :: xyz(1:3)     ! location of bank particle
        real(8) :: E            ! Incident energy for CE
        integer :: G            ! Incident energy group if in MG mode (PCQS시 delayed neutron energy group)
		integer :: idx			! material or isotope index
		real(8) :: time         
		real(8) :: beta(1:8)    
		real(8) :: lambda(1:8)  
		! --- PCQS RELATED
		INTEGER :: prec_G       ! Energy group for precursor itself   (모핵종의 에너지 그룹; npg)
		REAL(8) :: lambda_d     ! PrecBank banking 할 때 사용된 beta 값
    end type PrecBank
	
	! (NOT USED) type :: VRCBank 
    ! (NOT USED)     real(8) :: wgt          ! weight of bank site
    ! (NOT USED)     real(8) :: xyz(1:3)     ! location of bank particle
    ! (NOT USED)     real(8) :: uvw(1:3)     ! location of bank particle
    ! (NOT USED)     real(8) :: E            ! Incident energy for CE
    ! (NOT USED)     integer :: G            ! Incident energy group if in MG mode.
	! (NOT USED) end type VRCBank
	
    ! Source and fission bank
    type(Bank), allocatable, target :: source_bank(:)
    type(Bank), allocatable, target :: fission_bank(:)
    type(Bank), allocatable         :: temp_bank(:)
    type(Bank)                      :: thread_bank(15000)
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
	
	! *************************************************************************************************************
	! +++ Source and fission bank (LONG) --- SHARES bank_idx, split_idx, init_idx
	! *************************************************************************************************************
    type(Bank_LONG), allocatable, target :: source_bank_LONG (:)
    type(Bank_LONG), allocatable, target :: fission_bank_LONG(:)
	type(Bank_LONG), allocatable         :: temp_bank_LONG(:)
    type(Bank_LONG)                      :: thread_bank_LONG (15000)
    !$OMP THREADPRIVATE(thread_bank_LONG)
	
	!> HEAP STACK APPROACH / COULD STORE [arr_PwTL] & [arr_code] as well
	TYPE(Bank_LONG), ALLOCATABLE :: mat_source_bank_LONG(:,:)
	! *************************************************************************************************************
	
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
	
	integer :: MPI_bank
	INTEGER :: MPI_precbank
	INTEGER :: MPI_bank_LONG
	
	contains 
    subroutine MPI_banktype() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:9) :: blocklength, displacement, oldtype 
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 3
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 1
		blocklength(7) = 1
		blocklength(8) = 8
		blocklength(9) = 8
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex +   intex
		displacement(6) = 8*realex + 2*intex
		displacement(7) = 8*realex + 2*intex + logicex
		displacement(8) = 8*realex + 2*intex + logicex + realex
		displacement(9) = 8*realex + 2*intex + logicex + realex + 8*realex
		
		oldtype(0:3) = MPI_double_precision
		oldtype(4:5) = MPI_INTEGER
		oldtype(6)   = MPI_logical
		oldtype(7:9) = MPI_double_precision
		
		call MPI_TYPE_STRUCT (10, blocklength, displacement, oldtype,MPI_bank, ierr)
		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	end subroutine 
	
    subroutine MPI_precbanktype() 
		integer :: ierr 
		integer :: realex, intex, logicex
		integer, dimension(0:9) :: blocklength, displacement, oldtype 
		
		blocklength(0) = 1       ! real(8)
		blocklength(1) = 3       ! real(8)
		blocklength(2) = 1       ! real(8)
		blocklength(3) = 1       ! integer
		blocklength(4) = 1       ! integer
		blocklength(5) = 1       ! real(8)
		blocklength(6) = 8       ! real(8)
		blocklength(7) = 8       ! real(8)
		! PCQS REALTED
		blocklength(8) = 1 ! INTEGER
		blocklength(9) = 1 ! REAL(8)
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr) 
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		
		displacement(0)  = 0
		displacement(1)  = realex
		displacement(2)  = 4*realex
		displacement(3)  = 5*realex
		displacement(4)  = 5*realex + intex 
		displacement(5)  = 5*realex + 2*intex
		displacement(6)  = 5*realex + 2*intex + realex
		displacement(7)  = 5*realex + 2*intex + 9*realex
        displacement(8)  = 5*realex + 2*intex + 9*realex + 8 * realex
        displacement(9)  = 5*realex + 2*intex + 9*realex + 8 * realex + intex
		
		oldtype(0:2)   = MPI_double_precision
		oldtype(3:4)   = MPI_INTEGER
		oldtype(5:7)   = MPI_double_precision
        oldtype(8)     = MPI_INTEGER
        oldtype(9)     = MPI_double_precision
				
		call MPI_TYPE_STRUCT (10, blocklength, displacement, oldtype,MPI_precbank, ierr)
		call MPI_TYPE_COMMIT (MPI_precbank, ierr)
	end subroutine MPI_precbanktype

    subroutine MPI_banktype_ifp() 
		integer :: ierr 
		integer :: realex, intex, logicex
        integer, dimension(0:13) :: blocklength, displacement, oldtype
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 3
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 1
		blocklength(7) = 1
		blocklength(8) = 8
		blocklength(9) = 8
        ! ADJOINT MPI_BANK
        blocklength(10) = latent
        blocklength(11) = latent
        blocklength(12) = latent
		! DUMMY
		blocklength(13) = 1
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex +   intex
		displacement(6) = 8*realex + 2*intex
		displacement(7) = 8*realex + 2*intex + logicex
		displacement(8) = 8*realex + 2*intex + logicex + realex
		displacement(9) = 8*realex + 2*intex + logicex + realex + 8*realex
        displacement(10)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex
        displacement(11)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex
        displacement(12)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex
		displacement(13)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex
		
		oldtype(0:3) = MPI_double_precision
		oldtype(4:5) = MPI_INTEGER
		oldtype(6)   = MPI_logical
		oldtype(7:9) = MPI_double_precision
        oldtype(10)  = MPI_INTEGER
        oldtype(11)  = MPI_double_precision
        oldtype(12)  = MPI_double_precision
		oldtype(13)  = MPI_double_precision
		
		call MPI_TYPE_STRUCT (14, blocklength, displacement, oldtype,MPI_bank, ierr)
		call MPI_TYPE_COMMIT (MPI_bank, ierr)
	end subroutine MPI_banktype_ifp
	
    subroutine MPI_banktype_ifp_LONG() 
		integer :: ierr 
		integer :: realex, intex, logicex
        integer, dimension(0:17) :: blocklength, displacement, oldtype
		
		blocklength(0) = 1
		blocklength(1) = 3
		blocklength(2) = 3
		blocklength(3) = 1
		blocklength(4) = 1
		blocklength(5) = 1
		blocklength(6) = 1
		blocklength(7) = 1
		blocklength(8) = 8
		blocklength(9) = 8
        ! ADJOINT MPI_BANK
        blocklength(10) = latent
        blocklength(11) = latent
        blocklength(12) = latent
		! DUMMY
		blocklength(13) = 1
		! IFP DISTRIBUTION TALLY
		blocklength(14) = 1
		blocklength(15) = 1
		blocklength(16) = nainfo_src
		blocklength(17) = nainfo_src
		
		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr)
		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr)
		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr)
		displacement(0) = 0
		displacement(1) =   realex
		displacement(2) = 4*realex
		displacement(3) = 7*realex
		displacement(4) = 8*realex
		displacement(5) = 8*realex +   intex
		displacement(6) = 8*realex + 2*intex
		displacement(7) = 8*realex + 2*intex + logicex
		displacement(8) = 8*realex + 2*intex + logicex + realex
		displacement(9) = 8*realex + 2*intex + logicex + realex + 8*realex
        displacement(10)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex
        displacement(11)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex
        displacement(12)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex
		displacement(13)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex
		displacement(14)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex + realex
		displacement(15)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex + realex +   intex
		displacement(16)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex + realex + 2*intex
		displacement(17)= 8*realex + 2*intex + logicex + realex + 8*realex + 8*realex + latent*intex + latent*realex + latent*realex + realex + 2*intex + nainfo_src*intex
		
		oldtype(0:3) = MPI_double_precision
		oldtype(4:5) = MPI_INTEGER
		oldtype(6)   = MPI_logical
		oldtype(7:9) = MPI_double_precision
        oldtype(10)  = MPI_INTEGER
        oldtype(11)  = MPI_double_precision
        oldtype(12)  = MPI_double_precision
		oldtype(13)  = MPI_double_precision
		oldtype(14)  = MPI_INTEGER
		oldtype(15)  = MPI_INTEGER
		oldtype(16)  = MPI_INTEGER
		oldtype(17)  = MPI_double_precision
		
		call MPI_TYPE_STRUCT (18, blocklength, displacement, oldtype,MPI_bank_LONG, ierr)
		call MPI_TYPE_COMMIT (MPI_bank_LONG, ierr)
	end subroutine MPI_banktype_ifp_LONG
	
    ! (MUTED)	subroutine MPI_precbanktype_ifp() 
	! (MUTED)		integer :: ierr 
	! (MUTED)		integer :: realex, intex, logicex
	! (MUTED)		integer, dimension(0:10) :: blocklength, displacement, oldtype 
	! (MUTED)		
	! (MUTED)		
	! (MUTED)		blocklength(0) = 1
	! (MUTED)		blocklength(1) = 3
	! (MUTED)		blocklength(2) = 1
	! (MUTED)		blocklength(3) = 1
	! (MUTED)		blocklength(4) = 1
	! (MUTED)		blocklength(5) = 1
	! (MUTED)		blocklength(6) = 8
	! (MUTED)		blocklength(7) = 8
    ! (MUTED)	    ! ADJOINT MPI_PRECBANK
    ! (MUTED)	    blocklength(8) = latent
    ! (MUTED)	    blocklength(9) = latent
    ! (MUTED)	    blocklength(10)= latent
	! (MUTED)		
	! (MUTED)		call MPI_TYPE_EXTENT(MPI_double_precision, realex, ierr) 
	! (MUTED)		call MPI_TYPE_EXTENT(MPI_INTEGER, intex, ierr) 
	! (MUTED)		call MPI_TYPE_EXTENT(MPI_logical, logicex, ierr) 
	! (MUTED)		displacement(0) = 0; displacement(1) = realex; displacement(2) = 4*realex
	! (MUTED)		displacement(3) = 5*realex
	! (MUTED)		displacement(4) = 5*realex + intex 
	! (MUTED)		displacement(5) = 5*realex + 2*intex
	! (MUTED)		displacement(6) = 5*realex + 2*intex + realex
	! (MUTED)		displacement(7) = 5*realex + 2*intex + 9*realex
    ! (MUTED)	    displacement(8) = 5*realex + 2*intex + 9*realex + 8 * realex
    ! (MUTED)	    displacement(9) = 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex
    ! (MUTED)	    displacement(10)= 5*realex + 2*intex + 9*realex + 8 * realex + latent*intex + latent*realex
	! (MUTED)		
	! (MUTED)		oldtype(0:2) = MPI_double_precision
	! (MUTED)		oldtype(3:4)   = MPI_INTEGER
	! (MUTED)		oldtype(5:7)   = MPI_double_precision
    ! (MUTED)	    oldtype(8)     = MPI_INTEGER
    ! (MUTED)	    oldtype(9)     = MPI_double_precision
    ! (MUTED)	    oldtype(10)    = MPI_double_precision
	! (MUTED)				
	! (MUTED)		call MPI_TYPE_STRUCT (11, blocklength, displacement, oldtype,MPI_precbank, ierr)
	! (MUTED)		call MPI_TYPE_COMMIT (MPI_precbank, ierr)
	! (MUTED)	end subroutine 	
	
end module
