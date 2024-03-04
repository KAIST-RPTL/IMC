module particle_header

    use constants
    use variables
    use geometry_header,     only: base_univ, universe
    use bank_header
    use ace_header,          only: n_unr

    implicit none

    private

!===============================================================================
! LOCALCOORD describes the location of a particle local to a single
! universe. When the geometry consists of nested universes, a particle will have
! a list of coordinates in each level
!===============================================================================

    type :: LocalCoord

        ! Indices in various arrays for this level
        integer :: cell      = NONE
        integer :: universe  = NONE
        integer :: lattice   = NONE
        integer :: lattice_x = NONE
        integer :: lattice_y = NONE
        integer :: lattice_z = NONE
        
        ! Particle position and direction for this level
        real(8) :: xyz(3)
        real(8) :: uvw(3)
        real(8) :: dist
    contains
        procedure :: reset => reset_coord
    
    end type LocalCoord

    type, public :: Particle
       
        ! Particle coordinates
        integer          :: n_coord          ! number of current coordinates
        integer          :: cell_instance    ! offset for distributed properties
        type(LocalCoord) :: coord(MAX_COORD) ! coordinates for all levels
        
        ! Particle coordinates before crossing a surface
        integer :: last_n_coord         ! number of current coordinates
        integer :: last_cell(MAX_COORD) ! coordinates for all levels
        
        ! Energy Data
        real(8)    :: E      ! post-collision energy
        real(8)    :: last_E ! pre-collision energy
        integer    :: g      ! post-collision energy group (MG only)
        integer    :: last_g ! pre-collision energy group (MG only)
        
        ! Other physical data
        real(8)    :: wgt           ! particle weight
        real(8)    :: mu            ! angle of scatter
        logical    :: alive         ! is particle alive?
        
        ! Pre-collision physical data
        real(8)    :: last_xyz(3)         ! previous coordinates
        real(8)    :: last_uvw(3)         ! previous direction coordinates
        real(8)    :: last_wgt            ! pre-collision particle weight
        
        ! Indices for various arrays
        integer    :: material      ! index for current material
        integer    :: last_material ! index for last material
        
        ! Temperature of the current cell
        real(8)    :: sqrtkT        ! sqrt(k_Boltzmann * temperature) in MeV
        real(8)    :: last_sqrtKT   ! last temperature
        real(8)    :: kT            ! temperature in MeV
        
        ! Statistical data
        integer    :: n_collision   ! # of collisions
        integer    :: n_cross       ! # of surface cross
        
        ! Tag for S(a,b)
        logical    :: yes_sab = .false.

        ! VRC trace
        logical :: vrc_traced = .false.
		
		! Tetrahedron 
		integer :: tet_face
		logical :: in_tet = .false. 
		integer :: tet
		integer :: tet_prev
		
		! OTHER INFORMATION...
		integer :: iso 
		real(8) :: time = 0d0
        real(8), allocatable :: urn(:)
		
		! (TSOH-IFP)	! IFP RELATED
        ! (TSOH-IFP)	integer, allocatable :: delayedarr(:)
        ! (TSOH-IFP)	real(8), allocatable :: delayedlam(:)
        ! (TSOH-IFP)	real(8), allocatable :: nlifearr(:)
        ! (TSOH-IFP)	real(8)              :: trvltime 
		
		! TSOH-IFP: IFP-REALTED INFORMATION
		INTEGER :: i_Parent     ! Index for parent (애비가 누구야?)
		INTEGER :: dIMT			! Corresponds to delayedarr(:)
		REAL(8) :: dlam         ! Corresponds to delayedlam(:)
		REAL(8) :: nlife		! Corresponds to nlifearr(:)
		REAL(8) :: trvltime     ! Traveled distance of the neutron from its born:  Modified to time
		REAL(8) :: trvlength    ! Traveled distance of the neutron from its born
    contains
        procedure :: clear
        procedure :: initialize
        procedure :: set => set_particle
    end type Particle

contains

!===============================================================================
! RESET_COORD clears data from a single coordinate level
!===============================================================================
    elemental subroutine reset_coord(this)
        class(LocalCoord), intent(inout) :: this
        this % cell = NONE
        this % universe = NONE
        this % lattice = NONE
        this % lattice_x = NONE
        this % lattice_y = NONE
        this % lattice_z = NONE
    end subroutine reset_coord
    
!===============================================================================
! INITIALIZE sets default attributes for a particle from the source bank
!===============================================================================  
    subroutine initialize(this)
  
        class(Particle) :: this
        
        ! Clear coordinate lists
        call this % clear()
        
        ! Set particle to neutron that's alive
        this % alive = .true.
        
        ! clear attributes
        this % material          = NONE
        this % last_material     = NONE
        this % wgt               = ONE
        this % last_wgt          = ONE
        this % sqrtkT            = 0
        this % kT                = 0
        this % n_collision       = 0
        this % n_cross           = 0
        this % g                 = 1
        
        ! Set up base level coordinates
        this % coord(1) % universe = base_univ
        this % n_coord = 1
        this % last_n_coord = 1
        
        this % vrc_traced = .false.
		this % time = 0
		
		this % in_tet = .false. 
		this % tet = 0 
		this % tet_prev = 0

		! TSOH-IFP: ADJOINT RELATED
		this % i_Parent  = 0
		this % dIMT      = 0
		this % dlam      = ZERO
		this % nlife     = ZERO
		this % trvltime  = ZERO
		this % trvlength = ZERO

        if(.not. allocated(this%urn)) then
            allocate(this % urn(1:n_unr)); this % urn = 0D0
        endif

    end subroutine initialize
  
!===============================================================================
! SET_PARTICLE sets the particle from the source bank
!===============================================================================
    subroutine SET_PARTICLE(this, source)
        class(Particle) :: this
        type(bank)        :: source
        integer :: zidx,ridx
        
        this % coord(1) % xyz(:) = source % xyz(:)
        this % coord(1) % uvw(:) = source % uvw(:)
        this % wgt               = source % wgt
        this % E                 = source % E
        this % G                 = source % G
        this % time              = source % time
		
        ! (TSOH-IFP)	if(do_ifp)then
        ! (TSOH-IFP)	    this % delayedarr  = source % delayedarr
        ! (TSOH-IFP)	    this % delayedlam  = source % delayedlam
        ! (TSOH-IFP)	    this % nlifearr    = source % nlifearr
        ! (TSOH-IFP)	    this % trvltime          = 0.D0
        ! (TSOH-IFP)	endif

		! TSOH-IFP: ADJOINT RELATED
		this % i_Parent  = source % i_Current
		this % dIMT      = source % dIMT
		this % dlam      = source % dlam
		this % nlife     = source % nlife
		this % trvltime  = ZERO
		this % trvlength = ZERO

        ! MSR 
        !if(source%delayed) print *, 'PREC', source%xyz(1:3), source%G, this%wgt
        if(do_fuel_mv .and. source % delayed .and. curr_cyc > n_inact ) then
            zidx = floor((this%coord(1)%xyz(3)-core_base)/(core_height/real(N_core_axial,8)))+1
            !ridx = floor((this%coord(1)%xyz(1)**2+this%coord(1)%xyz(2)**2)/core_radius**2*real(n_core_radial,8))+1
            zidx = max(1,min(n_core_axial, zidx))
            !print *, 'prec', this%coord(1)%xyz(3)-core_base, zidx, source%G
            core_prec(source%G,zidx,1) = core_prec(source%G,zidx,1) + this % wgt
        endif
		
    end subroutine SET_PARTICLE
    
    
!===============================================================================
! CLEAR_PARTICLE resets all coordinate levels for the particle
!===============================================================================

    subroutine clear(this)
        class(Particle) :: this
    
        integer :: i
    
        ! remove any coordinate levels
        do i = 1, MAX_COORD
            call this % coord(i) % reset()
        end do
    end subroutine clear

end module 
