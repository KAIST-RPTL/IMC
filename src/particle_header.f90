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
        !integer    :: ep     ! (ONLY in UEG): Egrid
        
        ! Other physical data
        real(8)    :: wgt           ! particle weight
        real(8)    :: mu            ! angle of scatter
        logical    :: alive         ! is particle alive?
        
        ! Pre-collision physical data
		REAL(8)    :: last_xyz(3)         ! previous location on base universe
        real(8)    :: last_uvw(3)         ! previous direction coordinates
        real(8)    :: last_wgt            ! pre-collision particle weight
        
        ! Indices for various arrays
        integer    :: material      ! index for current material
        integer    :: last_material ! index for last material
        
        ! Temperature of the current cell
        real(8)    :: sqrtkT        ! sqrt(k_Boltzmann * temperature) in MeV
        real(8)    :: last_sqrtKT   ! last temperature
        real(8)    :: kT            ! temperature in MeV

        ! Density fraction
        real(8)    :: dens          ! Density fraction
        
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
		
		integer :: iso 
        
		real(8) :: time = 0d0
		
        ! Secondary particles created
        real(8), allocatable :: urn(:)

		! IFP RELATED
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr  (1:latent)
        real(8) :: trvltime              ! Traveled distance of the neutron from its born:  Modified to time
		REAL(8) :: trvlength    		 ! Traveled distance of the neutron from its born
		
		! IFP ADJOINT FLUX RELATED (TSOH-IFP)
		INTEGER :: i_Parent               ! Index for parent (애비가 누구야?)
		INTEGER :: arr_code(1:nainfo_src) ! --- (X,Y,Z,G) for neutron absorption (ADJOINT related)
		REAL(8) :: arr_PwTL(1:nainfo_src) ! --- WEIGHT    for neutron absorption (ADJOINT related)
		
		! IFP ADJOINT FLUX RELATED (We select nainfo_src out of the store information when exceeds the length) (TSOH-IFP)
		INTEGER, ALLOCATABLE :: ptc_code(:) ! MESH WISE (FOR SPATIAL DISTRIBUTION) / ENERGY BIN WISE (FOR SPECTRUM)
		REAL(8), ALLOCATABLE :: ptc_PwTL(:) ! MESH WISE (FOR SPATIAL DISTRIBUTION) / ENERGY BIN WISE (FOR SPECTRUM)
		
		! IFP ADJOINT FLUX RELATED (PROGENITOR) (TSOH-IFP)
		REAL(8) :: ptc_wgt0	                       ! Initial weight of the PTC
		INTEGER :: n_prog   		               ! Number of providing progenitor indexing (START BY ZERO / +1 FOR EVERY COLLISION)

    contains
        procedure :: clear
        procedure :: initialize
        procedure :: set => set_particle
		PROCEDURE :: set_LONG   => set_particle_LONG
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
        !this % rotated = .false.
    
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
        this % g = 1
        
        ! Set up base level coordinates
        this % coord(1) % universe = base_univ
        this % n_coord = 1
        this % last_n_coord = 1
        
        this % vrc_traced = .false.
		this % time = 0
		
		this % in_tet = .false. 
		this % tet = 0 
		this % tet_prev = 0

        this % dens = 1d0


		! IFP RELATED
		this % delayedarr(1:latent) = 0
		this % delayedlam(1:latent) = ZERO
		this % nlifearr(1:latent)   = ZERO
		this % trvltime             = ZERO
		this % trvlength            = ZERO

        if(.not. allocated(this%urn)) then
            allocate(this % urn(1:n_unr)); this % urn = 0D0
        endif

		! IFP BASED ADJOINT FLUX RELATED
		this % arr_code(1:nainfo_src) = 0
		this % arr_PwTL(1:nainfo_src) = 0.d0
		IF(ALLOCATED(this % ptc_code)) DEALLOCATE(this % ptc_code)
		IF(ALLOCATED(this % ptc_PwTL)) DEALLOCATE(this % ptc_PwTL)
		ALLOCATE(this % ptc_code(999))
		ALLOCATE(this % ptc_PwTL(999))
		this % ptc_wgt0 = 0.d0
		this % n_prog   = 0

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

		! IFP RELATED
        this % delayedarr  = source % delayedarr
        this % delayedlam  = source % delayedlam
        this % nlifearr    = source % nlifearr
        this % trvltime    = ZERO

        ! MSR
        !if(source%delayed) print *, 'PREC', source%xyz(1:3), source%G, this%wgt
        if(do_fuel_mv .and. source % delayed .and. curr_cyc > n_inact ) then
            zidx = floor((this%coord(1)%xyz(3)-core_base)/(core_height/real(N_core_axial,8)))+1
            zidx = max(1,min(n_core_axial, zidx))
            ridx = floor((this%coord(1)%xyz(1)**2+this%coord(1)%xyz(2)**2)/core_radius**2*real(n_core_radial,8))+1
            ridx = max(1,min(n_core_radial, ridx))
            !print *, 'prec', zidx, ridx, source%G, this % wgt
            core_prec(zidx, ridx, source%G) = core_prec(zidx, ridx, source%G) + this % wgt
            !print *, icore, 'coreprec', core_prec(:, ridx, 1)
        endif
    end subroutine SET_PARTICLE
    
!===============================================================================
! SET_PARTICLE sets the particle from the source bank (LONG) (TSOH-IFP)
!===============================================================================
    subroutine SET_PARTICLE_LONG(this, source)
		use randoms, only: rang
		USE FMFD_header
		implicit none 
        class(Particle) :: this
        type(bank_LONG) :: source
        integer :: zidx
		integer :: node_xyz(3)   
		INTEGER :: code_xyz
		
        this % coord(1) % xyz(:) = source % xyz(:)
        this % coord(1) % uvw(:) = source % uvw(:)
        this % wgt               = source % wgt
        this % E                 = source % E
        this % G                 = source % G
        this % time              = source % time
		
		! IFP RELATED
        this % delayedarr  = source % delayedarr
        this % delayedlam  = source % delayedlam
        this % nlifearr    = source % nlifearr
        this % trvltime    = ZERO
		
		!> IFP BASED ADJOINT FLUX RELATED
		this % i_Parent  = source % i_Current
		this % trvlength = ZERO
		this % arr_code = source % arr_code
		this % arr_PwTL = source % arr_PwTL

		!> INITIAL WEIGHT FROM THE PREVIOUS BANK
		this % ptc_wgt0 = source % wgt
		
		! determine z mesh idx 
		if (do_fuel_mv .and. source%delayed) then 
			zidx = floor((this%coord(1)%xyz(3)-core_base)/(core_height/real(N_core_axial,8)))+1
			core_prec(source%G,zidx,1) = core_prec(source%G,zidx,1) + this % wgt
		endif 
    end subroutine SET_PARTICLE_LONG
    
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
