module particle_header

    use constants
    use variables
    use geometry_header,     only: base_univ, universe
    use ace_header, only: n_unr
    use bank_header

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
        ! Basic data
        !integer(8) :: id            ! Unique ID
        
        ! Particle coordinates
        integer          :: n_coord          ! number of current coordinates
        integer          :: cell_instance    ! offset for distributed properties
        type(LocalCoord) :: coord(MAX_COORD) ! coordinates for all levels
        
        
        !class(universe), pointer :: univ     !> universe pointer 
        
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
        
        !! Pre-collision physical data
        !real(8)    :: last_xyz_current(3) ! coordinates of the last collision or
                                            !  reflective/periodic surface crossing
                                            !  for current tallies
        !real(8)    :: last_xyz(3)         ! previous coordinates
        real(8)    :: last_uvw(3)         ! previous direction coordinates
        real(8)    :: last_wgt            ! pre-collision particle weight
        !real(8)    :: absorb_wgt          ! weight absorbed for survival biasing
        
        !! What event last took place
        !logical    :: fission       ! did the particle cause implicit fission
        !integer    :: event         ! scatter, absorption
        !integer    :: event_nuclide ! index in nuclides array
        !integer    :: event_MT      ! reaction MT
        !integer    :: delayed_group ! delayed group
        
        !! Post-collision physical data
        !integer    :: n_bank        ! number of fission sites banked
        !real(8)    :: wgt_bank      ! weight of fission sites banked
        !integer    :: n_delayed_bank(MAX_DELAYED_GROUPS) ! number of delayed fission
                                                         ! sites banked
        
        !! Indices for various arrays
        !integer    :: surface       ! index for surface particle is on
        !integer    :: cell_born     ! index for cell particle was born in
        integer    :: material      ! index for current material
        integer    :: last_material ! index for last material
        
        ! Temperature of the current cell
        real(8)    :: sqrtkT        ! sqrt(k_Boltzmann * temperature) in MeV
        real(8)    :: last_sqrtKT   ! last temperature
        real(8)    :: kT            ! temperature in MeV
        
        ! Statistical data
        integer    :: n_collision   ! # of collisions
        integer    :: n_cross       ! # of surface cross
        
        !! Track output
        !logical    :: write_track = .false.

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
        !integer(8) :: n_secondary = 0
        !type(Bank) :: secondary_bank(MAX_SECONDARY)
        real(8), allocatable :: urn(:)
        
        ! ADJOINT : IFP related
        integer :: delayedarr(1:latent)
        real(8) :: delayedlam(1:latent)
        real(8) :: nlifearr(1:latent)
        real(8) :: trvltime ! Traveled distance of the neutron from its born:  Modified to time

        
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
        !this % surface           = NONE
        !this % cell_born         = NONE
        this % material          = NONE
        this % last_material     = NONE
        !this % last_sqrtkT       = NONE
        this % wgt               = ONE
        this % last_wgt          = ONE
        !this % absorb_wgt        = ZERO
        !this % n_bank            = 0
        !this % wgt_bank          = ZERO
        this % sqrtkT            = 0
        this % kT                = 0
        this % n_collision       = 0
        this % n_cross           = 0
        !this % fission           = .false.
        !this % delayed_group     = 0
        !this % n_delayed_bank(:) = 0
        this % g = NONE
        
        ! Set up base level coordinates
        this % coord(1) % universe = base_univ
        this % n_coord = 1
        this % last_n_coord = 1
        
        this % vrc_traced = .false.
		this % time = 0
		
		this % in_tet = .false. 
		this % tet = 0 
		this % tet_prev = 0

        ! ADJOINT related
        this % delayedarr(1:latent) = 0
        this % delayedlam(1:latent) = ZERO
        this % nlifearr(1:latent)   = ZERO
        this % trvltime             = ZERO

        if(.not. allocated(this%urn)) allocate(this%urn(1:n_unr))
		
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

        ! ADJOINT related
        this % delayedarr  = source % delayedarr
        this % delayedlam  = source % delayedlam
        this % nlifearr    = source % nlifearr
        this % trvltime          = 0.D0

        ! MSR
        !if(source%delayed) print *, 'PREC', source%xyz(1:3), source%G, this%wgt
        if(do_fuel_mv .and. source % delayed) then
            zidx = floor((this%coord(1)%xyz(3)-core_base)/(core_height/real(N_core_axial,8)))+1
            !ridx = floor((this%coord(1)%xyz(1)**2+this%coord(1)%xyz(2)**2)/core_radius**2*real(n_core_radial,8))+1
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
