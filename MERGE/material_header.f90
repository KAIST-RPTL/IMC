module material_header 
    
    implicit none 

    type Material_CE
        character(len=20)    :: mat_name       ! User-defined name
        integer              :: mat_type       ! 1 fuel 2 clad 3 coolant (cool)
        integer              :: n_iso = 0      ! number of isotopes (nuclides)
        integer, allocatable :: ace_idx(:)     ! index in nuclides array
        integer, allocatable :: ace_idx0(:)    ! index storage for pre-cor
        integer, allocatable :: ace_idx1(:)    ! index storage for pre-cor
        real(8), allocatable :: numden(:)      ! nuclide atom density (#/b-cm)
        real(8), allocatable :: full_numden(:) ! AD for all inventory (#/b-cm)
        real(8), allocatable :: full_numden0(:)! storage for predictor-corrector
        real(8), allocatable :: full_numden1(:)! storage for predictor-corrector
        real(8)              :: temp           ! temperature (MeV)
        real(8)              :: density_gpcc   ! total density in g/cm^3
        real(8)              :: vol            ! material volume (cm^3)
        real(8)              :: rgb(3)       ! Color for PLOT option
        integer, allocatable :: iso_idx(:)     ! Replacing iso_idx in depletion
        integer, allocatable :: zaid(:)
        
        ! Does this material contain fissionable nuclides? Is it depletable?
        logical :: fissionable = .false.
        logical :: depletable = .false.

        ! DUPLICABLE
        logical :: duplicable = .false.
        integer :: geom_count = 0

        real(8) :: flux = 0.0d0
        real(8) :: flux0, flux1, kappa
        real(8) :: pwr  = 0.0d0
        real(8), allocatable :: ogxs(:,:), ogxs0(:,:), ogxs1(:,:)

        integer, allocatable :: dtmc(:)

        !(21/10/12) eflux testing...
        real(8), allocatable :: eflux(:), e2flux(:)
        real(8), allocatable :: eflux0(:), eflux1(:)
        real(8), allocatable :: e2flux0(:), e2flux1(:)

        !(21/11/23) materialwise-fratio
        real(8), allocatable :: fratio(:,:)
        
        ! Isotopes for S(a,b)?
        logical :: sab = .false.
        logical :: therm = .false.
        integer, allocatable :: sablist(:)

        ! Doppler broadening
        logical :: db = .false.

        ! UEG treatment
        real(8), allocatable :: macro_ueg(:,:)
        logical, allocatable :: ures(:)
        integer, allocatable :: uresidx(:)

        real(8) :: ace_temp
        
    end type Material_CE
    
    integer :: n_materials ! # of materials
    
    type(Material_CE), allocatable, target :: materials(:), materials_temp(:)
    type(Material_CE), pointer :: CE_mat_ptr

    
    contains
    
    function find_CE_mat_idx (this, mat_id) result (idx) 
        type(Material_CE) :: this(:) 
        character(*) :: mat_id 
        integer :: i, idx, n 
        
        n = size(this)
        do i = 1, n
            if (trim(this(i)%mat_name) == mat_id) then 
                idx = i 
                return 
            endif
        enddo 
        print *, "no such CE mat id : ", mat_id 
        stop 
        
    end function 

    
    
    
    
    
end module 
