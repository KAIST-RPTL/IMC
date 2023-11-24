module variables 

implicit none

! steady-state parameters 
    integer :: n_totcyc, n_inact, n_act, n_batch
    integer :: curr_bat
    integer :: curr_cyc, curr_act
    integer :: n_history
    integer :: ngen         ! neutron generation size
    integer :: b_inact      ! # of inactive cycles for batch 0 calcultion
    integer :: t_inact      ! true inactive cycles
    integer :: t_totcyc     ! true total cycles
    
    integer :: E_mode !> 0 for MG // 1 for CE
    
    
    integer :: tally_switch !> o for off // 1 for on
    real(8) :: keff, k_col, k_tl
    real(8) :: k_vrc, fiss_vrc, loss_vrc, fiss_last, keff_vrc
    real(8) :: DMC_loss, DMC_prod, DMC_keff
	
    real(8), allocatable:: kprt(:)
    real(8) :: Nominal_Power
    real(8) :: cyc_power, avg_power, cyc_power0, w_tot, w_totprev
    real(8), allocatable :: cyc_p_arr(:)
    !FLUX TEST
    
    ! Considered for NFY interpolation
    real(8) :: fis_thermal = 0.d0
    real(8) :: fis_epi = 0.d0
    real(8) :: fis_fast = 0.d0
    real(8) :: thermal_ub = 112.d-6 ! sqrt(0.0251eV*500keV)
    real(8) :: fast_lb = 2.6458 ! sqrt(500keV*14MeV)

    !FPCUT
    real(8) :: fpcut = 0d0

    ! Xenon
    logical :: Xe_search = .false.

	real(8) :: wgt_min_dyn   = 1.0d-1
	real(8) :: wgt_split_dyn = 2.0d0 
	
    ! depletaion parameters
    logical :: chain_reduction = .false.
    logical :: do_burn = .false.
    
	! Tetrahedral parameters 
	logical :: do_gmsh = .false.
	logical :: do_gmsh_VRC = .false.
	
	! Transient parameters 
	logical :: do_transient = .false. 
	logical :: do_PCQS = .false. 
	logical :: do_DMC = .false. 
	real(8) :: k_steady = 1.0d0 
	
    ! Unresolved Resonance -- Need to be fixed later
    logical :: do_ures = .false.
    real(8) :: ures_cut = 1E-9

    ! Multigroup Tally
    logical :: do_mgtally   = .true.
    integer :: n_mg
    real(8), allocatable :: Ebin(:), micro_flux(:), micro_fis(:)

    ! 1G XS tally
    real(8) :: ogflx, ogtot, ogcap, ogfis, ogabs
    real(8) :: u238cap, u238flx 

    ! UNIONIZED GRID
    logical :: do_ueg = .false.

    ! OGXS DIRECT TALLY
    logical :: do_rx_tally = .false.

    ! IFP ADJOINT
    logical :: do_ifp = .false.
    integer, parameter :: latent = 10
    real(8) :: betaeff, gentime
    real(8) :: denom,  gen_numer, gen_prompt, beta_numer(8), lam_denom(8),denom_prompt
    real(8), allocatable :: betaarr(:,:), genarr(:), alphaarr(:), lamarr(:,:), betad(:,:)

    ! MSR parameter
    real(8) :: fuel_speed = 0.d0
    logical :: do_fuel_mv = .false.
    real(8) :: t_rc       = 0.d0
    real(8) :: core_height
    real(8) :: core_radius
    real(8) :: core_base
    integer :: n_core_axial, n_core_radial
    real(8), allocatable :: core_prec(:,:,:)
    integer :: MSR_leak, MSR_leak_kcol

    ! MODIFIED ( Oct. 29 2023 )
    ! Modified for RZ

    integer :: nr, nz
    real(8) :: axial_axis(2)
    real(8), allocatable :: velocity_r(:,:), velocity_z(:,:) !> nr, nz
    real(8), allocatable :: &
        active_r(:), active_z(:) !> nr+1, nz+1 sized
    real(8) :: riser_r
    real(8) :: t_recirc
    
    integer :: flowtype = 0
	
    !==============================================================================
    ! MPI parameters 
	integer ::    & 
    & ncore,      & !Total number of cores
    & icore,      & !My core id
    & score = 0 , & !Server rank
    & ierr,       & !Error
    & core,        & !Global communicator
    & iscore = .false.


    ! Error message
    character(len=80) :: err_msg !> error message when exception handler is called
    
    ! (MODIFYING...) ACE library path
    character(len=80):: acelib
    character(len=80), allocatable :: libname(:), libpath(:)
    logical, allocatable :: acerecord(:)

    character(len=80) :: directory
    character(len=20) :: title

contains

    subroutine para_range(n1, n2, nprocs, irank, ista, iend)
        integer :: iwork1, iwork2 
        integer, intent(in) :: n1, n2, nprocs, irank 
        integer, intent(inout) :: ista, iend
        
        iwork1 = (n2 - n1 + 1) / nprocs
        iwork2 = MOD(n2 - n1 + 1, nprocs)
        ista = irank * iwork1 + n1 + MIN(irank, iwork2)
        iend = ista + iwork1 - 1
        if (iwork2 > irank) iend = iend + 1
    
    end subroutine 

end module 
