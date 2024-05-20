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
    logical :: do_iso_ueg = .false.

    ! OGXS DIRECT TALLY
    logical :: do_rx_tally = .false.

    ! IFP ADJOINT
    logical :: do_ifp = .false.
    integer, parameter :: latent = 10
    real(8) :: betaeff, gentime
    real(8) :: denom,  gen_numer, gen_prompt, beta_numer(8), lam_denom(8),denom_prompt
    real(8), allocatable :: betaarr(:,:), genarr(:), alphaarr(:), lamarr(:,:), betad(:,:)

	! *** TSOH (IFP-BASED ADJOINT DISTRIBUTION CALCULATION) *** !

	! IFP ADJOINT DISTRIBUTION / SPATIAL & ENERGY SEPARATELY TREATED
	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! (*) FOR SPATIAL ADJOINT DISTRIBUTION (APPLICABLE ONLY FOR x,y,z geometry)
	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! -> FOR MG: GROUP-WISE ADJOINT DISTRIBUTION
	! -> FOR CE: NOT SUPPORTED AT THE MOMENT...
	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! (*) FOR ADJOINT SPECTRUM
	! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	! -> FOR MG: GROUP-WISE ADJOINT SPECTRUM
	! -> FOR CE: USER-DEFINED ENERGY GROUP EMPLOYED
	LOGICAL :: do_IFP_LONG     = .FALSE.
	LOGICAL :: tally_adj_flux  = .FALSE.
	LOGICAL :: meshon_adjflux  = .FALSE.
	LOGICAL :: tally_for_spect = .FALSE.
	LOGICAL :: tally_adj_spect = .FALSE.
	INTEGER :: latent_long = 20
	INTEGER, PARAMETER :: nainfo_src = 20
	
	!> RELATED TO ADJOINT FLUX DISTRIBUTION TALLY
	INTEGER :: num_adj_group  = 1   ! # OF ADJ GROUP FOR SPATIAL TALLYING
    real(8) :: fm0_adj(3)           ! (x0,y0,z0)
    real(8) :: fm1_adj(3)           ! (x1,y1,z1)
    real(8) :: fm2_adj(3)           ! (x1-x0,y1-y0,z1-z0)
    integer :: nfm_adj(3) = 0       ! (nx,ny,nz)
    real(8) :: dfm_adj(3)           ! (dx,dy,dz)
	
	REAL(8), ALLOCATABLE :: adj_phi_cyc_VEC(:)       ! VECTORIZED adj_phi_cyc (x,y,z,g) order
	REAL(8), ALLOCATABLE :: FOR_phi_cyc_VEC(:)       ! VECTORIZED FOR_phi_cyc (x,y,z,g) order
	REAL(8), ALLOCATABLE :: adj_phi_cyc    (:,:,:,:) !       #group (adj), node_x, node_y, node_z
	REAL(8), ALLOCATABLE :: FOR_phi_cyc    (:,:,:,:) !       #group (FOR), node_x, node_y, node_z
 	REAL(8), ALLOCATABLE :: adj_phi      (:,:,:,:,:) ! #cyc, #group (adj), node_x, node_y, node_z
	REAL(8), ALLOCATABLE :: FOR_phi      (:,:,:,:,:) ! #cyc, #group (FOR), node_x, node_y, node_z
	
	REAL(8), ALLOCATABLE :: adj_phi_STORE(:,:,:,:)   !       #group (adj), node_x, node_y, node_z
	
	LOGICAL :: zigzag_adjflux = .FALSE.
	INTEGER :: n_zz_adj       = 0
	INTEGER :: zz_div_adj     = 0
	INTEGER, ALLOCATABLE :: zzf0_adj(:)
	INTEGER, ALLOCATABLE :: zzf1_adj(:)
	INTEGER, ALLOCATABLE :: zzf2_adj(:)

	!> RELATED TO ADJOINT SPECTRUM TALLY
	INTEGER :: idx_egroup     = 1                 ! INDEX FOR ENERGY GROUP STRUCTURE (1: SERPENT 70group / 2: LANL 30group / 3: HELIOS 70group)
	INTEGER :: n_egroup_spect                     ! NUMBER OF ENERGY GROUP STRUCTURE
	REAL(8), ALLOCATABLE :: egroup_spect    (:)   ! ENERGY GROUP STRUCTURE                  [SIZE: n_egroup_spect]
	REAL(8), ALLOCATABLE :: egroup_spect_bin(:)   ! ENERGY GROUP STRUCTURE (MID VALUE; BIN) [SIZE: n_egroup_spect-1]
	REAL(8), ALLOCATABLE :: FOR_phi_ene_cyc (:)   ! ENERGY SPECTRUM (ENERGY GROUP BIN-WISE) [SIZE: n_egroup_spect-1]
	REAL(8), ALLOCATABLE :: FOR_phi_ene     (:,:) ! ENERGY SPECTRUM (ENERGY GROUP BIN-WISE) [SIZE: n_egroup_spect-1]
	
	REAL(8), ALLOCATABLE :: adj_phi_ene_cyc(:) ! ENERGY SPECTRUM (ENERGY GROUP BIN-WISE) [SIZE: n_egroup_spect-1]
	REAL(8), ALLOCATABLE :: adj_phi_ene  (:,:) ! ENERGY SPECTRUM (ENERGY GROUP BIN-WISE) [SIZE: n_egroup_spect-1]

    ! MSR parameter
    real(8), allocatable :: fuel_speed(:), active_mesh(:), fuel_stay_time(:)
    real(8) :: fuel_bulk_speed
    integer :: n_mesh_axial = 0
    logical :: do_fuel_mv = .false.
    real(8) :: t_rc       = 0.d0
    real(8) :: core_height
    real(8) :: core_radius
    real(8) :: core_base
    integer :: n_core_axial, n_core_radial
    real(8), allocatable :: core_prec(:,:,:)
    integer :: MSR_leak, MSR_leak0

    ! MODIFIED ( Oct. 29 2023 ) / Modified for RZ
    integer :: nr, nz
    real(8) :: axial_axis(2)
    real(8), allocatable :: velocity_r(:,:), velocity_z(:,:) !> nr, nz
    real(8), allocatable :: &
        active_r(:), active_z(:) !> nr+1, nz+1 sized
    real(8) :: riser_r
    real(8) :: t_recirc
    
    integer :: flowtype = 0

    ! KNF CR movement
    logical :: do_surf_mv = .false.
    
    ! KNF Tgrid
    logical :: do_temp_grid = .false.
	
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
    real(8), allocatable :: libtemp(:)
    logical, allocatable :: acerecord(:)

    character(len=80) :: directory
    character(len=20) :: title

    integer :: ncell, nsurf, nlatt, nuniv, npcell

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
