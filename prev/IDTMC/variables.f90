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

    !****
    integer:: latidx(3)
    !$OMP THREADPRIVATE(latidx)
    
    
    integer :: tally_switch !> o for off // 1 for on
    real(8) :: keff, k_col, k_tl
    real(8) :: k_vrc, fiss_vrc, loss_vrc, fiss_last, keff_vrc
    real(8) :: DMC_loss, DMC_prod, DMC_keff
    
    real(8), allocatable:: kprt(:)
    real(8) :: Nominal_Power
    real(8) :: cyc_power, avg_power, cyc_power0, w_tot, w_totprev

    real(8) :: wgt_min_dyn   = 1.0d-1
    real(8) :: wgt_split_dyn = 2.0d0 
    
    ! depletaion parameters
    logical :: do_burn = .false.
    
    ! Tetrahedral parameters 
	logical :: do_gmsh = .false.
	logical :: do_gmsh_VRC = .false.
    
    ! Transient parameters 
	logical :: do_transient = .false. 
	logical :: do_PCQS = .false. 
	logical :: do_DMC = .false. 
	real(8) :: k_steady = 1.0d0 



    !==============================================================================
    ! MPI parameters 
    integer ::    & 
    & ncore,      & !Total number of cores
    & icore,      & !My core id
    & score = 0 , & !Server rank
    & ierr,       & !Error
    & core,       & !Global communicator
    & iscore = .false.

    ! Error message
    character(len=80) :: err_msg !> error message when exception handler is called

end module 
