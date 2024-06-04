subroutine premc
    use constants
    use variables,  only : E_mode, do_burn
    use input_reader
    use ace_xs, only : setugrid
    use bank_header, only: source_bank
    use simulation, only: bank_initialize, draw_geometry
    use FMFD, only: FMFD_allocation, fmfdon, fake_MC, nfm
    use ENTROPY
    use DEPLETION_MODULE
    use transient, only: set_transient
    use TH_HEADER, only: th_on
    use TEMPERATURE, only: TH_INITIAL
    use TALLY, only: SET_MC_TALLY
    !use TALLY, only: p_MC, e_MC

    
    implicit none
    
    integer :: i, j
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call read_ctrl
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        call read_inventory
        call read_CE_mat
    endif
    call read_geom 
    if(tally_switch > 0) call read_tally
    if ( th_on ) then
        call READ_TH
        call TH_INITIAL
    end if

    ! =========================================================================
    ! Set tally parameters
    call SET_MC_TALLY
    !print *, 'MCTALLY_DONE'
    
    !===========================================================================
    !Set lethargy grid for hash-based energy look-up
    call setugrid
    !=======21/12/02 MOD: use UGRID in DEPLETION ==========
    call read_depletion
    if (do_burn) then 
        if(icore==score)  print '(A28)', '    Reading Depletion Lib...' 
        !call getdepletionlibrary
        call getENDFdepletionlibrary
    endif 
    !==============================================================================
    !Set material library for burnup and equilibrium Xe135 search 
    call setmat
    ! ==========================================================================
    if ( icore == score ) call ENTRP_INIT
    if ( mprupon ) call SET_PRUP

    ! ==========================================================================
    ! FMFD calculation
    if ( fmfdon ) call FMFD_allocation
    
    !===========================================================================
    !Source bank initialize
    ! ADJOINT : LATENT setup
    ! later swap to input parameter
    allocate(betaarr(n_act-latent,0:8)); allocate(genarr(n_act-latent));
    allocate(alphaarr(n_act-latent)); allocate(lamarr(n_act-latent,0:8))
    allocate(betad(0:8,n_act)); betad = 0.d0
    allocate(source_bank(ngen))
    call bank_initialize(source_bank)
	call MPI_banktype()
	call MPI_precbanktype()
	!call MPI_vrcbanktype()
	
    !===========================================================================
    !Set transient variables
	call set_transient()

    !===========================================================================
    !Draw geometry 
	call draw_geometry()
    if(icore==score) print *, '    PreMC Done!'
end subroutine
