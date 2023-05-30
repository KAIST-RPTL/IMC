subroutine premc
    use constants
    use variables,  only : E_mode, do_burn
    use input_reader
    use ace_xs, only : setugrid
    use bank_header, only: source_bank
    use simulation, only: bank_initialize, PARA_RANGE, PARA_RANGE2,&
                          draw_geometry
    use FMFD, only: FMFD_allocation, fmfdon, fake_MC, nfm
    use ENTROPY
    use DEPLETION_MODULE
    use transient, only: set_transient
    use TH_HEADER, only: th_on
    use TEMPERATURE, only: TH_INITIAL
    use TALLY, only: SET_MC_TALLY
    use SIMULATION_HEADER, only: source_read
    use FMFD_HEADER, only: p_dep_mc, p_dep_dt, p_dep_dt_pert
    use PERTURBATION, only: perton
    use COSAMPLING, only: n_pert
    !use TALLY, only: p_MC, e_MC
	
	use hex_fmfd, only: hf_allocate
    
    implicit none
    integer:: ista, iend
    real(8):: tt0, tt1
    integer:: iwork1, iwork2, mmm
    
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call READ_CTRL
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        call read_inventory
        call read_CE_mat
    endif
    call READ_GEOM
    if(tally_switch > 0) call read_tally
    if ( th_on ) then
        call READ_TH
        call TH_INITIAL
    end if

    call READ_DEPLETION
    if (do_burn) then 
        if(icore==score)  print '(A28)', '    Reading Depletion Lib...' 
        call getdepletionlibrary

        if ( iscore ) then
        if ( fmfdon ) then
            allocate(p_dep_dt(n_act,nfm(1),nfm(2),nfm(3)))
            ! Power distribution
            if(perton) allocate(p_dep_dt_pert(n_act,n_pert,nfm(1),nfm(2),nfm(3)))
        else
            allocate(p_dep_mc(n_act,nfm(1),nfm(2),nfm(3)))
        end if
        end if

        ! material indice for parallel calculation
        if(icore==score) print *, 'NMAT', n_materials
        allocate(mp1(0:ncore-1),mp2(0:ncore-1))
        iwork1 = n_materials / ncore
        iwork2 = mod(n_materials,ncore)
        do mmm = 0, ncore-1
        mp1(mmm) = mmm*iwork1+1+min(mmm,iwork2)
        mp2(mmm) = mp1(mmm)+iwork1-1
        if ( iwork2 > mmm ) mp2(mmm) = mp2(mmm) + 1
        end do
    endif 

    ! =========================================================================
    ! Set tally parameters
    call SET_MC_TALLY
    
    !==============================================================================
    !Set material library for burnup and equilibrium Xe135 search 
    call setmat
    
    !===========================================================================
    !Set lethargy grid for hash-based energy look-up
    call setugrid

    ! ==========================================================================
    if ( entrpon ) call ENTRP_INIT
    if ( bprupon ) call SET_PRUP

    ! ==========================================================================
    ! FMFD calculation
    if ( fmfdon ) then ! LINKPOINT
	    if (dduct < 0.0) then
		    call FMFD_ALLOCATION
		else
		    call hf_allocate
		end if
	end if
    
    !===========================================================================
    !Source bank initialize
    !call para_range(1,ngen,ncore,icore,ista,iend)
    call para_range2(1,ngen,ncore,icore,bank_size)
    allocate(source_bank(bank_size))
    call bank_initialize(source_bank)
    call MPI_banktype()
    call MPI_precbanktype()

    ! =========================================================================
    ! fission source reading
    if ( source_read ) call SOURCE_READING(bank_size)
    
    !===========================================================================
    !Set transient variables
    call set_transient()
    
    !===========================================================================
    !Draw geometry
    call draw_geometry()

end subroutine
