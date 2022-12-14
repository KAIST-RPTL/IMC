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
    integer :: iwork1, iwork2, mm
    
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call read_ctrl
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        !call read_inventory
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
    
    
    !===========================================================================
    !Set lethargy grid for hash-based energy look-up
    !=======21/12/02 MOD: use UGRID in DEPLETION ==========
    call read_depletion
    if (do_burn) then 
        if(icore==score)  print '(A28)', '    Reading Depletion Lib...' 
        call getENDFdepletionlibrary

        ngeom = 0; allocate(tmpgeom(1:n_materials))

        do mm = 1, n_materials
            !print *, mm, n_materials, materials(mm)%depletable, ngeom, icore
            if(materials(mm)%depletable) then
                ngeom = ngeom + 1
                tmpgeom(ngeom) = mm
            endif
        enddo
        
        if(ngeom == 0 .and. do_burn) then
            print *, '      No Depletable Region Defined'
            stop
        endif

        tmpgeom = tmpgeom(1:ngeom)
            
        iwork1 = ngeom / ncore
        if(mod(ngeom,ncore)/=0) iwork1 = iwork1 + 1

        allocate(mpigeom(iwork1, 0:ncore-1)); mpigeom = 0
        
        do mm = 1, ngeom
            if(mod(mm,ncore)>0) then
                mpigeom(1+mm/ncore,mod(mm,ncore)-1) = tmpgeom(mm)
            else
                mpigeom(mm/ncore,ncore-1) = tmpgeom(mm)
            endif
        enddo
        totgeom = ngeom
        deallocate(tmpgeom); ngeom = iwork1
        if(icore==score)then
            print *, 'NGEOM', ngeom
            do mm = 0, ncore-1
           !     print *, mm, mpigeom(1:ngeom,mm)
            enddo
        endif

    endif 

    call read_mgtally
    call setugrid

    if ( do_ueg ) then
        call setuegrid
!        iwork1 = num_iso / ncore
!        if(mod(num_iso, ncore)/=0) iwork1 = iwork1 + 1
!        allocate(mpiace(iwork1, 0:ncore-1)); mpiace = 0
!        do mm = 1, num_iso
!            if(mod(mm, ncore) > 0) then
!                mpiace(1+mm/ncore,mod(mm,ncore)-1) = mm
!            else
!                mpiace(mm/ncore,ncore-1) = mm
!            endif
!        enddo
!        nace = iwork1
!        if(icore==score) then
!            print *, 'ISO DISTRIBUTED...'
!            do mm = 0, ncore-1
!                print *, mm, mpiace(1:nace,mm)
!            enddo
!        endif
    endif

    ! UNIONIZED GRID TREATMENT
    call setueg

    ! UNRESOLVED RESONANCE UEG
    if(do_ures) call setures

    ! NEED TO BE DISCARDED LATER
    !call setogxs
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
	
	
end subroutine
