subroutine premc
    use constants
    use variables
    use input_reader
    use ace_xs, only : setugrid, getMicroXS, GET_OTF_DB_MIC
    use bank_header, only: source_bank
    use simulation, only: bank_initialize, draw_geometry
    use FMFD, only: FMFD_allocation, fmfdon, fake_MC, nfm
    use ENTROPY
    use DEPLETION_MODULE
    use transient, only: set_transient
    use TH_HEADER, only: th_on
    use TEMPERATURE, only: TH_INITIAL
    use TALLY, only: SET_MC_TALLY, SET_MC_TALLY_ADJ
    !use TALLY, only: p_MC, e_MC
    use FMFD_HEADER, only: p_dep_mc, p_dep_dt, p_dep_dt_pert
    use COSAMPLING, only: n_pert
    use PERTURBATION, only: perton
    
    implicit none
    integer :: iwork1, iwork2, mm, i, j, zaid, iso, iso_, rx, imat, niso, cnt
    real(8) :: ttt0, ttt1
    real(8) :: xs1(6), urn(n_unr), xs2
    type(AceFormat), pointer :: ac
    logical :: found
    real(8), allocatable :: numden(:)
    integer, allocatable :: sablist(:), aceidx(:)
    
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call read_ctrl
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        ttt0 = omp_get_wtime() 
        call read_CE_mat
        ttt1 = omp_get_wtime()
        if(icore == score) print *, '    Time to read CE_mat.inp [s]:', ttt1-ttt0
    endif
    ttt0 = omp_get_wtime() 
    call read_geom('geom.inp', 0)
    !call read_geom
    ttt1 = omp_get_wtime()
    if(icore == score) print *, '    Time to read geom.inp [s]:', ttt1-ttt0
    if(tally_switch > 0) call read_tally
    if ( th_on ) then
        call READ_TH
        call TH_INITIAL
    end if
	
    ! =========================================================================
    ! Set tally parameters
    call SET_MC_TALLY
	
	! =========================================================================
	! SET ADJOINT INFO TALLY
	CALL SET_MC_TALLY_ADJ
	
    !===========================================================================
    !Set lethargy grid for hash-based energy look-up
    !=======21/12/02 MOD: use UGRID in DEPLETION ==========
    call read_depletion
    if (do_burn) then 
        if(icore==score)  print '(A28)', '    Reading Depletion Lib...' 
        call getENDFdepletionlibrary
        if ( iscore ) then
        if ( fmfdon ) then
            allocate(p_dep_dt(n_act,nfm(1),nfm(2),nfm(3)))
            ! Power distribution
            if(perton) allocate(p_dep_dt_pert(n_act,n_pert,nfm(1),nfm(2),nfm(3)))
        else
            allocate(p_dep_mc(n_act,nfm(1),nfm(2),nfm(3)))
        end if
        end if

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
            print *, 'NMAT', n_materials
            do mm = 0, ncore-1
           !     print *, mm, mpigeom(1:ngeom,mm)
            enddo
        endif
    endif 

	! *** PERTAINS FOR CE CALCULATION ONLY
	IF(E_mode == 1) THEN
		ace = ace(1:num_iso)
		if(sab_iso > 0) sab = sab(1:sab_iso)
		if(therm_iso > 0) therm = therm(1:therm_iso)
		call read_mgtally
		call setugrid
		call setDBPP(.false.)
	END IF
	IF(E_mode == 0) do_mgtally = .FALSE.

    ! Shrink material
    if( do_burn ) then
        do imat = 1, n_materials
            if ( .not. materials(imat) % depletable ) cycle
            niso = count( materials(imat) % numden > 0d0 )
    
            materials(imat) % n_iso = niso
            allocate( numden (niso) )
            allocate( sablist(niso) )
            allocate( aceidx (niso) )
    
            cnt = 0
            do i = 1, materials(imat) % n_iso
                if ( materials(imat) % numden(i) > 0d0 ) then
                    cnt = cnt + 1
                    numden(cnt) = materials(imat) % numden(i)
                    sablist(cnt)= materials(imat) % sablist(i)
                    aceidx(cnt) = materials(imat) % ace_idx(i)
                endif
            enddo
    
            if ( cnt /= niso ) print *, 'WTF?', trim(materials(imat)%mat_name), cnt, niso
    
            call move_alloc ( numden, materials(imat) % numden ) 
            call move_alloc ( sablist, materials(imat) % sablist ) 
            call move_alloc ( aceidx, materials(imat) % ace_idx ) 
        enddo
    endif

    IF(ALLOCATED(ugrid)) deallocate(ugrid)
    call setugrid

    if ( do_iso_ueg .or. do_ueg) then
        call setuegrid
    endif

    ! UNIONIZED GRID TREATMENT
    call setueg

    ! UNRESOLVED RESONANCE UEG
    if(do_ures) call setures

    ! NEED TO BE DISCARDED LATER
    call setogxs
    !==============================================================================
    ! Set material library for burnup and equilibrium Xe135 search 
    call setmat
	
    ! ==========================================================================
    if ( icore == score ) call ENTRP_INIT
    if ( mprupon ) call SET_PRUP

    ! ==========================================================================
    ! FMFD calculation
    if ( fmfdon ) call FMFD_allocation
    
    ! ===========================================================================
    ! Source bank initialize / INCLUDES adjoint calculation
	! ===========================================================================
    allocate(source_bank(ngen))
    call bank_initialize(source_bank)
    ! --- IFP BASED CALCULATION
    allocate(betaarr(n_act-latent,0:8)); allocate(genarr(n_act-latent));
    allocate(alphaarr(n_act-latent)); allocate(lamarr(n_act-latent,0:8))
    allocate(betad(0:8,n_act)); betad = 0.d0
	! --- TALLY ADJFLUX OPTION CHECKED DURING SET_MC_TALLY.F90
	IF(tally_adj_flux) THEN
		ALLOCATE(adj_phi    (n_act-latent_LONG,num_adj_group,nfm_adj(1),nfm_adj(2),nfm_adj(3)))
		ALLOCATE(FOR_phi    (n_act            ,num_adj_group,nfm_adj(1),nfm_adj(2),nfm_adj(3)))
		ALLOCATE(adj_phi_cyc(                  num_adj_group,nfm_adj(1),nfm_adj(2),nfm_adj(3)))
		ALLOCATE(FOR_phi_cyc(                  num_adj_group,nfm_adj(1),nfm_adj(2),nfm_adj(3)))
		ALLOCATE(adj_phi_cyc_VEC(num_adj_group*nfm_adj(1)*nfm_adj(2)*nfm_adj(3)))
		ALLOCATE(FOR_phi_cyc_VEC(num_adj_group*nfm_adj(1)*nfm_adj(2)*nfm_adj(3)))
		adj_phi     = 0.d0
		FOR_phi     = 0.d0
		adj_phi_cyc = 0.d0
		FOR_phi_cyc = 0.d0
		adj_phi_cyc_VEC = 0.d0
		FOR_phi_cyc_VEC = 0.d0
	END IF
	! --- TALLY ENERGY SPECTRUM (TSOH-IFP)
	IF(tally_adj_spect) THEN
		IF(NOT(tally_for_spect)) THEN
			IF(icore == score) PRINT*,' [WARNING] TALLY ADJOINT SPECTRUM REQUIRES TALLYING FORWARD SPECTRUM'
			tally_for_spect = .TRUE.
		END IF
	END IF
	IF(tally_for_spect) THEN
		! +++ FOR MG SIMULATION, TALLY FOR EACH GROUP
		IF(E_mode == 0) THEN
			ALLOCATE(FOR_phi_ene_cyc(      n_group)); FOR_phi_ene_cyc = 0.d0
			ALLOCATE(FOR_phi_ene    (n_act,n_group)); FOR_phi_ene     = 0.d0
			IF(tally_adj_spect) THEN
				ALLOCATE(adj_phi_ene_cyc(             n_group)); adj_phi_ene_cyc = 0.d0
				ALLOCATE(adj_phi_ene    (n_act-latent,n_group)); adj_phi_ene     = 0.d0				
			END IF
		! +++ FOR CE SIMULATION, TALLY FOR USER PRESCRIBED ENERGY STRUCTURE (AT THE MOMENT, DEFAULT GROUP OF SERPENT USED)
		ELSE
			SELECT CASE(idx_egroup)
				CASE(1)
					n_egroup_spect = n_egroup_spect1; ALLOCATE(egroup_spect(n_egroup_spect))
					  egroup_spect =   egroup_spect1
				CASE(2)
					n_egroup_spect = n_egroup_spect2; ALLOCATE(egroup_spect(n_egroup_spect))
					  egroup_spect =   egroup_spect2 
				CASE(3)
					n_egroup_spect = n_egroup_spect3; ALLOCATE(egroup_spect(n_egroup_spect))
					  egroup_spect =   egroup_spect3
			END SELECT
			ALLOCATE(egroup_spect_bin(n_egroup_spect-1))
			DO i = 1,n_egroup_spect-1
				egroup_spect_bin(i) = 0.5d0 * (egroup_spect(i) + egroup_spect(i+1))
			END DO
			ALLOCATE(FOR_phi_ene_cyc(      n_egroup_spect-1)); FOR_phi_ene_cyc = 0.d0
			ALLOCATE(FOR_phi_ene    (n_act,n_egroup_spect-1)); FOR_phi_ene     = 0.d0
			IF(tally_adj_spect) THEN
				ALLOCATE(adj_phi_ene_cyc(             n_egroup_spect-1)); adj_phi_ene_cyc = 0.d0
				ALLOCATE(adj_phi_ene    (n_act-latent,n_egroup_spect-1)); adj_phi_ene     = 0.d0				
			END IF
		END IF
	END IF
	
	! +++ INITIALIZE THE BANK ACCORDINGLY
	call MPI_banktype_ifp()
	call MPI_precbanktype()
	IF(do_IFP_LONG) THEN
		call MPI_banktype_ifp_LONG()
	END IF	
	
    !===========================================================================
    ! Set transient variables
	call set_transient()
	
    !===========================================================================
    ! Draw geometry 
    if( icore == score ) call draw_geometry()
    call MPI_BARRIER(core, ierr)

    ! 23/12/01 : DBRC application
    if(n_iso0K > 0) then
        do i = 1, n_iso0K
            ACE0KLOOP: do j = 1, num_iso
                if( ace0K(i)%zaid == ace(j)%zaid ) then
                    if(icore==score) print *, '    0K Data Connected to ', trim(ace(j)%xslib)
                    ace(j) % resonant = i
                endif
            enddo ACE0KLOOP
        enddo
    endif
	
end subroutine premc
