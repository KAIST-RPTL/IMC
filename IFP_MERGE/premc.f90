subroutine premc
    use constants
    use variables,  only : E_mode, do_burn
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
    use TALLY, only: SET_MC_TALLY
    !use TALLY, only: p_MC, e_MC
    use FMFD_HEADER, only: p_dep_mc, p_dep_dt, p_dep_dt_pert
    use COSAMPLING, only: n_pert
    use PERTURBATION, only: perton
    
    implicit none
    integer :: iwork1, iwork2, mm, i, j, zaid, iso, iso_, rx
    real(8) :: xs1(6), urn(n_unr), xs2
    type(AceFormat), pointer :: ac
    logical :: found
    
    !===========================================================================
    !Read input geometry and cross section data
    call init_var
    call read_ctrl
    if (E_mode == 0) then 
        if(icore==score) print '(A30)', '    Multi-group Energy Mode...' 
        call read_MG_XS
    elseif (E_mode == 1) then 
        if(icore==score) print '(A29)', '    Continuous Energy Mode...' 
        call read_CE_mat
    endif
    call read_geom('geom.inp', 0)
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


    ! 23/12/04 : Preprocessor
    do i = 1, n_materials
        if ( .not. materials(i) % db ) cycle
        do iso = 1, materials(i) % n_iso
            if( abs(materials(i) % temp - ace(materials(i)%ace_idx(iso)) % temp) > 1E-3*K_B ) then 
                found = .false.
                ISO_LOOP: do iso_ = 1, num_iso
                    if( abs(materials(i) % temp - ace(iso_) % temp) < 1E-3 * K_B .and. &
                        ace(materials(i)%ace_idx(iso)) % zaid == ace(iso_) % zaid) then
                        materials(i) % ace_idx(iso) = iso_
                        found = .true.
                    exit ISO_LOOP
                    endif
                enddo ISO_LOOP
                    
                if( .not. found ) then
                    num_iso = num_iso + 1
                    ace(num_iso) = ace(materials(i)%ace_idx(iso))
                    ac => ace(num_iso)
                    ac % temp = materials(i) % temp
                    do j = 1, ac % NXS(3)
                        if ( ac % E(j) < 1d0 ) exit
                        call GET_OTF_DB_MIC(materials(i)%temp, materials(i)%ace_idx(iso), ac % E(j), xs1)
                        ac % sigt(j) = xs1(1)
                        ac % sigel(j) = xs1(2)
                        ac % sigd(j) = xs1(3)-xs1(4)
                        ac % sigf(j) = xs1(4)
                        do rx = 1, ac % NXS(5)
                            call GET_OTF_DB_MT(materials(i)%temp, materials(i)%ace_idx(iso), ac % E(j), rx, xs2)
                            ac % sig_MT(rx) % cx(j) = xs2
                        enddo
                    end do
                    nullify(ac)
                    materials(i) % ace_idx(iso) = num_iso
                    if(icore==score) print *, trim(materials(i)%mat_name), ': Adjusted XS for ', trim(ace(num_iso) % xslib), ' to', ace(num_iso) % temp/K_B
                else
                    if(icore==score) print *, trim(materials(i)%mat_name), ': Linked XS to ', trim(ace(materials(i)%ace_idx(iso)) % xslib), ' with T:', ace(materials(i)%ace_idx(iso)) % temp / K_B
                endif
            elseif( abs(materials(i) % temp - ace(materials(i)% ace_idx(iso)) % temp) < 1E-3 * K_B) then
                if(icore==score) print *, 'WARNING: Invalid Temperature for ', trim(materials(i)%mat_name), materials(i)%temp/K_B, ace(materials(i)%ace_idx(iso))%temp/K_B
            else
                if(icore==score) print *, trim(materials(i)%mat_name), ': no adjust required for ', trim(ace(materials(i)%ace_idx(iso))%xslib)
            endif
        enddo
    enddo



    ace = ace(1:num_iso)
    if(sab_iso > 0) sab = sab(1:sab_iso)
    if(therm_iso > 0) therm = therm(1:therm_iso)
    call read_mgtally
    call setugrid

    if ( do_ueg ) then
        call setuegrid
    endif

    ! UNIONIZED GRID TREATMENT
    call setueg

    ! UNRESOLVED RESONANCE UEG
    if(do_ures) call setures

    ! NEED TO BE DISCARDED LATER
    call setogxs
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
	! --- IFP BASED CALCULATION
    ALLOCATE(betaarr (n_act-latent,0:8))
	ALLOCATE(genarr  (n_act-latent))
    ALLOCATE(alphaarr(n_act-latent))
	ALLOCATE(lamarr  (n_act-latent,0:8))
	! --- CONVENTIONAL BETA CALCULATION
    ALLOCATE(betad   (0:8,n_act)); betad = 0.d0

    allocate(source_bank(ngen))
    call bank_initialize(source_bank)
	
	! (TSOH-IFP): MODIFIED BANK STRUCTURE
	CALL MPI_banktype()
	CALL MPI_precbanktype()
	! call MPI_banktype_ifp()
	! call MPI_precbanktype_ifp()
	! call MPI_vrcbanktype()
	
    !===========================================================================
    !Set transient variables
	call set_transient()
	
    !===========================================================================
    !Draw geometry 
	call draw_geometry()
	
!    if(icore==score) then
!        do iso = 1, size(materials)
!            do i = -44,0
!                !call GET_OTF_DB_MIC(600d0*K_B,iso,sqrt(1d1)**dble(i),xs1)
!                !xs2 = getMicroXS(iso, qrt(1d1)**(dble(i)))
!                xs1 = getMacroXS(materials(iso), (10d0)**(dble(i)/4d0), 600d0*K_B, urn)
!                print '(A,A,E10.2,2E18.10)', 'XSTEST: ', trim(materials(iso)%mat_name), (10d0)**(dble(i)/4d0), xs1(1), xs1(2)
!            enddo
!            print *, 'TEMP:', materials(iso)%temp/K_B
!        enddo
!    endif


    ! 23/12/04 : Doppler Broadening Preprocessor...        

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

end subroutine
