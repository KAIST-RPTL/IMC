module simulation
    use omp_lib
    use mpi
    use constants
    use tracking,           only : transport, transport_DT, transport_VRC, transport_dynamic, transport_pcqs, transport_PCQS_init
    use variables
    use particle_header
    use randoms
    use simulation_header
    use bank_header
    use geometry,           only : find_cell, find_cell_xyz, cell_contains
    use geometry_header,    only : cells, surfaces, universes, sgrid, lattices, lattice_coord
    use surface_header,     only : surface
    use XS_header
    use material_header
    use DEPLETION_MODULE,   only : inline_xenon
    use ENTROPY
    use MPRUP,              only : GENSIZE, MPRUP_DIST, CYCLECHANGE
    use ace_header,         only : ace
    use tally,              only : TallyCoord, TallyFlux, TallyPower,tally_buf, tally1, & 
                                   tally2, meshon, FM_ID, tallyon, &
                                   PROCESS_TALLY, NORM_TALLY, &
                                   TALLY_THREAD_INITIAL, FindTallyBin, &
								   meshon_tet_vrc, mesh_power, &
                                   MC_tally
    use FMFD,               only : FMFD_initialize_thread, FMFD_initialize, &
                                   FMFD_solve, n_skip, n_acc, fsd, &
                                   process_FMFD, NORM_FMFD, fmfdon, &
                                   nfm, k_fmfd, fake_MC, MCBU
    use TEMPERATURE,        only : TEMP_SOLVE, NORM_TH, PROCESS_TH, &
                                   TEMP_CONVERGE, TEMP_DISTRIBUTE
    use TH_HEADER,          only : th_on
    use tetrahedral, 		only : find_tet, transport_tet, tet, find_tet_old
	use transient 
	use DMC, 				only : set_dynamic_bank
	use PCQS!, 				only : PCQS_keff, set_PCQS_source, n_pcqs_totcyc,n_pcqs_act, PKE_init, PKE_keff_tally, PKE_keff0!
	use rgbimage_m
	
    implicit none 
    
	interface gatherBank
		module procedure  gatherSourceBank
		module procedure  gatherPrecBank
		module procedure  gatherVRCBank
	end interface 
	
	integer, parameter :: bank_max = 1000000
	
    contains 
    
subroutine simulate_history(bat,cyc)
    use FMFD_HEADER, only: n_fake
    implicit none
    integer, intent(in):: bat
    integer, intent(in):: cyc
    integer :: i, j, k, isize, i_surf
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: tid, my_id
    integer :: ista, iend
    real(8) :: Jtemp
    real(8), allocatable :: shape(:), rcvbuflong(:)
    integer :: id(3)
    real(8) :: rcv_buf, rcv_buf_long(8)
    real(8), allocatable :: rcv_msh(:)
    integer :: realex, intex, restype, ndata, idata
    integer, dimension(0:4) :: blocklength, displacement, oldtype 
    integer, allocatable :: ircnt(:), idisp(:) 
    real(8) :: time1, time2
	integer :: i_bin(4)
    real(8) :: adj_sum, totwgt
    real(8) :: t_fm1, t_fm2

    MSR_leak = 0d0
    if (allocated(fission_bank)) call move_alloc(fission_bank, source_bank)
    if ( icore == score ) then
        call SHENTROPY(source_bank)
        call FET_CALC(source_bank)
        if ( mprupon ) call GENSIZE(bat,cyc)
    end if
    isize = size(source_bank)
    if ( mprupon .or. ( .not. mprupon .and. genup ) ) &
        call MPRUP_DIST(isize,source_bank(:))
    if ( .not. mprupon .and. genup ) call CYCLECHANGE(cyc)
    allocate(fission_bank(0))
	
    !> Distribute source_bank to slave nodes 
    call para_range(1, isize, ncore, icore, ista, iend)        
    k_col = 0; k_tl = 0; k_vrc = 0; fiss_vrc = 0; loss_vrc = 0;
    if ( fmfdon ) call FMFD_initialize()
    if ( do_burn ) MC_tally(bat,:,:,:,:,:,:) = 0
    cyc_power = 0;
    !if(.not. allocated(cyc_p_arr)) allocate(cyc_p_arr(0:ncore-1))
    !cyc_p_arr = 0;
    
    ! ADJOINT related
    denom = 0; beta_numer = 0; gen_numer = 0; lam_denom = 0; !gen_prompt = 0

	n_col = 0; n_cross = 0 
	
	if (allocated(prec_bank)) deallocate(prec_bank)
	allocate(prec_bank(0))
	if (allocated(prompt_bank)) deallocate(prompt_bank)
	allocate(prompt_bank(0))

	
    !$omp parallel private(p) shared(source_bank, fission_bank, temp_bank, prec_bank, ista, iend)
      thread_bank(:)%wgt = 0; bank_idx = 0; prec_idx = 0 ; init_idx = 0
      if (tallyon .and. .not. fmfdon) call TALLY_THREAD_INITIAL(cyc)
      if ( fmfdon ) call FMFD_initialize_thread()
      !$omp do reduction(+: k_col, k_tl) 
        do i= ista, iend 
            call p%initialize()
            call p%set(source_bank(i))
			
			! initialize p%tet value (TODO)
			if (do_gmsh .and. curr_cyc > n_inact) then 
				i_bin = FindTallyBin(p)
				if ( i_bin(1) > 0 ) then 
					p%tet = find_tet(p%coord(p%n_coord)%xyz)
					p%in_tet = .true.
				else 
					p%in_tet = .false.
				endif 
				!if (p%tet .le. 0) p%in_tet = .false.
			endif 
            do while (p%alive == .true.)
				call transport(p)
            enddo 
		
            !if buffer is almost full -> add to the fission bank
            !if (bank_idx > int(size(thread_bank)*0.01*(80-OMP_GET_THREAD_NUM()))) then 
            if ( bank_idx > 7500 ) then 
              !$omp critical
				call gatherBank(fission_bank, thread_bank, bank_idx)
              !$omp end critical
                bank_idx = 0
            endif
			
            if ( prec_idx > 10000 ) then 
			  !$omp critical
				call gatherBank(prec_bank, prec_thread, prec_idx)
			  !$omp end critical
				prec_idx = 0
            endif
			
            if ( init_idx > 10000 ) then 
			  !$omp critical
				call gatherBank(prompt_bank, thread_bank_init, init_idx)
			  !$omp end critical
				init_idx = 0
            endif
        enddo
      !$omp end do

      ! print *, icore, 'CYC done'
		
		
      !$omp critical
		call gatherBank(fission_bank, thread_bank, bank_idx)
        !> normalize thread tally parameters (can be done outside critical)
        if ( tallyon .and. .not. fmfdon ) call NORM_TALLY(bat,cyc)
        if ( fmfdon ) call NORM_FMFD(cyc)
        if ( th_on .and. .not. fmfdon ) call NORM_TH()
        
      !$omp end critical
	  
      !$omp critical
		call gatherBank(prec_bank, prec_thread, prec_idx)
      !$omp end critical
	  
	  !$omp critical
		call gatherBank(prompt_bank, thread_bank_init, init_idx)
	  !$omp end critical
	  
    !$omp end parallel
	
	
	
    !call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    !> Process tallied FMFD parameters ==========================================
    if ( tallyon .and. .not. fmfdon ) call PROCESS_TALLY(bat,cyc)
    if ( fmfdon .and. cyc > n_skip ) call PROCESS_FMFD(bat,cyc)
    if ( th_on .and. .not. fmfdon ) call PROCESS_TH()
    
    !> Gather keff from the slave nodes =========================================
    call MPI_REDUCE(k_col,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    k_col = rcv_buf
	
    call MPI_REDUCE(k_tl,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    k_tl = rcv_buf
	

    call MPI_REDUCE(cyc_power,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    cyc_power = rcv_buf
    
	if (icore == score) avg_power = avg_power + cyc_power
	
    call MPI_REDUCE(n_col,i,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    n_col = i
	n_col_avg = n_col_avg + n_col 
	
	
	
    if(do_mgtally) then
        allocate(rcvbuflong(n_mg))
        call MPI_REDUCE(micro_flux, rcvbuflong, n_mg, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
        micro_flux(1:n_mg) = rcvbuflong(1:n_mg)
        deallocate(rcvbuflong)
!        allocate(rcvbuflong(n_mg))
!        call MPI_REDUCE(micro_fis, rcvbuflong, n_mg, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
!        micro_fis(1:n_mg) = rcvbuflong(1:n_mg)
!        deallocate(rcvbuflong)
!        call MPI_REDUCE(ogflx, rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!        ogflx = rcv_buf
!        call MPI_REDUCE(ogtot, rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!        ogtot = rcv_buf
!        call MPI_REDUCE(ogcap, rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!        ogcap = rcv_buf
!        call MPI_REDUCE(ogfis, rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!        ogfis = rcv_buf
    endif
	
    if(do_fuel_mv) then
        if( .not. allocated( rcv_msh ) ) allocate(rcv_msh(n_core_axial))
        do k = 1,8
            do j = 1,n_core_radial
            call MPI_REDUCE(core_prec(:,j,k),rcv_msh,k,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
            core_prec(:,j,k) = rcv_msh
            enddo
        enddo
        if(curr_cyc>n_inact .and. icore==score) then
            do i = 1,8
                do j = 1,n_core_radial
                    write(prt_fuel_mv,'(<N_CORE_AXIAL>e15.6)') core_prec(1:n_core_axial,j,i)
                enddo
            enddo
        endif
        core_prec = 0.d0
        !print *, 'leak', MSR_leak
        call MPI_ALLREDUCE(MSR_leak,ndata,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        MSR_leak = ndata
        !if(icore==score) print *, 'leak_red', MSR_leak
    endif
        
        !if(icore == score) print *, curr_cyc, isize, n_col	
    if(do_ifp .and. icore==score .and. curr_cyc>n_inact) then 
        betad(:,curr_cyc-n_inact) = 0.d0
        isize = size(fission_bank)
        totwgt = 0.d0
        do i = 1,isize
            totwgt = totwgt+fission_bank(i)%wgt
            if(fission_bank(i)%delayed .and. fission_bank(i)%G>0) &
            betad(fission_bank(i)%G,curr_cyc-n_inact) = &
            betad(fission_bank(i)%G,curr_cyc-n_inact) + fission_bank(i)%wgt
        enddo
        betad(0,curr_cyc-n_inact) = sum(betad(1:8,curr_cyc-n_inact))
        betad(:,curr_cyc-n_inact) = betad(:,curr_cyc-n_inact)/totwgt
        !print *, curr_cyc-n_inact, totwgt, betad(0:6,curr_cyc-n_inact)
    endif
	
	
    !> Calculate k_eff ==========================================================
    k_col = k_col / real(ngen,8)
    k_tl  = k_tl  / real(ngen,8) 
    keff  = (k_tl + k_col) / 2.0d0
    !keff = k_col
    
    if (icore == score) write(prt_keff,*) keff, k_col, k_tl
    
    call MPI_BCAST(keff, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
    
	
    if (do_ifp .and. curr_cyc > n_inact+latent .and. icore == score) then
        do i = 1,8
            betaarr(curr_cyc-n_inact-latent,i) = beta_numer(i)/denom
            if(lam_denom(i)>0) then
                lamarr(curr_cyc-n_inact-latent,i)  = beta_numer(i)/lam_denom(i)
            else
                lamarr(curr_cyc-n_inact-latent,i)  = 0.d0
            endif
        enddo
        
        betaarr(curr_cyc-n_inact-latent,0) = sum(beta_numer)/denom
        if(sum(lam_denom)>0) then
            lamarr(curr_cyc-n_inact-latent,0)  = sum(beta_numer)/sum(lam_denom)
        else
            lamarr(curr_cyc-n_inact-latent,0)  = 0.d0
        endif

        genarr(curr_cyc-n_inact-latent)  = gen_numer /denom
        !alphaarr(curr_cyc-n_inact-latent)= -(denom*sum(beta_numer(1:8)))/gen_numer 
        alphaarr(curr_cyc-n_inact-latent) = &
                -betaarr(curr_cyc-n_inact-latent,0)/genarr(curr_cyc-n_inact-latent)
        !alphaarr(curr_cyc-n_inact-latent)= -sum(betaarr(curr_cyc-n_inact-latent,1:8))/gen_prompt*denom_prompt
        write(prt_adjoint,*) 'GENAlpha', genarr(curr_cyc-n_inact-latent), alphaarr(curr_cyc-n_inact-latent), gen_numer, denom
        write(prt_adjoint,*) 'TOTAL', (betaarr(curr_cyc-n_inact-latent,0)), (lamarr(curr_cyc-n_inact-latent,0))
        do i = 1,8
            write(prt_adjoint,*) i,betaarr(curr_cyc-n_inact-latent,i), lamarr(curr_cyc-n_inact-latent,i)
        enddo
        writE(prt_adjoint,*) ' '

    endif
    !> Gather fission_bank from the slave nodes =================================        
    ndata = size(fission_bank)
    allocate(ircnt(1:ncore))
    allocate(idisp(1:ncore))
    do i = 1, ncore
        idata = ndata
        call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
        ircnt(i) = idata
    enddo 
    idisp(1) = 0
    do i = 2, ncore 
        idisp(i) = idisp(i-1) + ircnt(i-1)
    enddo
    allocate(temp_bank(1:sum(ircnt)))
    call MPI_ALLGATHERV(fission_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
    deallocate(ircnt); deallocate(idisp); deallocate(fission_bank) 
    call move_alloc(temp_bank, fission_bank)
	
    
	if (do_DMC .and. curr_cyc > n_inact) then 
		!> Gather prec_bank from the slave nodes =================================        
		ndata = size(prec_bank)
		allocate(ircnt(1:ncore))
		allocate(idisp(1:ncore))
		do i = 1, ncore
			idata = ndata
			call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
			ircnt(i) = idata
		enddo 
		idisp(1) = 0
		do i = 2, ncore 
			idisp(i) = idisp(i-1) + ircnt(i-1)
		enddo
		allocate(prec_bank_temp(1:sum(ircnt)))
		call MPI_ALLGATHERV(prec_bank,ndata,MPI_precbank,prec_bank_temp,ircnt,idisp,MPI_precbank,MPI_COMM_WORLD,ierr)
		deallocate(ircnt); deallocate(idisp); deallocate(prec_bank) 
		call move_alloc(prec_bank_temp, prec_bank)
	
	
		!> Gather prompt_bank from the slave nodes ======================   
		ndata = size(prompt_bank)
		allocate(ircnt(1:ncore))
		allocate(idisp(1:ncore))
		do i = 1, ncore
			idata = ndata
			call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
			ircnt(i) = idata
		enddo 
		idisp(1) = 0
		do i = 2, ncore 
			idisp(i) = idisp(i-1) + ircnt(i-1)
		enddo
		allocate(temp_bank(1:sum(ircnt)))
		call MPI_ALLGATHERV(prompt_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
		deallocate(ircnt); deallocate(idisp); deallocate(prompt_bank) 
		call move_alloc(temp_bank, prompt_bank)
		
		
		
		
	endif 

	
    !> Normalize source weight  =================================================
    isize = size(fission_bank)
    fission_bank(:)%wgt = real(ngen,8)/real(isize,8)
    if(do_fuel_mv) fission_bank(:)%wgt = real(ngen,8) / real(isize+MSR_leak,8)	
	fission_bank(:)%time = 0 
	
    !> Normalize tally (flux & power) ===========================================
    if (tally_switch>0 .and. curr_act > 0 .and. do_transient == .false.) then 
		isize = size(TallyFlux) 
        call MPI_REDUCE(TallyFlux, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
		TallyFlux = tally_buf
        call MPI_REDUCE(TallyPower, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
		TallyPower = tally_buf
		
		
        if (icore == score) then 
			!print *, sum(TallyPower), cyc_power, sum(TallyPower)/cyc_power, Nominal_Power
            TallyFlux(:) = TallyFlux(:) * dble(isize) / sum(TallyFlux)
            !TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / sum(TallyPower)
            TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / cyc_power
			
			
			if (do_gmsh) then
				do i = 1, isize
					TallyFlux(i)  = TallyFlux(i)/tet(i)%vol
					TallyPower(i) = TallyPower(i)/tet(i)%vol
				enddo
			else 
				do i = 1, isize
					TallyFlux(i)  = TallyFlux(i)/TallyCoord(i)%vol
					TallyPower(i) = TallyPower(i)/TallyCoord(i)%vol
				enddo
			endif 
			
            !> Print TallyFlux
            write(prt_flux, 100) TallyFlux(:)
            write(prt_powr, 100) TallyPower(:)
        100 format(<isize>ES15.7)
		

        endif 
		
		TallyFlux(:) =0
		TallyPower(:)=0
    endif
	
    !> Solve FMFD and apply FSD shape feedback ==================================
    if ( fmfdon .and. cyc > n_skip .and. .not. MCBU ) then 
            t_fm1 = MPI_WTIME()
            call FMFD_SOLVE(bat,cyc)
            t_fm2 = MPI_WTIME()
        if ( icore == score .and. bat > 0 ) t_det(bat,cyc) = t_fm2-t_fm1
        if ( fake_MC .and. cyc == n_fake ) fmfdon = .false.
        !if ( cyc == n_inact ) fmfdon = .false.
    endif 
	
    ! =========================================================================
    ! temperature distribution
    if ( th_on ) then
        call TEMP_SOLVE
        call TEMP_DISTRIBUTE
    end if
    if ( th_on .and. cyc >= n_act  ) call TEMP_CONVERGE


    !> Inline Equilibrium Xe-135 
    !if(Xe_search .and. cyc > n_inact)then
    !    call inline_xenon()
    !endif
    
    !> initialize the global tally parameters ==================================
    k_col = 0.0d0; k_tl = 0.0d0; 
    

end subroutine 



!===============================================================================
! DINAMIC_MC - DMC simulation subroutine.
!===============================================================================
subroutine dynamic_MC()
    implicit none
    integer :: i, j, k, isize
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: ista, iend
    real(8) :: rcv_buf, temp1, temp2 
	integer :: ndata, idata
    integer, allocatable :: ircnt(:), idisp(:) 
	integer :: size_bank
	real(8) :: temp, beta, w_nav, t0, t1
	real(8) :: w_d, w_p
	real(8) :: time1, time2, time3, time4 
	character(100) :: filename
	integer :: i_bin(4)
	
	if (.not. do_DMC) return 
	time1 = omp_get_wtime()
! 0. Initialize keff tally 
	DMC_loss = 0.0d0
	DMC_prod = 0.0d0
	DMC_keff = 0.0d0
	n_col = 0 
	n_cross = 0
	
! 1. source는 dynamic MC를 위한 dynamic_bank에서 가져옴
	call set_dynamic_bank()
	w_totprev = w_tot
	w_tot = sum(dynamic_bank(:)%wgt) 
	
	
! 2. tracking 하면서 새로운 dynamic_bank 채움
    if (allocated(dynamic_bank)) call move_alloc(dynamic_bank, source_bank)
    isize = size(source_bank)
    allocate(dynamic_bank(0))
	allocate(prec_bank_local(0)) 
	
    !> Distribute source_bank to slave nodes 
    call para_range(1, isize, ncore, icore, ista, iend)
    k_col = 0; k_tl = 0; loss_vrc = 0 
	cyc_power = 0 ; mesh_power = 0
	
	time3 = omp_get_wtime()
	
	k = 0 
	EMP_SRC: do 
	k = k + 1
	
	!print '(a,2i)', 'start', k, icore
	if (allocated(split_bank_temp)) deallocate(split_bank_temp)
	allocate(split_bank_temp(0))
	!write(prt_dynamic,*) "EMP_SRC CYCLE NO.", icore, k 
		
    !$omp parallel private(p) shared(source_bank, dynamic_bank, temp_bank,split_bank_temp)
      thread_bank(:)%wgt = 0; bank_idx = 0; split_idx = 0; prec_idx = 0; init_idx = 0; vrc_idx = 0 
		!write(prt_dynamic,*) "start end : ", ista, iend
      !$omp do reduction(+:DMC_loss, DMC_prod)
        do i= ista, iend 
			!write(prt_dynamic,*) "Particle ", i , source_bank(i)%wgt
            call p%initialize()
            call p%set(source_bank(i))
			
			! initialize p%tet value
			if (do_gmsh .and. curr_cyc > n_inact) then 
				p%tet = find_tet(p%coord(p%n_coord)%xyz)
			endif 
			
			
            do while (p%alive == .true.)
                call transport_dynamic(p)
            enddo 
			
			if ( prec_idx > 14000 ) then 
				!$omp critical
				call gatherBank(prec_bank_local, prec_thread, prec_idx)
				!$omp end critical
				prec_idx = 0
			endif
        enddo
      !$omp end do 
	  
		 
      !$omp critical
        !> gather thread fission bank
		call gatherBank(dynamic_bank, thread_bank, bank_idx)
      !$omp end critical
	  
      !$omp critical
        !> gather thread split bank
		call gatherBank(split_bank_temp, split_thread, split_idx)
      !$omp end critical
	  
      !$omp critical
        !> gather thread precursor bank
		call gatherBank(prec_bank_local, prec_thread, prec_idx)
      !$omp end critical
	  
	  
    !$omp end parallel
	
	  ista = 1; iend = size(split_bank_temp)
		
	  !wgt_min_dyn   = exp(DMC_loss/DMC_prod)*wgt_min_dyn
	  !wgt_split_dyn = exp(DMC_loss/DMC_prod)*wgt_split_dyn
		
	  !print *, icore, k, exp(1-DMC_loss/DMC_prod), iend
		
	  if (iend < 1) exit EMP_SRC
	  deallocate(source_bank)
	  call move_alloc(split_bank_temp, source_bank) 
	  !write (prt_dynamic,*) "REWIND SOURCE BANK "
	enddo EMP_SRC
	
    
	time4 = omp_get_wtime()
	
    !> Gather keff from the slave nodes =========================================
    call MPI_REDUCE(DMC_prod, rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    DMC_prod = rcv_buf
    
    call MPI_REDUCE(DMC_loss,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    DMC_loss = rcv_buf
    
    call MPI_ALLREDUCE(cyc_power,rcv_buf, 1,MPI_REAL8, MPI_SUM,MPI_COMM_WORLD,ierr)
    cyc_power = rcv_buf
	
	
    !> Calculate k_eff ==========================================================
    DMC_keff = DMC_prod / DMC_loss
    
	if (icore == score) cyc_power0 = cyc_power0 + cyc_power
	
	temp = n_col 
	!print *, icore, n_col
	
	!> Quantify computation burden 
    call MPI_REDUCE(temp,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    temp = rcv_buf
	n_col = int(temp)
	n_col_avg = n_col_avg + n_col 
	!if (icore==score) print *, n_col 
	
	if (icore == score) point_power(curr_timestep, curr_cyc-n_inact) = cyc_power / (del_t)
	 
	
	
    !> Gather dynamic_bank from the slave nodes =================================        
    ndata = size(dynamic_bank)
    allocate(ircnt(1:ncore))
    allocate(idisp(1:ncore))
    do i = 1, ncore
        idata = ndata
        call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
        ircnt(i) = idata
    enddo 
    idisp(1) = 0
    do i = 2, ncore 
        idisp(i) = idisp(i-1) + ircnt(i-1)
    enddo
    allocate(temp_bank(1:sum(ircnt)))
    call MPI_ALLGATHERV(dynamic_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
    deallocate(ircnt); deallocate(idisp); deallocate(dynamic_bank) 
    call move_alloc(temp_bank, dynamic_bank)
    
	
    !> Gather prec_bank_local from the slave nodes =================================        
    ndata = size(prec_bank_local)
    allocate(ircnt(1:ncore))
    allocate(idisp(1:ncore))
    do i = 1, ncore
        idata = ndata
        call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
        ircnt(i) = idata
    enddo 
    idisp(1) = 0
    do i = 2, ncore 
        idisp(i) = idisp(i-1) + ircnt(i-1)
    enddo
    allocate(prec_bank_temp(1:sum(ircnt)))
    call MPI_ALLGATHERV(prec_bank_local,ndata,MPI_precbank,prec_bank_temp,ircnt,idisp,MPI_precbank,MPI_COMM_WORLD,ierr)
    deallocate(ircnt); deallocate(idisp); deallocate(prec_bank_local) 
    call move_alloc(prec_bank_temp, prec_bank_local)
	
	
	call move_alloc(prec_bank, prec_bank_temp)
	i = size(prec_bank_temp)
	j = size(prec_bank_local)
	allocate(prec_bank(i+j))
	prec_bank(1:i)   = prec_bank_temp(1:i)
	prec_bank(i+1:i+j) = prec_bank_local(1:j)
	deallocate(prec_bank_local, prec_bank_temp)
	
	
	
    !> Normalize tally (flux & power) ===========================================
    if (tally_switch > 0 .and. do_transient == .true.) then 
		isize = size(TallyFlux) 
        call MPI_REDUCE(TallyFlux, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
		TallyFlux = tally_buf
        call MPI_REDUCE(TallyPower, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
		TallyPower = tally_buf
		
		
        if (icore == score) then 
            TallyFlux(:) = TallyFlux(:) * dble(isize) / sum(TallyFlux)
            !TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / sum(TallyPower)
            
			!if (meshon_tet_vrc) then 
			!	TallyPower(:) = TallyPower(:) * Nominal_Power * (mesh_power / cyc_power) * 1.0d6 / mesh_power
			!else 
				TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / cyc_power
			!endif 
			
			if (do_gmsh) then
				do i = 1, isize
					TallyFlux(i)  = TallyFlux(i)/tet(i)%vol
					TallyPower(i) = TallyPower(i)/tet(i)%vol
				enddo
			else 
				do i = 1, isize
					TallyFlux(i)  = TallyFlux(i)/TallyCoord(i)%vol
					TallyPower(i) = TallyPower(i)/TallyCoord(i)%vol
				enddo
			endif 
			
			
			
			
			!> Open print file 
			if (curr_timestep < 10) then 
				write (filename, "(A22,I1,A4)") "./DMC_data/power", curr_timestep, ".out"
			elseif (curr_timestep < 100) then                           
				write (filename, "(A22,I2,A4)") "./DMC_data/power", curr_timestep, ".out"
			elseif (curr_timestep < 1000) then                          
				write (filename, "(A22,I3,A4)") "./DMC_data/power", curr_timestep, ".out"
			endif 
			open(prt_powr, file=trim(filename), action="write",status="old", position="append")
			
            !> Print TallyPower
            !write(prt_flux, 100) TallyFlux(:)
            write(prt_powr, 100) TallyPower(:)
        100 format(<isize>ES15.7)
		
		
			close(prt_powr)
		

        endif 
		
		TallyFlux(:) =0
		TallyPower(:)=0
    endif	
	
	time2 = omp_get_wtime()

	if (icore==score ) write(prt_wgt,'(I4, 5E14.5)') curr_timestep, temp, time2-time1, time3-time1, time4-time3, time2-time4
	
	
	
    !> initialize the global tally parameters ==================================
    k_col = 0.0d0; k_tl = 0.0d0; 
    
! 3. dynamic_bank에서 prompt bank로 옮김 
	call move_alloc(dynamic_bank, prompt_bank)
	
	
end subroutine





!===============================================================================
! PCQS_MC - PCQS MC simulation subroutine.
!===============================================================================
subroutine PCQS_MC()
	use PCQS, only : PCQS_keff_cyc
    implicit none
    integer :: i, j, k, isize
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: ista, iend
    real(8) :: rcv_buf
	integer :: ndata, idata
    integer, allocatable :: ircnt(:), idisp(:) 
	integer :: size_bank
	real(8) :: temp, beta, w_nav, t0, t1
	real(8) :: w_d, w_p
	real(8) :: time1, time2, time3, time4, time5, time6, time7, time8
	character(100) :: filename
	real(8) :: avg, std
	real(8), allocatable :: PKE_tally_buf(:)
	
	if (.not. do_PCQS) return 
	time1 = omp_get_wtime()
	
	
	call set_PCQS_source() 
	
	
! 0. Initialize keff tally 
	PCQS_keff = 0.0d0
	PCQS_abs = 0.0d0 
	PCQS_prod = 0.0d0 
	PCQS_leak = 0.0d0 
	cyc_power = 0 ; mesh_power = 0
	n_col = 0; n_cross = 0 
	
	if (curr_cyc == n_pcqs_totcyc )  then 
		deallocate(prompt_bank) 
		allocate(prompt_bank(0)) 
	endif 
	
	
	if (curr_timestep == 0) PKE_keff0 = 1.0d0 
	

! 1. source는 
	w_totprev = w_tot
	w_tot = sum(fission_bank(:)%wgt)
	
! 2. tracking 하면서 새로운 dynamic_bank 채움
    if (allocated(fission_bank)) call move_alloc(fission_bank, source_bank)
    isize = size(source_bank)
	
	
    !> Distribute source_bank to slave nodes 
	!call para_range(1, isize, ncore, icore, ista, iend)
	ista = 1 
	iend = isize
	
	
	if (allocated(delayed_bank)) deallocate(delayed_bank)
	allocate(delayed_bank(0))
	
	if (allocated(vrc_bank)) deallocate(vrc_bank)
	allocate(vrc_bank(0))
	
	time3 = omp_get_wtime()
	
    !$omp parallel private(p) shared(source_bank, dynamic_bank, temp_bank)
      thread_bank(:)%wgt = 0; bank_idx = 0; vrc_idx = 0; init_idx = 0 
      !$omp do reduction(+:PCQS_keff)
        do i= ista, iend 
            call p%initialize()
            call p%set(source_bank(i))
			if (E_mode == 1 .and. p%E == 0)	p%alive = .false.
            do while (p%alive == .true.)
                call transport_PCQS(p)
				!print '(I,5E15.5)', icore, p%wgt, p%E, p%coord(1)%xyz!, p%coord(1)%uvw  
            enddo 
			
			if ( bank_idx > 13000 ) then 
				!$omp critical
				call gatherBank(delayed_bank, thread_bank, bank_idx)
				!$omp end critical
				bank_idx = 0
			endif
			
			!if (vrc_idx > 14000) then 
			!	!$omp critical
			!	call gatherBank(vrc_bank, vrc_thread, vrc_idx)
			!	!$omp end critical
			!	vrc_idx = 0
			!	print *, "vrc_idx exceed its limit", vrc_idx 
			!endif 
			
			
			if ( init_idx > 10000 .and. curr_cyc == n_pcqs_totcyc) then 
				!$omp critical
				call gatherBank(prompt_bank, thread_bank_init, init_idx)
				!$omp end critical
				init_idx = 0
			endif
			
        enddo
      !$omp end do 
	  
      !$omp critical
		call gatherBank(delayed_bank, thread_bank, bank_idx)
      !$omp end critical
	  
      !$omp critical
        !> gather thread vrc bank
		call gatherBank(vrc_bank, vrc_thread, vrc_idx)
      !$omp end critical
	  if (curr_cyc == n_pcqs_totcyc) then 
		  !$omp critical
			call gatherBank(prompt_bank, thread_bank_init, init_idx)
		  !$omp end critical
	  endif 
	  
    !$omp end parallel
	
	
	time4 = omp_get_wtime()
	
	
    !> Gather keff from the slave nodes =========================================
    
    call MPI_REDUCE(PCQS_keff,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    PCQS_keff = rcv_buf
    
    call MPI_ALLREDUCE(cyc_power,rcv_buf, 1,MPI_REAL8, MPI_SUM,MPI_COMM_WORLD,ierr)
    cyc_power = rcv_buf
	
    call MPI_ALLREDUCE(mesh_power,rcv_buf, 1,MPI_REAL8, MPI_SUM,MPI_COMM_WORLD,ierr)
    mesh_power = rcv_buf
	
    call MPI_REDUCE(PCQS_abs,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    PCQS_abs = rcv_buf
    call MPI_REDUCE(PCQS_leak,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    PCQS_leak = rcv_buf
    call MPI_REDUCE(PCQS_prod,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    PCQS_prod = rcv_buf

	!> Gather PKE parameters 
	allocate(PKE_tally_buf(N_PKE_tally)) 
    call MPI_REDUCE(PKE_tally,PKE_tally_buf,N_PKE_tally,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    PKE_tally = PKE_tally_buf
	deallocate(PKE_tally_buf) 
	
	if (curr_cyc > n_pcqs_inact) then 
		do i = 1, npg 
			PKE_beta_d(i) = PKE_beta_d(i) + PKE_tally(i) / PKE_tally(2*npg+1)
			PKE_lambda(i) = PKE_lambda(i) + PKE_tally(i) / PKE_tally(npg+i)
		!	PCQS_beta_cyc  (curr_cyc - n_pcqs_inact,i) = PKE_tally(i) / PKE_tally(2*npg+1)
		!	PCQS_lambda_cyc(curr_cyc - n_pcqs_inact,i) = PKE_tally(i) / PKE_tally(npg+i)
		enddo 
		PKE_gen = PKE_gen + PKE_tally(2*npg+2) / PKE_tally(2*npg+1)
		PKE_Z_tally1 = PKE_Z_tally1 + PKE_tally(2*npg+2)
		!PCQS_gen_cyc(curr_cyc - n_pcqs_inact) = PKE_tally(2*npg+2) / PKE_tally(2*npg+1)
		
	endif 
	PKE_tally(:) = 0 
	
	
	!> Quantify computation burden 
    call MPI_REDUCE(n_col,i,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    n_col = i
	n_col_avg = n_col_avg + n_col 
	
    call MPI_REDUCE(n_cross,i,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    n_cross = i
	n_cross_avg = n_cross_avg + n_cross 
	
	
    call MPI_REDUCE(isize,i,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
	
	
	if (curr_cyc > n_pcqs_inact) cyc_power0 = cyc_power0 + cyc_power
	
	
	if (icore==score) keff = PCQS_prod/(PCQS_abs+PCQS_leak)
	
	
    call MPI_BCAST(keff, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr)
	
	time2 = omp_get_wtime()
	
	if (icore==score .and. ( curr_cyc < n_pcqs_totcyc)) write(prt_wgt,'(4I,5E14.5)')  &
	!if (icore==score) write(prt_wgt,'(4I,5E14.5)')  &
		curr_timestep, curr_cyc, n_col, n_cross, time2-time1, time3-time1, time4-time3, time2-time4, cyc_power

	
	
	if (curr_cyc .le. n_pcqs_inact) then 
		if (icore==score) print '(I, F13.6, F10.2, a)', curr_cyc, keff, time2-time1, ' sec'
	else 
		PKE_keff_tally = PKE_keff_tally + keff
		if (icore==score) then 
			ndata = curr_cyc - n_pcqs_inact
			PCQS_keff_cyc(ndata) = keff
			PCQS_power_cyc(ndata) = cyc_power
			avg = sum(PCQS_keff_cyc(1:ndata))/ndata
			std = sqrt(dot_product((PCQS_keff_cyc(1:ndata)-avg),(PCQS_keff_cyc(1:ndata)-avg))/(ndata*(ndata-1)))
			if ( isnan(std) ) std = 0
			print '(I, F13.6, F10.2, a, F12.5, a, F10.1)', curr_cyc, keff, time2-time1, ' sec | avg', avg, ' std', std*1.0d5
			
			!> Average timestep power
			PCQS_power(curr_timestep,1) = cyc_power0 / real(n_pcqs_act,8) 
			PCQS_power(curr_timestep,2) = sqrt(dot_product((PCQS_power_cyc(1:ndata)-sum(PCQS_power_cyc(1:ndata))/ndata), &
											(PCQS_power_cyc(1:ndata)-sum(PCQS_power_cyc(1:ndata))/ndata))/(ndata*(ndata-1)))
		endif 
	endif 
	
	
	!if (curr_cyc == n_pcqs_totcyc .and. icore==score) then 
	!
	!	ndata = n_pcqs_act 
	!	do i = 1, npg 
	!	avg = sum(PCQS_beta_cyc(:,i))/ndata
	!	std = sqrt(dot_product((PCQS_beta_cyc(:,i)-avg),(PCQS_beta_cyc(:,i)-avg))/(ndata*(ndata-1)))
	!	print *, 'beta', i, avg, std 
	!	enddo 
	!	do i = 1, npg 
	!	avg = sum(PCQS_lambda_cyc(:,i))/ndata
	!	std = sqrt(dot_product((PCQS_lambda_cyc(:,i)-avg),(PCQS_lambda_cyc(:,i)-avg))/(ndata*(ndata-1)))
	!	print *, 'lambda', i, avg, std 
	!	enddo 
	!	avg = sum(PCQS_gen_cyc(:))/ndata
	!	std = sqrt(dot_product((PCQS_gen_cyc(:)-avg),(PCQS_gen_cyc(:)-avg))/(ndata*(ndata-1)))
	!	print *, 'gen', avg, std 
	!	
	!endif 
	
	
	
	
	
	
	
	!> Update PKE_gamma 
	if (curr_cyc == n_pcqs_totcyc ) PKE_gamma = log(PKE_n / PKE_n0) / del_t
	
	
	
    !> Gather vrc_bank from the slave nodes =================================        
	!ndata = size(vrc_bank)
    !allocate(ircnt(1:ncore))
    !allocate(idisp(1:ncore))
    !do i = 1, ncore
    !    idata = ndata
    !    call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
    !    ircnt(i) = idata
    !enddo 
    !idisp(1) = 0
    !do i = 2, ncore 
    !    idisp(i) = idisp(i-1) + ircnt(i-1)
    !enddo
    !allocate(vrc_bank_temp(1:sum(ircnt)))
    !call MPI_ALLGATHERV(vrc_bank,ndata,MPI_bank,vrc_bank_temp,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
    !deallocate(ircnt); deallocate(idisp); deallocate(vrc_bank) 
    !call move_alloc(vrc_bank_temp, vrc_bank)
	
	
	!if (do_gmsh_vrc .and. icore==score) then 
	!	if (curr_timestep < 10) then 
	!		write (filename, "(A22,I1,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
	!	elseif (i < 100) then 
	!		write (filename, "(A22,I2,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
	!	elseif (i < 1000) then 
	!		write (filename, "(A22,I3,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
	!	endif 
	!
	!	open(prt_tet_vrc, file=trim(filename), status="old", position="append", action="write")
	!
	!	ndata = size(vrc_bank)
	!	write(prt_tet_vrc,'(2i, e16.7)') ndata, curr_act, mesh_power
	!	do i = 1, ndata
	!		write(prt_tet_vrc,'(8e16.7,I,L)') vrc_bank(i)%xyz, vrc_bank(i)%uvw, vrc_bank(i)%E &
	!										 ,vrc_bank(i)%wgt, vrc_bank(i)%G, vrc_bank(i)%delayed
	!	enddo
	!endif 
	!close(prt_tet_vrc)
	
	
    !!> Normalize tally (flux & power) ===========================================
    !if (tally_switch > 0 .and. do_transient == .true.) then 
	!	isize = size(TallyFlux) 
    !    call MPI_REDUCE(TallyFlux, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
	!	TallyFlux = tally_buf
    !    call MPI_REDUCE(TallyPower, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
	!	TallyPower = tally_buf
	!	
	!	
    !    if (icore == score) then 
    !        TallyFlux(:) = TallyFlux(:) * dble(isize) / sum(TallyFlux)
	!		TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / cyc_power
	!		
	!		if (do_gmsh) then
	!			do i = 1, isize
	!				TallyFlux(i)  = TallyFlux(i)/tet(i)%vol
	!				TallyPower(i) = TallyPower(i)/tet(i)%vol
	!			enddo
	!		else 
	!			do i = 1, isize
	!				TallyFlux(i)  = TallyFlux(i)/TallyCoord(i)%vol
	!				TallyPower(i) = TallyPower(i)/TallyCoord(i)%vol
	!			enddo
	!		endif 
	!		
    !        !> Print TallyFlux
    !        write(prt_flux, 100) TallyFlux(:)
    !        write(prt_powr, 100) TallyPower(:)
    !    100 format(<isize>ES15.7)
	!	
    !
    !    endif 
	!	
	!	TallyFlux(:) =0
	!	TallyPower(:)=0
    !endif	
	
	
end subroutine PCQS_MC


!===============================================================================
! PCQS_INIT - Sample initial PCQS bank from steady-state 
!===============================================================================
subroutine PCQS_INIT()
	use PCQS, only : PKE_Z_tally2
    implicit none
    integer :: i, j, k, isize
    integer :: a,b,c, n
    type(particle) :: p
    logical :: found 
    integer :: ista, iend
    real(8) :: rcv_buf
	integer :: ndata, idata
    integer, allocatable :: ircnt(:), idisp(:) 
	integer :: size_bank
	real(8) :: temp, beta, w_nav, t0, t1
	real(8) :: w_d, w_p
	real(8) :: time1, time2
	character(100) :: filename
	
	
	if (.not. do_PCQS) return 
	
	PCQS_keff = 0.0d0 
	w_tot = sum(fission_bank(:)%wgt) 
! 2. tracking 하면서 새로운 dynamic_bank 채움
    if (allocated(fission_bank)) call move_alloc(fission_bank, source_bank)
    isize = size(source_bank)
	
	if (allocated(prompt_bank)) deallocate(prompt_bank)
	if (allocated(delayed_bank)) deallocate(delayed_bank)
	if(.not.allocated(prompt_bank)) allocate(prompt_bank(0))
	if(.not.allocated(delayed_bank)) allocate(delayed_bank(0))
	allocate(fission_bank(0))
	
	
    !> Distribute source_bank to slave nodes 
    call para_range(1, isize, ncore, icore, ista, iend)
	
	time1 = omp_get_wtime()
    !$omp parallel private(p) 
      thread_bank(:)%wgt = 0; thread_bank_init(:)%wgt = 0; split_thread(:)%wgt = 0; 
	  init_idx = 0; bank_idx = 0; split_idx = 0 
      !$omp do reduction(+:PCQS_keff)
        do i= ista, iend 
            call p%initialize()
            call p%set(source_bank(i))
            do while (p%alive == .true.)
                call transport_PCQS_init(p)
            enddo 
			
			if ( bank_idx > 14000 ) then 
				!$omp critical
				call gatherBank(fission_bank, thread_bank, bank_idx)
				!$omp end critical
				bank_idx = 0
			endif
			if ( init_idx > 13000 ) then 
				!$omp critical
				call gatherBank(prompt_bank, thread_bank_init, init_idx)
				!$omp end critical
				init_idx = 0
			endif
			if ( split_idx > 13000 ) then 
				!$omp critical
				call gatherBank(delayed_bank, split_thread, split_idx)
				!$omp end critical
				split_idx = 0
			endif
			
        enddo
      !$omp end do 
	  
	  !$omp critical
		call gatherBank(fission_bank, thread_bank, bank_idx)
	  !$omp end critical
	  
	  !$omp critical
		call gatherBank(prompt_bank, thread_bank_init, init_idx)
	  !$omp end critical
	  
	  !$omp critical
		call gatherBank(delayed_bank, split_thread, split_idx)
	  !$omp end critical
	  
    !$omp end parallel
	
	
	
	time2 = omp_get_wtime()
	
	
	
    !> Gather fission_bank from the slave nodes =================================        
    ndata = size(fission_bank)
    allocate(ircnt(1:ncore))
    allocate(idisp(1:ncore))
    do i = 1, ncore
        idata = ndata
        call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
        ircnt(i) = idata
    enddo 
    idisp(1) = 0
    do i = 2, ncore 
        idisp(i) = idisp(i-1) + ircnt(i-1)
    enddo
    allocate(temp_bank(1:sum(ircnt)))
    call MPI_ALLGATHERV(fission_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
    deallocate(ircnt); deallocate(idisp); deallocate(fission_bank) 
    call move_alloc(temp_bank, fission_bank)	
	
    isize = size(fission_bank)
    fission_bank(:)%wgt = real(ngen,8)/real(isize,8) 
	
	
	
	!!> Gather prompt_bank from the slave nodes =================================        
	!ndata = size(prompt_bank)
	!
	!if(ndata > bank_max) then 
	!	call move_alloc(prompt_bank, temp_bank)
	!	call combing_source(temp_bank, prompt_bank, bank_max) 
	!	deallocate(temp_bank) 
	!	ndata = size(prompt_bank)
	!endif 
	!
	!allocate(ircnt(1:ncore))
	!allocate(idisp(1:ncore))
	!do i = 1, ncore
	!	idata = ndata
	!	call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
	!	ircnt(i) = idata
	!enddo 
	!idisp(1) = 0
	!do i = 2, ncore 
	!	idisp(i) = idisp(i-1) + ircnt(i-1)
	!enddo
	!allocate(temp_bank(1:sum(ircnt)))
	!call MPI_ALLGATHERV(prompt_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
	!deallocate(ircnt); deallocate(idisp); deallocate(prompt_bank) 
	!call move_alloc(temp_bank, prompt_bank)
	!
	!
	!!> Gather delayed_bank from the slave nodes =================================        
	!ndata = size(delayed_bank)
	!allocate(ircnt(1:ncore))
	!allocate(idisp(1:ncore))
	!do i = 1, ncore
	!	idata = ndata
	!	call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
	!	ircnt(i) = idata
	!enddo 
	!idisp(1) = 0
	!do i = 2, ncore 
	!	idisp(i) = idisp(i-1) + ircnt(i-1)
	!enddo
	!allocate(temp_bank(1:sum(ircnt)))
	!call MPI_ALLGATHERV(delayed_bank,ndata,MPI_bank,temp_bank,ircnt,idisp,MPI_bank,MPI_COMM_WORLD,ierr)
	!deallocate(ircnt); deallocate(idisp); deallocate(delayed_bank) 
	!call move_alloc(temp_bank, delayed_bank)
	
    !call MPI_REDUCE(PCQS_keff,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    !PCQS_keff = rcv_buf / w_tot
	!keff = PCQS_keff 
    !call MPI_BCAST(keff, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
	
	!if (icore==score) then 
	!	rcv_buf = (sum(prompt_bank(:)%wgt) + sum(delayed_bank(:)%wgt))/w_tot
	!	print '(5E15.5)', rcv_buf, PCQS_keff, rcv_buf - PCQS_keff, (sum(delayed_bank(:)%wgt))/w_tot, (sum(delayed_bank(:)%wgt))/w_tot-PCQS_keff
	!endif 
	!
	!deallocate(prompt_bank, delayed_bank)
	
end subroutine PCQS_INIT

!===============================================================================
! BANK_INITIALIZE - Initialize bank for the very first simulation 
!===============================================================================
subroutine bank_initialize(this)
    use FMFD_HEADER, only: fcr, zigzagon, zz_div, zzc0, zzc1, zzc2, dfm, ncm
    use MATERIAL_HEADER, only: Material_CE
    use GEOMETRY_HEADER, only: universes
    use GEOMETRY, only: cell_xyz
    use TALLY, only: CM_ID
    use PCMFD, only: OUT_OF_ZZ
    implicit none
    class(bank)    :: this(:)
    type(surface), pointer :: surfptr
    integer        :: i, j, cell_idx, mat_idx, univ_idx, iso
    integer        :: n_hist, i_hist, i_cell
    real(8)        :: min(3), max(3), xyz(3), e1, e2
    logical        :: found
    integer        :: xy(2), m, n
    integer        :: univ_id, cell_id, mat_id
    integer        :: id(3), m0, n0
    
    
    if(icore==score) print *, '   Initializing Source...'
    
	
    if (.not. allocated(sgrid)) then 
        print *, 'not allocated'
        do i = 1, size(cells(i_cell)%pos_surf_idx)
            surfptr => surfaces(cells(i_cell)%pos_surf_idx(i))
            if (surfptr%surf_type == pz) then 
                max(3) = surfptr%parmtrs(1)
            elseif (surfptr%surf_type == sqcz) then 
                min(1) = -surfptr%parmtrs(3)
                min(2) = -surfptr%parmtrs(3)
                max(1) =  surfptr%parmtrs(3)
                max(2) =  surfptr%parmtrs(3)
            elseif (surfptr%surf_type == cylz) then
                min(1) = -surfptr%parmtrs(3)
                min(2) = -surfptr%parmtrs(3)
                max(1) =  surfptr%parmtrs(3)
                max(2) =  surfptr%parmtrs(3)
            endif 
        enddo 
    else 
        min(1:3) = sgrid(1:3)
        max(1:3) = sgrid(4:6)
    endif 
    
    univ_idx = 0
    do i = 1, size(this) 
        this(i) % wgt = 1
        found         = .false.
        this(i) % uvw = rand_vec()

!        ! Initialize for Latent
!        if(do_ifp .and. latent > 0) then
!            allocate( this(i) % delayedarr(1:latent) )
!            allocate( this(i) % delayedlam(1:latent) )
!            allocate( this(i) %   nlifearr(1:latent) )
!        endif
        
        ! multigroup MC
        if (E_mode == 0) then
            this(i) % G = 1
            search_MG: do while (found == .false.) 
                do j = 1, 3
                    this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
                enddo
                univ_idx = 0; xyz = this(i)%xyz
                call find_cell_xyz(xyz, univ_idx, cell_idx)
                if (cells(cell_idx)%mat_idx == 0) then 
                    found = .false. 
                    cycle search_MG
                endif
                mat_idx = cells(cell_idx)%mat_idx
                !print '(I3, I3, A20)', i, cell_idx, XS_MG(mat_idx)%mat_id
                do j = 1, size(XS_MG(mat_idx)%sig_fis)  
                    if(XS_MG(mat_idx)%sig_fis(j) > 0.0001) then 
                        found = .true.  
                        exit search_MG
                    else 
                        found = .false.
                        cycle search_MG
                    endif 
                enddo 
                
            enddo search_MG
        

        ! continuous energy MC
        elseif (E_mode == 1) then
                    
            e1 = rang()*32.d0 + 1.d0
            e2 = rang()
            this(i)%E = fes(e1) + e2*(fes(e1+1)-fes(e1))
            
			
			if (do_gmsh) then 
				search_CE_tet: do while (found == .false.) 
					do j = 1, 3
						this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
					enddo
					xyz = this(i)%xyz
					found = .true.
				enddo search_CE_tet
				
			else
                if(fmfdon) then
                search_CMFD: do while ( found == .false. ) 
                    do j = 1, 3
                        this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
                    enddo

                    if ( zigzagon ) then
                    id = CM_ID(this(i)%xyz(:))
                    if ( .not. OUT_OF_ZZ(id(1),id(2)) &
                        .and. id(1) > 0 .and. id(1) <= ncm(1) ) exit search_CMFD
                    else
                        exit search_CMFD
                    end if
                    
                enddo search_CMFD
                else
				search_CE: do while ( found == .false.) 
					do j = 1, 3
						this(i) % xyz(j) = rang()*(max(j)-min(j)) + min(j)
					enddo
					univ_idx = 0; xyz = this(i)%xyz
                    !print *, i, xyz(:)
					call find_cell_xyz(xyz, univ_idx, cell_idx)
					! within a material region
					if (cells(cell_idx)%mat_idx == 0) then 
						found = .false. 
						cycle search_CE
					endif
					mat_idx = cells(cell_idx)%mat_idx
					! within a fissionable material region
                    found = materials(mat_idx)%fissionable
					
				enddo search_CE 
                !print *, 'DONE', icore, i, '/', size(this), this(i)%xyz(1:3)
                endif
            endif
            
        endif
        
    enddo 
    
    if(icore==score) print *, '   Source Initializing Complete...'
    
	 
end subroutine


	
	
  
	subroutine draw_geometry() 
		integer :: i, j, i_plot
		integer :: univ_idx, cell_idx
		real(8) :: x, y, z, xyz(3) 
		integer :: mat_rgb(3, n_materials), rgb_blk(3) 
        !integer :: mat_rgb(n_materials, 3), rgb_blk(3) 
		real(8) :: rn 
		type(rgbimage) :: im
		character(24) :: plottitle
		
		
		if (.not. plotgeom) return 
		if (icore /=score) return 
		
		!> make material color map 
		do i = 1, n_materials
            if(sum(materials(i) % rgb) < 0) then
    			do j = 1, 3
    				call random_number(rn) 
    				mat_rgb(j,i) = floor(rn*255)
    			enddo 
            else ! RGB assigned
                mat_rgb(1:3,i) = materials(i) % rgb(1:3)
            endif
		enddo 
		rgb_blk(1:3) = 0
		
		
		
		do i_plot = 1, n_plot
		!> make bitmap 
		select case (plottype(i_plot)) 
		case(1) 
			call im%init(plt_nx(i_plot), plt_ny(i_plot))
			do i = 1, plt_nx(i_plot)
				do j = 1, plt_ny(i_plot)
					xyz(1) = plt_x0(i_plot) + (i-0.5)*plt_dx(i_plot)
					xyz(2) = plt_y0(i_plot) + (j-0.5)*plt_dy(i_plot)
					xyz(3) = plt_z0(i_plot)
					univ_idx = 0
					call find_cell_xyz(xyz, univ_idx, cell_idx)
                    !print *, xyz(1),xyz(2),univ_idx,cell_idx
					if (cells(cell_idx)%mat_idx < 1) then 
						call im%set_pixel(i, j, rgb_blk)
					else 
						call im%set_pixel(i, j, mat_rgb(1:3, cells(cell_idx)%mat_idx))
					endif 
				enddo 
			enddo 
			
		case(2) 
			call im%init(plt_ny(i_plot), plt_nz(i_plot))
			do i = 1, plt_ny(i_plot)
				do j = 1, plt_nz(i_plot)
					xyz(1) = plt_x0(i_plot)
					xyz(2) = plt_y0(i_plot) + (i-0.5)*plt_dy(i_plot)
					xyz(3) = plt_z0(i_plot) + (j-0.5)*plt_dz(i_plot)
					univ_idx = 0
					call find_cell_xyz(xyz, univ_idx, cell_idx)
					if (cells(cell_idx)%mat_idx < 1) then 
						call im%set_pixel(i, j, rgb_blk)
					else 
						call im%set_pixel(i, j, mat_rgb(1:3, cells(cell_idx)%mat_idx))
					endif 
				enddo 
			enddo 
		case(3) 
			call im%init(plt_nx(i_plot), plt_nz(i_plot))
			do i = 1, plt_nx(i_plot)
				do j = 1, plt_nz(i_plot)
					xyz(1) = plt_x0(i_plot) + (i-0.5)*plt_dx(i_plot)
					xyz(2) = plt_y0(i_plot) 
					xyz(3) = plt_z0(i_plot) + (j-0.5)*plt_dz(i_plot)
					univ_idx = 0
					call find_cell_xyz(xyz, univ_idx, cell_idx)
					if (cells(cell_idx)%mat_idx < 1) then 
						call im%set_pixel(i, j, rgb_blk)
					else 
						call im%set_pixel(i, j, mat_rgb(1:3, cells(cell_idx)%mat_idx))
					endif 
				enddo 
			enddo 
		
		end select 
		
		plottitle = adjustl(trim(plotlist(i_plot))) // ".ppm"
		
		write(*,'(a30,i1,a3,i1,a1)', advance='no')  "    Making ppm image file... (" , i_plot, " / ", n_plot, ")"
		! output image into file 
		call im%write(plottitle)
		!call system("python /home/guest/HyeonTae/src_temp/15_Test/src_merge2/convert_ppm_to_jpg.py")
		!call system("mv fig.jpg "//trim(plotlist(i_plot))//".jpg") 
		write(*,*) " - ", trim(plotlist(i_plot)), ".jpg"
		call system("python /home/guest/Inyup/convert_ppm_to_jpg.py "//trim(plottitle))

		enddo 
		!call system("python /home/guest/Inyup/convert_ppm_to_jpg.py")
		
!		deallocate(plt_x0, plt_y0, plt_z0, plt_x1, plt_y1, plt_z1, &
!					plt_dx, plt_dy, plt_dz, plt_nx, plt_ny, plt_nz, &
!					plotlist, plottype) 
					
					
	end subroutine
  
  	


    
!subroutine para_range(n1, n2, nprocs, irank, ista, iend)
!    integer :: iwork1, iwork2 
!    integer, intent(in) :: n1, n2, nprocs, irank 
!    integer, intent(inout) :: ista, iend
!    
!    iwork1 = (n2 - n1 + 1) / nprocs
!    iwork2 = MOD(n2 - n1 + 1, nprocs)
!    ista = irank * iwork1 + n1 + MIN(irank, iwork2)
!    iend = ista + iwork1 - 1
!    if (iwork2 > irank) iend = iend + 1
!
!end subroutine 

subroutine gatherSourceBank(bank_main, bank_thread, bank_idx)
	type(bank), allocatable, intent(inout) :: bank_main(:)
	type(bank), intent(inout) :: bank_thread(:) 
	type(bank), allocatable :: bank_temp(:) 
	integer, intent(inout) :: bank_idx
	integer :: isize
	real(8) :: time1, time2, time3, time4 
	
	time1 = omp_get_wtime()
	isize = size(bank_main)
	if(allocated(bank_temp)) deallocate(bank_temp)
	allocate(bank_temp(1:isize+bank_idx)) 
	time2 = omp_get_wtime()
	if (isize>0) bank_temp(1:isize) = bank_main(:)
	deallocate(bank_main)
	time3 = omp_get_wtime()
	bank_temp(isize+1:isize+bank_idx) = bank_thread(1:bank_idx)
	call move_alloc(bank_temp, bank_main)
	time4 = omp_get_wtime()
	
	!if (icore==score) print '(a, 3F12.5)', 'test', time2-time1, time3-time2, time4-time3
	
	
end subroutine

subroutine gatherPrecBank(bank_main, bank_thread, bank_idx)
	type(PrecBank), allocatable, intent(inout) :: bank_main(:)
	type(PrecBank), intent(inout) :: bank_thread(:) 
	type(PrecBank), allocatable :: bank_temp(:) 
	integer, intent(inout) :: bank_idx
	integer :: isize
	
	isize = size(bank_main)
	if(allocated(bank_temp)) deallocate(bank_temp)
	allocate(bank_temp(1:isize+bank_idx)) 
	if (isize>0) bank_temp(1:isize) = bank_main(:)
	deallocate(bank_main)
	bank_temp(isize+1:isize+bank_idx) = bank_thread(1:bank_idx)
	call move_alloc(bank_temp, bank_main)
end subroutine

subroutine gatherVRCBank(bank_main, bank_thread, bank_idx)
	type(VRCBank), allocatable, intent(inout) :: bank_main(:)
	type(VRCBank), intent(inout) :: bank_thread(:) 
	type(VRCBank), allocatable :: bank_temp(:) 
	integer, intent(inout) :: bank_idx
	integer :: isize
	
	isize = size(bank_main)
	if(allocated(bank_temp)) deallocate(bank_temp)
	allocate(bank_temp(1:isize+bank_idx)) 
	if (isize>0) bank_temp(1:isize) = bank_main(:)
	deallocate(bank_main)
	bank_temp(isize+1:isize+bank_idx) = bank_thread(1:bank_idx)
	call move_alloc(bank_temp, bank_main)
end subroutine



end module
