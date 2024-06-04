module simulation
    use FMFD,               only : MCBU
    use tracking,           only : transport, transport_dynamic, &
                                   transport_pcqs, transport_PCQS_init
    use simulation_header, only: t_det
    use GEOMETRY, only : find_cell_xyz
    use geometry_header, only : sgrid
    use surface_header, only : surface
    use ENTROPY, only : entrpon, mprupon, genup, up_sign, rampup, SHENTROPY
    use MPRUP, only : GENSIZE, CYCLECHANGE
    use tally, only : tallyon, findtallybin, tallyflux, tallypower, &
                    tallycoord, tally_buf, mesh_power, PROCESS_TALLY, &
                    NORM_TALLY, TALLY_THREAD_INITIAL, MC_tally
    use FMFD, only : fmfdon, n_skip, fake_MC, PROCESS_FMFD, FMFD_SOLVE, &
                    FMFD_INITIALIZE, FMFD_INITIALIZE_THREAD, NORM_FMFD
    use TEMPERATURE, only : TEMP_SOLVE, NORM_TH, PROCESS_TH, &
                                   TEMP_CONVERGE, TEMP_DISTRIBUTE
    use TH_HEADER, only : th_on
    use tetrahedral, only : find_tet, transport_tet, tet, find_tet_old
    use DMC, only : SET_DYNAMIC_BANK
    use BANK_HEADER
    use CONSTANTS, only: fes, prt_flux, prt_powr, prt_wgt, pz, sqcz, cylz
    use PCQS
    use PARTICLE_HEADER, only: particle
    use rgbimage_m
   
    implicit none 
    
    interface gatherBank
        module procedure :: gatherSourceBank
        module procedure :: gatherPrecBank
        module procedure :: gatherVRCBank
    end interface 
    
    contains 
    
subroutine simulate_history(bat,cyc)
    use FMFD_HEADER, only: n_fake, dual_fmfd, fwgt
    implicit none
    integer, intent(in):: bat
    integer, intent(in):: cyc
    integer :: i, j, k, isize, i_surf
    type(particle) :: p
    integer :: tid, my_id
    integer :: ista, iend
    real(8) :: Jtemp
    real(8), allocatable :: shape(:)
    real(8) :: rcv_buf
    integer :: realex, intex, restype, ndata, idata
    integer, dimension(0:4) :: blocklength, displacement, oldtype 
    integer, allocatable :: ircnt(:), idisp(:), bsize(:)
    integer:: bank_all, dsize, bank_size0
    integer:: i_bin(4)
    integer:: icore0(1), icore1(1)
    real(8) :: t_fm1, t_fm2

    if (allocated(fission_bank)) call move_alloc(fission_bank, source_bank)
    if ( entrpon ) call SHENTROPY(source_bank,cyc) ! need to be revised
        !call FET_CALC(source_bank)
    if ( iscore .and. mprupon ) call GENSIZE(bat,cyc)
    if ( mprupon .or. ( .not. mprupon .and. genup ) ) then
        call MPI_BCAST(genup,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(mprupon,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(up_sign,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
    end if
    if ( curr_cyc /= 1 .and. up_sign .and. mprupon ) then
        call MPI_BCAST(ngen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        source_bank(:)%wgt = real(ngen,8)/real(ngen-rampup,8)
    end if
    if ( .not. mprupon .and. genup ) call CYCLECHANGE(cyc)
    allocate(fission_bank(0))

    !> Distribute source_bank to slave nodes 
    k_col = 0; k_tl = 0; k_vrc = 0; fiss_vrc = 0; loss_vrc = 0;
    if ( fmfdon ) call FMFD_initialize()
    if ( do_burn ) MC_tally(bat,:,:,:,:,:,:) = 0
    cyc_power = 0;

    n_col = 0; n_cross = 0
    
    if (allocated(prec_bank)) deallocate(prec_bank)
    allocate(prec_bank(0))
    if (allocated(prompt_bank)) deallocate(prompt_bank)
    allocate(prompt_bank(0))

    !$omp parallel private(p) shared(source_bank, fission_bank, temp_bank, prec_bank)
      thread_bank(:)%wgt = 0; bank_idx = 0; prec_idx = 0 ; init_idx = 0
      if ( tallyon .and. .not. fmfdon ) call TALLY_THREAD_INITIAL(cyc)
      if ( fmfdon ) call FMFD_initialize_thread()
      !$omp do reduction(+:k_col, k_tl)
        do i= 1, bank_size
            call p%initialize()
            call p%set(source_bank(i))
            
            ! initialize p%tet value (TODO)
            if (do_gmsh .and. curr_cyc > n_inact) then 
                i_bin = FindTallyBin(p)
                if ( i_bin(1) > 0 ) then
                    p%tet = find_tet(p%coord(1)%xyz)
                    p%in_tet = .true.
                else
                    p%in_tet = .false.
                end if
                if (p%tet .le. 0) p%in_tet = .false.
            endif 
            
            do while (p%alive == .true.)
                call transport(p)
            enddo 
            
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
    
    if ( iscore ) avg_power = avg_power + cyc_power


    call MPI_REDUCE(n_col,i,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
    n_col = i
    n_col_avg = n_col_avg + n_col

    !> Calculate k_eff ==========================================================
    k_col = k_col / real(ngen,8)
    k_tl  = k_tl  / real(ngen,8) 
    keff  = (k_tl + k_col) / 2.0d0
    !keff = k_col
    
    !if (icore == score) write(prt_keff,*) keff, k_col, k_tl
    
    call MPI_BCAST(keff, 1, MPI_DOUBLE_PRECISION, score, MPI_COMM_WORLD, ierr) 
    
    !> fission neutron source sharing
    bank_size = size(fission_bank)
    if ( curr_cyc /= 1 ) then
    allocate(bsize(ncore))
    do i = 1, ncore
        bank_size0 = bank_size
        call MPI_BCAST(bank_size0,1,MPI_INTEGER,i-1,MPI_COMM_WORLD,ierr)
        bsize(i) = bank_size0
    end do

    if ( icore == score ) then
        icore0 = minloc(bsize(1:ncore))
        icore1 = maxloc(bsize(1:ncore))
        dsize = (bsize(icore1(1))-bsize(icore0(1)))/2
    end if

    call MPI_BCAST(icore0,1,MPI_INTEGER,score,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(icore1,1,MPI_INTEGER,score,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(dsize,1,MPI_INTEGER,score,MPI_COMM_WORLD,ierr)

    ! storing extra sources
    if ( dsize /= 0 ) then
    allocate(temp_bank(bsize(icore1(1))))
    if ( icore == icore1(1)-1 ) temp_bank = fission_bank
    call MPI_BCAST(temp_bank(bsize(icore1(1))-dsize+1:),dsize, &
        MPI_BANK,icore1(1)-1,MPI_COMM_WORLD,ierr)

    ! largest fission source
    if ( icore == icore1(1)-1 ) then
        deallocate(fission_bank)
        bank_size = bank_size-dsize
        allocate(fission_bank(1:bank_size))
        fission_bank(:) = temp_bank(1:bank_size)
    end if

    ! smallest fission source
    if ( icore == icore0(1)-1 ) then
        temp_bank(1:bank_size) = fission_bank(1:bank_size)
        deallocate(fission_bank)
        bank_size = bank_size+dsize
        allocate(fission_bank(1:bank_size))
        fission_bank(1:bsize(icore0(1))) = temp_bank(1:bsize(icore0(1)))
        fission_bank(bsize(icore0(1))+1:) = temp_bank(bsize(icore1(1))-dsize+1:)
    end if
    deallocate(temp_bank,bsize)
    end if
    end if

    ! DMC
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

    ! neutron weight normalization
    call MPI_ALLREDUCE(bank_size,bank_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    fission_bank(:)%wgt  = real(ngen,8)/real(bank_all,8)
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
    !call inline_xenon()
    
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
            !    TallyPower(:) = TallyPower(:) * Nominal_Power * (mesh_power / cyc_power) * 1.0d6 / mesh_power
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
            if (E_mode == 1 .and. p%E == 0)    p%alive = .false.
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
            !    !$omp critical
            !    call gatherBank(vrc_bank, vrc_thread, vrc_idx)
            !    !$omp end critical
            !    vrc_idx = 0
            !    print *, "vrc_idx exceed its limit", vrc_idx 
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
        !    PCQS_beta_cyc  (curr_cyc - n_pcqs_inact,i) = PKE_tally(i) / PKE_tally(2*npg+1)
        !    PCQS_lambda_cyc(curr_cyc - n_pcqs_inact,i) = PKE_tally(i) / PKE_tally(npg+i)
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
    !    ndata = n_pcqs_act 
    !    do i = 1, npg 
    !    avg = sum(PCQS_beta_cyc(:,i))/ndata
    !    std = sqrt(dot_product((PCQS_beta_cyc(:,i)-avg),(PCQS_beta_cyc(:,i)-avg))/(ndata*(ndata-1)))
    !    print *, 'beta', i, avg, std 
    !    enddo 
    !    do i = 1, npg 
    !    avg = sum(PCQS_lambda_cyc(:,i))/ndata
    !    std = sqrt(dot_product((PCQS_lambda_cyc(:,i)-avg),(PCQS_lambda_cyc(:,i)-avg))/(ndata*(ndata-1)))
    !    print *, 'lambda', i, avg, std 
    !    enddo 
    !    avg = sum(PCQS_gen_cyc(:))/ndata
    !    std = sqrt(dot_product((PCQS_gen_cyc(:)-avg),(PCQS_gen_cyc(:)-avg))/(ndata*(ndata-1)))
    !    print *, 'gen', avg, std 
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
    !    if (curr_timestep < 10) then 
    !        write (filename, "(A22,I1,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
    !    elseif (i < 100) then 
    !        write (filename, "(A22,I2,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
    !    elseif (i < 1000) then 
    !        write (filename, "(A22,I3,A4)") "./tet_vrc/data/tet_vrc", curr_timestep, ".dat"
    !    endif 
    !
    !    open(prt_tet_vrc, file=trim(filename), status="old", position="append", action="write")
    !
    !    ndata = size(vrc_bank)
    !    write(prt_tet_vrc,'(2i, e16.7)') ndata, curr_act, mesh_power
    !    do i = 1, ndata
    !        write(prt_tet_vrc,'(8e16.7,I,L)') vrc_bank(i)%xyz, vrc_bank(i)%uvw, vrc_bank(i)%E &
    !                                         ,vrc_bank(i)%wgt, vrc_bank(i)%G, vrc_bank(i)%delayed
    !    enddo
    !endif 
    !close(prt_tet_vrc)
    
    
    !!> Normalize tally (flux & power) ===========================================
    !if (tally_switch > 0 .and. do_transient == .true.) then 
    !    isize = size(TallyFlux) 
    !    call MPI_REDUCE(TallyFlux, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
    !    TallyFlux = tally_buf
    !    call MPI_REDUCE(TallyPower, tally_buf, isize, MPI_DOUBLE_PRECISION, MPI_SUM, score, MPI_COMM_WORLD, ierr)
    !    TallyPower = tally_buf
    !    
    !    
    !    if (icore == score) then 
    !        TallyFlux(:) = TallyFlux(:) * dble(isize) / sum(TallyFlux)
    !        TallyPower(:) = TallyPower(:) * Nominal_Power * 1.0d6 / cyc_power
    !        
    !        if (do_gmsh) then
    !            do i = 1, isize
    !                TallyFlux(i)  = TallyFlux(i)/tet(i)%vol
    !                TallyPower(i) = TallyPower(i)/tet(i)%vol
    !            enddo
    !        else 
    !            do i = 1, isize
    !                TallyFlux(i)  = TallyFlux(i)/TallyCoord(i)%vol
    !                TallyPower(i) = TallyPower(i)/TallyCoord(i)%vol
    !            enddo
    !        endif 
    !        
    !        !> Print TallyFlux
    !        write(prt_flux, 100) TallyFlux(:)
    !        write(prt_powr, 100) TallyPower(:)
    !    100 format(<isize>ES15.7)
    !    
    !
    !    endif 
    !    
    !    TallyFlux(:) =0
    !    TallyPower(:)=0
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
    !    call move_alloc(prompt_bank, temp_bank)
    !    call combing_source(temp_bank, prompt_bank, bank_max) 
    !    deallocate(temp_bank) 
    !    ndata = size(prompt_bank)
    !endif 
    !
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
    !    idata = ndata
    !    call MPI_BCAST(idata, 1, MPI_INTEGER, i-1, MPI_COMM_WORLD, ierr) 
    !    ircnt(i) = idata
    !enddo 
    !idisp(1) = 0
    !do i = 2, ncore 
    !    idisp(i) = idisp(i-1) + ircnt(i-1)
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
    !    rcv_buf = (sum(prompt_bank(:)%wgt) + sum(delayed_bank(:)%wgt))/w_tot
    !    print '(5E15.5)', rcv_buf, PCQS_keff, rcv_buf - PCQS_keff, (sum(delayed_bank(:)%wgt))/w_tot, (sum(delayed_bank(:)%wgt))/w_tot-PCQS_keff
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
    integer        :: i0, j0, cell_idx, mat_idx, univ_idx, iso
    integer        :: n_hist, i_hist, i_cell
    real(8)        :: min(3), max(3), xyz(3), e1, e2
    logical        :: found
    integer        :: xy(2), m, n
    integer        :: univ_id, cell_id, mat_id, n_univ
    integer        :: id(3), m0, n0
    
    
    if(icore==score) print *, '   Initializing Source...'
    
    if (.not. allocated(sgrid)) then 
        print *, 'not allocated'
        do i0 = 1, size(cells(i_cell)%pos_surf_idx)
            surfptr => surfaces(cells(i_cell)%pos_surf_idx(i0))
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
    n_univ = universes(univ_idx)%ncell
    do i0 = 1, bank_size
        this(i0) % wgt = 1
        found         = .false.
        this(i0) % uvw = rand_vec()
        
        ! multigroup MC
        if (E_mode == 0) then
            this(i0) % G = 1
            search_MG: do while (found == .false.) 
                do j0 = 1, 3
                    this(i0) % xyz(j0) = rang()*(max(j0)-min(j0)) + min(j0)
                enddo
                call find_cell_xyz(this(i0)%xyz, univ_idx, cell_idx)
                if (cells(cell_idx)%mat_idx == 0) then 
                    found = .false. 
                    cycle search_MG
                endif
                mat_idx = cells(cell_idx)%mat_idx
                !print '(I3, I3, A20)', i, cell_idx, XS_MG(mat_idx)%mat_id
                do j0 = 1, size(XS_MG(mat_idx)%sig_fis)  
                    if(XS_MG(mat_idx)%sig_fis(j0) > 0.0001) then 
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
            this(i0)%E = fes(e1) + e2*(fes(e1+1)-fes(e1))
            
            
            if (do_gmsh) then 
                search_CE_tet: do while (found == .false.) 
                    do j0 = 1, 3
                        this(i0) % xyz(j0) = rang()*(max(j0)-min(j0)) + min(j0)
                    enddo
                    xyz = this(i0)%xyz
                    found = .true.
                enddo search_CE_tet
                
            else 
                if ( fmfdon ) then
                search_CMFD: do while ( found == .false. ) 
                    do j0 = 1, 3
                        this(i0) % xyz(j0) = rang()*(max(j0)-min(j0)) + min(j0)
                    enddo

                    if ( zigzagon ) then
                    id = CM_ID(this(i0)%xyz(:))
                    !print'(2i5,3x,2f15.7)', id(1:2), this(i0)%xyz(1:2)
                    !if ( id(1) == 0 ) print*, OUT_OF_ZZ(id(1),id(2))
                    if ( .not. OUT_OF_ZZ(id(1),id(2)) &
                        .and. id(1) > 0 .and. id(1) <= ncm(1) ) exit search_CMFD
                    else
                    exit search_CMFD
                    end if
                    
                enddo search_CMFD
                if ( this(i0)%xyz(1) < -74 .and. this(i0)%xyz(2) < -74 ) then
                    id = CM_ID(this(i0)%xyz(:))
                    print*, "id", id(1:2)
                    print*, this(i0)%xyz(1:2), OUT_OF_ZZ(id(1),id(2))
                end if
                if ( this(i0)%xyz(1) >  74 .and. this(i0)%xyz(2) >  74 ) then
                    id = CM_ID(this(i0)%xyz(:))
                    print*, "id", id(1:2)
                    print*, this(i0)%xyz(1:2), OUT_OF_ZZ(id(1),id(2))
                end if


                else
                search_CE: do while ( found == .false. ) 
                    do j0 = 1, 3
                        this(i0) % xyz(j0) = rang()*(max(j0)-min(j0)) + min(j0)
                    enddo

                    ! within a universe
                    do j0 = 1, n_univ
                        cell_idx = universes(univ_idx)%cell(j0)
                        if ( cell_xyz(cells(cell_idx),this(i0)%xyz) ) then
                        if ( cells(cell_idx)%mat_idx /= 0 ) then
                            found = .true.
                            exit
                        end if
                        end if
                    end do

!                    call find_cell_xyz(this(i0)%xyz, univ_idx, cell_idx)
!                    ! within a material region
!                    if (cells(cell_idx)%mat_idx == 0) then 
!                        found = .false. 
!                        cycle search_CE
!                    endif
!                    mat_idx = cells(cell_idx)%mat_idx
!                    ! within a fissionable material region
!                    do j0 = 1, materials(mat_idx)%n_iso
!                        iso = materials(mat_idx)%ace_idx(j0)
!                        if(ace(iso)%jxs(2) /= 0) then 
!                            found = .true.  
!                        endif 
!                    enddo 
                    
                enddo search_CE 
                end if
            endif
            
        endif
        
    enddo 
    
    if(icore==score) print *, '   Source Initializing Complete...'

end subroutine


  
subroutine draw_geometry() 
    integer :: i, j, i_plot
    integer :: univ_idx, cell_idx
    real(8) :: x, y, z, xyz(3) 
    integer :: mat_rgb(n_materials, 3), rgb_blk(3) 
    real(8) :: rn 
    type(rgbimage) :: im
    character(24) :: plottitle
    
    
    if (.not. plotgeom) return 
    if (icore /=score) return 
    
    !> make material color map 
    do i = 1, n_materials
        do j = 1, 3
            call random_number(rn) 
            mat_rgb(i,j) = floor(rn*255)
        enddo 
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
                if (cells(cell_idx)%mat_idx < 1) then 
                    call im%set_pixel(i, j, rgb_blk)
                else 
                    call im%set_pixel(i, j, mat_rgb(cells(cell_idx)%mat_idx,:))
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
                    call im%set_pixel(i, j, mat_rgb(cells(cell_idx)%mat_idx,:))
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
                    call im%set_pixel(i, j, mat_rgb(cells(cell_idx)%mat_idx,:))
                endif 
            enddo 
        enddo 
    
    end select 
    
    plottitle = plotlist(i_plot) // ".ppm"
    
    write(*,'(a30,i1,a3,i1,a1)', advance='no')  "    Making ppm image file... (" , i_plot, " / ", n_plot, ")"
    ! output image into file 
    call im%write("fig.ppm")
    call system("python /home/guest/HyeonTae/src_temp/15_Test/src_merge2/convert_ppm_to_jpg.py")
    !call system("python ./convert_ppm_to_jpg.py")
    call system("mv fig.jpg "//trim(plotlist(i_plot))//".jpg") 
    write(*,*) " - ", trim(plotlist(i_plot)), ".jpg"

    enddo 
    
    deallocate(plt_x0, plt_y0, plt_z0, plt_x1, plt_y1, plt_z1, &
                plt_dx, plt_dy, plt_dz, plt_nx, plt_ny, plt_nz, &
                plotlist, plottype) 
                
                
end subroutine
  
      


    
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
    
subroutine para_range2(n1, n2, nprocs, irank, bsize)
    integer, intent(in) :: n1, n2, nprocs, irank 
    integer, intent(inout) :: bsize
    integer :: iwork
    
    bsize = (n2 - n1 + 1) / nprocs
    iwork = MOD(n2 - n1 + 1, nprocs)
    if (iwork > irank) bsize = bsize + 1

end subroutine 

subroutine gatherSourceBank(bank_main, bank_thread, bank_idx)
    type(bank), allocatable, intent(inout) :: bank_main(:)
    type(bank), intent(inout) :: bank_thread(:) 
    type(bank), allocatable :: bank_temp(:) 
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
