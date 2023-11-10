program main
use CONSTANTS, only: prt_dynamic, prt_tet_vrc, prt_prec, prt_delayed, prt_keff
use ENTROPY,    only : mprupon, bprupon, entrp0, ENTRP_INIT2
use TH_HEADER
use simulation 
use DEPLETION_MODULE, only: nstep_burnup, MPI_REDUCE_BURNUP, INIT_BURNUP, &
                        DEPLETION
use mpi
use transient
use TEMPERATURE, only: TEMP_DISTRIBUTE, POWER_NORM_START, POWER_IMC_TO_START
use TALLY, only: k_eff, TallyFlux
use STATISTICS
use PRINTER
use FMFD, only: DET_POWER, INTRA_PIN_DTMC
use FMFD_HEADER, only: acc_skip, n_skip
use communication

!use ace_header, only: udelta, Emin, ugrid, nugrid

implicit none

integer :: i, j
integer :: provide 
real(8) :: k_sum 
logical :: isopened
real(8),allocatable :: tally_val(:,:)
integer :: nsize
real(8), allocatable :: ttemp(:,:,:)
character(100) :: filename
real(8) :: kavg, kstd
real(8) :: time1, time2, time3, time4
real(8) :: err
integer :: ix, iy, iz, tmp

!> Preparation for parallelization ===============================================
!call omp_set_num_threads(13)
call MPI_Init_thread(MPI_THREAD_SINGLE, provide, ierr)
core = MPI_COMM_WORLD
call MPI_COMM_RANK(core,icore,ierr)
call MPI_COMM_SIZE(core,ncore,ierr)

!if ( icore == score ) call INPUT_GEN

!> PreMC : Read input / Initialize / Set Random Seed etc. ========================
call premc
if(do_child) then
    call INIT_CHILD

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(t_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(t_comm_fuel, nth(1)*nth(2)*nth(3), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rho_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
    print *, 'BCAST WELL', sum(t_comm_cool)/dble(n_channels*(nth(3)+1)), sum(t_comm_fuel)/dble(nth(1)*nth(2)*nth(3)), sum(rho_comm_cool), icore

    call TH_ASSIGN_GRID

end if

if ( icore == score ) call TIME_MEASURE

!> Stead-state Simlulation Start =================================================
call START_MSG

curr_bat = 0
BATCH : do

    if ( n_batch == 1 ) curr_bat = 1
    if ( n_batch > 1 ) call BATCH_MSG(curr_bat)


BURNUP : do
    
    if ( do_burn ) then
        if ( iscore ) call BURNUP_MSG
        call INIT_BURNUP
        !open(990,file='COS.out',action='write',status='replace')
    end if

    ! for transient calculation
    if ( tally_switch > 0 .and. iscore .and. .not. do_transient ) then
        open(prt_flux,file='flux.out',action="write",status="replace")
        open(prt_powr,file='power.out',action="write",status="replace")
    end if

    ! for DMC calculation
    if ( do_DMC ) then
        open(prt_dynamic,file="dynamicMC.out",action="write",status="replace")
        open(prt_wgt,file="tetstop.out",action="write",status="replace")
    end if
    
    ! for GMESH calculation
    if (do_gmsh_vrc) then 
        do i = 1, n_timestep
            if (i < 10) then 
                write (filename, "(A22,I1,A4)") "./tet_vrc/data/tet_vrc", i, ".dat"
            elseif (i < 100) then 
                write (filename, "(A22,I2,A4)") "./tet_vrc/data/tet_vrc", i, ".dat"
            elseif (i < 1000) then 
                write (filename, "(A22,I3,A4)") "./tet_vrc/data/tet_vrc", i, ".dat"
            else 
                print *, "ERROR :: TOO MANY ACTIVE CYCLES FOR GMSH_VRC - ", n_timestep
                stop
            endif 
            open(prt_tet_vrc, file=trim(filename), action="write",status="replace")
            close(prt_tet_vrc)
        enddo 
    endif 
    
    ! for DMC calculation
    if (do_DMC .and. tally_switch > 0 .and. iscore) then 
        call system("mkdir -p ./DMC_data/") 
        do i = 1, n_timestep
            if (i < 10) then 
                write (filename, "(A22,I1,A4)") "./DMC_data/power", i, ".out"
            elseif (i < 100) then                           
                write (filename, "(A22,I2,A4)") "./DMC_data/power", i, ".out"
            elseif (i < 1000) then                          
                write (filename, "(A22,I3,A4)") "./DMC_data/power", i, ".out"
            else 
                print *, "ERROR :: TOO MANY timesteps for DMC TET - ", n_timestep
                stop
            endif 
            open(prt_powr, file=trim(filename), action="write",status="replace")
            close(prt_powr)
        enddo 
    endif

TH : do
    ! steady-state calculation
    if ( bprupon ) call ENTRP_INIT2
    time3 = omp_get_wtime()
    curr_cyc = 0
    CYC: Do
        curr_cyc = curr_cyc + 1
        curr_act = curr_cyc - n_inact
        !> history wise transport simulation
        time1 = MPI_WTIME();
        call simulate_history(curr_bat,curr_cyc)
        time2 = MPI_WTIME();
        if ( iscore .and. curr_bat /= 0 ) t_tot(curr_bat,curr_cyc) = time2-time1
        call RUN_MSG(curr_bat,time1,time2)

        if ( curr_bat == 0 ) then
        if ( curr_cyc == n_inact ) exit CYC
        else
        if ( curr_cyc <= n_inact ) cycle CYC
        end if

        ! transient calculation
        if ( do_DMC ) then
        curr_time = 0 
        call normalizeInitialSource()
        cyc_power0 = 0
        w_tot = sum(fission_bank(:)%wgt) 
        call save_condition()
        !call save_MG_XS()  ! ---- save MG XS for perturbation
        Do curr_timestep = 1, n_timestep
            !write (prt_dynamic, *) "========== NEW TIME STEP ",curr_timestep," =========="
            !call adjust_MG_XS() ! -------------------------- Adjust MG XS accordingly
            !call adjust_CE_MAT() ! CE mat number density change
            call condition_change()
            time1 = omp_get_wtime()
            call Dynamic_MC()
            curr_time = curr_time + del_t
            time2 = omp_get_wtime()
            if ( iscore ) then
                print '(I3,F9.5,E13.4,A18,F9.1,A)', curr_timestep, DMC_keff, &
                    cyc_power, ' | elapsted time : ', time2-time1, " sec"
                write(prt_dynamic,*),curr_timestep, DMC_keff, cyc_power
            end if
        Enddo 
        !call restore_MG_XS() ! ----- restore MG XS for the next cycle
        call restore_condition()
        call finalize_src()
        end if

        if ( curr_cyc == n_totcyc ) exit

    Enddo CYC
    
    if (do_DMC) then 
        close(prt_dynamic)
        !close(prt_wgt)
    endif 
    
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    
    call MPI_BCAST(t_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(t_comm_fuel, nth(1)*nth(2)*nth(3), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(rho_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)

    !> PCQS Transient MC =============================================================================
    if (.not. do_PCQS) goto 99
    !open(prt_dynamic,file="PCQS.out",action="write",status="replace")
    open(prt_wgt,file="PCQS_power.out",action="write",status="replace")
    open(prt_prec,file="PCQS_wgt.out",action="write",status="replace")
    open(prt_delayed,file="PCQS_corrector.out",action="write",status="replace")
    
    curr_time = 0
    call PKE_init() 
    call PCQS_INIT()
    !PCQS_keff = AVG(k_eff(curr_bat,n_inact+1:n_totcyc))
    call save_condition() 
    curr_timestep = 0 
    if ( icore == score ) print *, 'timestep ', curr_timestep
    Do curr_cyc = 1, n_pcqs_totcyc
        call PCQS_MC() 
    Enddo 
    call solve_PKE() 
    
    Do curr_timestep = 1, n_timestep
        
        if ( icore == score ) then
            print *, '' 
            print '(a10,I)', ' Timestep ', curr_timestep
        endif
        call condition_change()
        cyc_power0 = 0
        
        !> Predictor 
        corrector = .false. 
        Do curr_cyc = 1, n_pcqs_totcyc
            call PCQS_MC() 
        Enddo 
        
        !> Solve PKE (Corrector) 
        call solve_PKE() 
        
        
        kavg = sum(PCQS_keff_cyc(:))/n_pcqs_act
        kstd = sqrt(dot_product((PCQS_keff_cyc(:)-kavg), &
            (PCQS_keff_cyc(:)-kavg))/(n_pcqs_act*(n_pcqs_act-1)))
        
        
        if (icore==score) write(prt_delayed,'(E11.3,2F16.6, F15.2, 3E15.5)') &
                        curr_time, PKE_keff, kavg, kstd*1d5, PKE_amp, &
                        PCQS_power(curr_timestep,1)*PKE_f, &
                        100d0*PCQS_power(curr_timestep,2)/PCQS_power(curr_timestep,1)
        
        curr_time = curr_time + del_t
        
        !> Zero tally bank 
        PKE_beta_tally1 = 0;  PKE_beta_tally2 = 0;
        PKE_lambda_tally1 = 0; PKE_lambda_tally2 = 0;
        PKE_prec_tally1 = 0;  PKE_prec_tally2 = 0;
        PKE_gen_tally1 = 0;  PKE_gen_tally2 = 0;
        PKE_keff_tally = 0; PKE_Z_tally1 = 0
        
    Enddo 
    
    call restore_condition()
    call finalize_src()
    !close(prt_dynamic)
    close(prt_wgt)
    close(prt_prec)
    close(prt_delayed)
    
    !> End of PCQS 
    ! ===============================================================================================
    
99  if (do_gmsh_vrc) close(prt_tet_vrc)

    call CYCLE_TALLY_MSG(curr_bat)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)    
    if ( th_on ) then
        if ( icore == score ) then
        allocate(ttemp(nfm(1),nfm(2),nfm(3)))
        do ii = 1, nfm(1)
        do jj = 1, nfm(2)
        do kk = 1, nfm(3)
            ttemp(ii,jj,kk) = AVG(p_fmfd(curr_bat,:,ii,jj,kk))
        end do
        end do
        end do
        call DET_POWER(ttemp(:,:,:))
        call TEMP_SOLVE
        call TEMP_DISTRIBUTE
        call TEMP_CONVERGE(err)
        deallocate(ttemp)
        end if

    elseif(do_child .and. th_iter < th_iter_max) then !START OPTION
        if(icore==score) print *, 'THITER', th_iter, th_iter_max
        th_iter = th_iter + 1
        if(icore==score) then
            
            t_save  = t_fuel
            ! ASSIGN POWER
    
            !if(err < th_cvg_crit) exit TH
    
            allocate(ttemp(nfm(1),nfm(2),nfm(3)))
            do ii = 1, nfm(1)
            do jj = 1, nfm(2)
            do kk = 1, nfm(3)
                ttemp(ii,jj,kk) = AVG(p_fmfd(curr_bat,:,ii,jj,kk))
            end do
            end do
            end do
            if(fmfdon) call DET_POWER(ttemp(:,:,:))
            
            call POWER_NORM_START
            if(icore==score) then
                do iz = 1, nth(3)
                    print *, 'POWERZ', iz, th_iter
                    do iy = 1, nth(2)
                        write(*,'(<nth(1)>F10.3)') (pp(ix, iy, iz), ix = 1,nth(1))
                    enddo
                enddo
            endif
            call POWER_IMC_TO_START
    
            pp = 0d0 ; pp_thread = 0d0
        
    
        endif
            call INIT_CHILD
            call TEMP_CONVERGE(err)
            call MPI_BCAST(err, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr)
            !if(icore==score) print *, 'ERROR', th_iter, err
        
        call MPI_BCAST(t_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(t_comm_fuel, nth(1)*nth(2)*nth(3), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(rho_comm_cool, n_channels * (nth(3)+1), MPI_REAL8, score, MPI_COMM_WORLD, ierr)
            

        call TH_ASSIGN_GRID
        if(allocated(ttemp)) deallocate(ttemp)
    else
        exit
    end if



end do TH

    !> Check burnup loop exit condition
    if ( do_burn ) then
        !> initialization of CMFD parameters
        if ( istep_burnup .and. fmfdon == 0 ) then
            n_skip = 0
            acc_skip = 1
        end if

        !> Gather Burnup Tallies
        call MPI_reduce_burnup()

        !> Intra-pin flux distribution for iDTMC
        if ( fmfdon ) call INTRA_PIN_DTMC

        !> Make & Solve depletion matrix
        call DEPLETION
    
        if ( istep_burnup > nstep_burnup ) exit BURNUP
    else
        time4 = omp_get_wtime()
        call END_MSG(curr_bat,time3,time4)
        exit BURNUP
    end if

    call MPI_BARRIER(core,ierr)

    time4 = omp_get_wtime()
    call END_MSG(curr_bat,time3,time4)

end do BURNUP

    if ( curr_bat == n_batch ) then
        call BATCH_TALLY_MSG
        exit BATCH
    end if
    curr_bat = curr_bat + 1
    if ( n_batch > 1 .and. curr_bat == 1 ) then
        n_inact  = t_inact
        n_totcyc = t_totcyc
    end if

end do BATCH

if (allocated(source_bank)) deallocate(source_bank)
inquire(unit=prt_flux, opened=isopened)
if ( isopened ) close(prt_flux)
inquire(unit=prt_powr, opened=isopened)
if ( isopened ) close(prt_powr)
close(prt_keff)


if ( tally_switch > 0 .and. icore == score .and. .not. do_transient) then
    nsize = size(TallyFlux)
    allocate(tally_val(1:nsize, 1:n_act))
    
    open(prt_flux,file="flux.out",action="read",status="old")
    do i = 1, n_act 
        read(prt_flux, '(<nsize>ES15.7)') tally_val(:,i)
    enddo 
    close(prt_flux)
    open(prt_flux,file="flux.out",action="write",status="replace")
    do i = 1, nsize
        write(prt_flux, '(2ES15.7)') sum(tally_val(i,:))/n_act, STD_M(tally_val(i,:))
    enddo 
    close(prt_flux) 
    
    
    open(prt_powr,file="power.out",action="read",status="old")
    do i = 1, n_act 
        read(prt_powr, '(<nsize>ES15.7)') tally_val(:,i)
    enddo 
    close(prt_powr)
    open(prt_powr,file="power.out",action="write",status="replace")
    do i = 1, nsize
        write(prt_powr, '(2ES15.7)') sum(tally_val(i,:))/n_act, STD_M(tally_val(i,:))
    enddo 
    close(prt_powr) 
    deallocate(tally_val)
endif 

if ( tally_switch > 0 .and. icore == score .and. do_transient) then


    nsize = size(TallyPower)
    allocate(tally_val(1:nsize, 1:n_act))
    
    do curr_timestep = 1, n_timestep
        if (curr_timestep < 10) then 
            write (filename, "(A22,I1,A4)") "./DMC_data/power", curr_timestep, ".out"
        elseif (curr_timestep < 100) then                           
            write (filename, "(A22,I2,A4)") "./DMC_data/power", curr_timestep, ".out"
        elseif (curr_timestep < 1000) then                          
            write (filename, "(A22,I3,A4)") "./DMC_data/power", curr_timestep, ".out"
        endif 
        open(prt_powr,file=trim(filename),action="read",status="old")
        do i = 1, n_act 
            read(prt_powr, '(<nsize>ES15.7)') tally_val(:,i)
        enddo 
        close(prt_powr)
        
        open(prt_powr,file=trim(filename),action="write",status="replace")
        do i = 1, nsize
            write(prt_powr, '(2ES15.7)') sum(tally_val(i,:))/n_act, STD_M(tally_val(i,:))
        enddo 
        close(prt_powr) 
    enddo 
    deallocate(tally_val)
    
endif 

call MPI_BCAST(err, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr)

print *, 'WTF...', icore

call MPI_FINALIZE(ierr)

contains


subroutine INPUT_GEN
    implicit none
    integer:: cnum
    real(8):: vol, density

    density = -10.4
    vol = 66.7559060


    cnum = 100

!    open(124,file='FA.inp')
!    do jj = 1, 16
!    do ii = 1, 16
!        if ( ( ii == 1  .and. jj == 1  ) .or. &
!             ( ii == 2  .and. jj == 1  ) .or. &
!             ( ii == 15 .and. jj == 1  ) .or. &
!             ( ii == 16 .and. jj == 1  ) .or. &
!             ( ii == 1  .and. jj == 2  ) .or. &
!             ( ii == 16 .and. jj == 2  ) .or. &
!             ( ii == 4  .and. jj == 3  ) .or. &
!             ( ii == 5  .and. jj == 3  ) .or. &
!             ( ii == 12 .and. jj == 3  ) .or. &
!             ( ii == 13 .and. jj == 3  ) .or. &
!             ( ii == 3  .and. jj == 4  ) .or. &
!             ( ii == 6  .and. jj == 4  ) .or. &
!             ( ii == 11 .and. jj == 4  ) .or. &
!             ( ii == 14 .and. jj == 4  ) .or. &
!             ( ii == 3  .and. jj == 5  ) .or. &
!             ( ii == 6  .and. jj == 5  ) .or. &
!             ( ii == 11 .and. jj == 5  ) .or. &
!             ( ii == 14 .and. jj == 5  ) .or. &
!             ( ii == 4  .and. jj == 6  ) .or. &
!             ( ii == 5  .and. jj == 6  ) .or. &
!             ( ii == 12 .and. jj == 6  ) .or. &
!             ( ii == 13 .and. jj == 6  ) .or. &
!             ( ii == 8  .and. jj == 7  ) .or. &
!             ( ii == 9  .and. jj == 7  ) .or. &
!             ( ii == 7  .and. jj == 8  ) .or. &
!             ( ii == 10 .and. jj == 8  ) .or. &
!             ( ii == 7  .and. jj == 9  ) .or. &
!             ( ii == 10 .and. jj == 9  ) .or. &
!             ( ii == 8  .and. jj == 10 ) .or. &
!             ( ii == 9  .and. jj == 10 ) .or. &
!             ( ii == 4  .and. jj == 11 ) .or. &
!             ( ii == 5  .and. jj == 11 ) .or. &
!             ( ii == 12 .and. jj == 11 ) .or. &
!             ( ii == 13 .and. jj == 11 ) .or. &
!             ( ii == 3  .and. jj == 12 ) .or. &
!             ( ii == 6  .and. jj == 12 ) .or. &
!             ( ii == 11 .and. jj == 12 ) .or. &
!             ( ii == 14 .and. jj == 12 ) .or. &
!             ( ii == 3  .and. jj == 13 ) .or. &
!             ( ii == 6  .and. jj == 13 ) .or. &
!             ( ii == 11 .and. jj == 13 ) .or. &
!             ( ii == 14 .and. jj == 13 ) .or. &
!             ( ii == 4  .and. jj == 14 ) .or. &
!             ( ii == 5  .and. jj == 14 ) .or. &
!             ( ii == 12 .and. jj == 14 ) .or. &
!             ( ii == 13 .and. jj == 14 ) .or. &
!             ( ii == 1  .and. jj == 15 ) .or. &
!             ( ii == 16 .and. jj == 15 ) .or. &
!             ( ii == 1  .and. jj == 16 ) .or. &
!             ( ii == 2  .and. jj == 16 ) .or. &
!             ( ii == 15 .and. jj == 16 ) .or. &
!             ( ii == 16 .and. jj == 16 ) ) then
!        write(124,1), "% =========", "(",ii,",",jj,")"
!        cnum = cnum+1
!        write(124,2), "cell", cnum, 2,ii,jj, "Z", ii, jj, 1, "&", -1
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 2,ii,jj, "Z", ii, jj, 2, "&",  1, -2
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 2,ii,jj, "Z", ii, jj, 3, "&",  2, -3
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 2,ii,jj, "gap", "&",  3, -4
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 2,ii,jj, "can", "&",  4, -5
!        cnum = cnum+1
!        write(124,5), "cell", cnum, 2,ii,jj, "moderator", "&",  5
!        write(124,*)
!
!    else if ( ( ii == 4  .and. jj == 2  ) .or. &
!              ( ii == 13 .and. jj == 2  ) .or. &
!              ( ii == 2  .and. jj == 4  ) .or. &
!              ( ii == 15 .and. jj == 4  ) .or. &
!              ( ii == 6  .and. jj == 6  ) .or. &
!              ( ii == 11 .and. jj == 6  ) .or. &
!              ( ii == 6  .and. jj == 11 ) .or. &
!              ( ii == 11 .and. jj == 11 ) .or. &
!              ( ii == 2  .and. jj == 13 ) .or. &
!              ( ii == 15 .and. jj == 13 ) .or. &
!              ( ii == 4  .and. jj == 15 ) .or. &
!              ( ii == 13 .and. jj == 15 ) ) then
!        write(124,1), "% =========", "(",ii,",",jj,")"
!        cnum = cnum+1
!        write(124,2), "cell", cnum, 3,ii,jj, "A", ii, jj, 1, "&", -1
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 3,ii,jj, "A", ii, jj, 2, "&",  1, -2
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 3,ii,jj, "A", ii, jj, 3, "&",  2, -3
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 3,ii,jj, "gap", "&",  3, -4
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 3,ii,jj, "can", "&",  4, -5
!        cnum = cnum+1
!        write(124,5), "cell", cnum, 3,ii,jj, "moderator", "&",  5
!        write(124,*)
!
!    else if ( ( ii == 4  .and. jj == 4  ) .or. &
!              ( ii == 4  .and. jj == 5  ) .or. &
!              ( ii == 4  .and. jj == 12 ) .or. &
!              ( ii == 4  .and. jj == 13 ) .or. &
!              ( ii == 5  .and. jj == 4  ) .or. &
!              ( ii == 5  .and. jj == 5  ) .or. &
!              ( ii == 5  .and. jj == 12 ) .or. &
!              ( ii == 5  .and. jj == 13 ) .or. &
!              ( ii == 12 .and. jj == 4  ) .or. &
!              ( ii == 12 .and. jj == 5  ) .or. &
!              ( ii == 12 .and. jj == 12 ) .or. &
!              ( ii == 12 .and. jj == 13 ) .or. &
!              ( ii == 13 .and. jj == 4  ) .or. &
!              ( ii == 13 .and. jj == 5  ) .or. &
!              ( ii == 13 .and. jj == 12 ) .or. &
!              ( ii == 13 .and. jj == 13 ) ) then
!        write(124,1), "% =========", "(",ii,",",jj,")"
!        cnum = cnum+1
!        write(124,2), "cell", cnum, 4,ii,jj, "G", ii, jj, 1, "&", -1
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 4,ii,jj, "G", ii, jj, 2, "&",  1, -2
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 4,ii,jj, "G", ii, jj, 3, "&",  2, -3
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 4,ii,jj, "gap", "&",  3, -4
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 4,ii,jj, "can", "&",  4, -5
!        cnum = cnum+1
!        write(124,5), "cell", cnum, 4,ii,jj, "moderator", "&",  5
!        write(124,*)
!
!    else if ( ( ii == 8  .and. jj == 8  ) .or. &
!              ( ii == 9  .and. jj == 8  ) .or. &
!              ( ii == 8  .and. jj == 9  ) .or. &
!              ( ii == 9  .and. jj == 9  ) ) then
!        write(124,1), "% =========", "(",ii,",",jj,")"
!        cnum = cnum+1
!        write(124,2), "cell", cnum, 5,ii,jj, "I", ii, jj, 1, "&", -1
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 5,ii,jj, "I", ii, jj, 2, "&",  1, -2
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 5,ii,jj, "I", ii, jj, 3, "&",  2, -3
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 5,ii,jj, "gap", "&",  3, -4
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 5,ii,jj, "can", "&",  4, -5
!        cnum = cnum+1
!        write(124,5), "cell", cnum, 5,ii,jj, "moderator", "&",  5
!        write(124,*)
!
!    else
!        write(124,1), "% =========", "(",ii,",",jj,")"
!        cnum = cnum+1
!        write(124,2), "cell", cnum, 1,ii,jj, "F", ii, jj, 1, "&", -1
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 1,ii,jj, "F", ii, jj, 2, "&",  1, -2
!        cnum = cnum+1
!        write(124,3), "cell", cnum, 1,ii,jj, "F", ii, jj, 3, "&",  2, -3
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 1,ii,jj, "gap", "&",  3, -4
!        cnum = cnum+1
!        write(124,4), "cell", cnum, 1,ii,jj, "can", "&",  4, -5
!        cnum = cnum+1
!        write(124,5), "cell", cnum, 1,ii,jj, "moderator", "&",  5
!        write(124,*)
!    end if
!    end do
!    end do
!    close(124)
!
!    1 format(a,2x,a,i2.2,a,i2.2,a)
!    2 format(a,2x,i4,2x,i1,i2.2,i2.2,2x,a,i2.2,i2.2,i1,5x,a,2x,i2)
!    3 format(a,2x,i4,2x,i1,i2.2,i2.2,2x,a,i2.2,i2.2,i1,5x,a,2(2x,i2))
!    4 format(a,2x,i4,2x,i1,i2.2,i2.2,2x,a,8x,a,2(2x,i2))
!    5 format(a,2x,i4,2x,i1,i2.2,i2.2,2x,a,2x,a,2x,i2)


    open(124,file='FA.inp')
    do jj = 1, 16
    do ii = 1, 16
        if ( ( ii == 1  .and. jj == 1  ) .or. &
             ( ii == 2  .and. jj == 1  ) .or. &
             ( ii == 15 .and. jj == 1  ) .or. &
             ( ii == 16 .and. jj == 1  ) .or. &
             ( ii == 1  .and. jj == 2  ) .or. &
             ( ii == 16 .and. jj == 2  ) .or. &
             ( ii == 4  .and. jj == 3  ) .or. &
             ( ii == 5  .and. jj == 3  ) .or. &
             ( ii == 12 .and. jj == 3  ) .or. &
             ( ii == 13 .and. jj == 3  ) .or. &
             ( ii == 3  .and. jj == 4  ) .or. &
             ( ii == 6  .and. jj == 4  ) .or. &
             ( ii == 11 .and. jj == 4  ) .or. &
             ( ii == 14 .and. jj == 4  ) .or. &
             ( ii == 3  .and. jj == 5  ) .or. &
             ( ii == 6  .and. jj == 5  ) .or. &
             ( ii == 11 .and. jj == 5  ) .or. &
             ( ii == 14 .and. jj == 5  ) .or. &
             ( ii == 4  .and. jj == 6  ) .or. &
             ( ii == 5  .and. jj == 6  ) .or. &
             ( ii == 12 .and. jj == 6  ) .or. &
             ( ii == 13 .and. jj == 6  ) .or. &
             ( ii == 8  .and. jj == 7  ) .or. &
             ( ii == 9  .and. jj == 7  ) .or. &
             ( ii == 7  .and. jj == 8  ) .or. &
             ( ii == 10 .and. jj == 8  ) .or. &
             ( ii == 7  .and. jj == 9  ) .or. &
             ( ii == 10 .and. jj == 9  ) .or. &
             ( ii == 8  .and. jj == 10 ) .or. &
             ( ii == 9  .and. jj == 10 ) .or. &
             ( ii == 4  .and. jj == 11 ) .or. &
             ( ii == 5  .and. jj == 11 ) .or. &
             ( ii == 12 .and. jj == 11 ) .or. &
             ( ii == 13 .and. jj == 11 ) .or. &
             ( ii == 3  .and. jj == 12 ) .or. &
             ( ii == 6  .and. jj == 12 ) .or. &
             ( ii == 11 .and. jj == 12 ) .or. &
             ( ii == 14 .and. jj == 12 ) .or. &
             ( ii == 3  .and. jj == 13 ) .or. &
             ( ii == 6  .and. jj == 13 ) .or. &
             ( ii == 11 .and. jj == 13 ) .or. &
             ( ii == 14 .and. jj == 13 ) .or. &
             ( ii == 4  .and. jj == 14 ) .or. &
             ( ii == 5  .and. jj == 14 ) .or. &
             ( ii == 12 .and. jj == 14 ) .or. &
             ( ii == 13 .and. jj == 14 ) .or. &
             ( ii == 1  .and. jj == 15 ) .or. &
             ( ii == 16 .and. jj == 15 ) .or. &
             ( ii == 1  .and. jj == 16 ) .or. &
             ( ii == 2  .and. jj == 16 ) .or. &
             ( ii == 15 .and. jj == 16 ) .or. &
             ( ii == 16 .and. jj == 16 ) ) then
        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "Z", ii, jj, 1
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 1
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  6.16512E+20"
        write(124,8),                        "92238.711nc  2.24490E+22"
        write(124,8),                        " 8016.711nc  4.61311E+22"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "Z", ii, jj, 2
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 2
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  6.16512E+20"
        write(124,8),                        "92238.711nc  2.24490E+22"
        write(124,8),                        " 8016.711nc  4.61311E+22"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "Z", ii, jj, 3
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 3
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  6.16512E+20"
        write(124,8),                        "92238.711nc  2.24490E+22"
        write(124,8),                        " 8016.711nc  4.61311E+22"
        write(124,1), "END_MAT"
        write(124,*)

    else if ( ( ii == 4  .and. jj == 2  ) .or. &
              ( ii == 13 .and. jj == 2  ) .or. &
              ( ii == 2  .and. jj == 4  ) .or. &
              ( ii == 15 .and. jj == 4  ) .or. &
              ( ii == 6  .and. jj == 6  ) .or. &
              ( ii == 11 .and. jj == 6  ) .or. &
              ( ii == 6  .and. jj == 11 ) .or. &
              ( ii == 11 .and. jj == 11 ) .or. &
              ( ii == 2  .and. jj == 13 ) .or. &
              ( ii == 15 .and. jj == 13 ) .or. &
              ( ii == 4  .and. jj == 15 ) .or. &
              ( ii == 13 .and. jj == 15 ) ) then
        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "A", ii, jj, 1
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 1
        write(124,9), "n_iso",          "=", 10
        write(124,5), "isotopes",       "=", "92235.711nc  4.3869971E+20"
        write(124,8),                        "92238.711nc  2.1225324E+22"
        write(124,8),                        " 8016.711nc  4.7536482E+22"
        write(124,8),                        "64152.711nc  5.8050553E+18"
        write(124,8),                        "64154.711nc  6.2453349E+19"
        write(124,8),                        "64155.711nc  4.2125976E+20"
        write(124,8),                        "64156.711nc  5.7891287E+20"
        write(124,8),                        "64157.711nc  4.3977916E+20"
        write(124,8),                        "64158.711nc  6.9360859E+20"
        write(124,8),                        "64160.711nc  6.0276792E+20"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "A", ii, jj, 2
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 2
        write(124,9), "n_iso",          "=", 10
        write(124,5), "isotopes",       "=", "92235.711nc  4.3869971E+20"
        write(124,8),                        "92238.711nc  2.1225324E+22"
        write(124,8),                        " 8016.711nc  4.7536482E+22"
        write(124,8),                        "64152.711nc  5.8050553E+18"
        write(124,8),                        "64154.711nc  6.2453349E+19"
        write(124,8),                        "64155.711nc  4.2125976E+20"
        write(124,8),                        "64156.711nc  5.7891287E+20"
        write(124,8),                        "64157.711nc  4.3977916E+20"
        write(124,8),                        "64158.711nc  6.9360859E+20"
        write(124,8),                        "64160.711nc  6.0276792E+20"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "A", ii, jj, 3
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 3
        write(124,9), "n_iso",          "=", 10
        write(124,5), "isotopes",       "=", "92235.711nc  4.3869971E+20"
        write(124,8),                        "92238.711nc  2.1225324E+22"
        write(124,8),                        " 8016.711nc  4.7536482E+22"
        write(124,8),                        "64152.711nc  5.8050553E+18"
        write(124,8),                        "64154.711nc  6.2453349E+19"
        write(124,8),                        "64155.711nc  4.2125976E+20"
        write(124,8),                        "64156.711nc  5.7891287E+20"
        write(124,8),                        "64157.711nc  4.3977916E+20"
        write(124,8),                        "64158.711nc  6.9360859E+20"
        write(124,8),                        "64160.711nc  6.0276792E+20"
        write(124,1), "END_MAT"
        write(124,*)

    else if ( ( ii == 4  .and. jj == 4  ) .or. &
              ( ii == 4  .and. jj == 5  ) .or. &
              ( ii == 4  .and. jj == 12 ) .or. &
              ( ii == 4  .and. jj == 13 ) .or. &
              ( ii == 5  .and. jj == 4  ) .or. &
              ( ii == 5  .and. jj == 5  ) .or. &
              ( ii == 5  .and. jj == 12 ) .or. &
              ( ii == 5  .and. jj == 13 ) .or. &
              ( ii == 12 .and. jj == 4  ) .or. &
              ( ii == 12 .and. jj == 5  ) .or. &
              ( ii == 12 .and. jj == 12 ) .or. &
              ( ii == 12 .and. jj == 13 ) .or. &
              ( ii == 13 .and. jj == 4  ) .or. &
              ( ii == 13 .and. jj == 5  ) .or. &
              ( ii == 13 .and. jj == 12 ) .or. &
              ( ii == 13 .and. jj == 13 ) ) then

    else if ( ( ii == 8  .and. jj == 8  ) .or. &
              ( ii == 9  .and. jj == 8  ) .or. &
              ( ii == 8  .and. jj == 9  ) .or. &
              ( ii == 9  .and. jj == 9  ) ) then

    else
        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "F", ii, jj, 1
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 1
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  7.33270E+20"
        write(124,8),                        "92238.711nc  2.23336E+22"
        write(124,8),                        " 8016.711nc  4.61337E+22"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "F", ii, jj, 2
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 2
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  7.33270E+20"
        write(124,8),                        "92238.711nc  2.23336E+22"
        write(124,8),                        " 8016.711nc  4.61337E+22"
        write(124,1), "END_MAT"
        write(124,*)

        write(124,1), "MAT"
        write(124,2), "mat_name",       "=", "F", ii, jj, 3
        write(124,3), "density_gpcc",   "=", density
        write(124,4), "vol",            "=", vol
        write(124,5), "fissionable",    "=", "T"
        write(124,5), "depletable",     "=", "T"
        write(124,7), "N_DTMC",         '=', ii, jj, 1, 3
        write(124,6), "n_iso",          "=", 3
        write(124,5), "isotopes",       "=", "92235.711nc  7.33270E+20"
        write(124,8),                        "92238.711nc  2.23336E+22"
        write(124,8),                        " 8016.711nc  4.61337E+22"
        write(124,1), "END_MAT"
        write(124,*)
    end if
    end do
    end do
    close(124)

    1 format(a)
    2 format(t4,a,t20,a,1x,a,2i2.2,i1)
    3 format(t4,a,t20,a,1x,f5.1)
    4 format(t4,a,t20,a,1x,f11.7)
    5 format(t4,a,t20,a,1x,a)
    6 format(t4,a,t20,a,1x,i1)
    8 format(t20,2x,a)
    9 format(t4,a,t20,a,1x,i2)
    7 format(t4,a,t20,a,4i4)

    stop

end subroutine


end program 
