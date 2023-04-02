module PRINTER
    use VARIABLES, only: icore, score, iscore, ncore, n_batch, n_inact, n_act, &
                        n_totcyc, ngen, t_inact, t_totcyc, curr_cyc, keff
    use FMFD_HEADER, only: p_fmfd, k_fmfd, fmfdon, cmfdon, nfm, mvec, &
                        dual_fmfd, n_acc, fake_MC, DTMCBU, MCBU
    use DEPLETION_MODULE, only: preco, porc, istep_burnup, burn_step, &
                        do_burn
    use ENTROPY, only: ccrt, scrt, rampup, crt1c, crt2c, entrp0
    use OMP_LIB, only: OMP_GET_MAX_THREADS
    use TALLY, only: tallyon, k_eff
    use STATISTICS
    use SIMULATION_HEADER, only: t_tot, t_mc, t_det
    implicit none
    integer:: ii, jj, kk
    character(80):: dfile, dfile1


    contains

! =============================================================================
! BATCH_MSG
! =============================================================================
subroutine BATCH_MSG(curr_bat)
    integer, intent(in):: curr_bat

    if ( icore /= score ) return

    write(*,10), '   =========================================='
    write(*,11), '    Batch calculation #', curr_bat
    write(*,10), '   =========================================='
    write(*,*)

    10 format(A)
    11 format(A,i3)

end subroutine

! =============================================================================
! START_MSG
! =============================================================================
subroutine START_MSG
    use FMFD_HEADER, only: n_fake, pfmfdon, fmfd2mc, fm0, fm1, nfm, &
                            zigzagon, zzf0, n_zz, fcr, fcz, pcmfdon, quarter, &
                            inactive_cmfd, acc_skip
    use ENTROPY, only: mprupon
    implicit none

if ( icore==score ) then  

    write(*,*)
    write(*,*)
    write(*,10), '  > No. of Threads per Node   ', omp_get_max_threads()
    write(*,10), '  > No. of MPI Nodes          ', ncore
    write(*,10), '  > No. of Batches            ', n_batch
    write(*,10), '  > No. of Histories per Cycle', ngen
    if ( n_batch > 1 .or. mprupon ) then
    write(*,11), '  > Skip Cycles:', t_inact , &
                 '  /  Active Cycles:', t_totcyc-t_inact
    else
    write(*,11), '  > Skip Cycles:', n_inact , &
                 '  /  Active Cycles:', n_totcyc-n_inact
    end if
    write(*,*)
!    if (tally_switch > 0) then 
!        write(*,*), ' > Tally is On :: See tally.inp'
!    else 
!        write(*,*), ' > Tally is OFF' 
!    endif 

    if ( tallyon ) then
        write(*,*), ' > Tally is on'
    end if

    if ( mprupon ) then
    if ( ccrt == 1 .and. scrt == 1 ) then
        write(*,14), '  > m-PRUP is on', rampup, "  (1)", crt1c, "  (2)", crt2c
    else if ( ccrt == 1 .and. scrt == 2 ) then
        write(*,15), '  > m-PRUP is on', rampup, "  (1)", crt1c, "  (2)", crt2c
    else if ( ccrt == 2 .and. scrt == 1 ) then
        write(*,16), '  > m-PRUP is on', rampup, "  (1)", crt1c, "  (2)", crt2c
    else
        write(*,17), '  > m-PRUP is on', rampup, "  (1)", crt1c, "  (2)", crt2c
    end if
    else
        write(*,*), ' > m-PRUP is OFF'
    end if

    if ( fake_MC ) then
        write(*,13), '  > Fake MC is On ( No. of skips : ', n_fake, ' )'
    end if

    if ( fmfdon ) then
        if ( fmfd2mc ) then
        if ( pfmfdon ) then
        if ( cmfdon ) then
        if ( pcmfdon ) then
        write(*,*), ' > pFMFD with pCMFD is On'
        else
        write(*,*), ' > pFMFD with CMFD is On'
        end if
        else
        write(*,*), ' > pFMFD is On'
        end if
        else
        if ( cmfdon ) then
        if ( pcmfdon ) then
        write(*,*), ' > pFMFD with pCMFD is On'
        else
        write(*,*), ' > pFMFD with CMFD is On'
        end if
        else
        write(*,*), ' > FMFD is On'
        end if
        end if
        else
        if ( pfmfdon ) then
        if ( cmfdon ) then
        if ( pcmfdon ) then
        write(*,*), ' > pFMFD with pCMFD is On (w/o feedback)'
        else
        write(*,*), ' > pFMFD with CMFD is On (w/o feedback)'
        end if
        else
        write(*,*), ' > pFMFD is On (w/o feedback)'
        end if
        else
        if ( cmfdon ) then
        write(*,*), ' > FMFD with CMFD is On (w/o feedback)'
        else
        write(*,*), ' > FMFD is On (w/o feedback)'
        end if
        end if
        end if
        write(*,18), '   >> accumulation length = ', n_acc
        write(*,19), '   >> mesh grid = ', fm0(1:3), fm1(1:3), nfm(1:3)
        if ( zigzagon ) &
        write(*,20), '   >> zigzag = ', zzf0(2:n_zz+1)
        if ( cmfdon ) &
        write(*,20), '   >> 1-node CMFD = ', fcr, fcz
        if ( quarter ) &
        write(*,20), '   >> quarter core is only considered'
        if ( dual_fmfd ) &
        write(*,20), '   >> dual FMDF system is on'
        if ( inactive_CMFD ) &
        write(*,20), '   >> Inactive CMFD & active FMFD'
        if ( acc_skip > 0 ) &
        write(*,18), '   >> skip accumulation cycles = ', acc_skip
    else
        write(*,*), ' > FMFD is OFF'
    end if

    write(*,*)
    write(*,*), '   Transport Simulation Starts...' 

endif

10 format(A30,I9)
11 format(A16,I5,A19,I5)
12 format(A16)
13 format(A,I2,A)
14 format(A,I10,2(A,ES10.2))
15 format(A,I10,A,ES10.2,A,F10.1)
16 format(A,I10,A,F8.1,A,ES10.2)
17 format(A,I10,A,F8.1,A,F14.1)
18 format(A,i5)
19 format(A,6F7.1,3I4)
20 format(A,10I4)

end subroutine

! =============================================================================
! BURNUP_MSG
! =============================================================================
subroutine BURNUP_MSG
    if ( preco == 0 ) then
    write(*,10), '   =========================================='
    write(*,11), '      Burnup step', istep_burnup
    write(*,12), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
    write(*,10), '   =========================================='

    else
    if ( porc == 1 ) then
    write(*,10), '   =========================================='
    write(*,11), '      Burnup step', istep_burnup, '(predictor)'
    write(*,12), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
    write(*,10), '   =========================================='

    elseif ( porc == 2 ) then
    write(*,10), '   =========================================='
    write(*,11), '      Burnup step', istep_burnup, '(corrector)'
    write(*,12), burn_step(istep_burnup)/86400.d0, ' CUMULATIVE DAYS'
    write(*,10), '   =========================================='
    end if
    end if

    10 format(A45)
    11 format(A17,I4,2x,a)
    12 format(F14.2,A16)

end subroutine


! =============================================================================
! RUN_MSG
! =============================================================================
subroutine RUN_MSG(bat,time1,time2)
use ENTROPY, only: up_sign
implicit none
integer:: bat
real(8):: time1, time2
    
if ( icore == score ) then
    if ( curr_cyc <= n_inact .or. bat == 0 ) then
    if ( up_sign ) then
    write(*,10), curr_cyc, time2-time1, "sec", entrp0, " | ", "keff", keff, &
                 "//", ngen
    up_sign = .false.
    else
    write(*,11), curr_cyc, time2-time1, "sec", entrp0, " | ", "keff", keff
    end if
    else
        k_eff(bat,curr_cyc) = keff
        write(*,12), curr_cyc, time2-time1, "sec", entrp0, " | ", &
            "keff", keff, &
            "avg", AVG(k_eff(bat,n_inact+1:curr_cyc)), &
            "SD", PCM(STD_M(k_eff(bat,n_inact+1:curr_cyc)))
    end if
end if

10 format(i8,f9.2,1x,a,f10.5,1x,a,1x,a,f9.5,2x,a,i10)
11 format(i8,f9.2,1x,a,f10.5,1x,a,1x,a,f9.5)
12 format(i8,f9.2,1x,a,f10.5,1x,a,1x,2(a,f9.5,3x),a,f9.3)
13 format(I5,10ES15.7)

end subroutine


! =============================================================================
! CYCLE_TALLY_MSG
! =============================================================================
subroutine CYCLE_TALLY_MSG(bat)
    use tally, only: MC_tally
    use TH_HEADER, only: th_on, t_fuel, t_bulk
    use FMFD_HEADER, only: p_dep_mc, p_dep_dt, k_real, p_dep_dt_pert
    use MATERIAL_HEADER, only: materials, n_materials
    use PERTURBATION, only: perton
    use COSAMPLING, only: n_pert
    implicit none
    integer, intent(in):: bat
    integer:: cc, mm, nn, rr, ci
    real(8), allocatable:: zavgf(:,:,:,:), zavgp(:,:,:,:)  ! (c,i,j,r)
    integer:: nsum
    real(8):: vsum
    logical:: yes

    if ( icore /= score ) return

    ! multiplication factor
    if ( fmfdon .and. bat /= 0 ) then
    write(*,*)
    write(*,*), "   DTMC keff"
    do ii = 1, n_inact
        write(*,10), ii, k_fmfd(bat,ii)
    end do
    do ii = n_inact+1, n_totcyc
        if ( ii == n_totcyc ) then
        if ( preco == 1 ) then
            ! predictor
            if ( porc == 1 ) then
            write(*,11), ii, k_fmfd(bat,ii), AVG(k_fmfd(bat,n_inact+1:ii)), &
                    PCM(STD_M(k_fmfd(bat,n_inact+1:ii))), "iDTMC (predictor)"

            ! corrector
            elseif ( porc == 2 ) then
            write(*,11), ii, k_fmfd(bat,ii), AVG(k_fmfd(bat,n_inact+1:ii)), &
                    PCM(STD_M(k_fmfd(bat,n_inact+1:ii))), "iDTMC (corrector)"

            end if
        else
            write(*,11), ii, k_fmfd(bat,ii), AVG(k_fmfd(bat,n_inact+1:ii)), &
                    PCM(STD_M(k_fmfd(bat,n_inact+1:ii))), "iDTMC"
        end if
        else
        write(*,12), ii, k_fmfd(bat,ii), AVG(k_fmfd(bat,n_inact+1:ii)), &
                    PCM(STD_M(k_fmfd(bat,n_inact+1:ii)))
        end if
    end do
    if(perton) then
        do ii = n_inact+1,n_totcyc
            write(*,'(A,I5,I5,F10.6,F10.2)') 'REAL', istep_burnup, ii-n_inact, &
                AVG(k_real(bat,ii,1:n_pert)), &
                PCM(STD_S(k_real(bat,ii,1:n_pert)))
        enddo
        write(*,*)
    endif
    end if

    if ( bat == n_batch ) then
    ! power distribution normalization
    if ( tallyon ) &
    call NORM_DIST(MC_tally(1:n_batch,1:n_act,1,1,:,:,:))
    if ( fmfdon ) &
    call NORM_DIST(p_fmfd(1:n_batch,1:n_act,:,:,:))
    end if

    ! computing time
    if ( bat == 1 .and. .not. do_burn ) then
    t_MC = t_tot - t_det
    write(*,*)
    write(*,*), "   Computing time"
    do ii = 1, n_totcyc
        write(*,15), ii, AVG(t_MC(1:,ii)), AVG(t_det(1:,ii)), AVG(t_tot(1:,ii))
    end do
    write(*,*)
    end if

!    ! temperature distribution
!    if ( th_on ) then
!    write(*,*), "   Temperature distribution"
!    write(*,14), t_fuel/k_b
!    write(*,*)
!    write(*,14), t_bulk/k_b
!    write(*,*)
!    end if
!    end if


    ! bunrup dependent pin power distribution
    if ( DO_BURN ) then
    if ( DTMCBU .and. .not. MCBU ) then
    if ( istep_burnup == 0 ) then
    ! find if the file exists
    nsum = 0
    dfile = 'dep_dt0.out'
    do
    inquire(file=trim(dfile),exist=yes)
    if ( yes ) then
        nsum = nsum + 1
        if ( nsum < 10 ) then
            write(dfile,'(a,i1,a)'), 'dep_dt',nsum,'.out'
        else
            write(dfile,'(a,i2,a)'), 'dep_dt',nsum,'.out'
        end if

    else
        exit
    end if
    end do
    ! open a new file
    open(46,file=trim(dfile))
    close(46)
    end if
    ! parameter generation
    do ii = 1, n_act
    ! --- average
    do jj = 1, nfm(1)
    do kk = 1, nfm(2)
        p_dep_dt(ii,jj,kk,1) = AVG(p_dep_dt(ii,jj,kk,1:nfm(3)))
        if(.not. perton) cycle
        do mm = 1,n_pert
            p_dep_dt_pert(ii,mm,jj,kk,1) = &
                AVG(p_dep_dt_pert(ii,mm,jj,kk,1:nfm(3)))
        enddo
    end do
    end do
    ! --- summation
    vsum = 0; nsum = 0
    do jj = 1, nfm(1)
    do kk = 1, nfm(2)
        if ( isnan(p_dep_dt(ii,jj,kk,1)) ) cycle
        nsum = nsum + 1
        vsum = vsum + p_dep_dt(ii,jj,kk,1)
    end do
    end do
    p_dep_dt(ii,:,:,1) = p_dep_dt(ii,:,:,1)*nsum/dble(vsum)
    
    if(perton) then
    vsum = 0; nsum = 0
    do jj = 1, nfm(1)
    do kk = 1, nfm(2)
        if ( isnan(sum(p_dep_dt_pert(ii,:,jj,kk,1))) ) cycle
        nsum = nsum + n_act
        vsum = vsum + sum(p_dep_dt_pert(ii,:,jj,kk,1))
    end do
    end do
    endif
    p_dep_dt_pert(ii,:,:,:,1) = p_dep_dt_pert(ii,:,:,:,1)*nsum/dble(vsum)

    end do
    open(46,file=trim(dfile),access='append',status='old')
!    write(46,*), " HERE : pin power", " | step : ", istep_burnup
!    do jj = nfm(2), 1, -1
!    write(46,1), (p_dep_dt(1,ii,jj,1), ii = 1, nfm(1))
!    end do
!    write(46,*)
!    do jj = nfm(2), 1, -1
!    write(46,1), (AVG(p_dep_dt(1:3,ii,jj,1)), ii = 1, nfm(1))
!    end do
!    write(46,*)
!    do jj = nfm(2), 1, -1
!    write(46,1), (AVG(p_dep_dt(1:5,ii,jj,1)), ii = 1, nfm(1))
!    end do
!    write(46,*)
    do jj = nfm(2), 1, -1
    write(46,1), (AVG(p_dep_dt(1:n_act,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(46,*)
    
    if (perton) then
    write(*,*) 'PERTURBED AVG'
    do jj = nfm(2), 1, -1
    write(46,1), (AVG(p_dep_dt_pert(n_act,1:n_pert,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(46,*)
    endif

!    write(46,*), " HERE : SD of pin power"
    write(*,*) 'APPARENT SD'
    do jj = nfm(2), 1, -1
    write(46,1), (STD_M(p_dep_dt(1:n_act,ii,jj,1)) &
        /AVG(p_dep_dt(1:n_act,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(46,*)
    
    if (perton) then
    write(*,*) 'PERTURBED SD'
    do jj = nfm(2), 1, -1
    write(46,1), (STD_S(p_dep_dt_pert(n_act,1:n_pert,ii,jj,1)) &
        /AVG(p_dep_dt_pert(n_act,1:n_pert,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(46,*)
    endif
    
    close(46)

    write(*,*) 'WRITTEN DEP_DT FILE: ', trim(dfile)
    else
    if ( istep_burnup == 0 ) then
    ! find if the file exists
    nsum = 1
    dfile = 'dep_mc1.out'
    do
    inquire(file=trim(dfile),exist=yes)
    if ( yes ) then
        nsum = nsum + 1
        if ( nsum < 10 ) then
            write(dfile,'(a,i1,a)'), 'dep_mc',nsum,'.out'
        else
            write(dfile,'(a,i2,a)'), 'dep_mc',nsum,'.out'
        end if

    else
        exit
    end if
    end do
    ! open a new file
    open(45,file=trim(dfile))
    close(45)
    end if
    ! parameter generation
    do ii = 1, n_act
    ! --- average
    do jj = 1, nfm(1)
    do kk = 1, nfm(2)
        p_dep_mc(ii,jj,kk,1) = AVG(p_dep_mc(ii,jj,kk,1:nfm(3)))
    end do
    end do
    ! --- summation
    vsum = 0; nsum = 0
    do jj = 1, nfm(1)
    do kk = 1, nfm(2)
        if ( isnan(p_dep_mc(ii,jj,kk,1)) ) cycle
        nsum = nsum + 1
        vsum = vsum + p_dep_mc(ii,jj,kk,1)
    end do
    end do
        p_dep_mc(ii,:,:,1) = p_dep_mc(ii,:,:,1)*nsum/dble(vsum)
    end do
    open(45,file=trim(dfile),access='append',status='old')
    !write(45,*), " HERE : pin power", " | step : ", istep_burnup
    do jj = nfm(2), 1, -1
    write(45,1), (p_dep_mc(1,ii,jj,1), ii = 1, nfm(1))
    end do
    write(45,*)
    do jj = nfm(2), 1, -1
    write(45,1), (AVG(p_dep_mc(1:3,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(45,*)
    do jj = nfm(2), 1, -1
    write(45,1), (AVG(p_dep_mc(1:5,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(45,*)
    do jj = nfm(2), 1, -1
    write(45,1), (AVG(p_dep_mc(1:n_act,ii,jj,1)), ii = 1, nfm(1))
    end do
    write(45,*)

!    write(45,*), " HERE : SD of pin power"
!    do jj = nfm(2), 1, -1
!    write(45,1), (STD_M(p_dep_mc(1:n_act,ii,jj,1)) &
!        /AVG(p_dep_mc(1:n_act,ii,jj,1)), ii = 1, nfm(1))
!    end do
!    write(45,*)
    close(45)
    end if
    end if

    1 format(1000ES15.7)
    2 format(2I4,1000ES15.7)
    10 format(1X,I5,F10.6)
    11 format(1X,I5,2F10.6,F10.2,2X,A)
    12 format(1X,I5,2F10.6,F10.2,2X)
    16 format(4X,2F10.6,F10.2,2x,a)
    14 format(<nfm(1)>ES15.7)
    15 format(4X,I4,3F12.3)



!    ! intra pin power distribution
!    if ( DO_BURN ) then
!    if ( DTMCBU .and. .not. MCBU ) then
!    if ( istep_burnup == 0 ) then
!    ! find if the file exists
!    nsum = 1
!    dfile1 = 'intra_dt1.out'
!    do
!    inquire(file=trim(dfile1),exist=yes)
!    if ( yes ) then
!        nsum = nsum + 1
!        if ( nsum < 10 ) then
!            write(dfile1,'(a,i1,a)'), 'intra_dt',nsum,'.out'
!        else
!            write(dfile1,'(a,i2,a)'), 'intra_dt',nsum,'.out'
!        end if
!
!    else
!        exit
!    end if
!    end do
!    ! open a new file
!    open(48,file=trim(dfile1))
!    close(48)
!    end if
!    ! parameter print
!    open(48,file=trim(dfile1),access='append',status='old')
!    do ii = 1, n_materials
!        if ( .not. materials(ii)%depletable ) cycle
!        write(48,*), materials(ii)%flux
!    end do
!    write(48,*)
!    close(48)
!
!    else
!    if ( istep_burnup == 0 ) then
!    ! find if the file exists
!    nsum = 1
!    dfile1 = 'intra_mc1.out'
!    do
!    inquire(file=trim(dfile1),exist=yes)
!    if ( yes ) then
!        nsum = nsum + 1
!        if ( nsum < 10 ) then
!            write(dfile1,'(a,i1,a)'), 'intra_mc',nsum,'.out'
!        else
!            write(dfile1,'(a,i2,a)'), 'intra_mc',nsum,'.out'
!        end if
!
!    else
!        exit
!    end if
!    end do
!    ! open a new file
!    open(49,file=trim(dfile1))
!    close(49)
!    end if
!    ! parameter generation
!    open(49,file=trim(dfile1),access='append',status='old')
!    do ii = 1, n_materials
!        if ( .not. materials(ii)%depletable ) cycle
!        write(49,*), materials(ii)%flux
!    end do
!    write(49,*)
!    close(49)
!    end if
!    end if

end subroutine


! =============================================================================
! END_MSG
! =============================================================================
subroutine END_MSG(bat,time3,time4)
    implicit none
    integer:: bat
    real(8):: time3, time4

if ( icore == score ) then
    if ( n_batch > 1 .and. bat == 0 ) then
        write(*,*)
        return
    end if
    write(*,*)
    write(*,*), '   Simulation of Burnup Step Terminated...'
    write(*,10), "    - Elapsed time    : ", &
        time4 - time3, 'sec', (time4-time3)/60, 'min'
    if ( preco == 1 ) then
        ! predictor
    if ( porc == 1 ) then
    write(*,12), "    - Step Final keff : ", &
        AVG(k_eff(bat,n_inact+1:n_totcyc)), "+/-", &
        PCM(STD_M(k_eff(bat,n_inact+1:n_totcyc))), "(predcitor)"

        ! corrector
    elseif ( porc == 2 ) then
    write(*,12), "    - Step Final keff : ", &
        AVG(k_eff(bat,n_inact+1:n_totcyc)), "+/-", &
        PCM(STD_M(k_eff(bat,n_inact+1:n_totcyc))), "(corrector)"

    end if

    else
    write(*,11), "    - Step Final keff : ", &
        AVG(k_eff(bat,n_inact+1:n_totcyc)), "+/-", &
        PCM(STD_M(k_eff(bat,n_inact+1:n_totcyc)))
    end if
    write(*,*)

    10 format(A,F10.3,A4,F8.2,A4)
    11 format(A,F10.6,A4,F8.3)
    12 format(A,F10.6,A4,F8.3,1X,A)

end if

end subroutine


! =============================================================================
! BATCH_TALLY_MSG
! =============================================================================
subroutine BATCH_TALLY_MSG
    use TALLY, only: n_type, ttally, MC_tally, ttally, MC_tally, n_tcycle, &
                    tcycle
    use FMFD_HEADER, only: wholecore, k_fmfd2, p_fmfd2, inactive_CMFD
    implicit none
    real(8), allocatable:: k_avg(:,:)
    real(8), allocatable:: t_avg(:,:,:,:,:)
    integer:: xx, yy, zz

    if ( icore /= score ) return

    allocate(k_avg(n_batch,n_act))
    ! apparent standard deviation of the multiplication factor
    write(*,11), '   =========================================='
    write(*,*), "   Apparent standard deviation of the multiplication factor"
    if ( .not. fmfdon ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = PCM(STD_M(k_eff(ii,n_inact+1:n_inact+jj)))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, AVG(k_avg(1:n_batch,ii))
    end do
    write(*,*)
    else
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = PCM(STD_M(k_fmfd(ii,n_inact+1:n_inact+jj)))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, AVG(k_avg(1:n_batch,ii))
    end do
    write(*,*)

    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = PCM(STD_M(k_fmfd2(ii,n_inact+1:n_inact+jj)))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, AVG(k_avg(1:n_batch,ii))
    end do
    write(*,*)
    end if
    end if

    ! real standard deviation of the multiplication factor
    if ( n_batch > 1 ) then
    write(*,11), '   =========================================='
    write(*,*), "   Real standard deviation of the multiplication factor"
    if ( .not. fmfdon ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = AVG(k_eff(ii,n_inact+1:n_inact+jj))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, PCM(STD_S(k_avg(1:n_batch,ii)))
    end do
    write(*,*)
    else
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = AVG(k_fmfd(ii,n_inact+1:n_inact+jj))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, PCM(STD_S(k_avg(1:n_batch,ii)))
    end do
    write(*,*)

    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        k_avg(ii,jj) = AVG(k_fmfd2(ii,n_inact+1:n_inact+jj))
    end do
    end do
    do ii = 1, n_act
        write(*,10), ii, PCM(STD_S(k_avg(1:n_batch,ii)))
    end do
    write(*,*)
    end if
    end if
    end if

    if ( wholecore ) then
        !call BATCH_TALLY_MSG_ACTIVE
        return
    end if

    allocate(t_avg(n_batch,n_act,nfm(1),nfm(2),nfm(3)))
    ! apparent standard deviation of the pin-wise information
    if ( tallyon .or. fmfdon ) then
    write(*,11), '   =========================================='
    write(*,*), "   Apparent standard deviation of the pin-wise information"
    if ( .not. fmfdon .or. inactive_cmfd ) then
    if ( .not. do_burn ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
        t_avg(ii,jj,1,1,1) = STD_P(MC_tally(ii,1:jj,1,1,:,:,:))
    end do
    end do
    do ii = 1, n_act
        write(*,12), ii, AVG(t_avg(1:n_batch,ii,1,1,1))
    end do
    write(*,*)
    end if
    end if
    if ( fmfdon ) then
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        t_avg(ii,jj,1,1,1) = STD_P(p_fmfd(ii,1:jj,:,:,:))
    end do
    end do
    do ii = 1, n_act
        write(*,12), ii, AVG(t_avg(1:n_batch,ii,1,1,1))
    end do
    write(*,*)

    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
        t_avg(ii,jj,1,1,1) = STD_P(p_fmfd2(ii,1:jj,:,:,:))
    end do
    end do
    do ii = 1, n_act
        write(*,12), ii, AVG(t_avg(1:n_batch,ii,1,1,1))
    end do
    write(*,*)
    end if
    end if
    end if

    if ( n_batch > 1 ) then
    ! real standard deviation of the pin-wise information
    if ( tallyon .or. fmfdon ) then
    write(*,11), '   =========================================='
    write(*,*), "   Real standard deviation of the pin-wise information"
    if ( .not. fmfdon .or. inactive_CMFD ) then
    if ( .not. do_burn ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(MC_tally(ii,1:jj,1,1,xx,yy,zz))
    end do
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_act
        write(*,12), ii, STD_PS(t_avg(1:n_batch,ii,:,:,:))
    end do
    write(*,*)
    end if
    end if
    if ( fmfdon ) then
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd(ii,1:jj,xx,yy,zz))
    end do
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_act
        write(*,12), ii, STD_PS(t_avg(1:n_batch,ii,:,:,:))
    end do
    write(*,*)
    
    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd2(ii,1:jj,xx,yy,zz))
    end do
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_act
        write(*,12), ii, STD_PS(t_avg(1:n_batch,ii,:,:,:))
    end do
    write(*,*)
    end if
    end if
    end if
    end if

!    ! 2D pin power distribution
!    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   2D pin power distribution"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    do jj = 1, n_act
!    do xx = 1, nfm(1)
!    do yy = 1, nfm(2)
!    do zz = 1, nfm(3)
!        t_avg(1,jj,xx,yy,zz) = AVG(MC_tally(1,1:jj,1,1,xx,yy,zz))
!    end do
!        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nfm(3))
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nfm(2), 1, -1
!    write(*,16), (t_avg(1,tcycle(ii),xx,yy,1), xx = 1, nfm(1))
!    end do
!    write(*,*)
!    end do
!    else
!    write(*,*), "   FMFD"
!    do jj = 1, n_act
!    do xx = 1, nfm(1)
!    do yy = 1, nfm(2)
!    do zz = 1, nfm(3)
!        t_avg(1,jj,xx,yy,zz) = AVG(p_fmfd(1,1:jj,xx,yy,zz))
!    end do
!        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nfm(3))
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nfm(2), 1, -1
!    write(*,16), (t_avg(1,tcycle(ii),xx,yy,1), xx = 1, nfm(1))
!    end do
!    write(*,*)
!    end do
!    end if
!    end if

    ! batch-averaged 2D pin power distribution
    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
    write(*,11), '   =========================================='
    write(*,*), "   2D pin power distribution"
    if ( .not. fmfdon .or. inactive_CMFD ) then
    if ( .not. do_burn ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(MC_tally(ii,1:jj,1,1,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if
    end if
    if ( fmfdon .and. .not. inactive_CMFD ) then
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd(ii,1:jj,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do

    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd2(ii,1:jj,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if
    end if

    if ( fmfdon .and. inactive_CMFD ) then
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = p_fmfd(ii,jj,xx,yy,zz)
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if
    end if

    ! 2D pin apparent standard deviation distribution
    if ( tallyon .or. fmfdon ) then
    write(*,11), '   =========================================='
    write(*,*), "   2D apparent standard deviation distribution"
    if ( .not. fmfdon .or. inactive_CMFD ) then
    if ( .not. do_burn ) then
    write(*,*), "   MC"
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    call NORM_DIST(t_avg(:,:,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (STD_M(t_avg(1,1:tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if 
    end if
    if ( fmfdon .and. .not. inactive_CMFD ) then
    write(*,*), "   FMFD"
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    call NORM_DIST(t_avg(:,:,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (STD_M(t_avg(1,1:tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if
    end if

    ! 2D pin real standard deviation distribution
    if ( n_batch > 1 ) then
    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
    write(*,11), '   =========================================='
    write(*,*), "   2D real standard deviation distribution"
    if ( .not. fmfdon .and. inactive_CMFD ) then
    if ( do_burn ) then
    write(*,*), "   MC"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(MC_tally(ii,1:jj,1,1,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (STD_S(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    end do
    end if
    end if
    if ( fmfdon ) then
    write(*,*), "   FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd(ii,1:jj,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (STD_S(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do

    if ( dual_fmfd ) then
    write(*,*), "   DUAL FMFD"
    do ii = 1, n_batch
    do jj = 1, n_act
    do xx = 1, nfm(1)
    do yy = 1, nfm(2)
    do zz = 1, nfm(3)
        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd2(ii,1:jj,xx,yy,zz))
    end do
        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nfm(3))
    end do
    end do
    end do
    end do
    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
    do ii = 1, n_tcycle
    write(*,*), " - cycle :", tcycle(ii)
    do yy = nfm(2), 1, -1
    write(*,16), (STD_S(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nfm(1))
    end do
    write(*,*)
    end do
    end if
    end if
    end if
    end if


!    ! batch-averaged 3D pin power distribution
!    if ( tallyon .or. fmfdon ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   3D pin power distribution"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nfm(1)
!    do yy = 1, nfm(2)
!    do zz = 1, nfm(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(MC_tally(ii,1:jj,1,1,xx,yy,zz))
!    end do
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    ii = n_tcycle
!    do zz = nfm(3), 1, -1
!    do yy = nfm(2), 1, -1
!    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,zz)), xx = 1, nfm(1))
!    end do
!    end do
!    write(*,*)
!    else
!    write(*,*), "   FMFD"
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nfm(1)
!    do yy = 1, nfm(2)
!    do zz = 1, nfm(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd(ii,1:jj,xx,yy,zz))
!    end do
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    ii = n_tcycle
!    do zz = nfm(3), 1, -1
!    do yy = nfm(2), 1, -1
!    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,zz)), xx = 1, nfm(1))
!    end do
!    end do
!    write(*,*)
!
!    if ( dual_fmfd ) then
!    write(*,*), "   DUAL FMFD"
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nfm(1)
!    do yy = 1, nfm(2)
!    do zz = 1, nfm(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(p_fmfd2(ii,1:jj,xx,yy,zz))
!    end do
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    ii = n_tcycle
!    do zz = nfm(3), 1, -1
!    do yy = nfm(2), 1, -1
!    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,zz)), xx = 1, nfm(1))
!    end do
!    end do
!    write(*,*)
!    end if
!    end if
!    end if

    10 format(4X,I4,F10.2)
    11 format(A)
    12 format(4X,I4,ES15.7)
    15 format(4X,3ES15.7)
    16 format(<nfm(1)>ES15.7)

end subroutine

!! =============================================================================
!! BATCH_TALLY_MSG
!! =============================================================================
!subroutine BATCH_TALLY_MSG_ACTIVE
!    use TALLY, only: n_type, ttally, MC_tally, ttally, MC_tally, n_tcycle, &
!                    tcycle
!    use FMFD_HEADER, only: nam, fcr, fcz
!    implicit none
!    real(8), allocatable:: t_avg(:,:,:,:,:)
!    integer:: xx, yy, zz
!
!    ! active core region
!    nam(1) = nfm(1)-2*fcr
!    nam(2) = nfm(2)-2*fcr
!    nam(3) = nfm(3)-2*fcz
!
!    allocate(t_avg(n_batch,n_act,nam(1),nam(2),nam(3)))
!    ! apparent standard deviation of the pin-wise information
!    if ( tallyon .or. fmfdon ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   Apparent standard deviation of the pin-wise information"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        MC_tally(:,:,1,1,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!        t_avg(ii,jj,1,1,1) = STD_P(t_avg(ii,1:jj,:,:,:))
!    end do
!    end do
!    do ii = 1, n_act
!        write(*,12), ii, AVG(t_avg(1:n_batch,ii,1,1,1))
!    end do
!    write(*,*)
!    else
!    write(*,*), "   FMFD"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        p_fmfd(:,:,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!        t_avg(ii,jj,1,1,1) = STD_P(t_avg(ii,1:jj,:,:,:))
!    end do
!    end do
!    do ii = 1, n_act
!        write(*,12), ii, AVG(t_avg(1:n_batch,ii,1,1,1))
!    end do
!    write(*,*)
!    end if
!    end if
!
!    ! real standard deviation of the pin-wise information
!    if ( n_batch > 1 ) then
!    if ( tallyon .or. fmfdon ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   Real standard deviation of the pin-wise information"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        MC_tally(:,:,1,1,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!    end do
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_act
!        write(*,12), ii, STD_PS(t_avg(1:n_batch,ii,:,:,:))
!    end do
!    write(*,*)
!    else
!    write(*,*), "   FMFD"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        p_fmfd(:,:,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!    end do
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_act
!        write(*,12), ii, STD_PS(t_avg(1:n_batch,ii,:,:,:))
!    end do
!    write(*,*)
!    end if
!    end if
!    end if
!
!    ! 2D pin power distribution
!    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   2D pin power distribution"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        MC_tally(:,:,1,1,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(1,jj,xx,yy,zz) = AVG(t_avg(1,1:jj,xx,yy,zz))
!    end do
!        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nam(3))
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nam(2), 1, -1
!    write(*,16), (t_avg(1,tcycle(ii),xx,yy,1), xx = 1, nam(1))
!    end do
!    write(*,*)
!    end do
!    else
!    write(*,*), "   FMFD"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        p_fmfd(:,:,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(1,jj,xx,yy,zz) = AVG(t_avg(1,1:jj,xx,yy,zz))
!    end do
!        t_avg(1,jj,xx,yy,:) = sum(t_avg(1,jj,xx,yy,:))/dble(nam(3))
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nam(2), 1, -1
!    write(*,16), (t_avg(1,tcycle(ii),xx,yy,1), xx = 1, nam(1))
!    end do
!    write(*,*)
!    end do
!    end if
!    end if
!
!!    ! batch-averaged 2D pin power distribution
!!    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
!!    write(*,11), '   =========================================='
!!    write(*,*), "   2D pin power distribution"
!!    if ( .not. fmfdon ) then
!!    write(*,*), "   MC"
!!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!!        MC_tally(:,:,1,1,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!!    do ii = 1, n_batch
!!    do jj = 1, n_act
!!    do xx = 1, nam(1)
!!    do yy = 1, nam(2)
!!    do zz = 1, nam(3)
!!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!!    end do
!!        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nam(3))
!!    end do
!!    end do
!!    end do
!!    end do
!!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!!    do ii = 1, n_tcycle
!!    write(*,*), " - cycle :", tcycle(ii)
!!    do yy = nam(2), 1, -1
!!    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nam(1))
!!    end do
!!    write(*,*)
!!    end do
!!    else
!!    write(*,*), "   FMFD"
!!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!!        p_fmfd(:,:,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!!    do ii = 1, n_batch
!!    do jj = 1, n_act
!!    do xx = 1, nam(1)
!!    do yy = 1, nam(2)
!!    do zz = 1, nam(3)
!!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!!    end do
!!        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nam(3))
!!    end do
!!    end do
!!    end do
!!    end do
!!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!!    do ii = 1, n_tcycle
!!    write(*,*), " - cycle :", tcycle(ii)
!!    do yy = nam(2), 1, -1
!!    write(*,16), (AVG(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nam(1))
!!    end do
!!    write(*,*)
!!    end do
!!    end if
!!    end if
!
!    ! 2D pin real standard deviation distribution
!    if ( n_batch > 1 ) then
!    if ( tallyon .or. fmfdon .and. n_tcycle >= 1 ) then
!    write(*,11), '   =========================================='
!    write(*,*), "   2D real standard deviation distribution"
!    if ( .not. fmfdon ) then
!    write(*,*), "   MC"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        MC_tally(:,:,1,1,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!    end do
!        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nam(3))
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nam(2), 1, -1
!    write(*,16), (STD_S(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nam(1))
!    end do
!    end do
!    else
!    write(*,*), "   FMFD"
!    t_avg(:,:,1:nam(1),1:nam(2),1:nam(3)) = &
!        p_fmfd(:,:,fcr+1:nfm(1)-fcr,fcr+1:nfm(2)-fcr,fcz+1:nfm(3)-fcz)
!    do ii = 1, n_batch
!    do jj = 1, n_act
!    do xx = 1, nam(1)
!    do yy = 1, nam(2)
!    do zz = 1, nam(3)
!        t_avg(ii,jj,xx,yy,zz) = AVG(t_avg(ii,1:jj,xx,yy,zz))
!    end do
!        t_avg(ii,jj,xx,yy,:) = sum(t_avg(ii,jj,xx,yy,:))/dble(nam(3))
!    end do
!    end do
!    end do
!    end do
!    call NORM_DIST(t_avg(1:n_batch,1:n_act,:,:,:))
!    do ii = 1, n_tcycle
!    write(*,*), " - cycle :", tcycle(ii)
!    do yy = nam(2), 1, -1
!    write(*,16), (STD_S(t_avg(1:n_batch,tcycle(ii),xx,yy,1)), xx = 1, nam(1))
!    end do
!    end do
!    end if
!    end if
!    end if
!
!    10 format(4X,I4,F10.2)
!    11 format(A)
!    12 format(4X,I4,ES15.7)
!    15 format(4X,3ES15.7)
!    16 format(<nam(1)>ES15.7)
!
!end subroutine
!
!
!subroutine HY
!    use TALLY, only: MC_tally, MC_stally, MC_scat, n_tgroup, n_type
!    use FMFD_HEADER, only: nfm, dfm
!    use INPUT_READER, only: directory, filename
!    implicit none
!    integer:: gg
!    integer:: ii, jj
!    real(8):: ptemp(nfm(1),nfm(2),nfm(3))
!
!    if ( icore /= score ) return
!    filename = trim(directory)//'group_constants.out'
!
!    open(1,file=trim(filename))
!    write(1,*), "keff : ", AVG(k_eff(1,n_inact+1:n_totcyc)), "+/-", &
!        PCM(STD_M(k_eff(1,n_inact+1:n_totcyc)))
!
!!    do ii = 1, nfm(1)
!!    do jj = 1, nfm(2)
!!    do kk = 1, nfm(3)
!!        ptemp(ii,jj,kk) = AVG(MC_tally(1,:,1,1,ii,jj,kk))
!!    end do
!!    end do
!!    end do
!!    write(1,3), ptemp
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!    do gg = 1, n_tgroup
!        write(1,1), ii, jj, dfm(1), dfm(2), &
!                    AVG(MC_tally(1,:,1,gg,ii,jj,1)), &
!                    AVG(MC_tally(1,:,2,gg,ii,jj,1)), &
!                    AVG(MC_tally(1,:,4,gg,ii,jj,1)), mod(dble(gg),2D0), &
!                    AVG(MC_tally(1,:,3,gg,ii,jj,1)), &
!                    AVG(MC_tally(1,:,5,gg,ii,jj,1)), &
!                    AVG(MC_scat(1,:,gg,1,ii,jj,1)), &
!                    AVG(MC_scat(1,:,gg,2,ii,jj,1))
!    end do
!        write(1,2), AVG(MC_tally(1,:,n_type,gg,ii,jj,1))
!    do gg = 1, n_tgroup
!        write(1,2), AVG(MC_stally(1,:,1,gg,ii,jj,1,2)), &
!                    AVG(MC_stally(1,:,1,gg,ii,jj,1,1)), &
!                    AVG(MC_stally(1,:,1,gg,ii,jj,1,4)), &
!                    AVG(MC_stally(1,:,1,gg,ii,jj,1,3))
!        write(1,2), AVG(MC_stally(1,:,2,gg,ii,jj,1,2)), &
!                    AVG(MC_stally(1,:,2,gg,ii,jj,1,1)), &
!                    AVG(MC_stally(1,:,2,gg,ii,jj,1,4)), &
!                    AVG(MC_stally(1,:,2,gg,ii,jj,1,3))
!    end do
!    end do
!    end do
!    1 format(2i4,2f7.3,100es15.7)
!    2 format(100es15.7)
!    3 format(<nfm(1)>ES15.7)
!    close(1)
!
!end subroutine

end module
