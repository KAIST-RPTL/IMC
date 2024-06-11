module PERTURBATION
    use FMFD_HEADER
    !###use FMFD, only: D_TILDA_CALCULATION, D_PHAT_CALCULATION, PFMFD_MATRIX, SET_ZERO_M, POWER
    use CMFD, only: L_DTILDA, ONE_NODE_CMFD, FM_TRANSPOSE
    use PCMFD, only: L_PDHAT, L_PMATRIX
    use STATISTICS, only: AVG, STD_S, STD_M, PCM, COV_M
    use RANDOMS, only: rang
    use COSAMPLING
    use omp_lib
    implicit none
    logical:: perton = .false.
    ! Excluded, as they are only used for unused subroutines
!    real(8), allocatable, dimension(:,:,:):: avg_sigt, avg_siga, avg_nufi
!    real(8), allocatable, dimension(:,:,:):: std_sigt, std_siga, std_nufi
    real(8), allocatable:: Mfm0(:,:,:,:), fm_s0(:,:,:)
    real(8), allocatable:: del_f(:,:,:), del_m(:,:,:,:)
    real(8), allocatable:: adjphi(:,:,:)

    contains

! =============================================================================
! 
! =============================================================================
subroutine PERTURBATION_KEFF(bat, k_eff,cyc)
    use VARIABLES, only: n_act, n_inact, n_totcyc, icore, score, iscore
    use CONSTANTS, only: pi
    use STATISTICS, only: AVG_P
    use MPI
    use FMFD_HEADER
    implicit none
    real(8), intent(in):: k_eff
    integer, intent(in):: cyc, bat
    real(8):: unphi(nfm(1),nfm(2),nfm(3))
    real(8):: inv_k
    integer:: xx, yy, ierr, it

    real(8) :: cos1, cos2, ktmp, ktmp_std

    real(8) :: t1, t2

    if (.not. allocated(fphi2) ) then
        sqrt2 = sqrt(2D0)
        sqrt2pi = sqrt2*pi
        allocate(fphi2(nfm(1),nfm(2),nfm(3)))
        allocate(del_m(nfm(1),nfm(2),nfm(3),7))
        allocate(del_f(nfm(1),nfm(2),nfm(3)))
        allocate(Mfm0(nfm(1),nfm(2),nfm(3),7))
        allocate(fm_s0(nfm(1),nfm(2),nfm(3)))
        allocate(acc_sigt(n_totcyc,nfm(1),nfm(2),nfm(3)))
        allocate(acc_siga(n_totcyc,nfm(1),nfm(2),nfm(3)))
        allocate(acc_nufi(n_totcyc,nfm(1),nfm(2),nfm(3)))
        allocate(acc_phi(n_totcyc,nfm(1),nfm(2),nfm(3)))
        allocate(acc_Jn(n_totcyc,nfm(1),nfm(2),nfm(3),6))
        allocate(acc_J0(n_totcyc,nfm(1),nfm(2),nfm(3),6))
        allocate(acc_J1(n_totcyc,nfm(1),nfm(2),nfm(3),6))
        allocate(adjphi(nfm(1),nfm(2),nfm(3)))
        allocate(acc_kap(n_totcyc,nfm(1),nfm(2),nfm(3)))
        k_pert  = 1
        k_pert2 = 1
    end if


    ! initialization of parameters
    if ( cyc <= acc_skip ) return
    if ( iscore ) then
    where ( fphi1 /= 0 ) 
    acc_phi(cyc,:,:,:)  = fm_avg(:,:,:)%phi
    acc_sigt(cyc,:,:,:) = fm_avg(:,:,:)%sig_t
    acc_siga(cyc,:,:,:) = fm_avg(:,:,:)%sig_a
    acc_nufi(cyc,:,:,:) = fm_avg(:,:,:)%nusig_f
    acc_kap(cyc,:,:,:)  = fm_avg(:,:,:)%kappa
    end where
    end if

    !call MPI_BCAST(acc_phi ,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(acc_sigt,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(acc_siga,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(acc_nufi,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    !call MPI_BCAST(acc_kap ,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,0,MPI_COMM_WORLD,ierr)

    if ( cyc <= n_inact ) return

    !if ( icore == score ) call RECALL2

    ! static initialization
    if ( iscore ) then
    inv_k = 1D0/k_eff
    where ( fphi1 /= 0 )
    unphi(:,:,:)   = fphi1(:,:,:)
    fm_t(:,:,:)    = fm_avg(:,:,:)%sig_t
    fm_a(:,:,:)    = fm_avg(:,:,:)%sig_a
    fm_nf(:,:,:)   = fm_avg(:,:,:)%nusig_f 
    fphi1(:,:,:)   = fm_avg(:,:,:)%phi
    elsewhere
    unphi = 0.d0;
    fm_t  = 0.d0;
    fm_a  = 0.d0;
    fm_nf = 0.d0;
    fphi1 = 0.d0;
    end where
    do ii = 1, 6
    where ( fphi1 /= 0 )
    fmJ0(:,:,:,ii) = fm_avg(:,:,:)%J0(ii)
    fmJ1(:,:,:,ii) = fm_avg(:,:,:)%J1(ii)
    elsewhere
    fmJ0(:,:,:,ii) = 0.d0
    fmJ1(:,:,:,ii) = 0.d0;
    end where
    end do
    fmJn = fmJ1-fmJ0
    fmD  = 1D0 / (3D0 * fm_t)
    where( fphi1 == 0 ) fmD = 0
    where( fphi1 == 0 .or. fphi1 < 1E-13 ) fm_nf = 0
    !k_pert2 = k_eff
    !k_pert3 = k_eff
    where ( fphi1 /= 0 )
    fphi2 = unphi*unphi
    fm_s0 = fm_nf*unphi !*unphi
    end where
    ! ---
    call D_TILDA_CALCULATION  ! &&&
    call D_PHAT_PERTURB(0)
    call PFMFD_MATRIX
    !write(*,*) 'CPR', 0, sum(fmD), sum(fmJ0), sum(fmJ1), sum(fmJn)
    !write(*,*) 'HATTIL', 0, sum(fmDt), sum(fmDh), sum(fm_a), sum(fphi1)
    where(isnan(fm_s0)) fm_s0 = 0.d0;
    !write(*,*) 'FMS0', sum(fm_s0), sum(fm_nf), sum(unphi)
    if ( zigzagon ) call SET_ZERO_M
    Mfm0 = Mfm
    ! --- adjoint
    if ( cyc == n_inact + 1 ) then
        t1 = omp_get_wtime()
        if(.not. allocated(ptphi)) allocate(ptphi  (n_pert,nfm(1),nfm(2),nfm(3)))
        allocate(tmpfphi1(nfm(1),nfm(2),nfm(3)))
        tmpfphi1 = fphi1
        call FM_TRANSPOSE   ! &&&
        t2 = omp_get_wtime()
        if(icore==score) print *, 'TRANSPOSE', t2 - t1
        t1 = omp_get_wtime()
        call POWER_TMP(k_pert2(1))
        t2 = omp_get_wtime()
        if(icore==score) print *, 'POWER_TMP', t2 - t1
        where ( fphi1 /= 0 )
            adjphi(:,:,:) = fphi1(:,:,:)
        elsewhere
            adjphi(:,:,:) = 0D0
        endwhere
        fphi1 = tmpfphi1
        deallocate(tmpfphi1)
    end if
    end if
    
    cos1 = omp_get_wtime()
    ! correlation
    call COSAMPLING_MAIN(cyc)
    
    cos2 = omp_get_wtime()
    if(iscore) write(*,*) 'TIME', cos2-cos1

    ! multiplication factor calculation
    do mm = 1, n_pert
    if ( iscore ) then
        ! parameters
        where ( fphi1 /= 0 )
        fm_t(:,:,:)  = coxs(:,:,:,mm,1,it)
        fm_a(:,:,:)  = coxs(:,:,:,mm,2,it)
        fm_nf(:,:,:) = coxs(:,:,:,mm,3,it)
        fphi1(:,:,:) = ptphi(mm,:,:,:)        ! flux sampling
        del_f(:,:,:) = (fm_nf(:,:,:)-fm_avg(:,:,:)%nusig_f) * inv_k
        elsewhere
        fm_t = 0.d0
        fm_a = 0.d0
        fm_nf= 0.d0
        del_f= 0.d0
        fphi1= 0.d0
        end where
        do ii = 1, 6
            where ( fphi1 /= 0 )
                fmJ0(:,:,:,ii) = fm_avg(:,:,:)%J0(ii)
                fmJ1(:,:,:,ii) = fm_avg(:,:,:)%J1(ii)
            elsewhere
                fmJ0(:,:,:,ii) = 0.d0
                fmJ1(:,:,:,ii) = 0.d0
            end where
        end do

        fmJn = fmJ1-fmJ0
        ! perturbed migration matrix
        !write(*,*) 'CPR', mm, sum(fmD), sum(fmJ0), sum(fmJ1), sum(fmJn)
        fmD = 1D0/(3D0*fm_t)
        where ( fm_t == 0 ) fmD = 0
        !where ( fphi1 == 0 .or. fphi1 < 1E-13 ) fm_nf = 0
        !fphi1 = unphi
        call D_TILDA_CALCULATION   ! &&&
        call D_PHAT_PERTURB(mm)
        call PFMFD_MATRIX
        if ( zigzagon ) call SET_ZERO_M
        del_m = Mfm-Mfm0
        k_pert2(mm) = sum(MAT_MUL(del_m,del_f,unphi)*adjphi)/sum(fm_s0*adjphi)
        k_pert2(mm) = 1D0/(inv_k-k_pert2(mm))
        write(*,'(I3,F12.5,F12.5,F8.3)') mm, k_pert2(mm), AVG(k_pert2(1:mm)), PCM(STD_S(k_pert2(1:mm)))
        !writE(*,*) 'TEST', sum(fm_t), sum(fm_a), sum(fm_nf)
        if(it == maxiter) then
        where(fm_avg%nusig_f /= 0)
            p_dep_dt_pert(cyc-n_inact,mm,:,:,:) = &
                fm_avg(:,:,:)%kappa / fm_avg(:,:,:)%nusig_f &
                * fm_nf(:,:,:) * fphi1
        elsewhere
            p_dep_dt_pert(cyc-n_inact,mm,:,:,:) = 0D0
        endwhere
        endif
    end if

    ! for power distribution
!    if ( cyc == n_totcyc ) then
!    if ( mm <= 60 ) then
!    !if ( cyc == n_totcyc ) then
!        call ONE_NODE_CMFD(k_pert3(mm),1D-6)
!       ! if ( icore == score ) p_pert(mm,:,:,:) = fm_nf(:,:,:)*fphi1(:,:,:)
!    end if
!    end if

    !write(*,*) AVG(k_pert2(1:n_pert)), PCM(STD_S(k_pert2(1:n_pert)))
    ! standard deviation
    if(iscore) then
        write(*,*) 'ACT', it, cyc-n_inact+1, AVG(k_pert2(:)), PCM(STD_S(k_pert2(:)))
        write(*,*) 'PERTK', AVG(k_pert2(:))
        k_real(bat,cyc,1:n_pert) = k_pert2(1:n_pert)
    endif
    enddo
end subroutine

! =============================================================================
subroutine RECALL2

    ! writing
!    open(6,file='recall_parameter')
!    write(6,10), acc_phi
!    write(6,10), acc_sigt
!    write(6,10), acc_siga
!    write(6,10), acc_nufi
!    write(6,10), acc_J0 
!    write(6,10), acc_J1 
!    write(6,10), fm_avg(:,:,:)%phi
!    write(6,10), fm_avg(:,:,:)%sig_t
!    write(6,10), fm_avg(:,:,:)%sig_a
!    write(6,10), fm_avg(:,:,:)%nusig_f
!    do ii = 1, 6
!    write(6,10), fm_avg(:,:,:)%J0(ii)
!    write(6,10), fm_avg(:,:,:)%J1(ii)
!    end do
!    write(6,10), fphi1
!    close(6)
!    stop


    ! reading
    open(6,file='recall_parameter')
    read(6,10), acc_phi
    read(6,10), acc_sigt
    read(6,10), acc_siga
    read(6,10), acc_nufi
    read(6,10), acc_J0
    read(6,10), acc_J1
    acc_Jn = acc_J1-acc_J0
    read(6,10), fm_avg(:,:,:)%phi
    read(6,10), fm_avg(:,:,:)%sig_t
    read(6,10), fm_avg(:,:,:)%sig_a
    read(6,10), fm_avg(:,:,:)%nusig_f
    do ii = 1, 6
    read(6,10), fm_avg(:,:,:)%J0(ii)
    read(6,10), fm_avg(:,:,:)%J1(ii)
    end do
    read(6,10), fphi1
    close(6)

    10 format(<nfm(1)>ES16.8)

end subroutine


! =============================================================================
! EXCLUDED, SINCE THEY ARE NOT IN USE
! =============================================================================
!subroutine STD_FMFD_PARAMETERS(cyc)
!    implicit none
!    integer, intent(in):: cyc
!    real(8):: avg1, avg0, std1, std0, cov
!
!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!        if ( fphi1(ii,jj,kk) == 0 ) cycle
!        ! flux
!        avg0 = AVG(acc_phi(acc_skip+1:cyc,ii,jj,kk))
!        std0 = STD_M(acc_phi(acc_skip+1:cyc,ii,jj,kk))/avg0
!
!        ! total cross section
!        avg1 = AVG(acc_sigt(acc_skip+1:cyc,ii,jj,kk))
!        std1 = STD_M(acc_sigt(acc_skip+1:cyc,ii,jj,kk))/avg1
!        cov  = 2D0*COV_M(acc_sigt(acc_skip+1:cyc,ii,jj,kk), &
!               acc_phi(acc_skip+1:cyc,ii,jj,kk))/(avg0*avg1)
!        avg_sigt(ii,jj,kk) = avg1/avg0
!        std_sigt(ii,jj,kk) = avg_sigt(ii,jj,kk)*sqrt(std1*std1+std0*std0-cov)
!
!        ! absorption cross section
!        avg1 = AVG(acc_siga(acc_skip+1:cyc,ii,jj,kk))
!        std1 = STD_M(acc_siga(acc_skip+1:cyc,ii,jj,kk))/avg1
!        cov  = 2D0*COV_M(acc_siga(acc_skip+1:cyc,ii,jj,kk), &
!               acc_phi(acc_skip+1:cyc,ii,jj,kk))/(avg0*avg1)
!        avg_siga(ii,jj,kk) = avg1/avg0
!        std_siga(ii,jj,kk) = avg_siga(ii,jj,kk)*sqrt(std1*std1+std0*std0-cov)
!
!        ! fission cross section
!        avg1 = AVG(acc_nufi(acc_skip+1:cyc,ii,jj,kk))
!        std1 = STD_M(acc_nufi(acc_skip+1:cyc,ii,jj,kk))/avg1
!        cov  = 2D0*COV_M(acc_nufi(acc_skip+1:cyc,ii,jj,kk), &
!               acc_phi(acc_skip+1:cyc,ii,jj,kk))/(avg0*avg1)
!        avg_nufi(ii,jj,kk) = avg1/avg0
!        std_nufi(ii,jj,kk) = avg_nufi(ii,jj,kk)*sqrt(std1*std1+std0*std0-cov)
!    end do
!    end do
!    end do
!
!end subroutine



! =============================================================================
! =============================================================================
!subroutine PERTURB_M
!    implicit none
!    include 'mkl_vml.f90'
!    real(8):: err(1), inverr(1)
!
!    ! perturbed cross section
!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!        if ( fphi1(ii,jj,kk) == 0 ) cycle
!        err(1) = 2d0*rang()-1d0
!        call vderfinv(1,err,inverr)
!        fm_t(ii,jj,kk) = std_sigt(ii,jj,kk)*sqrt2*inverr(1)+avg_sigt(ii,jj,kk)
!        err(1) = 2d0*rang()-1d0
!        call vderfinv(1,err,inverr)
!        fm_a(ii,jj,kk) = std_siga(ii,jj,kk)*sqrt2*inverr(1)+avg_siga(ii,jj,kk)
!    end do
!    end do
!    end do
!    fmD = 1D0/(3D0*fm_t)
!    where( fphi1 == 0 ) fmD = 0
!
!    ! perturbed migration matrix
!!    call D_TILDA_CALCULATION   ! &&&
!!    call D_PHAT_CALCULATION
!!    call PFMFD_MATRIX
!!    if ( zigzagon ) call SET_ZERO_M
!!    del_m = Mfm-Mfm0
!
!end subroutine
!
!
!! =============================================================================
!! 
!! =============================================================================
!subroutine PERTURB_F(keff)
!    implicit none
!    include 'mkl_vml.f90'
!    real(8), intent(in):: keff
!    real(8):: err(1), inverr(1)
!
!    ! perturbed cross section
!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!        if ( fphi1(ii,jj,kk) == 0 .or. avg_nufi(ii,jj,kk) == 0 ) cycle
!        err(1) = 2d0*rang()-1d0
!        call vderfinv(1,err,inverr)
!        fm_nf(ii,jj,kk) = std_nufi(ii,jj,kk)*sqrt2*inverr(1)+avg_nufi(ii,jj,kk)
!        del_f(ii,jj,kk) = fm_nf(ii,jj,kk) - fm_avg(ii,jj,kk)%nusig_f
!    end do
!    end do
!    end do
!    where( fphi1 == 0 .or. fphi1 < 1E-13 ) fm_nf = 0
!
!end subroutine


! =============================================================================
! 
! =============================================================================
function MAT_MUL(aa,bb,cc)
    real(8):: MAT_MUL(nfm(1),nfm(2),nfm(3))
    real(8), intent(in):: aa(:,:,:,:), bb(:,:,:), cc(:,:,:)

    MAT_MUL = 0
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    if ( cc(ii,jj,kk) == 0 ) cycle
    MAT_MUL(ii,jj,kk) = (bb(ii,jj,kk)-aa(ii,jj,kk,4))*cc(ii,jj,kk)
    if ( ii /= 1 )      MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,3)*cc(ii-1,jj,kk)
    if ( ii /= nfm(1) ) MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,5)*cc(ii+1,jj,kk)
    if ( jj /= 1 )      MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,2)*cc(ii,jj-1,kk)
    if ( jj /= nfm(2) ) MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,6)*cc(ii,jj+1,kk)
    if ( kk /= 1 )      MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,1)*cc(ii,jj,kk-1)
    if ( kk /= nfm(3) ) MAT_MUL(ii,jj,kk) = MAT_MUL(ii,jj,kk) &
                                          - aa(ii,jj,kk,7)*cc(ii,jj,kk+1)
    !write(*,*) 'MARK', ii, jj, kk, MAT_MUL(ii,jj,kk)
    !write(*,*) cc(ii,jj,kk), bb(ii,jj,kk)
    !write(*,*) aa(ii,jj,kk,1:7)
    if(isnan(MAT_MUL(ii,jj,kk))) MAT_MUL = 0.d0 ! TEMPORARY
    end do
    end do
    end do

end function


! =============================================================================
! 
! =============================================================================
function MAT_MUL2(aa,bb,cc)
    real(8):: MAT_MUL2(nfm(1),nfm(2),nfm(3))
    real(8), intent(in):: aa(:,:,:,:), bb(:,:,:), cc(:,:,:)

    MAT_MUL2 = 0
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    if ( cc(ii,jj,kk) == 0 ) cycle
    MAT_MUL2(ii,jj,kk) = (bb(ii,jj,kk)-aa(ii,jj,kk,4))*cc(ii,jj,kk)*cc(ii,jj,kk)
    if ( ii /= 1 )      MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,3)*cc(ii-1,jj,kk)*cc(ii-1,jj,kk)
    if ( ii /= nfm(1) ) MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,5)*cc(ii+1,jj,kk)*cc(ii+1,jj,kk)
    if ( jj /= 1 )      MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,2)*cc(ii,jj-1,kk)*cc(ii,jj-1,kk)
    if ( jj /= nfm(2) ) MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,6)*cc(ii,jj+1,kk)*cc(ii,jj+1,kk)
    if ( kk /= 1 )      MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,1)*cc(ii,jj,kk-1)*cc(ii,jj,kk-1)
    if ( kk /= nfm(3) ) MAT_MUL2(ii,jj,kk) = MAT_MUL2(ii,jj,kk) &
                                          - aa(ii,jj,kk,7)*cc(ii,jj,kk+1)*cc(ii,jj,kk+1)
    end do
    end do
    end do

end function

! =============================================================================
! 
! =============================================================================
function PT(aa,bb,cc)
    real(8):: PT
    real(8), intent(in):: aa(:,:,:,:), bb(:,:,:), cc(:,:,:)

    PT = 0
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    if ( cc(ii,jj,kk) == 0 ) cycle
    PT = PT + (bb(ii,jj,kk)-aa(ii,jj,kk,4))*cc(ii,jj,kk)
    if ( ii /= 1 )      PT = PT - aa(ii,jj,kk,3)*cc(ii-1,jj,kk)
    if ( ii /= nfm(1) ) PT = PT - aa(ii,jj,kk,5)*cc(ii+1,jj,kk)
    if ( jj /= 1 )      PT = PT - aa(ii,jj,kk,2)*cc(ii,jj-1,kk)
    if ( jj /= nfm(2) ) PT = PT - aa(ii,jj,kk,6)*cc(ii,jj+1,kk)
    if ( kk /= 1 )      PT = PT - aa(ii,jj,kk,1)*cc(ii,jj,kk-1)
    if ( kk /= nfm(3) ) PT = PT - aa(ii,jj,kk,7)*cc(ii,jj,kk+1)
    end do
    end do
    end do

end function


! =============================================================================
! 
! =============================================================================
subroutine SAMPLING_FROM_PDF(cyc)
    use RANDOMS, only: RANG
    implicit none
    integer, intent(in):: cyc
    real(8):: v_min, v_max, v_diff
    real(8):: val(0:10)
    real(8):: cdf(0:10)
    real(8):: urn
    integer:: steps = 1D1
    

    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    if ( fphi1(ii,jj,kk) == 0 ) cycle

    ! ===== total cross section =====
    ! minimum & maximum
    v_min = minval(acc_sigt(acc_skip+1:cyc,ii,jj,kk))
    v_max = maxval(acc_sigt(acc_skip+1:cyc,ii,jj,kk))
    v_diff = (v_max-v_min)/real(steps-1)
    val(0) = v_min-v_diff*5D-1
    do mm = 1, steps
        val(mm) = val(mm-1)+v_diff
    end do

    ! cumulative density function (CDF)
    cdf = 0
    do mm = acc_skip+1, cyc
    do nn = 1, steps
    if ( acc_sigt(mm,ii,jj,kk) < val(nn) ) then
        cdf(nn:) = cdf(nn:) + 1
        exit
    end if
    end do
    end do
    cdf = cdf / cdf(steps)

    ! sampling
    urn = rang()
    do mm = 1, steps
    if ( urn < cdf(mm) ) then
        fm_t(ii,jj,kk) = val(mm-1)+v_diff*(urn-cdf(mm-1))/(cdf(mm)-cdf(mm-1))
        exit
    end if
    end do


    ! ===== absorption cross section =====
    ! minimum & maximum
    v_min = minval(acc_siga(acc_skip+1:cyc,ii,jj,kk))
    v_max = maxval(acc_siga(acc_skip+1:cyc,ii,jj,kk))
    v_diff = (v_max-v_min)/real(steps-1)
    val(0) = v_min-v_diff*5D-1
    do mm = 1, steps
        val(mm) = val(mm-1)+v_diff
    end do

    ! cumulative density function (CDF)
    cdf = 0
    do mm = acc_skip+1, cyc
    do nn = 1, steps
    if ( acc_siga(mm,ii,jj,kk) < val(nn) ) then
        cdf(nn:) = cdf(nn:) + 1
        exit
    end if
    end do
    end do
    cdf = cdf / cdf(steps)

    ! sampling
    urn = rang()
    do mm = 1, steps
    if ( urn < cdf(mm) ) then
        fm_a(ii,jj,kk) = val(mm-1)+v_diff*(urn-cdf(mm-1))/(cdf(mm)-cdf(mm-1))
        exit
    end if
    end do


    ! ===== nu fission cross section =====
    ! minimum & maximum
    v_min = minval(acc_nufi(acc_skip+1:cyc,ii,jj,kk))
    v_max = maxval(acc_nufi(acc_skip+1:cyc,ii,jj,kk))
    v_diff = (v_max-v_min)/real(steps-1)
    val(0) = v_min-v_diff*5D-1
    do mm = 1, steps
        val(mm) = val(mm-1)+v_diff
    end do

    ! cumulative density function (CDF)
    cdf = 0
    do mm = acc_skip+1, cyc
    do nn = 1, steps
    if ( acc_nufi(mm,ii,jj,kk) < val(nn) ) then
        cdf(nn:) = cdf(nn:) + 1
        exit
    end if
    end do
    end do
    cdf = cdf / cdf(steps)

    ! sampling
    urn = rang()
    do mm = 1, steps
    if ( urn < cdf(mm) ) then
        fm_nf(ii,jj,kk) = val(mm-1)+v_diff*(urn-cdf(mm-1))/(cdf(mm)-cdf(mm-1))
        del_f(ii,jj,kk) = fm_nf(ii,jj,kk)-fm_avg(ii,jj,kk)%nusig_f
        exit
    end if
    end do


!    write(*,1), "====="
!    write(*,1), "vvv", v_min, v_max, v_diff
!    write(*,1), "raw", acc_nufi(acc_skip:cyc,ii,jj,kk)
!    write(*,1), "val", val
!    write(*,1), "cdf", cdf
!    write(*,1), "urn", urn
!    write(*,2), fm_avg(ii,jj,kk)%nusig_f, fm_nf(ii,jj,kk)
!    1 format(A,20es15.7)
!    2 format(20ES15.7)
!
!    if ( jj == 3 ) pause


!    print*, ii, jj, kk, "==============="
!    print*, acc_sigt(acc_skip+1:cyc,ii,jj,kk)
!    print*, acc_siga(acc_skip+1:cyc,ii,jj,kk)
!    print*, acc_nufi(acc_skip+1:cyc,ii,jj,kk)
!    print*, fm_t(ii,jj,kk), fm_a(ii,jj,kk), fm_nf(ii,jj,kk)
!    print*, "===================="
!    print*
!    if ( ii == 3 ) pause

    end do
    end do
    end do
    fmD = 1D0/(3D0*fm_t)
    where( fphi1 == 0 ) fmD = 0

    ! perturbed migration matrix
!    call D_TILDA_CALCULATION   ! &&&
!    call D_PHAT_CALCULATION
!    call PFMFD_MATRIX
!    if ( zigzagon ) call SET_ZERO_M
!    del_m = Mfm-Mfm0

end subroutine



! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_PHAT_PERTURB(pt)
    integer, intent(in):: pt
    integer :: i, j, k
    fmDh = 0.d0
    ! inner region (outgoing direction)
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        fmDh(i,j,k,1) = (fmJ0(i,j,k,1)-5D-1*fmDt(i,j,k,1) &
                    *(fphi1(i,j,k)-fphi1(i-1,j,k)))/fphi1(i,j,k)
        if ( i /= nfm(1) ) & ! x1
        fmDh(i,j,k,2) = (fmJ1(i,j,k,2)+5D-1*fmDt(i,j,k,2) &
                    *(fphi1(i+1,j,k)-fphi1(i,j,k)))/fphi1(i,j,k)
        if ( j /= 1 ) &      ! y0
        fmDh(i,j,k,3) = (fmJ0(i,j,k,3)-5D-1*fmDt(i,j,k,3) &
                    *(fphi1(i,j,k)-fphi1(i,j-1,k)))/fphi1(i,j,k)
        if ( j /= nfm(2) ) & ! y1
        fmDh(i,j,k,4) = (fmJ1(i,j,k,4)+5D-1*fmDt(i,j,k,4) &
                    *(fphi1(i,j+1,k)-fphi1(i,j,k)))/fphi1(i,j,k)
        if ( k /= 1 ) &      ! y0
        fmDh(i,j,k,5) = (fmJ0(i,j,k,5)-5D-1*fmDt(i,j,k,5) &
                    *(fphi1(i,j,k)-fphi1(i,j,k-1)))/fphi1(i,j,k)
        if ( k /= nfm(3) ) & ! y1
        fmDh(i,j,k,6) = (fmJ1(i,j,k,6)+5D-1*fmDt(i,j,k,6) &
                    *(fphi1(i,j,k+1)-fphi1(i,j,k)))/fphi1(i,j,k)
        !write(*,*) 'FFF', i,j,k, fphi1(i,j,k), fmDh(i,j,k,1:6)
    end do
    end do
    end do

    !write(*,*) 'HAT1', sum(fmDh), sum(fphi1), sum(fmDt)

    ! boundary
    if ( .not. zigzagon ) then
    i = 1;      fmDh(i,:,:,1) = -fmJn(i,:,:,1)/fphi1(i,:,:)
    i = nfm(1); fmDh(i,:,:,2) = +fmJn(i,:,:,2)/fphi1(i,:,:)
    j = 1;      fmDh(:,j,:,3) = -fmJn(:,j,:,3)/fphi1(:,j,:)
    j = nfm(2); fmDh(:,j,:,4) = +fmJn(:,j,:,4)/fphi1(:,j,:)
    else
    do i = 1, zz_div
        !write(*,*) 'ZZ', zzf0(i), zzf1(i), zzf2(i)
        fmDh(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) = &
            -fmJn(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) &
            /fphi1(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:)
        fmDh(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) = &
            +fmJn(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) &
            /fphi1(zzf2(i),zzf0(i)+1:zzf0(i+1),:)
        fmDh(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) = &
            -fmJn(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) &
            /fphi1(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:)
        fmDh(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) = &
            +fmJn(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) &
            /fphi1(zzf0(i)+1:zzf0(i+1),zzf2(i),:)
    end do
    end if
    k = 1;      fmDh(:,:,k,5) = -fmJn(:,:,k,5)/fphi1(:,:,k)
    k = nfm(3); fmDh(:,:,k,6) = +fmJn(:,:,k,6)/fphi1(:,:,k)

    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if(fphi1(i,j,k)==0) fmDh(i,j,k,1:6) = 0.d0
    enddo
    enddo
    enddo

end subroutine

subroutine POWER_TMP (k_eff,phi2)
    use SOLVERS, only: BiCGSTAB_PRE, BiCGSTAB_ILU
    use PRECONDITIONER, only: FINE_ILU_INITIAL, FINE_ILU
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in), optional:: phi2(:,:,:)
    real(8), parameter:: ONE = 1D0
    real(8):: phi0(nfm(1),nfm(2),nfm(3))
    real(8):: sorphi(nfm(1),nfm(2),nfm(3))
    integer:: iter, iter_max = 1D5
    real(8):: err
    real(8):: tt0, tt1, idx
    real(8):: kpre
    real(8):: ww = 1.7

    if ( .not. allocated(fn) ) call FINE_ILU_INITIAL
    call FINE_ILU

    err = ONE
    iter = 1
    tt0 = omp_get_wtime()
!    idx = 1
    if ( dual_fmfd .and. present(phi2) ) fphi1 = phi2
    do while ( ( err > 1D-9 ) .and. ( iter < iter_max ) )
        !call CPU_TIME(tt0)
        iter = iter + 1
        !sorphi = ww*fphi1+(1D0-ww)*phi0
        phi0 = fphi1
        kpre = k_eff
        fm_s = fm_nf(:,:,:)*fphi1(:,:,:)/k_eff
        !fm_s = fm_nf(:,:,:)*sorphi(:,:,:)/k_eff
        !fphi1 = BiCGSTAB_PRE(Mfm(:,:,:,:),fm_s(:,:,:))
        fphi1 = BiCGSTAB_ILU(mvec3(:,:),fm_s(:,:,:))
        k_eff = k_eff*sum(fm_nf*fphi1*fm_nf*fphi1) &
               / sum(fm_nf*fphi1*fm_nf*phi0)
        !err = abs(k_eff-kpre)/k_eff
        if ( err < 1D-3 ) then
            tt1 = omp_get_wtime()
            print*, 'CURR', tt1-tt0
            exit
        end if
        err = sum(abs(fphi1-phi0))/sum(fphi1)
        !print*, iter, k_eff, err
        !call CPU_TIME(tt1)
    enddo
    tt1 = omp_get_wtime()
    print*, "time", tt1-tt0

    
end subroutine 
end module
