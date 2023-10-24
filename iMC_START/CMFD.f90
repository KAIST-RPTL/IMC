module CMFD
    use VARIABLES   ! ***
    use FMFD_HEADER
    !use VARIABLES, only: icore, score, iscore
    use PCMFD, only: OUT_OF_ZZ, OUT_OF_ZZ1, IN_ZZ
    use MPI, only: MPI_WTIME
    implicit none
    integer:: ix0, ix1, iy0, iy1, iz0, iz1

    contains


! =============================================================================
! ONE_NODE_CMFD
! =============================================================================
subroutine ONE_NODE_CMFD(k_eff,ecvg,fphi2)
    use SOLVERS ! ***
    !use SOLVERS, only: BICG_G, BiCG_L, BICG_G_ILU, BICG_L_ILU!, SORL, SORG
    use PCMFD, only: L_PDHAT, L_PBC, L_PMATRIX, L_POUTJ, L_PSOURCE, L_PREFJ
    use MPI, only: MPI_COMM_WORLD, MPI_REAL8
    use PRECONDITIONER
    use VARIABLES, only: curr_cyc
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in):: ecvg  ! eigenvalue convergence
    real(8), intent(in), optional:: fphi2(:,:,:)
    real(8) :: error, k_eff0, k_eff1, mpie
    integer :: global, local
    integer :: iter, iter_max = 3D2 ! 2D2
    logical :: cvg  ! converged

    real(8):: jsweep, eigen, gpara, fixed, source, &
              module1, module2, homo, lpara
    ! current sweeping
    ! global calculation
    ! global parameters + ILU
    ! local calculation
    ! local source
    ! modulation 1
    ! modulation 2
    ! Reference J + homogenization
    ! local parameters + ILU
!    jsweep = 0
!    eigen = 0
!    gpara = 0
!    fixed = 0
!    source = 0
!    module1 = 0
!    module2 = 0
!    homo = 0
!    lpara = 0


    if ( pcmfdon ) then
    if ( icore == score ) then
    call L2G
    call L_DTILDA
    call L_PDHAT
    call D_BC
    call L_PBC
    call L_PMATRIX
    where( fphi1 == 0 ) Mfm(:,:,:,4) = 0
    cvg = .false.
    if ( isnan(sum(Mfm)) ) then
        k_eff = 3D0
        cvg = .true.
    end if
    if ( dual_fmfd .and. present(fphi2) ) fphi1 = fphi2
    end if
    call MPI_BCAST(cvg,1,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    if ( cvg ) return
    call MPI_BCAST(Mfm,n_nodes*7,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    call ILU_DECOMPOSE(Mfm)

    iter  = 1; error = 1; cvg = .false.
    do
    if ( icore == score ) then
    ! ------------------------------- GLOBAL
    call G_PDHAT
    call G_PMATRIX
    call GLOBAL_ILU
    k_eff1 = k_eff


    !do 
    do global = 1, 5
    k_eff0 = k_eff
    cphi0 = cphi1
    cm_s = cm_nf*cphi0/k_eff
    cphi1 = BiCG_G_ILU(mvec2,cm_s)
    k_eff = k_eff*sum(cm_nf*cphi1*cm_nf*cphi1) &
           / sum(cm_nf*cphi0*cm_nf*cphi1)
    !if ( abs(k_eff0-k_eff)/k_eff < 1D-9 ) exit
    if ( abs(k_eff0-k_eff)/k_eff < 1D-12 ) exit
    end do
    error = abs(k_eff-k_eff1)/k_eff
    !error = sum(abs(fphi1-fphi0))/sum(fphi1)
    !print*, iter, k_eff, error
    end if

    ! convergence test
    !if ( icore == score .and. ( error < 1D-9 .or. isnan(k_eff)  &
    if ( icore == score .and. ( error < ecvg .or. isnan(k_eff)  &
        .or. iter > iter_max ) ) cvg = .true.
    call MPI_BCAST(cvg,1,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    if ( cvg ) exit

    ! ------------------------------- LOCAL
    if ( icore == score ) call G_PINJ
    
    do local = 1, 2
    if ( icore == score ) then
    call G2L
    fphi0 = fphi1
    call L_PSOURCE(k_eff)
    end if

    call MPI_BCAST(fm_s,n_nodes,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    call LINEATION2(fm_s(:,:,:))
    fphi1(:,:,:) = BICG_L_ILU(mvec1(:,:,:),svec1(:,:))
    if ( icore == score ) call L_POUTJ
    end do

    if ( icore == score ) then
    call L_PREFJ
    call G_XS
    iter = iter + 1
    end if
    end do

    return
    end if

end subroutine


! =============================================================================
! ONE_NODE_CMFD
! =============================================================================
subroutine ONE_NODE_ADJOINT(k_eff,fphi2)
    use SOLVERS, only: BICG_G, BiCG_L, BICG_G_ILU, BICG_L_ILU!, SORL, SORG
    use PCMFD, only: L_PDHAT, L_PBC, L_PMATRIX, L_POUTJ, L_PSOURCE, L_PREFJ
    use MPI, only: MPI_COMM_WORLD, MPI_REAL8
    use PRECONDITIONER
    use VARIABLES, only: curr_cyc
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in), optional:: fphi2(:,:,:)
    real(8) :: error, k_eff0, k_eff1, mpie
    integer :: global, local
    integer :: iter, iter_max = 3D2 ! 2D2
    real(8) :: tt0, tt1
    logical :: cvg  ! converged


    if ( pcmfdon ) then
    if ( icore == score ) then
    call L2G
    call L_DTILDA
    call L_PDHAT
    call D_BC
    call L_PBC
    call L_PMATRIX
    call FM_TRANSPOSE
    where( fphi1 == 0 ) Mfm(:,:,:,4) = 0
    cvg = .false.
    if ( isnan(sum(Mfm)) ) then
        k_eff = 3D0
        cvg = .true.
    end if
    if ( dual_fmfd .and. present(fphi2) ) fphi1 = fphi2
    end if
    call MPI_BCAST(cvg,1,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    if ( cvg ) return
    call MPI_BCAST(Mfm,n_nodes*7,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    call ILU_DECOMPOSE(Mfm)

    iter  = 1; error = 1; cvg = .false.
    do
    if ( icore == score ) then
    ! ------------------------------- GLOBAL
    call G_PDHAT
    call G_PMATRIX
    call CM_TRANSPOSE
    call GLOBAL_ILU
    k_eff1 = k_eff


    do 
    k_eff0 = k_eff
    cphi0 = cphi1
    cm_s = cm_nf*cphi0/k_eff
    cphi1 = BiCG_G_ILU(mvec2,cm_s)
    k_eff = k_eff*sum(cm_nf*cphi1*cm_nf*cphi1) &
           / sum(cm_nf*cphi0*cm_nf*cphi1)
    if ( abs(k_eff0-k_eff)/k_eff < 1D-9 ) exit
    end do
    error = abs(k_eff-k_eff1)/k_eff
    end if

    ! convergence test
    if ( icore == score .and. ( error < 1D-9 .or. isnan(k_eff)  &
        .or. iter > iter_max ) ) cvg = .true.
    call MPI_BCAST(cvg,1,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    if ( cvg ) exit

    ! ------------------------------- LOCAL
    if ( icore == score ) call G_PINJ



    
    do local = 1, 1
    if ( icore == score ) then
    call G2L
    fphi0 = fphi1
    call L_PSOURCE(k_eff)
    end if
    call MPI_BCAST(fm_s,n_nodes,MPI_REAL8,score,MPI_COMM_WORLD,mpie)
    call LINEATION2(fm_s(:,:,:))
    fphi1(:,:,:) = BICG_L_ILU(mvec1(:,:,:),svec1(:,:))
    if ( icore == score ) call L_POUTJ
    end do

    if ( icore == score ) then
    call L_PREFJ
    call G_XS
    iter = iter + 1
    end if
    end do
    return
    end if

end subroutine


subroutine FM_TRANSPOSE
    implicit none
    real(8), allocatable:: Mfm1(:,:), imn(:,:,:)
    integer:: iid

    if ( .not. allocated(Mfm1) ) then
        allocate(Mfm1(nfm(1)*nfm(2)*nfm(3),7))
        allocate(imn(nfm(1),nfm(2),nfm(3)))

        iid  = 0
        imn  = 0
        do kk = 1, nfm(3)
        do jj = 1, nfm(2)
        do ii = 1, nfm(1)
            if ( fphi1(ii,jj,kk) == 0 ) cycle
            iid = iid + 1
            imn(ii,jj,kk) = iid
    
            ! data transfer
            Mfm1(iid,1) = Mfm(ii,jj,kk,1)
            Mfm1(iid,2) = Mfm(ii,jj,kk,2)
            Mfm1(iid,3) = Mfm(ii,jj,kk,3)
            Mfm1(iid,5) = Mfm(ii,jj,kk,5)
            Mfm1(iid,6) = Mfm(ii,jj,kk,6)
            Mfm1(iid,7) = Mfm(ii,jj,kk,7)
    
        end do
        end do
        end do
    end if

    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        if ( imn(ii,jj,kk) == 0 ) cycle
        if ( kk /= 1 ) then
        if ( imn(ii,jj,kk-1) /= 0 ) then
            Mfm(ii,jj,kk,1) = Mfm1(imn(ii,jj,kk-1),7)
        end if
        end if
        if ( jj /= 1 ) then
        if ( imn(ii,jj-1,kk) /= 0 ) then
            Mfm(ii,jj,kk,2) = Mfm1(imn(ii,jj-1,kk),6)
        end if
        end if
        if ( ii /= 1 ) then
        if ( imn(ii-1,jj,kk) /= 0 ) then
            Mfm(ii,jj,kk,3) = Mfm1(imn(ii-1,jj,kk),5)
        end if
        end if
        if ( ii /= nfm(1) ) then
        if ( imn(ii+1,jj,kk) /= 0 ) then
            Mfm(ii,jj,kk,5) = Mfm1(imn(ii+1,jj,kk),3)
        end if
        end if
        if ( jj /= nfm(2) ) then
        if ( imn(ii,jj+1,kk) /= 0 ) then
            Mfm(ii,jj,kk,6) = Mfm1(imn(ii,jj+1,kk),2)
        end if
        end if
        if ( kk /= nfm(3) ) then
        if ( imn(ii,jj,kk+1) /= 0 ) then
            Mfm(ii,jj,kk,7) = Mfm1(imn(ii,jj,kk+1),1)
        end if
        end if
    end do
    end do
    end do

end subroutine



subroutine CM_TRANSPOSE
    implicit none
    real(8), allocatable:: Mcm1(:,:), imn(:,:,:)
    integer:: iid

    if ( .not. allocated(Mcm1) ) then
        allocate(Mcm1(ncm(1)*ncm(2)*ncm(3),7))
        allocate(imn(ncm(1),ncm(2),ncm(3)))

        iid  = 0
        imn  = 0
        do kk = 1, ncm(3)
        do jj = 1, ncm(2)
        do ii = 1, ncm(1)
            if ( fphi1(ii,jj,kk) == 0 ) cycle
            iid = iid + 1
            imn(ii,jj,kk) = iid
    
            ! data transfer
            Mcm1(iid,1) = Mcm(ii,jj,kk,1)
            Mcm1(iid,2) = Mcm(ii,jj,kk,2)
            Mcm1(iid,3) = Mcm(ii,jj,kk,3)
            Mcm1(iid,5) = Mcm(ii,jj,kk,5)
            Mcm1(iid,6) = Mcm(ii,jj,kk,6)
            Mcm1(iid,7) = Mcm(ii,jj,kk,7)
    
        end do
        end do
        end do
    end if

    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        if ( imn(ii,jj,kk) == 0 ) cycle
        if ( kk /= 1 ) then
        if ( imn(ii,jj,kk-1) /= 0 ) then
            Mcm(ii,jj,kk,1) = Mcm1(imn(ii,jj,kk-1),7)
        end if
        end if
        if ( jj /= 1 ) then
        if ( imn(ii,jj-1,kk) /= 0 ) then
            Mcm(ii,jj,kk,2) = Mcm1(imn(ii,jj-1,kk),6)
        end if
        end if
        if ( ii /= 1 ) then
        if ( imn(ii-1,jj,kk) /= 0 ) then
            Mcm(ii,jj,kk,3) = Mcm1(imn(ii-1,jj,kk),5)
        end if
        end if
        if ( ii /= ncm(1) ) then
        if ( imn(ii+1,jj,kk) /= 0 ) then
            Mcm(ii,jj,kk,5) = Mcm1(imn(ii+1,jj,kk),3)
        end if
        end if
        if ( jj /= ncm(2) ) then
        if ( imn(ii,jj+1,kk) /= 0 ) then
            Mcm(ii,jj,kk,6) = Mcm1(imn(ii,jj+1,kk),2)
        end if
        end if
        if ( kk /= ncm(3) ) then
        if ( imn(ii,jj,kk+1) /= 0 ) then
            Mcm(ii,jj,kk,7) = Mcm1(imn(ii,jj,kk+1),1)
        end if
        end if
    end do
    end do
    end do

end subroutine



! =============================================================================
! ONE_NODE_CMFD
! =============================================================================
subroutine CMFD_CALCULATION(k_eff)
    use SOLVERS, only: BICG_G_ILU, BICG_L_ILU, BICGSTAB_PRE
    use PCMFD, only: L_PDHAT, L_PBC, L_PMATRIX, L_POUTJ, L_PSOURCE, L_PREFJ
    use MPI, only: MPI_COMM_WORLD, MPI_REAL8
    use PRECONDITIONER
    implicit none
    real(8), intent(inout):: k_eff
    real(8) :: error, k_eff0, k_eff1, mpie
    integer :: global, local
    integer :: iter, iter_max = 3D2 ! 2D2
    real(8) :: tt0, tt1
    logical :: cvg  ! converged

    if ( pcmfdon ) then
    if ( icore == score ) then
    call L2G
    call L_DTILDA
    call L_PDHAT
    call G_PDHAT
    call G_PMATRIX
    call GLOBAL_ILU

    iter = 1
    do
    k_eff0 = k_eff
    cphi0 = cphi1
    cm_s = cm_nf*cphi1/k_eff
    cphi1 = BiCG_G_ILU(mvec2,cm_s)
    k_eff = k_eff*sum(cm_nf*cphi1*cm_nf*cphi1) &
           / sum(cm_nf*cphi0*cm_nf*cphi1)
    if ( abs(k_eff0-k_eff)/k_eff < 1D-10 .or. iter == iter_max ) exit
    !print*, iter, k_eff, abs(k_eff0-k_eff)/k_eff
    iter = iter + 1
    end do
    !call G2L2
    end if
    end if

end subroutine




subroutine MATRIX_OUTPUT(phi1,M)
    real(8), intent(in):: phi1(:,:,:)
    real(8), intent(in):: M(:,:,:,:)
    integer:: mn(0:nfm(1)+1,0:nfm(2)+1,0:nfm(3)+1)
    integer:: num

    mn = 0
    num = 0
    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        if ( phi1(ii,jj,kk) == 0 ) cycle
        num = num + 1
        mn(ii,jj,kk) = num
    end do
    end do
    end do

    open(12,file='matrix_apr14.out')
    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        if ( mn(ii,jj,kk) == 0 ) cycle
        if ( kk /= 1 .and. mn(ii,jj,kk-1) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii,jj,kk-1), M(ii,jj,kk,1)
        if ( jj /= 1 .and. mn(ii,jj-1,kk) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii,jj-1,kk), M(ii,jj,kk,2)
        if ( ii /= 1 .and. mn(ii-1,jj,kk) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii-1,jj,kk), M(ii,jj,kk,3)
        write(12,1), mn(ii,jj,kk), mn(ii,jj,kk),   M(ii,jj,kk,4)
        if ( ii /= nfm(1) .and. mn(ii+1,jj,kk) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii+1,jj,kk), M(ii,jj,kk,5)
        if ( jj /= nfm(2) .and. mn(ii,jj+1,kk) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii,jj+1,kk), M(ii,jj,kk,6)
        if ( kk /= nfm(3) .and. mn(ii,jj,kk+1) /= 0 ) &
        write(12,1), mn(ii,jj,kk), mn(ii,jj,kk+1), M(ii,jj,kk,7)
    end do
    end do
    end do
    close(12)
    1 format(2i9,es16.8)

end subroutine


! =============================================================================
! L2G homogenizes the reactor parameters from local to global (general)
! =============================================================================
subroutine L2G
    implicit none
    real(8):: ssum(1:3)

    ! -------------------------------------------------------------------------
    ! homogenization
    !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        cphi1(ii,jj,kk) = sum(fphi1(ix0:ix1,iy0:iy1,iz0:iz1))
        if ( OUT_OF_ZZ(ii,jj) ) cycle
        cm_t(ii,jj,kk) = sum(fm_t(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
        cm_a(ii,jj,kk) = sum(fm_a(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
        cm_nf(ii,jj,kk) = sum(fm_nf(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
    end do
    end do
    end do
    !$omp end parallel do
    cmD = 1D0 / (3D0 * cm_t)
    where ( cphi1 == 0 ) 
        cmD = 0
        cm_nf = 0
    end where
    cphi1 = cphi1 / (fcr*fcr*fcz)

    ! interface diffusion coefficient
    if ( .not. pcmfdon ) then
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        cmDt(ii,jj,kk,1) = 2D0*cmD(ii,jj,kk)/dcm(1)
        cmDt(ii,jj,kk,2) = 2D0*cmD(ii,jj,kk)/dcm(1)
        cmDt(ii,jj,kk,3) = 2D0*cmD(ii,jj,kk)/dcm(2)
        cmDt(ii,jj,kk,4) = 2D0*cmD(ii,jj,kk)/dcm(2)
        cmDt(ii,jj,kk,5) = 2D0*cmD(ii,jj,kk)/dcm(3)
        cmDt(ii,jj,kk,6) = 2D0*cmD(ii,jj,kk)/dcm(3)
    end do
    end do
    end do
    else
    deltc0 = 0
    do ii = 1, ncm(1)
        if ( ii /= 1 ) &
        deltc0(ii,:,:,1) = 2D0*cmD(ii,:,:)*cmD(ii-1,:,:) / &
                          ((cmD(ii,:,:)+cmD(ii-1,:,:))*dcm(1))
        if ( ii /= ncm(1) ) &
        deltc0(ii,:,:,2) = 2D0*cmD(ii+1,:,:)*cmD(ii,:,:) / &
                          ((cmD(ii+1,:,:)+cmD(ii,:,:))*dcm(1))
    end do
    do jj = 1, ncm(2)
        if ( jj /= 1 ) &
        deltc0(:,jj,:,3) = 2D0*cmD(:,jj,:)*cmD(:,jj-1,:) / &
                          ((cmD(:,jj,:)+cmD(:,jj-1,:))*dcm(2))
        if ( jj /= ncm(2) ) &
        deltc0(:,jj,:,4) = 2D0*cmD(:,jj+1,:)*cmD(:,jj,:) / &
                          ((cmD(:,jj,:)+cmD(:,jj+1,:))*dcm(2))
    end do
    do kk = 1, ncm(3)
        if ( kk /= 1 ) &
        deltc0(:,:,kk,5) = 2D0*cmD(:,:,kk)*cmD(:,:,kk-1) / &
                          ((cmD(:,:,kk)+cmD(:,:,kk-1))*dcm(3))
        if ( kk /= ncm(3) ) &
        deltc0(:,:,kk,6) = 2D0*cmD(:,:,kk+1)*cmD(:,:,kk) / &
                          ((cmD(:,:,kk+1)+cmD(:,:,kk))*dcm(3))
    end do
    where ( isnan(deltc0) ) deltc0 = 0
    end if


    ! -------------------------------------------------------------------------
    ! surface average
    if ( .not. pcmfdon ) then
        do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
            if ( OUT_OF_ZZ1(ii,jj) ) cycle
            ! x-direction
            if ( ii /= 1 ) then
            ssum(1:2) = 0;  id(1) = id0(1)+1
            do oo = 1, fcz; id(3) = id0(3)+oo
            do nn = 1, fcr; id(2) = id0(2)+nn
                ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),1)
                ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
            end do
            end do
            cmJn(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
            cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
            cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
            cmF(ii-1,jj,kk,2)  = cmF(ii,jj,kk,1)
            end if
            ! y-direction
            if ( jj /= 1 ) then
            ssum(1:2) = 0;  id(2) = id0(2)+1
            do oo = 1, fcz; id(3) = id0(3)+oo
            do mm = 1, fcr; id(1) = id0(1)+mm
                ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),3)
                ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
            end do
            end do
            cmJn(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
            cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
            cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
            cmF(ii,jj-1,kk,4)  = cmF(ii,jj,kk,3)
            end if
            ! z-direction
            if ( kk /= 1 ) then
            ssum(1:2) = 0;  id(3) = id0(3)+1
            do mm = 1, fcr; id(1) = id0(1)+mm
            do nn = 1, fcr; id(2) = id0(2)+nn
                ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),5)
                ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
            end do
            end do
            cmJn(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
            cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
            cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
            cmF(ii,jj,kk-1,6)  = cmF(ii,jj,kk,5)
            end if
        end do
        end do
        end do

        ! boundary surfaces
        do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        !   x-direction
        do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1:2) = 0
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(1,id(2),id(3),1)
            ssum(2) = ssum(2) + fmJn(nfm(1),id(2),id(3),2)
        end do
        end do
        cmJn(1,jj,kk,1)      = ssum(1) / (fcr*fcz)
        cmJn(ncm(1),jj,kk,2) = ssum(2) / (fcr*fcz)
        end do
        !   y-direction
        do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(1:2) = 0
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),1,id(3),3)
            ssum(2) = ssum(2) + fmJn(id(1),nfm(2),id(3),4)
        end do
        end do
        cmJn(ii,1,kk,3)      = ssum(1) / (fcr*fcz)
        cmJn(ii,ncm(2),kk,4) = ssum(2) / (fcr*fcz)
        end do
        end do
        !   z-direction
        do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1:2) = 0
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),1,5)
            ssum(2) = ssum(2) + fmJn(id(1),id(2),nfm(3),6)
        end do
        end do
        cmJn(ii,jj,1,5)      = ssum(1) / (fcr*fcr)
        cmJn(ii,jj,ncm(3),6) = ssum(2) / (fcr*fcr)
        end do
        end do

    else
        if ( zigzagon ) then
        !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
        do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1

            ! inner surfaces
            cmJ0(ii,jj,kk,1) = sum(fmJ0(ix0,iy0:iy1,iz0:iz1,1))/fc1
            cmJ1(ii,jj,kk,1) = sum(fmJ1(ix0,iy0:iy1,iz0:iz1,1))/fc1

            cmJ0(ii,jj,kk,2) = sum(fmJ0(ix1,iy0:iy1,iz0:iz1,2))/fc1
            cmJ1(ii,jj,kk,2) = sum(fmJ1(ix1,iy0:iy1,iz0:iz1,2))/fc1

            cmJ0(ii,jj,kk,3) = sum(fmJ0(ix0:ix1,iy0,iz0:iz1,3))/fc1
            cmJ1(ii,jj,kk,3) = sum(fmJ1(ix0:ix1,iy0,iz0:iz1,3))/fc1

            cmJ0(ii,jj,kk,4) = sum(fmJ0(ix0:ix1,iy1,iz0:iz1,4))/fc1
            cmJ1(ii,jj,kk,4) = sum(fmJ1(ix0:ix1,iy1,iz0:iz1,4))/fc1

            cmJ0(ii,jj,kk,5) = sum(fmJ0(ix0:ix1,iy0:iy1,iz0,5))/fc2
            cmJ1(ii,jj,kk,5) = sum(fmJ1(ix0:ix1,iy0:iy1,iz0,5))/fc2

            cmJ0(ii,jj,kk,6) = sum(fmJ0(ix0:ix1,iy0:iy1,iz1,6))/fc2
            cmJ1(ii,jj,kk,6) = sum(fmJ1(ix0:ix1,iy0:iy1,iz1,6))/fc2
    
            ! boundary surfaces
            cmJn(ii,jj,kk,1) = sum(fmJn(ix0,iy0:iy1,iz0:iz1,1))/fc1
            cmJn(ii,jj,kk,2) = sum(fmJn(ix1,iy0:iy1,iz0:iz1,2))/fc1
            cmJn(ii,jj,kk,3) = sum(fmJn(ix0:ix1,iy0,iz0:iz1,3))/fc1
            cmJn(ii,jj,kk,4) = sum(fmJn(ix0:ix1,iy1,iz0:iz1,4))/fc1
            cmJn(ii,jj,kk,5) = sum(fmJn(ix0:ix1,iy0:iy1,iz0,5))/fc2
            cmJn(ii,jj,kk,6) = sum(fmJn(ix0:ix1,iy0:iy1,iz1,6))/fc2

        end do
        end do
        end do
        !$omp end parallel do
        else
        !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
        do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
            ! inner surfaces
            cmJ0(ii,jj,kk,1) = sum(fmJ0(ix0,iy0:iy1,iz0:iz1,1))/fc1
            cmJ1(ii,jj,kk,1) = sum(fmJ1(ix0,iy0:iy1,iz0:iz1,1))/fc1

            cmJ0(ii,jj,kk,2) = sum(fmJ0(ix1,iy0:iy1,iz0:iz1,2))/fc1
            cmJ1(ii,jj,kk,2) = sum(fmJ1(ix1,iy0:iy1,iz0:iz1,2))/fc1

            cmJ0(ii,jj,kk,3) = sum(fmJ0(ix0:ix1,iy0,iz0:iz1,3))/fc1
            cmJ1(ii,jj,kk,3) = sum(fmJ1(ix0:ix1,iy0,iz0:iz1,3))/fc1

            cmJ0(ii,jj,kk,4) = sum(fmJ0(ix0:ix1,iy1,iz0:iz1,4))/fc1
            cmJ1(ii,jj,kk,4) = sum(fmJ1(ix0:ix1,iy1,iz0:iz1,4))/fc1

            cmJ0(ii,jj,kk,5) = sum(fmJ0(ix0:ix1,iy0:iy1,iz0,5))/fc2
            cmJ1(ii,jj,kk,5) = sum(fmJ1(ix0:ix1,iy0:iy1,iz0,5))/fc2

            cmJ0(ii,jj,kk,6) = sum(fmJ0(ix0:ix1,iy0:iy1,iz1,6))/fc2
            cmJ1(ii,jj,kk,6) = sum(fmJ1(ix0:ix1,iy0:iy1,iz1,6))/fc2

            ! boundary surfaces
            if ( ii == 1 ) &
            cmJn(ii,jj,kk,1) = sum(fmJn(ix0,iy0:iy1,iz0:iz1,1))/fc1
            if ( ii == ncm(1) ) &
            cmJn(ii,jj,kk,2) = sum(fmJn(ix1,iy0:iy1,iz0:iz1,2))/fc1
            if ( jj == 1 ) &
            cmJn(ii,jj,kk,3) = sum(fmJn(ix0:ix1,iy0,iz0:iz1,3))/fc1
            if ( jj == ncm(2) ) &
            cmJn(ii,jj,kk,4) = sum(fmJn(ix0:ix1,iy1,iz0:iz1,4))/fc1
            if ( kk == 1 ) &
            cmJn(ii,jj,kk,5) = sum(fmJn(ix0:ix1,iy0:iy1,iz0,5))/fc2
            if ( kk == ncm(3) ) &
            cmJn(ii,jj,kk,6) = sum(fmJn(ix0:ix1,iy0:iy1,iz1,6))/fc2
        end do
        end do
        end do
        !$omp end parallel do
        end if
    end if

end subroutine


! =============================================================================
! L_DTILDA (general)
! =============================================================================
subroutine L_DTILDA
    implicit none
    
    ! inner region
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 ) then         ! x0
            fmDt(ii,jj,kk,1) = 2D0*fmD(ii,jj,kk)*fmD(ii-1,jj,kk) &
                /((fmD(ii,jj,kk)+fmD(ii-1,jj,kk))*dfm(1))
        end if
        if ( ii /= nfm(1) ) then    ! x1
            fmDt(ii,jj,kk,2) = 2D0*fmD(ii+1,jj,kk)*fmD(ii,jj,kk) &
                /((fmD(ii+1,jj,kk)+fmD(ii,jj,kk))*dfm(1))
        end if
        if ( jj /= 1 ) then         ! y0
            fmDt(ii,jj,kk,3) = 2D0*fmD(ii,jj,kk)*fmD(ii,jj-1,kk) &
                /((fmD(ii,jj,kk)+fmD(ii,jj-1,kk))*dfm(2))
        end if
        if ( jj /= nfm(2) ) then    ! y1
            fmDt(ii,jj,kk,4) = 2D0*fmD(ii,jj+1,kk)*fmD(ii,jj,kk) &
                /((fmD(ii,jj+1,kk)+fmD(ii,jj,kk))*dfm(2))
        end if
        if ( kk /= 1 ) then         ! z0
            fmDt(ii,jj,kk,5) = 2D0*fmD(ii,jj,kk)*fmD(ii,jj,kk-1) &
                /((fmD(ii,jj,kk)+fmD(ii,jj,kk-1))*dfm(3))
        end if
        if ( kk /= nfm(3) ) then    ! z1
            fmDt(ii,jj,kk,6) = 2D0*fmD(ii,jj,kk+1)*fmD(ii,jj,kk) &
                /((fmD(ii,jj,kk+1)+fmD(ii,jj,kk))*dfm(3))
        end if
    end do
    end do
    end do

    deltf0 = fmDt

end subroutine

! =============================================================================
! L_DHAT
! =============================================================================
subroutine L_DHAT(Dt,phi,Jn,Dh)
    implicit none
    real(8):: Dt(:,:,:,:), phi(:,:,:), Jn(:,:,:,:), Dh(:,:,:,:)

    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        if ( ii /= 1 )      Dh(ii,jj,kk,1) = (Jn(ii,jj,kk,1)+Dt(ii,jj,kk,1) &
            *(phi(ii,jj,kk)-phi(ii-1,jj,kk)))/(phi(ii,jj,kk)+phi(ii-1,jj,kk))
        if ( ii /= nfm(1) ) Dh(ii,jj,kk,2) = (Jn(ii,jj,kk,2)+Dt(ii,jj,kk,2) &
            *(phi(ii+1,jj,kk)-phi(ii,jj,kk)))/(phi(ii+1,jj,kk)+phi(ii,jj,kk))
        if ( jj /= 1 )      Dh(ii,jj,kk,3) = (Jn(ii,jj,kk,3)+Dt(ii,jj,kk,3) &
            *(phi(ii,jj,kk)-phi(ii,jj-1,kk)))/(phi(ii,jj,kk)+phi(ii,jj-1,kk))
        if ( jj /= nfm(2) ) Dh(ii,jj,kk,4) = (Jn(ii,jj,kk,4)+Dt(ii,jj,kk,4) &
            *(phi(ii,jj+1,kk)-phi(ii,jj,kk)))/(phi(ii,jj+1,kk)+phi(ii,jj,kk))
        if ( kk /= 1 )      Dh(ii,jj,kk,5) = (Jn(ii,jj,kk,5)+Dt(ii,jj,kk,5) &
            *(phi(ii,jj,kk)-phi(ii,jj,kk-1)))/(phi(ii,jj,kk)+phi(ii,jj,kk-1))
        if ( kk /= nfm(3) ) Dh(ii,jj,kk,6) = (Jn(ii,jj,kk,6)+Dt(ii,jj,kk,6) &
            *(phi(ii,jj,kk+1)-phi(ii,jj,kk)))/(phi(ii,jj,kk+1)+phi(ii,jj,kk))
    end do
    end do
    end do

    ! Boundary condition
    ii = 1;      Dh(ii,:,:,1) = Jn(ii,:,:,1)/phi(ii,:,:)
    ii = nfm(1); Dh(ii,:,:,2) = Jn(ii,:,:,2)/phi(ii,:,:)
    jj = 1;      Dh(:,jj,:,3) = Jn(:,jj,:,3)/phi(:,jj,:)
    jj = nfm(2); Dh(:,jj,:,4) = Jn(:,jj,:,4)/phi(:,jj,:)
    kk = 1;      Dh(:,:,kk,5) = Jn(:,:,kk,5)/phi(:,:,kk)
    kk = nfm(3); Dh(:,:,kk,6) = Jn(:,:,kk,6)/phi(:,:,kk)

end subroutine

! =============================================================================
! D_BC (general)
! =============================================================================
subroutine D_BC
    implicit none

    ! diffusion coefficient at boundary
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        fmDt(ix0,:,:,1) = 2D0*fmD(ix0,:,:)/dfm(1); deltf0(ix0,:,:,1) = 0D0
        fmDt(ix1,:,:,2) = 2D0*fmD(ix1,:,:)/dfm(1); deltf0(ix1,:,:,2) = 0D0
    end do
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        fmDt(:,iy0,:,3) = 2D0*fmD(:,iy0,:)/dfm(2); deltf0(:,iy0,:,3) = 0D0
        fmDt(:,iy1,:,4) = 2D0*fmD(:,iy1,:)/dfm(2); deltf0(:,iy1,:,4) = 0D0
    end do
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        fmDt(:,:,iz0,5) = 2D0*fmD(:,:,iz0)/dfm(3); deltf0(:,:,iz0,5) = 0D0
        fmDt(:,:,iz1,6) = 2D0*fmD(:,:,iz1)/dfm(3); deltf0(:,:,iz1,6) = 0D0
    end do

    if ( zigzagon ) then
    do ii = 1, zz_div
        fmDt(zzf1(ii)+1,zzf0(ii)+1:zzf0(ii+1),:,1) = 0
        fmDt(zzf2(ii),zzf0(ii)+1:zzf0(ii+1),:,2)   = 0
        fmDt(zzf0(ii)+1:zzf0(ii+1),zzf1(ii)+1,:,3) = 0
        fmDt(zzf0(ii)+1:zzf0(ii+1),zzf2(ii),:,4)   = 0
    end do
    end if

end subroutine


! =============================================================================
! L_BC
! =============================================================================
subroutine L_BC(Dt,Dh)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:)

    deltf1 = Dh

    ! interface boundary
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( OUT_OF_ZZ(ii,jj) ) cycle
        if ( ii /= 1 ) then;        id(1) = id0(1)+1    ! x0
        deltf1(id(1),:,:,1) = -(Dt(id(1),:,:,1)-Dh(id(1),:,:,1)) &
                            / (1D0+2D0*Dt(id(1),:,:,1))
        end if
        if ( ii /= ncm(1) ) then;   id(1) = id0(1)+fcr  ! x1
        deltf1(id(1),:,:,2) = (Dt(id(1),:,:,2)+Dh(id(1),:,:,2)) &
                            / (1D0+2D0*Dt(id(1),:,:,2))
        end if
        if ( jj /= 1 ) then;        id(2) = id0(2)+1    ! y0
        deltf1(:,id(2),:,3) = -(Dt(:,id(2),:,3)-Dh(:,id(2),:,3)) &
                            / (1D0+2D0*Dt(:,id(2),:,3))
        end if
        if ( jj /= ncm(2) ) then;   id(2) = id0(2)+fcr  ! y1
        deltf1(:,id(2),:,4) = (Dt(:,id(2),:,4)+Dh(:,id(2),:,4)) &
                            / (1D0+2D0*Dt(:,id(2),:,4))
        end if
        if ( kk /= 1 ) then;        id(3) = id0(3)+1    ! z0
        deltf1(:,:,id(3),5) = -(Dt(:,:,id(3),5)-Dh(:,:,id(3),5)) &
                            / (1D0+2D0*Dt(:,:,id(3),5))
        end if
        if ( kk /= ncm(3) ) then;   id(3) = id0(3)+fcz  ! z1
        deltf1(:,:,id(3),6) = (Dt(:,:,id(3),6)+Dh(:,:,id(3),6)) &
                            / (1D0+2D0*Dt(:,:,id(3),6))
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! L_MATRIX
! =============================================================================
subroutine L_MATRIX(Dt,Dh,abso)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), abso(:,:,:)
    real(8), allocatable:: deno(:,:,:)  ! denominator

    ! Matrix formulation

    ! -------------------------------------------------------------------------
    !   migration term
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)

        if ( kk /= 1   )    Mfm(ii,jj,kk,1) = &
                -(deltf0(ii,jj,kk,5)+deltf1(ii,jj,kk,5))/dfm(3)
        if ( jj /= 1   )    Mfm(ii,jj,kk,2) = &
                -(deltf0(ii,jj,kk,3)+deltf1(ii,jj,kk,3))/dfm(2)
        if ( ii /= 1   )    Mfm(ii,jj,kk,3) = &
                -(deltf0(ii,jj,kk,1)+deltf1(ii,jj,kk,1))/dfm(1)
        if ( ii /= nfm(1) ) Mfm(ii,jj,kk,5) = &
                -(deltf0(ii,jj,kk,2)-deltf1(ii,jj,kk,2))/dfm(1)
        if ( jj /= nfm(2) ) Mfm(ii,jj,kk,6) = &
                -(deltf0(ii,jj,kk,4)-deltf1(ii,jj,kk,4))/dfm(2)
        if ( kk /= nfm(3) ) Mfm(ii,jj,kk,7) = &
                -(deltf0(ii,jj,kk,6)-deltf1(ii,jj,kk,6))/dfm(3)
        
        Mfm(ii,jj,kk,4) = &
            +(deltf0(ii,jj,kk,1)-deltf1(ii,jj,kk,1))/dfm(1) &
            +(deltf0(ii,jj,kk,2)+deltf1(ii,jj,kk,2))/dfm(1) &
            +(deltf0(ii,jj,kk,3)-deltf1(ii,jj,kk,3))/dfm(2) &
            +(deltf0(ii,jj,kk,4)+deltf1(ii,jj,kk,4))/dfm(2) &
            +(deltf0(ii,jj,kk,5)-deltf1(ii,jj,kk,5))/dfm(3) &
            +(deltf0(ii,jj,kk,6)+deltf1(ii,jj,kk,6))/dfm(3) &
            +abso(ii,jj,kk)

    end do
    end do
    end do


    ! -------------------------------------------------------------------------
    !   source term
    allocate(deno(nfm(1),nfm(2),nfm(3)))
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( OUT_OF_ZZ(ii,jj) ) cycle
        if ( ii /= 1 ) then;        id(1) = id0(1)+1    ! x0
            deno(id(1),:,:) = 1D0+2D0*Dt(id(1),:,:,1)
            jsrc(id(1),:,:,1) = 4D0*Dt(id(1),:,:,1)/deno(id(1),:,:)
            fsrc(id(1),:,:,1) = Dh(id(1),:,:,1)/deno(id(1),:,:)
        end if
        if ( ii /= ncm(1) ) then;   id(1) = id0(1)+fcr  ! x1
            deno(id(1),:,:) = 1D0+2D0*Dt(id(1),:,:,2)
            jsrc(id(1),:,:,2) = 4D0*Dt(id(1),:,:,2)/deno(id(1),:,:)
            fsrc(id(1),:,:,2) = Dh(id(1),:,:,2)/deno(id(1),:,:)
        end if
        if ( jj /= 1 ) then;        id(2) = id0(2)+1    ! y0
            deno(:,id(2),:) = 1D0+2D0*Dt(:,id(2),:,3)
            jsrc(:,id(2),:,3) = 4D0*Dt(:,id(2),:,3)/deno(:,id(2),:)
            fsrc(:,id(2),:,3) = Dh(:,id(2),:,3)/deno(:,id(2),:)
        end if
        if ( jj /= ncm(2) ) then;   id(2) = id0(2)+fcr  ! y1
            deno(:,id(2),:) = 1D0+2D0*Dt(:,id(2),:,4)
            jsrc(:,id(2),:,4) = 4D0*Dt(:,id(2),:,4)/deno(:,id(2),:)
            fsrc(:,id(2),:,4) = Dh(:,id(2),:,4)/deno(:,id(2),:)
        end if
        if ( kk /= 1 ) then;        id(3) = id0(3)+1    ! z0
            deno(:,:,id(3)) = 1D0+2D0*Dt(:,:,id(3),5)
            jsrc(:,:,id(3),5) = 4D0*Dt(:,:,id(3),5)/deno(:,:,id(3))
            fsrc(:,:,id(3),5) = Dh(:,:,id(3),5)/deno(:,:,id(3))
        end if
        if ( kk /= ncm(3) ) then;   id(3) = id0(3)+fcz  ! z1
            deno(:,:,id(3)) = 1D0+2D0*Dt(:,:,id(3),6)
            jsrc(:,:,id(3),6) = 4D0*Dt(:,:,id(3),6)/deno(:,:,id(3))
            fsrc(:,:,id(3),6) = Dh(:,:,id(3),6)/deno(:,:,id(3))
        end if
    end do
    end do
    end do
    deallocate(deno)

end subroutine

! =============================================================================
! G_DHAT
! =============================================================================
subroutine G_DHAT(Jn,Dt,vphi,sphi,Dh)
    implicit none
    real(8), intent(in) :: Jn(:,:,:,:), Dt(:,:,:,:), vphi(:,:,:), sphi(:,:,:,:)
    real(8), intent(out):: Dh(:,:,:,:)

    ! x0 +
    Dh(:,:,:,1) = (Jn(:,:,:,1)+Dt(:,:,:,1) &
        *(vphi(:,:,:)-sphi(:,:,:,1)))/(vphi(:,:,:)+sphi(:,:,:,1))
    ! x1 -
    Dh(:,:,:,2) = (Jn(:,:,:,2)+Dt(:,:,:,2) &
        *(sphi(:,:,:,2)-vphi(:,:,:)))/(sphi(:,:,:,2)+vphi(:,:,:))
    ! y0 +
    Dh(:,:,:,3) = (Jn(:,:,:,3)+Dt(:,:,:,3) &
        *(vphi(:,:,:)-sphi(:,:,:,3)))/(vphi(:,:,:)+sphi(:,:,:,3))
    ! y1 -
    Dh(:,:,:,4) = (Jn(:,:,:,4)+Dt(:,:,:,4) &
        *(sphi(:,:,:,4)-vphi(:,:,:)))/(sphi(:,:,:,4)+vphi(:,:,:))
    ! z0 +
    Dh(:,:,:,5) = (Jn(:,:,:,5)+Dt(:,:,:,5) &
        *(vphi(:,:,:)-sphi(:,:,:,5)))/(vphi(:,:,:)+sphi(:,:,:,5))
    ! z1 -
    Dh(:,:,:,6) = (Jn(:,:,:,6)+Dt(:,:,:,6) &
        *(sphi(:,:,:,6)-vphi(:,:,:)))/(sphi(:,:,:,6)+vphi(:,:,:))

end subroutine


! =============================================================================
! G_PDHAT
! =============================================================================
subroutine G_PDHAT
    implicit none

    !$omp parallel default(shared) private(ii,jj,kk)
    !$omp do
    do ii = 2, ncm(1)
        deltc1(ii,:,:,1) = (cmJ0(ii,:,:,1)-5D-1 &
            *deltc0(ii,:,:,1)*(cphi1(ii,:,:)-cphi1(ii-1,:,:)))/cphi1(ii,:,:)
    end do
    !$omp end do
    !$omp do
    do ii = 1, ncm(1)-1
        deltc1(ii,:,:,2) = (cmJ1(ii,:,:,2)+5D-1 &
            *deltc0(ii,:,:,2)*(cphi1(ii+1,:,:)-cphi1(ii,:,:)))/cphi1(ii,:,:)
    end do
    !$omp end do
    !$omp do
    do jj = 2, ncm(2)
        deltc1(:,jj,:,3) = (cmJ0(:,jj,:,3)-5D-1 &
            *deltc0(:,jj,:,3)*(cphi1(:,jj,:)-cphi1(:,jj-1,:)))/cphi1(:,jj,:)
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2)-1
        deltc1(:,jj,:,4) = (cmJ1(:,jj,:,4)+5D-1 &
            *deltc0(:,jj,:,4)*(cphi1(:,jj+1,:)-cphi1(:,jj,:)))/cphi1(:,jj,:)
    end do
    !$omp end do
    !$omp do
    do kk = 2, ncm(3)
        deltc1(:,:,kk,5) = (cmJ0(:,:,kk,5)-5D-1 &
            *deltc0(:,:,kk,5)*(cphi1(:,:,kk)-cphi1(:,:,kk-1)))/cphi1(:,:,kk)
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3)-1
        deltc1(:,:,kk,6) = (cmJ1(:,:,kk,6)+5D-1 &
            *deltc0(:,:,kk,6)*(cphi1(:,:,kk+1)-cphi1(:,:,kk)))/cphi1(:,:,kk)
    end do
    !$omp end do
    !$omp end parallel


    if ( .not. zigzagon ) then
    ii = 1;      deltc1(ii,:,:,1) = -cmJn(ii,:,:,1)/cphi1(ii,:,:)
    ii = ncm(1); deltc1(ii,:,:,2) = +cmJn(ii,:,:,2)/cphi1(ii,:,:)
    jj = 1;      deltc1(:,jj,:,3) = -cmJn(:,jj,:,3)/cphi1(:,jj,:)
    jj = ncm(2); deltc1(:,jj,:,4) = +cmJn(:,jj,:,4)/cphi1(:,jj,:)
    else
    do ii = 1, zz_div
        deltc1(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,1) = &
            -cmJn(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,1) &
            /cphi1(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:)
        deltc1(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,2) = &
            +cmJn(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,2) &
            /cphi1(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:)
        deltc1(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,3) = &
            -cmJn(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,3) &
            /cphi1(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:)
        deltc1(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,4) = &
            +cmJn(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,4) &
            /cphi1(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:)
    end do
    end if
    kk = 1;      deltc1(:,:,kk,5) = -cmJn(:,:,kk,5)/cphi1(:,:,kk)
    kk = ncm(3); deltc1(:,:,kk,6) = +cmJn(:,:,kk,6)/cphi1(:,:,kk)

end subroutine


! =============================================================================
! G_MATRIX
! =============================================================================
subroutine G_MATRIX(Dt,Dh,Jn,phi)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), Jn(:,:,:,:), phi(:,:,:)
    real(8):: deno    ! denominator of the parameter

   ! diffusion coefficient
   do ii = 1, ncm(1)
   do jj = 1, ncm(2)
   do kk = 1, ncm(3)
       if ( OUT_OF_ZZ(ii,jj) ) cycle
       if ( ii /= 1 ) then         ! x0
       deno = Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2)
       deltc0(ii,jj,kk,1) = (Dt(ii-1,jj,kk,2)*Dt(ii,jj,kk,1) &
           +Dh(ii,jj,kk,1)*Dh(ii-1,jj,kk,2))/deno
       deltc1(ii,jj,kk,1) = (Dt(ii-1,jj,kk,2)*Dh(ii,jj,kk,1) &
           +Dt(ii,jj,kk,1)*Dh(ii-1,jj,kk,2))/deno
       end if
       if ( ii /= ncm(1) ) then    ! x1
       deno = Dt(ii,jj,kk,2)+Dt(ii+1,jj,kk,1)+Dh(ii+1,jj,kk,1)-Dh(ii,jj,kk,2)
       deltc0(ii,jj,kk,2) = (Dt(ii,jj,kk,2)*Dt(ii+1,jj,kk,1) &
           +Dh(ii+1,jj,kk,1)*Dh(ii,jj,kk,2))/deno
       deltc1(ii,jj,kk,2) = (Dt(ii,jj,kk,2)*Dh(ii+1,jj,kk,1) &
           +Dt(ii+1,jj,kk,1)*Dh(ii,jj,kk,2))/deno
       end if
       if ( jj /= 1 ) then         ! y0
       deno = Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4)
       deltc0(ii,jj,kk,3) = (Dt(ii,jj-1,kk,4)*Dt(ii,jj,kk,3) &
           +Dh(ii,jj,kk,3)*Dh(ii,jj-1,kk,4))/deno
       deltc1(ii,jj,kk,3) = (Dt(ii,jj-1,kk,4)*Dh(ii,jj,kk,3) &
           +Dt(ii,jj,kk,3)*Dh(ii,jj-1,kk,4))/deno
       end if
       if ( jj /= ncm(2) ) then    ! y1
       deno = Dt(ii,jj,kk,4)+Dt(ii,jj+1,kk,3)+Dh(ii,jj+1,kk,3)-Dh(ii,jj,kk,4)
       deltc0(ii,jj,kk,4) = (Dt(ii,jj,kk,4)*Dt(ii,jj+1,kk,3) &
           +Dh(ii,jj+1,kk,3)*Dh(ii,jj,kk,4))/deno
       deltc1(ii,jj,kk,4) = (Dt(ii,jj,kk,4)*Dh(ii,jj+1,kk,3) &
           +Dt(ii,jj+1,kk,3)*Dh(ii,jj,kk,4))/deno
       end if
       if ( kk /= 1 ) then         ! z0
       deno = Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6)
       deltc0(ii,jj,kk,5) = (Dt(ii,jj,kk-1,6)*Dt(ii,jj,kk,5) &
           +Dh(ii,jj,kk,5)*Dh(ii,jj,kk-1,6))/deno
       deltc1(ii,jj,kk,5) = (Dt(ii,jj,kk-1,6)*Dh(ii,jj,kk,5) &
           +Dt(ii,jj,kk,5)*Dh(ii,jj,kk-1,6))/deno
       end if
       if ( kk /= ncm(3) ) then    ! z1
       deno = Dt(ii,jj,kk,6)+Dt(ii,jj,kk+1,5)+Dh(ii,jj,kk+1,5)-Dh(ii,jj,kk,6)
       deltc0(ii,jj,kk,6) = (Dt(ii,jj,kk,6)*Dt(ii,jj,kk+1,5) &
           +Dh(ii,jj,kk+1,5)*Dh(ii,jj,kk,6))/deno
       deltc1(ii,jj,kk,6) = (Dt(ii,jj,kk,6)*Dh(ii,jj,kk+1,5) &
           +Dt(ii,jj,kk+1,5)*Dh(ii,jj,kk,6))/deno
       end if
   end do
   end do
   end do


   ! boundary condition
   ! - square boundary
   if ( .not. zigzagon ) then
   ii = 1;      deltc1(ii,:,:,1) = Jn(ii,:,:,1) / phi(ii,:,:)
   ii = ncm(1); deltc1(ii,:,:,2) = Jn(ii,:,:,2) / phi(ii,:,:)
   jj = 1;      deltc1(:,jj,:,3) = Jn(:,jj,:,3) / phi(:,jj,:)
   jj = ncm(2); deltc1(:,jj,:,4) = Jn(:,jj,:,4) / phi(:,jj,:)

   ! - zigzag boundary
   else
   do ii = 1, zz_div
       ! deltc0
       deltc0(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,1) = 0
       deltc0(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,2)   = 0
       deltc0(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,3) = 0
       deltc0(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,4)   = 0
       ! deltc1
       deltc1(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,1) = &
           Jn(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,1) &
           / phi(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:)
       deltc1(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,2) = &
           Jn(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,2) &
           / phi(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:)
       deltc1(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,3) = &
           Jn(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,3) &
           / phi(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:)
       deltc1(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,4) = &
           Jn(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,4) &
           / phi(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:)
   end do
   end if
   kk = 1;      deltc1(:,:,kk,5) = Jn(:,:,kk,5) / phi(:,:,kk)
   kk = ncm(3); deltc1(:,:,kk,6) = Jn(:,:,kk,6) / phi(:,:,kk)


   ! cell components
   do kk = 1, ncm(3)
   do jj = 1, ncm(2)
   do ii = 1, ncm(1)
       if ( OUT_OF_ZZ(ii,jj) ) cycle
       ! conventional FDM
       if ( kk /= 1 )      Mcm(ii,jj,kk,1) = &
           -(deltc0(ii,jj,kk,5)+deltc1(ii,jj,kk,5))/(dfm(3)*fcz)
       if ( jj /= 1 )      Mcm(ii,jj,kk,2) = &
           -(deltc0(ii,jj,kk,3)+deltc1(ii,jj,kk,3))/(dfm(2)*fcr)
       if ( ii /= 1 )      Mcm(ii,jj,kk,3) = &
           -(deltc0(ii,jj,kk,1)+deltc1(ii,jj,kk,1))/(dfm(1)*fcr)
       if ( ii /= ncm(1) ) Mcm(ii,jj,kk,5) = &
           -(deltc0(ii,jj,kk,2)-deltc1(ii,jj,kk,2))/(dfm(1)*fcr)
       if ( jj /= ncm(2) ) Mcm(ii,jj,kk,6) = &
           -(deltc0(ii,jj,kk,4)-deltc1(ii,jj,kk,4))/(dfm(2)*fcr)
       if ( kk /= ncm(3) ) Mcm(ii,jj,kk,7) = &
           -(deltc0(ii,jj,kk,6)-deltc1(ii,jj,kk,6))/(dfm(3)*fcz)
       
       Mcm(ii,jj,kk,4)= &
           +(deltc0(ii,jj,kk,1)-deltc1(ii,jj,kk,1))/(dfm(1)*fcr) &
           +(deltc0(ii,jj,kk,2)+deltc1(ii,jj,kk,2))/(dfm(1)*fcr) &
           +(deltc0(ii,jj,kk,3)-deltc1(ii,jj,kk,3))/(dfm(2)*fcr) &
           +(deltc0(ii,jj,kk,4)+deltc1(ii,jj,kk,4))/(dfm(2)*fcr) &
           +(deltc0(ii,jj,kk,5)-deltc1(ii,jj,kk,5))/(dfm(3)*fcz) &
           +(deltc0(ii,jj,kk,6)+deltc1(ii,jj,kk,6))/(dfm(3)*fcz) &
           +cm_a(ii,jj,kk)

   end do
   end do
   end do

   ! zigzag boundary
   if ( zigzagon ) then
   where ( phi(:,:,:) == 0 ) 
       Mcm(:,:,:,1) = 0
       Mcm(:,:,:,2) = 0
       Mcm(:,:,:,3) = 0
       Mcm(:,:,:,4) = 1
       Mcm(:,:,:,5) = 0
       Mcm(:,:,:,6) = 0
       Mcm(:,:,:,7) = 0
   end where
   end if

end subroutine


! =============================================================================
! G_PMATRIX
! =============================================================================
subroutine G_PMATRIX
    implicit none

    do ii = 1, ncm(1)
        if ( ii /= 1 ) &
        Mcm(ii,:,:,3) = &
            -(deltc0(ii,:,:,1)+deltc1(ii-1,:,:,2))/dcm(1)
        if ( ii /= ncm(1) ) &
        Mcm(ii,:,:,5) = &
            -(deltc0(ii,:,:,2)+deltc1(ii+1,:,:,1))/dcm(1)
    end do
    do jj = 1, ncm(2)
        if ( jj /= 1 ) &
        Mcm(:,jj,:,2) = &
            -(deltc0(:,jj,:,3)+deltc1(:,jj-1,:,4))/dcm(2)
        if ( jj /= ncm(2) ) &
        Mcm(:,jj,:,6) = &
            -(deltc0(:,jj,:,4)+deltc1(:,jj+1,:,3))/dcm(2)
    end do
    do kk = 1, ncm(3)
        if ( kk /= 1 ) &
        Mcm(:,:,kk,1) = &
            -(deltc0(:,:,kk,5)+deltc1(:,:,kk-1,6))/dcm(3)
        if ( kk /= ncm(3) ) &
        Mcm(:,:,kk,7) = &
            -(deltc0(:,:,kk,6)+deltc1(:,:,kk+1,5))/dcm(3)
    end do

    Mcm(:,:,:,4) = &
    +(deltc0(:,:,:,1)+deltc1(:,:,:,1)+(deltc0(:,:,:,2)+deltc1(:,:,:,2)))/dcm(1) &
    +(deltc0(:,:,:,3)+deltc1(:,:,:,3)+(deltc0(:,:,:,4)+deltc1(:,:,:,4)))/dcm(2) &
    +(deltc0(:,:,:,5)+deltc1(:,:,:,5)+(deltc0(:,:,:,6)+deltc1(:,:,:,6)))/dcm(3) &
    +cm_a(:,:,:)

    ! zigzag boundary
    if ( zigzagon ) then
    where ( cphi1(:,:,:) == 0 ) 
        Mcm(:,:,:,1) = 0
        Mcm(:,:,:,2) = 0
        Mcm(:,:,:,3) = 0
        Mcm(:,:,:,4) = 0
        Mcm(:,:,:,5) = 0
        Mcm(:,:,:,6) = 0
        Mcm(:,:,:,7) = 0
    end where
    do ii = 1, zz_div
        Mcm(zzc1(ii)+1,zzc0(ii)+1:zzc0(ii+1),:,3) = 0
        Mcm(zzc2(ii),zzc0(ii)+1:zzc0(ii+1),:,5) = 0
        Mcm(zzc0(ii)+1:zzc0(ii+1),zzc1(ii)+1,:,2) = 0
        Mcm(zzc0(ii)+1:zzc0(ii+1),zzc2(ii),:,6) = 0
    end do
    end if

end subroutine

! =============================================================================
! G_INJ
! =============================================================================
subroutine G_INJ(Dt,Dh,phi)
    implicit none
    real(8), intent(in):: Dt(:,:,:,:), Dh(:,:,:,:), phi(:,:,:)
    real(8):: sflux, netj

    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        if ( OUT_OF_ZZ(ii,jj) ) cycle
        ! x-direction ---------------------------------------------------------
        if ( ii /= 1 ) then
        sflux = ((Dt(ii,jj,kk,1)-Dh(ii,jj,kk,1))*phi(ii,jj,kk) &
            +(Dt(ii-1,jj,kk,2)+Dh(ii-1,jj,kk,2))*phi(ii-1,jj,kk)) &
            /(Dt(ii-1,jj,kk,2)+Dt(ii,jj,kk,1)+Dh(ii,jj,kk,1)-Dh(ii-1,jj,kk,2))
        netj = -Dt(ii,jj,kk,1)*(phi(ii,jj,kk)-sflux) &
               +Dh(ii,jj,kk,1)*(phi(ii,jj,kk)+sflux)

        cmJ0(ii-1,jj,kk,2) = 25D-2*sflux-5D-1*netj
        cmJ1(ii,jj,kk,1)   = 25D-2*sflux+5D-1*netj
        end if
        ! y-direction ---------------------------------------------------------
        if ( jj /= 1 ) then
        sflux = ((Dt(ii,jj,kk,3)-Dh(ii,jj,kk,3))*phi(ii,jj,kk) &
            +(Dt(ii,jj-1,kk,4)+Dh(ii,jj-1,kk,4))*phi(ii,jj-1,kk)) &
            /(Dt(ii,jj-1,kk,4)+Dt(ii,jj,kk,3)+Dh(ii,jj,kk,3)-Dh(ii,jj-1,kk,4))
        netj = -Dt(ii,jj,kk,3)*(phi(ii,jj,kk)-sflux) &
               +Dh(ii,jj,kk,3)*(phi(ii,jj,kk)+sflux)

        cmJ0(ii,jj-1,kk,4) = 25D-2*sflux-5D-1*netj
        cmJ1(ii,jj,kk,3)   = 25D-2*sflux+5D-1*netj
        end if
        ! z-direction ---------------------------------------------------------
        if ( kk /= 1 ) then
        sflux = ((Dt(ii,jj,kk,5)-Dh(ii,jj,kk,5))*phi(ii,jj,kk) &
            +(Dt(ii,jj,kk-1,6)+Dh(ii,jj,kk-1,6))*phi(ii,jj,kk-1)) &
            /(Dt(ii,jj,kk-1,6)+Dt(ii,jj,kk,5)+Dh(ii,jj,kk,5)-Dh(ii,jj,kk-1,6))
        netj = -Dt(ii,jj,kk,5)*(phi(ii,jj,kk)-sflux) &
               +Dh(ii,jj,kk,5)*(phi(ii,jj,kk)+sflux)

        cmJ0(ii,jj,kk-1,6) = 25D-2*sflux-5D-1*netj
        cmJ1(ii,jj,kk,5)   = 25D-2*sflux+5D-1*netj
        end if
    end do
    end do
    end do

end subroutine

! =============================================================================
! G_PINJ
! =============================================================================
subroutine G_PINJ
    implicit none

    !$omp parallel default(shared) private(ii,jj,kk)
    !$omp do
    do ii = 2, ncm(1)
        cmJ0(ii-1,:,:,2) = +5D-1*deltc0(ii-1,:,:,2) &
            *(cphi1(ii,:,:)-cphi1(ii-1,:,:))+deltc1(ii,:,:,1)*cphi1(ii,:,:)
        cmJ1(ii,:,:,1)   = -5D-1*deltc0(ii,:,:,1) &
            *(cphi1(ii,:,:)-cphi1(ii-1,:,:))+deltc1(ii-1,:,:,2)*cphi1(ii-1,:,:)
    end do
    !$omp end do

    !$omp do
    do jj = 2, ncm(2)
        cmJ0(:,jj-1,:,4) = +5D-1*deltc0(:,jj-1,:,4) &
            *(cphi1(:,jj,:)-cphi1(:,jj-1,:))+deltc1(:,jj,:,3)*cphi1(:,jj,:)
        cmJ1(:,jj,:,3)   = -5D-1*deltc0(:,jj,:,3) &
            *(cphi1(:,jj,:)-cphi1(:,jj-1,:))+deltc1(:,jj-1,:,4)*cphi1(:,jj-1,:)
    end do
    !$omp end do

    !$omp do
    do kk = 2, ncm(3)
        cmJ0(:,:,kk-1,6) = +5D-1*deltc0(:,:,kk-1,6) &
            *(cphi1(:,:,kk)-cphi1(:,:,kk-1))+deltc1(:,:,kk,5)*cphi1(:,:,kk)
        cmJ1(:,:,kk,5)   = -5D-1*deltc0(:,:,kk,5) &
            *(cphi1(:,:,kk)-cphi1(:,:,kk-1))+deltc1(:,:,kk-1,6)*cphi1(:,:,kk-1)
    end do
    !$omp end do
    !$omp end parallel

    where ( isnan(cmJ0) ) cmJ0 = 0
    where ( isnan(cmJ1) ) cmJ1 = 0

end subroutine

! =============================================================================
! G2L carries out the flux and current modulation (from GLOBAL to LOCAL)
! =============================================================================
subroutine G2L
    implicit none

    !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1

    if ( cphi1(ii,jj,kk) == 0 ) cycle

    ! flux modulation
    fphi1(ix0:ix1,iy0:iy1,iz0:iz1) = fphi1(ix0:ix1,iy0:iy1,iz0:iz1) &
    /sum(fphi1(ix0:ix1,iy0:iy1,iz0:iz1))*n_lnodes*cphi1(ii,jj,kk)

    ! incoming partial current modulation
    if ( ii /= 1 ) &
        fmJ1(ix0,iy0:iy1,iz0:iz1,1) = fmJ1(ix0,iy0:iy1,iz0:iz1,1) &
        /sum(fmJ1(ix0,iy0:iy1,iz0:iz1,1))*fc1*cmJ1(ii,jj,kk,1)
    if ( ii /= ncm(1) ) &
        fmJ0(ix1,iy0:iy1,iz0:iz1,2) = fmJ0(ix1,iy0:iy1,iz0:iz1,2) &
        /sum(fmJ0(ix1,iy0:iy1,iz0:iz1,2))*fc1*cmJ0(ii,jj,kk,2)
    if ( jj /= 1 ) &
        fmJ1(ix0:ix1,iy0,iz0:iz1,3) = fmJ1(ix0:ix1,iy0,iz0:iz1,3) &
        /sum(fmJ1(ix0:ix1,iy0,iz0:iz1,3))*fc1*cmJ1(ii,jj,kk,3)
    if ( jj /= ncm(2) ) &
        fmJ0(ix0:ix1,iy1,iz0:iz1,4) = fmJ0(ix0:ix1,iy1,iz0:iz1,4) &
        /sum(fmJ0(ix0:ix1,iy1,iz0:iz1,4))*fc1*cmJ0(ii,jj,kk,4)
    if ( kk /= 1 ) &
        fmJ1(ix0:ix1,iy0:iy1,iz0,5) = fmJ1(ix0:ix1,iy0:iy1,iz0,5) &
        /sum(fmJ1(ix0:ix1,iy0:iy1,iz0,5))*fc2*cmJ1(ii,jj,kk,5)
    if ( kk /= ncm(3) ) &
        fmJ0(ix0:ix1,iy0:iy1,iz1,6) = fmJ0(ix0:ix1,iy0:iy1,iz1,6) &
        /sum(fmJ0(ix0:ix1,iy0:iy1,iz1,6))*fc2*cmJ0(ii,jj,kk,6)

    end do
    end do
    end do
    !$omp end parallel do

end subroutine



! =============================================================================
! G2L carries out the flux and current modulation (from GLOBAL to LOCAL)
! =============================================================================
subroutine G2L2
    implicit none

    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1

    if ( cphi1(ii,jj,kk) == 0 ) cycle
    fphi1(ix0:ix1,iy0:iy1,iz0:iz1) = cphi1(ii,jj,kk)
    fm_nf(ix0:ix1,iy0:iy1,iz0:iz1) = cm_nf(ii,jj,kk)

    end do
    end do
    end do

end subroutine


! =============================================================================
! L_SOURCE
! =============================================================================
subroutine L_SOURCE(phi1,k_eff,fm_nf,fm_s,fmJ0,fmJ1)
    implicit none
    real(8), intent(inout):: fm_s(:,:,:)
    real(8), intent(in):: phi1(:,:,:), fm_nf(:,:,:), k_eff
    real(8), intent(in):: fmJ0(:,:,:,:), fmJ1(:,:,:,:)

    ! neutron fission source
    fm_s(:,:,:) = fm_nf(:,:,:)*phi1(:,:,:)/k_eff

    ! interface BC
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then; id(1) = id0(1)+1           ! x0
            fm_s(id(1),:,:) = fm_s(id(1),:,:) &
                +(jsrc(id(1),:,:,1)*fmJ1(id(1),:,:,1) &
                +fsrc(id(1),:,:,1)*phi1(id(1)-1,:,:))/dfm(1)
        end if
        if ( ii /= ncm(1) ) then; id(1) = id0(1)+fcr    ! x1
            fm_s(id(1),:,:) = fm_s(id(1),:,:) & 
                +(jsrc(id(1),:,:,2)*fmJ0(id(1),:,:,2) &
                -fsrc(id(1),:,:,2)*phi1(id(1)+1,:,:))/dfm(1)
        end if
    end do
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        if ( jj /= 1 ) then; id(2) = id0(2)+1           ! y0
            fm_s(:,id(2),:) = fm_s(:,id(2),:) &
                +(jsrc(:,id(2),:,3)*fmJ1(:,id(2),:,3) &
                +fsrc(:,id(2),:,3)*phi1(:,id(2)-1,:))/dfm(2)
        end if
        if ( jj /= ncm(2) ) then; id(2) = id0(2)+fcr    ! y1
            fm_s(:,id(2),:) = fm_s(:,id(2),:) &
                +(jsrc(:,id(2),:,4)*fmJ0(:,id(2),:,4) &
                -fsrc(:,id(2),:,4)*phi1(:,id(2)+1,:))/dfm(2)
        end if
    end do
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( kk /= 1 ) then; id(3) = id0(3)+1           ! z0
            fm_s(:,:,id(3)) = fm_s(:,:,id(3)) &
                +(jsrc(:,:,id(3),5)*fmJ1(:,:,id(3),5) &
                +fsrc(:,:,id(3),5)*phi1(:,:,id(3)-1))/dfm(3)
        end if
        if ( kk /= ncm(3) ) then; id(3) = id0(3)+fcz    ! z1
            fm_s(:,:,id(3)) = fm_s(:,:,id(3)) &
                +(jsrc(:,:,id(3),6)*fmJ0(:,:,id(3),6) &
                -fsrc(:,:,id(3),6)*phi1(:,:,id(3)+1))/dfm(3)
        end if
    end do

end subroutine

! =============================================================================
! LINEATION converts 3D matrix to 1D array
! =============================================================================
subroutine LINEATION1(Mfm)
    implicit none
    real(8), intent(in):: Mfm(:,:,:,:)
    integer:: ij

    !$omp parallel do default(shared) private(ij)
    do ij = 1, anode
        mvec(ij,1:fcr,1:fcr,1:fcz,1:7) = &
            Mfm(ax(ij)+1:ax(ij)+fcr,ay(ij)+1:ay(ij)+fcr,az(ij)+1:az(ij)+fcz,1:7)
    end do
    !$omp end parallel do

end subroutine

subroutine LINEATION2(fm_s)
    implicit none
    real(8), intent(in):: fm_s(:,:,:)
    integer:: ij

    !$omp parallel do default(shared) private(ij)
    do ij = i_para0, i_para1
        !svec1(ij,:) = reshape(fm_s(ax(ij)+1:ax(ij)+fcr,ay(ij)+1:ay(ij)+fcr,az(ij)+1:az(ij)+fcz),[1])
        svec1(ij,:) = reshape(fm_s(ax(ij)+1:ax(ij)+fcr,ay(ij)+1:ay(ij)+fcr,az(ij)+1:az(ij)+fcz),(/n_lnodes/))
    end do
    !$omp end parallel do

end subroutine

! =============================================================================
! PRECONDITIONING
! =============================================================================
subroutine PRECONDITIONING(mv,pcd)
    real(8), intent(in) :: mv(:,:,:,:)
    real(8), intent(out):: pcd(:,:,:,:)
    
    pcd(:,:,:,:) = 1D0/mv(:,:,:,:)

end subroutine


! =============================================================================
! L_OUTJ
! =============================================================================
subroutine L_OUTJ(phi0,phi1,fmF,fmJ0,fmJ1,fmJn)
    implicit none
    real(8), intent(inout):: phi0(:,:,:), phi1(:,:,:)
    real(8), intent(inout):: fmF(:,:,:,:), fmJ0(:,:,:,:)
    real(8), intent(inout):: fmJ1(:,:,:,:), fmJn(:,:,:,:)

    ! outgoing partial current
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then; id(1) = id0(1)+1           ! x0
            fmJn(id(1),:,:,1) = +jsrc(id(1),:,:,1)*fmJ1(id(1),:,:,1) &
                +deltf1(id(1),:,:,1)*phi1(id(1),:,:) &
                +fsrc(id(1),:,:,1)*phi0(id(1)-1,:,:)
            fmF(id(1),:,:,1) = 4D0*fmJ1(id(1),:,:,1)-2D0*fmJn(id(1),:,:,1)
            fmJ0(id(1),:,:,1) = fmJ1(id(1),:,:,1)-fmJn(id(1),:,:,1)
        end if
        if ( ii /= ncm(1) ) then; id(1) = id0(1)+fcr    ! x1
            fmJn(id(1),:,:,2) = -jsrc(id(1),:,:,2)*fmJ0(id(1),:,:,2) &
                +deltf1(id(1),:,:,2)*phi1(id(1),:,:) &
                +fsrc(id(1),:,:,2)*phi0(id(1)+1,:,:)
            fmF(id(1),:,:,2) = 4D0*fmJ0(id(1),:,:,2)+2D0*fmJn(id(1),:,:,2)
            fmJ1(id(1),:,:,2) = fmJ0(id(1),:,:,2)+fmJn(id(1),:,:,2)
        end if
    end do
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        if ( jj /= 1 ) then; id(2) = id0(2)+1           ! y0
            fmJn(:,id(2),:,3) = +jsrc(:,id(2),:,3)*fmJ1(:,id(2),:,3) &
                +deltf1(:,id(2),:,3)*phi1(:,id(2),:) &
                +fsrc(:,id(2),:,3)*phi0(:,id(2)-1,:)
            fmF(:,id(2),:,3) = 4D0*fmJ1(:,id(2),:,3)-2D0*fmJn(:,id(2),:,3)
            fmJ0(:,id(2),:,3) = fmJ1(:,id(2),:,3)-fmJn(:,id(2),:,3)
        end if
        if ( jj /= ncm(2) ) then; id(2) = id0(2)+fcr    ! y1
            fmJn(:,id(2),:,4) = -jsrc(:,id(2),:,4)*fmJ0(:,id(2),:,4) &
                +deltf1(:,id(2),:,4)*phi1(:,id(2),:) &
                +fsrc(:,id(2),:,4)*phi0(:,id(2)+1,:)
            fmF(:,id(2),:,4) = 4D0*fmJ0(:,id(2),:,4)+2D0*fmJn(:,id(2),:,4)
            fmJ1(:,id(2),:,4) = fmJ0(:,id(2),:,4)+fmJn(:,id(2),:,4)
        end if
    end do
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( kk /= 1 ) then; id(3) = id0(3)+1           ! z0
            fmJn(:,:,id(3),5) = +jsrc(:,:,id(3),5)*fmJ1(:,:,id(3),5) &
                +deltf1(:,:,id(3),5)*phi1(:,:,id(3)) &
                +fsrc(:,:,id(3),5)*phi0(:,:,id(3)-1)
            fmF(:,:,id(3),5) = 4D0*fmJ1(:,:,id(3),5)-2D0*fmJn(:,:,id(3),5)
            fmJ0(:,:,id(3),5) = fmJ1(:,:,id(3),5)-fmJn(:,:,id(3),5)
        end if
        if ( kk /= ncm(3) ) then; id(3) = id0(3)+fcz    ! z1
            fmJn(:,:,id(3),6) = -jsrc(:,:,id(3),6)*fmJ0(:,:,id(3),6) &
                +deltf1(:,:,id(3),6)*phi1(:,:,id(3)) &
                +fsrc(:,:,id(3),6)*phi0(:,:,id(3)+1)
            fmF(:,:,id(3),6) = 4D0*fmJ0(:,:,id(3),6)+2D0*fmJn(:,:,id(3),6)
            fmJ1(:,:,id(3),6) = fmJ0(:,:,id(3),6)+fmJn(:,:,id(3),6)
        end if
    end do


    ! data sweeping for updating next-iteration incoming partial currents
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( ii /= 1 ) then; id(1) = id0(1)+1
            fmJ0(id(1)-1,:,:,2) = fmJ0(id(1),:,:,1)
        end if
        if ( ii /= ncm(1) ) then; id(1) = id0(1)+fcr
            fmJ1(id(1)+1,:,:,1) = fmJ1(id(1),:,:,2)
        end if
    end do
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
        if ( jj /= 1 ) then; id(2) = id0(2)+1
            fmJ0(:,id(2)-1,:,4) = fmJ0(:,id(2),:,3)
        end if
        if ( jj /= ncm(2) ) then; id(2) = id0(2)+fcr
            fmJ1(:,id(2)+1,:,3) = fmJ1(:,id(2),:,4)
        end if
    end do
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
        if ( kk /= 1 ) then; id(3) = id0(3)+1
            fmJ0(:,:,id(3)-1,6) = fmJ0(:,:,id(3),5)
        end if
        if ( kk /= ncm(3) ) then; id(3) = id0(3)+fcz
            fmJ1(:,:,id(3)+1,5) = fmJ1(:,:,id(3),6)
        end if
    end do

end subroutine


! =============================================================================
! L_REFJ
! =============================================================================
subroutine L_REFJ(fmF,fmJ0,fmJ1,fmJn,phi1)
    implicit none
    real(8), intent(in), dimension(:,:,:,:):: fmF, fmJ0, fmJ1
    real(8), intent(inout), dimension(:,:,:,:):: fmJn
    real(8), intent(in):: phi1(:,:,:)
    real(8):: ssum(2)

    ! boundary surface
    ii = 1;      fmJn(ii,:,:,1) = deltf1(ii,:,:,1)*phi1(ii,:,:)
    ii = nfm(1); fmJn(ii,:,:,2) = deltf1(ii,:,:,2)*phi1(ii,:,:)
    jj = 1;      fmJn(:,jj,:,3) = deltf1(:,jj,:,3)*phi1(:,jj,:)
    jj = nfm(2); fmJn(:,jj,:,4) = deltf1(:,jj,:,4)*phi1(:,jj,:)
    kk = 1;      fmJn(:,:,kk,5) = deltf1(:,:,kk,5)*phi1(:,:,kk)
    kk = nfm(3); fmJn(:,:,kk,6) = deltf1(:,:,kk,6)*phi1(:,:,kk)

    ! net current & surface average
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
        if ( OUT_OF_ZZ1(ii,jj) ) cycle
        ! x0
        if ( ii /= 1 ) then
        ssum = 0;       id(1) = id0(1)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),1)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),1)
        end do
        end do
        cmJn(ii,jj,kk,1) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,1)  = ssum(2) / (fcr*fcz)
        end if
        ! x1
        if ( ii /= ncm(1) ) then
        ssum = 0;       id(1) = id0(1)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),2)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),2)
        end do
        end do
        cmJn(ii,jj,kk,2) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,2)  = ssum(2) / (fcr*fcz)
        end if
        ! y0
        if ( jj /= 1 ) then
        ssum = 0;       id(2) = id0(2)+1
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),3)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),3)
        end do
        end do
        cmJn(ii,jj,kk,3) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,3)  = ssum(2) / (fcr*fcz)
        end if
        ! y1
        if ( jj /= ncm(2) ) then
        ssum = 0;       id(2) = id0(2)+fcr
        do oo = 1, fcz; id(3) = id0(3)+oo
        do mm = 1, fcr; id(1) = id0(1)+mm
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),4)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),4)
        end do
        end do
        cmJn(ii,jj,kk,4) = ssum(1) / (fcr*fcz)
        cmF(ii,jj,kk,4)  = ssum(2) / (fcr*fcz)
        end if
        ! z0
        if ( kk /= 1 ) then
        ssum = 0;       id(3) = id0(3)+1
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),5)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),5)
        end do
        end do
        cmJn(ii,jj,kk,5) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,5)  = ssum(2) / (fcr*fcr)
        end if
        ! z1
        if ( kk /= ncm(3) ) then
        ssum = 0;       id(3) = id0(3)+fcz
        do mm = 1, fcr; id(1) = id0(1)+mm
        do nn = 1, fcr; id(2) = id0(2)+nn
            ssum(1) = ssum(1) + fmJn(id(1),id(2),id(3),6)
            ssum(2) = ssum(2) + fmF(id(1),id(2),id(3),6)
        end do
        end do
        cmJn(ii,jj,kk,6) = ssum(1) / (fcr*fcr)
        cmF(ii,jj,kk,6)  = ssum(2) / (fcr*fcr)
        end if
    end do
    end do
    end do

    ! interface surface
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        ! x-direction
        if ( ii /= 1 ) then
            cmJn(ii,jj,kk,1) = (cmJn(ii,jj,kk,1)+cmJn(ii-1,jj,kk,2))/2D0
            cmJn(ii-1,jj,kk,2) = cmJn(ii,jj,kk,1)
            cmF(ii,jj,kk,1) = (cmF(ii,jj,kk,1)+cmF(ii-1,jj,kk,2))/2D0
            cmF(ii-1,jj,kk,2) = cmF(ii,jj,kk,1)
        end if
        ! y-direction
        if ( jj /= 1 ) then
            cmJn(ii,jj,kk,3) = (cmJn(ii,jj,kk,3)+cmJn(ii,jj-1,kk,4))/2D0
            cmJn(ii,jj-1,kk,4) = cmJn(ii,jj,kk,3)
            cmF(ii,jj,kk,3) = (cmF(ii,jj,kk,3)+cmF(ii,jj-1,kk,4))/2D0
            cmF(ii,jj-1,kk,4) = cmF(ii,jj,kk,3)
        end if
        ! z-direction
        if ( kk /= 1 ) then
            cmJn(ii,jj,kk,5) = (cmJn(ii,jj,kk,5)+cmJn(ii,jj,kk-1,6))/2D0
            cmJn(ii,jj,kk-1,6) = cmJn(ii,jj,kk,5)
            cmF(ii,jj,kk,5) = (cmF(ii,jj,kk,5)+cmF(ii,jj,kk-1,6))/2D0
            cmF(ii,jj,kk-1,6) = cmF(ii,jj,kk,5)
        end if
    end do
    end do
    end do

    ! boundary surfaces
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz
    !   x-direction
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1:2) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(1,id(2),id(3),1)
        ssum(2) = ssum(2) + fmJn(nfm(1),id(2),id(3),2)
    end do
    end do
    cmJn(1,jj,kk,1)      = ssum(1) / (fcr*fcz)
    cmJn(ncm(1),jj,kk,2) = ssum(2) / (fcr*fcz)
    end do
    !   y-direction
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; ssum(1:2) = 0
    do oo = 1, fcz; id(3) = id0(3)+oo
    do mm = 1, fcr; id(1) = id0(1)+mm
        ssum(1) = ssum(1) + fmJn(id(1),1,id(3),3)
        ssum(2) = ssum(2) + fmJn(id(1),nfm(2),id(3),4)
    end do
    end do
    cmJn(ii,1,kk,3)      = ssum(1) / (fcr*fcz)
    cmJn(ii,ncm(2),kk,4) = ssum(2) / (fcr*fcz)
    end do
    end do
    !   z-direction
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; ssum(1:2) = 0
    do mm = 1, fcr; id(1) = id0(1)+mm
    do nn = 1, fcr; id(2) = id0(2)+nn
        ssum(1) = ssum(1) + fmJn(id(1),id(2),1,5)
        ssum(2) = ssum(2) + fmJn(id(1),id(2),nfm(3),6)
    end do
    end do
    cmJn(ii,jj,1,5)      = ssum(1) / (fcr*fcr)
    cmJn(ii,jj,ncm(3),6) = ssum(2) / (fcr*fcr)
    end do
    end do


end subroutine
    
    
! =============================================================================
! G_XS produces the flux-volume-weight group constants
! ============================================================================= 
subroutine G_XS
    implicit none

    ! homogenization
    !$omp parallel do default(shared) private(ii,jj,kk,ix0,ix1,iy0,iy1,iz0,iz1)
    do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
    do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
    do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
        cphi1(ii,jj,kk) = sum(fphi1(ix0:ix1,iy0:iy1,iz0:iz1))
        cm_t(ii,jj,kk) = sum(fm_t(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
        cm_a(ii,jj,kk) = sum(fm_a(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
        cm_nf(ii,jj,kk) = sum(fm_nf(ix0:ix1,iy0:iy1,iz0:iz1) &
            *fphi1(ix0:ix1,iy0:iy1,iz0:iz1))/cphi1(ii,jj,kk)
    end do
    end do
    end do
    !$omp end parallel do
    cmD = 1D0 / (3D0 * cm_t)
    where ( cphi1 == 0 ) 
        cmD = 0
        cm_nf = 0
    end where
    cphi1 = cphi1 / (fcr*fcr*fcz)

    ! interface diffusion coefficient
    if ( .not. pcmfdon ) then
    do ii = 1, ncm(1)
    do jj = 1, ncm(2)
    do kk = 1, ncm(3)
        cmDt(ii,jj,kk,1) = 2D0*cmD(ii,jj,kk)/(dcm(1))
        cmDt(ii,jj,kk,2) = 2D0*cmD(ii,jj,kk)/(dcm(1))
        cmDt(ii,jj,kk,3) = 2D0*cmD(ii,jj,kk)/(dcm(2))
        cmDt(ii,jj,kk,4) = 2D0*cmD(ii,jj,kk)/(dcm(2))
        cmDt(ii,jj,kk,5) = 2D0*cmD(ii,jj,kk)/(dcm(3))
        cmDt(ii,jj,kk,6) = 2D0*cmD(ii,jj,kk)/(dcm(3))
    end do
    end do
    end do
    else
    deltc0 = 0
    !$omp parallel default(shared) private(ii,jj,kk)
    !$omp do
    do ii = 1, ncm(1)
        if ( ii /= 1 ) &
        deltc0(ii,:,:,1) = 2D0*cmD(ii,:,:)*cmD(ii-1,:,:) / &
                          ((cmD(ii,:,:)+cmD(ii-1,:,:))*dcm(1))
        if ( ii /= ncm(1) ) &
        deltc0(ii,:,:,2) = 2D0*cmD(ii+1,:,:)*cmD(ii,:,:) / &
                          ((cmD(ii+1,:,:)+cmD(ii,:,:))*dcm(1))
    end do
    !$omp end do
    !$omp do
    do jj = 1, ncm(2)
        if ( jj /= 1 ) &
        deltc0(:,jj,:,3) = 2D0*cmD(:,jj,:)*cmD(:,jj-1,:) / &
                          ((cmD(:,jj,:)+cmD(:,jj-1,:))*dcm(2))
        if ( jj /= ncm(2) ) &
        deltc0(:,jj,:,4) = 2D0*cmD(:,jj+1,:)*cmD(:,jj,:) / &
                          ((cmD(:,jj,:)+cmD(:,jj+1,:))*dcm(2))
    end do
    !$omp end do
    !$omp do
    do kk = 1, ncm(3)
        if ( kk /= 1 ) &
        deltc0(:,:,kk,5) = 2D0*cmD(:,:,kk)*cmD(:,:,kk-1) / &
                          ((cmD(:,:,kk)+cmD(:,:,kk-1))*dcm(3))
        if ( kk /= ncm(3) ) &
        deltc0(:,:,kk,6) = 2D0*cmD(:,:,kk+1)*cmD(:,:,kk) / &
                          ((cmD(:,:,kk+1)+cmD(:,:,kk))*dcm(3))
    end do
    !$omp end do
    !$omp end parallel
    end if

!    if ( curr_cyc > n_inact ) then
!        write(*,*), "QQQQQQQQQQQQQQQQQQ FLUX QQQQQQQQQQQQQQQQQQQQQ"
!        write(*,1), cphi1
!        1 format(<ncm(1)>es15.7)
!    end if


end subroutine

end module
