module FMFD_HEADER
     implicit none

    ! ==== surface numbering ====
    !  1    2    3    4    5    6
    ! x0 / x1 / y0 / y1 / z0 / z1
    ! ===========================

    ! indice
    integer:: id(3), id0(3), id1(3)
    integer:: mm, nn, oo    ! find mesh
    integer:: ii, jj, kk    ! coarse mesh

    ! =========================================================================
    ! FMFD & DTMC Calculation
    real(8):: k_mprup = 1D0
    real(8), allocatable :: k_fmfd(:,:)         ! (batches,cycles)
    real(8), allocatable :: k_fmfd2(:,:)        ! (batches,cycles)
    real(8), allocatable :: p_fmfd(:,:,:,:,:)   ! (batches,cycles,x,y,z)
    real(8), allocatable :: p_fmfd2(:,:,:,:,:)  ! (batches,cycles,x,y,z)
    real(8), allocatable :: p_dtmc(:,:,:,:)     ! (i,j,k,ring) power
    real(8), allocatable :: f_dtmc(:,:,:,:)     ! (i,j,k,ring) flux
    !! PERTURBATION
    real(8), allocatable :: k_real(:,:,:) ! (batches, cycles, n_pert)
    type :: FMFD_parameters
        real(8) :: phi
        real(8) :: sig_t 
        real(8) :: sig_a 
        real(8) :: nusig_f 
        real(8) :: kappa
        real(8) :: J0(8) ! LINKPOINT
        real(8) :: J1(8)
    end type
    logical :: fmfdon = .false.
    logical :: pfmfdon = .false.
    logical :: fake_MC = .false.
    logical :: fmfd2mc = .false.
    logical :: quarter = .false.
    logical :: inactive_cmfd = .false.
    integer :: n_fake
    integer :: n_skip, n_acc, acc_skip
    integer :: FMFD_type        ! 1 FMFD / 2 p-FMFD / 3 1-node CMFD
    real(8) :: a_fm(6), v_fm
    
    type :: FMFD_accumulation
        type(FMFD_parameters), allocatable :: fm(:,:,:)
    endtype 

    ! fine mesh parameters
    real(8), allocatable, dimension(:,:,:):: &
        fm_t, &     ! total XS
        fmD, &      ! diffusion coefficient
        fm_a, &     ! absorption XS
        fm_nf, &    ! nu x fission XS
        kappa, &    ! kappa x fission XS
        fm_s, &     ! neutrons source
        fphi0, &    ! fine mesh flux 0
        fphi1, &    ! fine mesh flux 1
        tmpfphi1

    real(8), allocatable, dimension(:,:,:,:):: &
        deltf0, &   ! delta tilda
        deltf1, &   ! delta hat
        fmJ0, &     ! partial current -
        fmJ1, &     ! partial current +
        fmJn, &     ! net current
        fmF, &      ! surface flux
        fmDt, &     ! D tilda
        fmDh, &     ! D hat
        Mfm         ! matrix elements for FM
    real(8), allocatable, dimension(:,:,:,:,:):: ptJn, ptJ0, ptJ1
    real(8), allocatable:: ptphi(:,:,:,:), ptkap(:,:,:,:)

    
    type(FMFD_parameters), allocatable :: fm(:,:,:)
    type(FMFD_parameters), allocatable :: fm_avg(:,:,:)
    type(FMFD_parameters), allocatable :: fm_thread(:,:,:)
    !$OMP THREADPRIVATE(fm_thread)
    type(FMFD_accumulation), allocatable, target :: acc(:)

    ! FMFD grid
    real(8):: fm0(3)    ! (x0,y0,z0)
    real(8):: fm1(3)    ! (x1,y1,z1)
    real(8):: fm2(3)    ! (x1-x0,y1-y0,z1-z0)
    integer:: nfm(3)    ! (nx,ny,nz)
    real(8):: dfm(3)    ! (dx,dy,dz)
    ! ---
    integer, allocatable:: lx0(:), lx1(:), ly0(:), ly1(:), lz0(:), lz1(:)
    integer, allocatable:: gx0(:), gx1(:), gy0(:), gy1(:), gz0(:), gz1(:)

    ! Zigzag
    logical:: zigzagon = .false.
    integer:: n_zz      ! # of zigzag points
    integer:: zz_div    ! # of divisions
    integer:: n_nodes
    ! - for fine mesh
    integer, allocatable:: zzf0(:) ! reference points
    integer, allocatable:: zzf1(:) ! lower points for zz0
    integer, allocatable:: zzf2(:) ! upper points for zz0
    ! - for coarse mesh
    integer, allocatable:: zzc0(:) ! reference points
    integer, allocatable:: zzc1(:) ! lower points for zz0
    integer, allocatable:: zzc2(:) ! upper points for zz0
    ! - active core
    logical:: wholecore = .false.
    integer:: nam(3)
    integer, allocatable:: zza0(:)  ! coarse mesh
    integer, allocatable:: zza1(:)
    integer, allocatable:: zza2(:)
    integer, allocatable:: zzb0(:)  ! fine mesh
    integer, allocatable:: zzb1(:)
    integer, allocatable:: zzb2(:)


    ! fission source distribution
    real(8), allocatable:: fsd_mc(:,:,:)
    real(8), allocatable:: fsd_fm(:,:,:)
    real(8), allocatable:: fsd(:,:,:)
    real(8), allocatable:: fsd2(:,:,:)

    ! =========================================================================
    ! One-Node CMFD Acceleration
    logical :: cmfdon = .false.
    logical :: pcmfdon = .false.

    type CMFD_PARAMETERS
        real(8):: phi
        real(8):: sig_t
        real(8):: sig_a
        real(8):: nusig_f
        real(8):: Jn(6)
        real(8):: J0(6)
        real(8):: J1(6)
    end type

    type(CMFD_PARAMETERS), allocatable:: cm(:,:,:)

    real(8), allocatable, dimension(:,:,:):: &
        cm_t, &     ! total XS
        cmD, &      ! diffusion coefficient
        cm_a, &     ! absorption XS
        cm_nf, &    ! nu x fission XS
        cm_s, &     ! neutrons source
        cphi0, &    ! coarse flux 0
        cphi1       ! coarse flux 1

    real(8), allocatable, dimension(:,:,:,:):: &
        deltc0, &   ! delta tilda
        deltc1, &   ! delta hat
        jsrc, &     ! current source for LOCAL
        fsrc, &     ! flux source for LOCAL
        cmJ0, &     ! partial current -
        cmJ1, &     ! partial current +
        cmJn, &     ! net current
        cmF, &      ! surface flux
        cmDt, &     ! D tilda
        cmDh, &     ! D hat
        Mcm         ! matrix elements for CM

    ! CMFD grid
    integer:: ncm(3)    ! (nx,ny,nz)
    real(8):: dcm(3)    ! size of coarse node
    integer:: fcr       ! fine > coarse in r-direction
    integer:: fcz       ! fine > coarse in z-direction
    integer:: fc1, fc2  ! # of fine cells

    ! Parallel calculation for one-node CMFD
    integer:: anode     ! # of active nodes
    integer:: n_lnodes  ! # of local nodes
    integer:: bs0, bs1  ! buffer size
    integer, allocatable:: ax(:), ay(:), az(:)
    real(8), allocatable:: mvec(:,:,:,:,:)  ! matrix vector
    real(8), allocatable:: mvec1(:,:,:)     ! matrix vector
    real(8), allocatable:: mvec2(:,:)       ! matrix vector
    real(8), allocatable:: svec(:,:,:,:)    ! source vector
    real(8), allocatable:: svec1(:,:)       ! source vector
    real(8), allocatable:: svec2(:)         ! source vector
    real(8), allocatable:: pcd(:,:,:,:)     ! preconditioner
    integer:: i_para0, i_para1
    integer, allocatable:: gcn(:,:,:)

    ! Dual FMFD system
    logical:: dual_fmfd = .false.
    real(8), allocatable:: fwgt(:)

    ! For depletion calculation
    integer:: nrings
    integer, allocatable:: buuniv(:,:,:)
    integer, allocatable:: buring(:,:,:)
    real(8), allocatable:: thflux(:,:,:,:)  ! thread-wise
    real(8), allocatable:: thkapa(:,:,:,:)
    !$OMP THREADPRIVATE(thflux,thkapa)
    real(8), allocatable:: tmflux(:,:,:,:)  ! temporary
    real(8), allocatable:: tmkapa(:,:,:,:)  
    real(8), allocatable:: buflux(:,:,:,:)  ! flux for burnup calculation
    real(8), allocatable:: bukapa(:,:,:,:)  ! kappa x sigma_f for burnup calc
    type CELL_FOR_DEPLETION
        integer:: univ
        integer:: ncell
        real(8), allocatable:: flux(:), kapa(:)
    end type
    type(CELL_FOR_DEPLETION), allocatable:: cede(:,:,:)
    type(CELL_FOR_DEPLETION), allocatable:: thde(:,:,:) ! OMP threads
    !$OMP THREADPRIVATE(thde)


    logical:: dep_dtmc
    real(8):: r_dep0, r_dep1, r_dep2  ! radius for pin depletion
    real(8), allocatable:: dep_sigt(:,:,:,:), &
                           dep_siga(:,:,:,:), &
                           dep_nufi(:,:,:,:), &
                           dep_J0(:,:,:,:,:), &
                           dep_J1(:,:,:,:,:), &
                           dep_phi1(:,:,:,:), &
                           dep_Jn(:,:,:,:,:)
    real(8), allocatable:: dth_sigt(:,:,:,:), &
                           dth_siga(:,:,:,:), &
                           dth_nufi(:,:,:,:), &
                           dth_J0(:,:,:,:,:), &
                           dth_J1(:,:,:,:,:), &
                           dth_phi1(:,:,:,:)

    logical:: DTMCBU = .false.
    logical:: MCBU = .false.
    logical:: BUGROUP = .false.
    integer:: n_rings
    real(8), allocatable:: rings(:), v_ring(:)
    integer, allocatable:: butype(:,:,:)        ! (i,j,k) burnup type
    real(8), allocatable:: intra_phi(:,:,:,:)   ! (i,j,k,ring) flux
    real(8), allocatable:: intra_pow(:,:,:,:)   ! (i,j,k,ring) power
    real(8), allocatable:: intra_kap(:,:,:,:,:) ! (cyc,i,j,k,ring) kappa
    real(8), allocatable:: p_dep_mc(:,:,:,:)    ! (cyc,i,j,k) power
    real(8), allocatable:: p_dep_dt(:,:,:,:)    ! (cyc,i,j,k) power
    real(8), allocatable:: p_dep_dt_pert(:,:,:,:,:) ! (cyc,n_pert,i,j,k) power
!    integer:: n_rtypes
!    type:: BU_CELL_FOR_RECONSTRUCTION
!        integer:: rtype
!        integer:: n_rings
!        real(8), allocatable:: rings(:)
!        real(8), allocatable:: v_ring(:)
!    end type
!    type(BU_CELL_FOR_RECONSTRUCTION), allocatable:: buce(:,:,:)
    ! --- for thread (OpenMP parallel)
    real(8), allocatable:: dth_phi(:,:,:,:)   ! flux
    real(8), allocatable:: dth_pow(:,:,:,:)   ! power
    !$OMP THREADPRIVATE(dth_phi,dth_pow)


    ! ILU decomposition for p-FMFD
    integer:: n_fnodes  ! # of filled nodes (no vacancy in zigzag region)
    integer, allocatable:: fx0(:), fx1(:), fy0(:), fy1(:), fz0(:), fz1(:)
    integer, allocatable:: fn(:,:,:)
    real(8), allocatable:: mvec3(:,:)
    !real(8) :: k_pert
    
    public :: SET_ZERO_M, PFMFD_MATRIX, D_PHAT_CALCULATION, D_HAT_CALCULATION, D_TILDA_CALCULATION  
    contains

! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_TILDA_CALCULATION
    implicit none
    
    fmDt = 0

    ! inner region
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 )      fmDt(ii,jj,kk,1) = 2D0*fmD(ii,jj,kk) &
            *fmD(ii-1,jj,kk)/((fmD(ii,jj,kk)+fmD(ii-1,jj,kk))*dfm(1))
        if ( ii /= nfm(1) ) fmDt(ii,jj,kk,2) = 2D0*fmD(ii+1,jj,kk) &
            *fmD(ii,jj,kk)/((fmD(ii+1,jj,kk)+fmD(ii,jj,kk))*dfm(1))
        if ( jj /= 1 )      fmDt(ii,jj,kk,3) = 2D0*fmD(ii,jj,kk) &
            *fmD(ii,jj-1,kk)/((fmD(ii,jj,kk)+fmD(ii,jj-1,kk))*dfm(2))
        if ( jj /= nfm(2) ) fmDt(ii,jj,kk,4) = 2D0*fmD(ii,jj+1,kk) &
            *fmD(ii,jj,kk)/((fmD(ii,jj+1,kk)+fmD(ii,jj,kk))*dfm(2))
        if ( kk /= 1 )      fmDt(ii,jj,kk,5) = 2D0*fmD(ii,jj,kk) &
            *fmD(ii,jj,kk-1)/((fmD(ii,jj,kk)+fmD(ii,jj,kk-1))*dfm(3))
        if ( kk /= nfm(3) ) fmDt(ii,jj,kk,6) = 2D0*fmD(ii,jj,kk+1) &
            *fmD(ii,jj,kk)/((fmD(ii,jj,kk+1)+fmD(ii,jj,kk))*dfm(3))
    end do
    end do
    end do

    where ( isnan(fmDt) ) fmDt = 0

end subroutine

subroutine D_HAT_CALCULATION
    integer :: i, j, k

    ! inner region
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        fmDh(i,j,k,1) = (fmJn(i,j,k,1)+fmDt(i,j,k,1) &
            *(fphi1(i,j,k)-fphi1(i-1,j,k)))/(fphi1(i,j,k)+fphi1(i-1,j,k))
        if ( i /= nfm(1) ) & ! x1
        fmDh(i,j,k,2) = (fmJn(i,j,k,2)+fmDt(i,j,k,2) &
            *(fphi1(i+1,j,k)-fphi1(i,j,k)))/(fphi1(i+1,j,k)+fphi1(i,j,k))
        if ( j /= 1 ) &      ! y0
        fmDh(i,j,k,3) = (fmJn(i,j,k,3)+fmDt(i,j,k,3) &
            *(fphi1(i,j,k)-fphi1(i,j-1,k)))/(fphi1(i,j,k)+fphi1(i,j-1,k))
        if ( j /= nfm(2) ) & ! y1
        fmDh(i,j,k,4) = (fmJn(i,j,k,4)+fmDt(i,j,k,4) &
            *(fphi1(i,j+1,k)-fphi1(i,j,k)))/(fphi1(i,j+1,k)+fphi1(i,j,k))
        if ( k /= 1 ) &      ! y0
        fmDh(i,j,k,5) = (fmJn(i,j,k,5)+fmDt(i,j,k,5) &
            *(fphi1(i,j,k)-fphi1(i,j,k-1)))/(fphi1(i,j,k)+fphi1(i,j,k-1))
        if ( k /= nfm(3) ) & ! y1
        fmDh(i,j,k,6) = (fmJn(i,j,k,6)+fmDt(i,j,k,6) &
            *(fphi1(i,j,k+1)-fphi1(i,j,k)))/(fphi1(i,j,k+1)+fphi1(i,j,k))
    end do
    end do
    end do

    ! boundary
    i = 1;      fmDh(i,:,:,1) = fmJn(i,:,:,1)/fphi1(i,:,:)
    i = nfm(1); fmDh(i,:,:,2) = fmJn(i,:,:,2)/fphi1(i,:,:)
    j = 1;      fmDh(:,j,:,3) = fmJn(:,j,:,3)/fphi1(:,j,:)
    j = nfm(2); fmDh(:,j,:,4) = fmJn(:,j,:,4)/fphi1(:,j,:)
    k = 1;      fmDh(:,:,k,5) = fmJn(:,:,k,5)/fphi1(:,:,k)
    k = nfm(3); fmDh(:,:,k,6) = fmJn(:,:,k,6)/fphi1(:,:,k)

end subroutine


! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_PHAT_CALCULATION
    integer :: i, j, k

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
        !write(*,*) 'TTT', i,j,k, fphi1(i,j,k), fmDh(i,j,k,1:6)
    end do
    end do
    end do

    ! boundary
    if ( .not. zigzagon ) then
    i = 1;      fmDh(i,:,:,1) = -fmJn(i,:,:,1)/fphi1(i,:,:)
    i = nfm(1); fmDh(i,:,:,2) = +fmJn(i,:,:,2)/fphi1(i,:,:)
    j = 1;      fmDh(:,j,:,3) = -fmJn(:,j,:,3)/fphi1(:,j,:)
    j = nfm(2); fmDh(:,j,:,4) = +fmJn(:,j,:,4)/fphi1(:,j,:)
    else
    do i = 1, zz_div
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

subroutine PFMFD_MATRIX
    integer :: i, j, k
    
    ! initialization
    Mfm = 0
    
    ! Mfm matrix set 
    do k = 1, nfm(3)
    do j = 1, nfm(2)
    do i = 1, nfm(1)
        if ( i /= 1 ) &         ! x0
            Mfm(i,j,k,3) = -(fmDt(i,j,k,1)+fmDh(i-1,j,k,2))/dfm(1)
        if ( i /= nfm(1) ) &    ! x1
            Mfm(i,j,k,5) = -(fmDt(i,j,k,2)+fmDh(i+1,j,k,1))/dfm(1)
        if ( j /= 1 ) &         ! y0
            Mfm(i,j,k,2) = -(fmDt(i,j,k,3)+fmDh(i,j-1,k,4))/dfm(2)
        if ( j /= nfm(2) ) &    ! y1
            Mfm(i,j,k,6) = -(fmDt(i,j,k,4)+fmDh(i,j+1,k,3))/dfm(2)
        if ( k /= 1 ) &         ! z0
            Mfm(i,j,k,1) = -(fmDt(i,j,k,5)+fmDh(i,j,k-1,6))/dfm(3)
        if ( k /= nfm(3) ) &    ! z1
            Mfm(i,j,k,7) = -(fmDt(i,j,k,6)+fmDh(i,j,k+1,5))/dfm(3)
        
        Mfm(i,j,k,4) = &
            +(fmDt(i,j,k,1)+fmDh(i,j,k,1)+fmDt(i,j,k,2)+fmDh(i,j,k,2))/dfm(1) &
            +(fmDt(i,j,k,3)+fmDh(i,j,k,3)+fmDt(i,j,k,4)+fmDh(i,j,k,4))/dfm(2) &
            +(fmDt(i,j,k,5)+fmDh(i,j,k,5)+fmDt(i,j,k,6)+fmDh(i,j,k,6))/dfm(3) &
            +fm_a(i,j,k)

    enddo
    enddo
    enddo

end subroutine PFMFD_MATRIX


subroutine SET_ZERO_M
    implicit none
    integer:: ij

    where( fphi1 == 0 ) 
    Mfm(:,:,:,1) = 0
    Mfm(:,:,:,2) = 0
    Mfm(:,:,:,3) = 0
    Mfm(:,:,:,4) = 0
    Mfm(:,:,:,5) = 0
    Mfm(:,:,:,6) = 0
    Mfm(:,:,:,7) = 0
    end where

    do ij = 1, zz_div
        Mfm(zzf1(ij)+1,zzf0(ij)+1:zzf0(ij+1),:,3) = 0
        Mfm(zzf2(ij),zzf0(ij)+1:zzf0(ij+1),:,5) = 0
        Mfm(zzf0(ij)+1:zzf0(ij+1),zzf1(ij)+1,:,2) = 0
        Mfm(zzf0(ij)+1:zzf0(ij+1),zzf2(ij),:,6) = 0
    end do

end subroutine

end module
