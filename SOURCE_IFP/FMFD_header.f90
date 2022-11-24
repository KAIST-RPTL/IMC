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
    real(8), allocatable :: k_fmfd(:,:)         ! (batches,cycles)
    real(8), allocatable :: p_fmfd(:,:,:,:,:)   ! (batches,cycles,x,y,z)
    real(8), allocatable :: e_fmfd(:,:,:)
    type :: FMFD_parameters
        real(8) :: phi
        real(8) :: sig_t 
        real(8) :: sig_a 
        real(8) :: nusig_f 
        real(8) :: Jn(6) 
        real(8) :: J0(6)
        real(8) :: J1(6)
    end type
    logical :: fmfdon = .false.
    logical :: pfmfdon = .false.
    logical :: fake_MC = .false.
    logical :: fmfd2mc = .true.
    integer :: n_fake
    integer :: n_skip, n_acc
    integer :: FMFD_type        ! 1 FMFD / 2 p-FMFD / 3 1-node CMFD
    real(8) :: a_fm(6), v_fm
    
    type :: FMFD_accumulation
        type(FMFD_parameters), allocatable :: fm(:,:,:)
    endtype 
    
    type(FMFD_parameters), allocatable :: fm(:,:,:)
    type(FMFD_parameters), allocatable :: fm_avg(:,:,:)
    type(FMFD_parameters), allocatable :: fm_thread(:,:,:)
    !$OMP THREADPRIVATE(fm_thread)
    type(FMFD_accumulation), allocatable :: acc(:)

    ! FMFD grid
    real(8):: fm0(3)    ! (x0,y0,z0)
    real(8):: fm1(3)    ! (x1,y1,z1)
    real(8):: fm2(3)    ! (x1-x0,y1-y0,z1-z0)
    integer:: nfm(3)    ! (nx,ny,nz)
    real(8):: dfm(3)    ! (dx,dy,dz)

    ! Zigzag
    logical:: zigzagon
    integer:: n_zz      ! # of zigzag points
    integer:: zz_div    ! # of divisions
    ! - for fine mesh
    integer, allocatable:: zzf0(:) ! reference points
    integer, allocatable:: zzf1(:) ! lower points for zz0
    integer, allocatable:: zzf2(:) ! upper points for zz0
    ! - for coarse mesh
    integer, allocatable:: zzc0(:) ! reference points
    integer, allocatable:: zzc1(:) ! lower points for zz0
    integer, allocatable:: zzc2(:) ! upper points for zz0

    ! fission source distribution
    real(8), allocatable:: fsd_mc(:,:,:)
    real(8), allocatable:: fsd_fm(:,:,:)
    real(8), allocatable:: fsd(:,:,:)


    ! =========================================================================
    ! One-Node CMFD Acceleration
    logical :: cmfdon = .false.

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
        cm_phi0, &  ! coarse flux 0
        cm_phi1     ! coarse flux 1

    real(8), allocatable, dimension(:,:,:,:):: &
        deltf0, &   ! delta tilda
        deltf1, &   ! delta hat
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
        Mcm, &      ! matrix elements for CM
        Mfm         ! matrix elements for FM

    ! CMFD grid
    integer:: ncm(3)    ! (nx,ny,nz)
    integer:: fcr       ! fine > coarse in r-direction
    integer:: fcz       ! fine > coarse in z-direction

    ! Parallel calculation for one-node CMFD
    integer:: anode     ! # of active nodes
    integer:: bs0, bs1  ! buffer size
    integer, allocatable:: ax(:), ay(:), az(:)
    real(8), allocatable:: mvec(:,:,:,:,:)  ! matrix vector
    real(8), allocatable:: svec(:,:,:,:)    ! source vector

end module
