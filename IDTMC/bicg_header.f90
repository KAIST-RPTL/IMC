
include 'mkl_rci.f90'
include 'mkl_sparse_handle.f90'
include 'mkl_spblas.f90'

module BICG_HEADER
    use MKL_SPBLAS
    use MKL_SPARSE_HANDLE
    use MKL_RCI_TYPE
    use MKL_RCI
    implicit none

!    !---------------------------------------------------------------------------
!    ! # of entries & matrix size of input matrix and output matrix
!    ! 행렬 크기와 성분 개수
!    !---------------------------------------------------------------------------
!    INTEGER:: N, total_entry
!    INTEGER:: L_entry, U_entry
!
!    !---------------------------------------------------------------------------
!    ! Define variables for constructing sparse matrix of original matrix
!    ! 인풋을 통해 만들 희소행렬 자료 정보
!    !---------------------------------------------------------------------------
!    INTEGER,dimension(:),allocatable::IA,JA, A_istart,A_iend
!    real(8),dimension(:),allocatable::A

    !---------------------------------------------------------------------------
    ! ILU Decomposition configuration
    !---------------------------------------------------------------------------
    integer::maxfill
    real(8)::tol
    integer size0
    parameter (size0=128)
    INTEGER ipar(SIZE0)
    DOUBLE PRECISION dpar(SIZE0)

!    !---------------------------------------------------------------------------
!    ! Define variables for entire sparse matrix obtained from ILU(0) decomposition
!    ! 인풋으로 ILU(0) 분해 후 얻은 전체 행렬 정보
!    !---------------------------------------------------------------------------
!    real(8),dimension(:),allocatable::bilu0
!
!    !---------------------------------------------------------------------------
!    ! Define variables for entire sparse matrix obtained from ILU(T) decomposition
!    ! 인풋으로 ILU(T) 분해 후 얻은 전체 행렬 정보
!    !---------------------------------------------------------------------------
!    real(8),dimension(:),allocatable::bilut
!    integer,dimension(:),allocatable::jbilut,ibilut
!
!
!    !---------------------------------------------------------------------------
!    ! Define variables for constructing sparse triangular matrix with original matrix
!    ! 인풋을 이용해 ILU 분해를 한 뒤 얻은 삼각행렬 정보
!    !---------------------------------------------------------------------------
!    real(8),dimension(:),allocatable::L_value
!    integer,dimension(:),allocatable::L_j,L_istart,L_iend
!
!    real(8),dimension(:),allocatable::U_value
!    integer,dimension(:),allocatable::U_j,U_istart,U_iend
!
!    !---------------------------------------------------------------------------
!    ! Define variables corresponding to Intel MKL's standard
!    ! - Contains sparse matrix information with only one variable (pointer)
!    ! 인텔 MKL 라이브러리 계산에 사용되기 위해 인텔 MKL 기준에 부합하는 형태로 저장해줌
!    !---------------------------------------------------------------------------
!    type(sparse_matrix_t)::L_handle, U_handle, A_handle
!
!    !---------------------------------------------------------------------------
!    ! Provide structural hints of L, U matrices to sparse triangular solver
!    ! 삼각행렬 solver에게 행렬 형태를 전달함으로써 효율적 계산
!    !---------------------------------------------------------------------------
!    type(matrix_descr)::L_descr, U_descr, A_descr
!
!    !---------------------------------------------------------------------------
!    ! Solution vector, RHS vector
!    !---------------------------------------------------------------------------
!    real(8),dimension(:),allocatable::global_rhs,global_sol
!
!    !---------------------------------------------------------------------------
!    ! BiCGstab variable
!    !---------------------------------------------------------------------------
!    real(8),dimension(:),allocatable::bicg_sol



    ! matrix element & entries
    ! --- local
    real(8), allocatable:: am(:)        ! A matrix & entries
    integer, allocatable:: ia(:), ja(:)
    real(8), allocatable:: bm(:,:)      ! B matrix & entries (ILU)
    integer, allocatable:: ib(:,:), jb(:,:)
    ! --- global
    real(8), allocatable:: gam(:)        ! A matrix & entries
    integer, allocatable:: gia(:), gja(:)
    real(8), allocatable:: gbm(:)        ! B matrix & entries (ILU)
    integer, allocatable:: gib(:), gjb(:)
    ! --- fine
    real(8), allocatable:: fam(:)        ! A matrix & entries
    integer, allocatable:: fia(:), fja(:)
    real(8), allocatable:: fbm(:)        ! B matrix & entries (ILU)
    integer, allocatable:: fib(:), fjb(:)


    ! --- local
    integer, allocatable:: ln(:,:), un(:,:)     ! # of lower & upper entries
    integer, allocatable:: li(:,:,:), ui(:,:,:) ! lower & upper entries
    ! --- global
    integer, allocatable:: gln(:), gun(:)
    integer, allocatable:: gli(:,:), gui(:,:)
    ! --- fine
    integer, allocatable:: fln(:), fun(:)
    integer, allocatable:: fli(:,:), fui(:,:)


end module
