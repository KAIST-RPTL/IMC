module TH_HEADER
    implicit none

    logical:: th_on = .false.   ! T/H feedback on?
    logical:: th_cycle = .true.

    integer :: th_iter_max = 5
    integer :: th_iter     = 0
    real(8) :: th_cvg_crit = 1E-2 ! TH Convergence

    ! material properties
    real(8), allocatable:: k_fuel(:,:), & ! thermal conductivity of fuel
                           k_clad(:,:), & ! thermal conductivity of cladding
                           k_cool(:,:), & ! thermal conductivity of coolant
                           h_cool(:,:), & ! enthalpy of coolant
                           u_cool(:,:), & ! viscosity of coolant
                           c_cool(:,:)    ! specific heat of coolant

    ! average temperature
    real(8), allocatable:: t_fuel(:,:,:)    ! fuel temperature
    real(8), allocatable:: t_clad(:,:,:)    ! cladding temperature
    real(8), allocatable:: t_bulk(:,:,:)    ! bulk coolant temperature
    real(8), allocatable:: t_save(:,:,:)    ! backup for convergence test

    ! AVERAGE DENSITY (UPDATED FROM START)
    real(8), allocatable:: rho_bulk(:,:,:)

    ! TH grid
    real(8):: th0(3), th1(3)    ! starting and end point of T/H grid
    real(8):: th2(3)            ! total size of mesh
    real(8):: dth(3)            ! size of each node
    integer:: nth(3)            ! # of nodes
    integer:: mth(2)            ! # of actual pin nodes
    integer:: n_channels        ! Number of channels: (x+1)*(y+1)

    ! parameters
    real(8), allocatable:: mt1(:), mt2(:), mt3(:), st(:)
    real(8):: mflow ! [kg]
    real(8):: t_in  = 568.95D0      ! inlet temperature [K]
    real(8):: p_in  = 14.91336D+6   ! inlet pressure [Pa]
    real(8):: mflux
    real(8):: t_out = 601.75D0      ! outlet temperature [K]
    integer:: ith(2)                ! index
    real(8), allocatable:: rth(:)   ! radial distance
    real(8):: dr0, dr1              ! node size
    real(8):: inv_dr0, inv_dr1      ! inverse of node size
    real(8), allocatable:: hh(:,:,:), pp(:,:,:), pp_thread(:,:,:)
    real(8):: h_gap = 1D4           ! heat transfer coefficient [W/m2K]
    real(8):: p_th                  ! pitch size
    real(8):: rr0, rr1              ! radius to fuel and clad
    integer:: npp                   ! No. of nodes that the power produces
    
    ! START
    logical :: do_th_fmfd
    real(8), allocatable :: dthx(:), dthy(:), dthz(:)
                                    ! Difference with dth: not uniform
    real(8) :: pitch(3)
    !real(8), allocatable :: fuel_th(:,:), rad_th(:,:), cld_th(:,:), gap_th(:,:)
    real(8), allocatable :: rad_th(:,:)
    real(8) :: fuel_th, cld_th, gap_th
    real(8), allocatable :: power_th(:,:)
    real(8) :: margin(3) ! Margin between boundary mesh and actual mesh
    integer :: n_rod
    
    ! REFLECTIVE BOUNDARY if 1
    integer :: refl_th_w = 0
    integer :: refl_th_e = 0
    integer :: refl_th_n = 0
    integer :: refl_th_s = 0

    ! RECEIVING FROM START
    real(8), allocatable :: t_comm_cool(:,:), t_comm_fuel(:,:)
    real(8), allocatable :: rho_comm_cool(:,:), rho_comm_fuel(:,:)

    ! INITIAL VALUE
    real(8) :: rho_init = 1E3
    
end module
