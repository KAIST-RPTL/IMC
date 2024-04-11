module TH_HEADER
    implicit none

    logical:: th_on = .false.   ! T/H feedback on?
    logical:: temp_grid_on = .false.
    logical:: th_cycle = .true.

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

    ! Coolant density [ fraction ]
    real(8), allocatable :: rho_bulk(:,:,:) ! bulk density fraction

    ! BU-dependent
    real(8), allocatable:: t_fuel_bu(:,:,:,:)    ! fuel temperature
    real(8), allocatable:: t_clad_bu(:,:,:,:)    ! cladding temperature
    real(8), allocatable:: t_bulk_bu(:,:,:,:)    ! bulk coolant temperature
    real(8), allocatable :: rho_bulk_bu(:,:,:,:) ! bulk density fraction

    logical :: do_temp_grid = .false.

    ! TH grid
    real(8):: th0(3), th1(3)    ! starting and end point of T/H grid
    real(8):: th2(3)            ! total size of mesh
    real(8):: dth(3)            ! size of each node
    integer:: nth(3)            ! # of nodes
    integer:: mth(2)            ! # of actual pin nodes

    ! parameters
    real(8), allocatable:: mt1(:), mt2(:), mt3(:), st(:)
    real(8):: mflow ! [kg]
    real(8):: t_in  = 568.95D0      ! inlet temperature [K]
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

    real(8):: cool_dens = 1d0
    integer:: symm = 1

end module
