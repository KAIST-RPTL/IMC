module TEMPERATURE
    use TH_HEADER
    use VARIABLES, only: Nominal_Power, icore, score
    use CONSTANTS, only: K_B, PI
    implicit none
    integer:: ii, jj, kk
    real(8):: tvol

    contains


!! =============================================================================
!! 
!! =============================================================================
!subroutine TH_DISTANCE
!    
!    xyz(:) = p%coord(1)%xyz(:)
!    uvw(:) = p%coord(1)%uvw(:)
!    j_xyz  = TH_ID(xyz(:))
!
!
!end subroutine

! =============================================================================
! 
! =============================================================================
subroutine TH_INITIAL

    ! parameter
    mth(1:2) = th2(1:2) / p_th

    ! index
    ith(:) = 30
    if ( rr0 /= rr1 ) then
        ith(2) = 50
    end if

    ! delta radius
    dr0 = rr0/3D1
    dr1 = (rr1-rr0)/2D1
    inv_dr0 = 1d0/dr0
    inv_dr1 = 1d0/dr1

    ! allocation & initialization
    allocate(mt1(0:ith(2)),mt2(0:ith(2)),mt3(0:ith(2)))
    allocate(st(0:ith(2)),rth(ith(2)))
    allocate(hh(nth(1),nth(2),0:nth(3)),pp(nth(1),nth(2),nth(3)))
    allocate(pp_thread(nth(1),nth(2),nth(3)))
    allocate(t_fuel(nth(1),nth(2),nth(3)))
    allocate(t_bulk(nth(1),nth(2),nth(3)))
    allocate(t_save(nth(1),nth(2),nth(3)))
    if ( rr0 /= rr1 ) allocate(t_clad(nth(1),nth(2),nth(3)))
    pp_thread = 0
    pp  = 0
    mt1 = 0
    mt2 = 0
    mt3 = 0
    st  = 0
    t_fuel = t_out*K_B
    t_bulk = t_fuel
    if ( allocated(t_clad) ) t_clad = t_fuel

    ! radial distance
    rth(1) = dr0/2D0
    do ii = 2, ith(1)
        rth(ii) = rth(ii-1) + dr0
    end do
    if ( rr0 == rr1 ) return
    rth(ith(1)+1) = rth(ith(1)) + (dr0+dr1)/2d0 ! + 9D-3 (gap thickness)
    do ii = ith(1)+2, ith(2), 1
        rth(ii) = rth(ii-1) + dr1
    end do

end subroutine


! =============================================================================
! 
! =============================================================================
subroutine TH_INSIDE(xyz,j_xyz,inside_th)
    real(8), intent(in):: xyz(3)
    integer, intent(out):: j_xyz(3)
    logical:: inside_th
    integer:: ij

    inside_th = .true.
    j_xyz  = TH_ID(xyz(:))
    do ij = 1, 3
        if ( j_xyz(ij) < 1 .or. nth(ij) < j_xyz(ij) ) then
            inside_th = .false.
            return
        end if
    end do

end subroutine

! =============================================================================
! 
! =============================================================================
function TH_ID(xyz) result(id)
    real(8):: id(3)
    real(8), intent(in):: xyz(:)

    id(:) = floor((xyz(:)-th0(:))/dth(:))+1

end function

! =============================================================================
! 
! =============================================================================
subroutine TH_COL(wgt,xst,xsnf,id)
    real(8), intent(in):: wgt, xst, xsnf
    integer, intent(in):: id(1:3)

    pp_thread(id(1),id(2),id(3)) = pp_thread(id(1),id(2),id(3)) + wgt/xst*xsnf

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine NORM_TH()
    implicit none

    pp(:,:,:) = pp(:,:,:) + pp_thread(:,:,:)

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine PROCESS_TH()
    use VARIABLES, only: Nominal_Power
    use MPI, only: MPI_SUM
    implicit none
    real(8), dimension(nth(1),nth(2),nth(3)):: pp_mpi
    integer:: dsize
    integer:: ierr

    ! MPI data gathering
    dsize = nth(1)*nth(2)*nth(3)
    pp_mpi = 0
    call MPI_REDUCE(pp,pp_mpi,dsize,15,MPI_SUM,score,0,ierr)
    pp = pp_mpi

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine TEMP_SOLVE
    implicit none
    real(8), dimension(0:ith(2)):: tt0, tt1
    real(8), parameter:: error = 1E-9
    real(8):: heat2cool
    real(8):: t_b
    integer:: ee

    if ( score /= icore ) return

    call POWER_NORM
    call T_BULK_CALC

    tt1(0:) = (t_in+t_out)/2D0
    do ii = 1, nth(1)
    do jj = 1, nth(2)
    do kk = 1, nth(3)
    t_b = t_bulk(ii,jj,kk)
    call HEAT_TRANSFER(t_b,heat2cool)
    do
        tt0(0:) = tt1(0:)
        call T_MATRIX(tt1(0:ith(2)),t_b,heat2cool,pp(ii,jj,kk))
        call T_SOLVE(mt1(0:ith(2)),mt2(0:ith(2)),mt3(0:ith(2)), &
                     tt1(0:ith(2)),st(0:ith(2)))
        if ( norm2(tt1-tt0)/norm2(tt1) < error ) exit
    end do
    call AVG_T_FUEL(tt1(0:ith(1)),t_fuel(ii,jj,kk))
    if ( rr0 /= rr1 ) t_clad(ii,jj,kk) = 5D-1*(tt1(ith(1)+1)+tt1(ith(2)))*K_B
    end do
    end do
    end do
    t_bulk = t_bulk*K_B


end subroutine

! =============================================================================
! 
! =============================================================================
subroutine POWER_NORM
    use FMFD_HEADER, only: FSD_FM, fmfdon
    implicit none
    integer:: npp_(mth(1),mth(2))

    ! # of nodes that power produces
    npp_ = 1
    if ( fmfdon ) then
    where ( fsd_FM(:,:,1) == 0 ) npp_ = 0
    else
    where ( pp(:,:,1) == 0 ) npp_ = 0
    end if
    npp = sum(npp_)

    ! normalization
    pp = Nominal_Power*1D6*pp/sum(pp) ! [W]
    tvol = rr0*rr0*PI*dth(3)*npp/dble(nth(1)*nth(2)) ! [cm3]
    pp = pp/tvol ! [W/cm3]

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine T_BULK_CALC
    implicit none
    real(8):: h_in
    real(8):: dtemp
    ! h : enthalpy [kJ/kg]
    ! m : mass flow rate [kg/sec]
    ! T : temperature [K]
    ! P : power [J/sec]

    ! mass flow rate
    h_in = T2H(t_in)    ! [J/kg]
    mflow = Nominal_Power*1D6 / ((T2H(t_out)-T2H(t_in))*nth(1)*nth(2)) ! [kg/s]

    ! enthalpy distribution
    hh(:,:,0) = h_in
    do kk = 1, nth(3)
        hh(:,:,kk) = hh(:,:,kk-1) + pp(:,:,kk)*tvol/mflow ! [J/kg]
    end do

    ! bulk temperature distribution
    do kk = 1, nth(3)
    do jj = 1, nth(2)
    do ii = 1, nth(1)
        t_bulk(ii,jj,kk) = (H2T(hh(ii,jj,kk))+H2T(hh(ii,jj,kk-1)))/2D0
    end do
    end do
    end do

    ! non-fission region
    dtemp = (t_out-t_in)/dble(nth(3))
    do ii = 1, nth(1)
    do jj = 1, nth(2)
    do kk = 1, nth(3)
        if ( pp(ii,jj,kk) == 0 ) t_bulk(ii,jj,kk) = t_in + dtemp*(kk-5D-1)
    end do
    end do
    end do

end subroutine


! =============================================================================
! 
! =============================================================================
subroutine HEAT_TRANSFER(t_b,hc)
    implicit none
    real(8), intent(in) :: t_b
    real(8), intent(out):: hc
    real(8):: De    ! equivalent diameter
    real(8):: Re    ! Reynolds number
    real(8):: Pr    ! Prandtl number
    real(8):: area  ! total area
    real(8):: aa    ! parameter
    real(8):: mu    ! viscosity
    real(8):: kk    ! conductivity

    aa = p_th*p_th-PI*rr1*rr1
    area = aa*npp/dble(nth(1)*nth(2))
    De = 4D0*aa/(2D0*PI*rr1)
    mu = T2U(t_b)
    kk = T2K(3,t_b)
    Re = De*mflow/(area*mu)*1D2 ! [cm][kg/sec]/[cm2]/[kg/m-sec]
    Pr = T2C(t_b)*mu/kk
    hc = 2.3D-2*kk/De*(Re**0.8D0)*(Pr**0.4D0)*1D2 ! [W/m2-K]

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine T_MATRIX(tt1,t_b,h_c,qq)
    real(8), intent(in):: tt1(0:ith(2))
    real(8), intent(in):: h_c, t_b, qq
    integer:: mm

    ! -------------------------------------------------------------------------
    ! only fuel
    if ( rr0 == rr1 ) then
        ! fuel element
        do mm = 1, ith(1)-1
        mt1(mm) = -rth(mm)  *T2K(1,tt1(mm-1))*inv_dr0   ! [W/m-K]
        mt3(mm) = -rth(mm+1)*T2K(1,tt1(mm))  *inv_dr0
        mt2(mm) = -mt1(mm)-mt3(mm)
        end do
        mt2(0)      = rth(1)*T2K(1,tt1(0))*inv_dr0
        mt3(0)      = -mt2(0)
        mt1(ith(1)) = -rth(ith(1))*T2K(1,tt1(ith(1)))*inv_dr0
        mt2(ith(1)) = -mt1(ith(1))+h_c*(rth(ith(1))+1.5D0*dr0)*1D-2
    
        ! source element
        do mm = 1, ith(1)-1
        st(mm) = 5D-1*(rth(mm+1)*rth(mm+1)-rth(mm)*rth(mm))*qq*1D2
        end do
        st(0) = 5D-1*rth(1)*rth(1)*qq*1D2
        st(ith(1)) = +h_c*(rth(ith(1))+1.5D0*dr0)*t_b*1D-2

    ! -------------------------------------------------------------------------
    ! fuel - gap - clad
    else
        ! fuel element
        do mm = 1, ith(1)-1
        mt1(mm) = -rth(mm)  *T2K(1,tt1(mm-1))*inv_dr0   ! [W/m-K]
        mt3(mm) = -rth(mm+1)*T2K(1,tt1(mm))  *inv_dr0   ! [W/m-K]
        mt2(mm) = -mt1(mm)-mt3(mm)
        end do
        mt2(0)      = rth(1)*T2K(1,tt1(0))*inv_dr0 ! [W/m-K]
        mt3(0)      = -mt2(0)
        mt1(ith(1)) = -rth(ith(1))*T2K(1,tt1(ith(1)))*inv_dr0
        mt3(ith(1)) = -h_gap*(rth(ith(1))+1.5D0*dr0)*1D-2
        mt2(ith(1)) = -mt1(ith(1))-mt3(ith(1))
    
        ! clad element
        do mm = ith(1)+2, ith(2)-1
        mt1(mm) = -rth(mm)  *T2K(2,tt1(mm-1))*inv_dr1   ! [W/m-K]
        mt3(mm) = -rth(mm+1)*T2K(2,tt1(mm))  *inv_dr1   ! [W/m-K]
        mt2(mm) = -mt1(mm)-mt3(mm)
        end do
        mt1(ith(1)+1) = -(rth(ith(1))+5D-1*dr0)*h_c*1D-2
        mt3(ith(1)+1) = -rth(ith(1)+1)*T2K(2,tt1(ith(1)+1))*inv_dr1
        mt2(ith(1)+1) = -mt1(ith(1)+1)-mt3(ith(1)+1)
        mt1(ith(2)) = -rth(ith(2))*T2K(2,tt1(ith(2)))*inv_dr1
        mt2(ith(2)) = -mt1(ith(2))+(rth(ith(2))+5D-1*dr1)*h_c*1D-2

        ! source element
        st(ith(1):) = 0D0
        do mm = 1, ith(1)-1
        st(mm) = 5D-1*(rth(mm+1)*rth(mm+1)-rth(mm)*rth(mm))*qq*1D2
        end do
        st(0) = 5D-1*rth(1)*rth(1)*qq*1D2
        st(ith(2)) = (rth(ith(2))+5D-1*dr1)*h_c*t_b*1D-2

    end if

end subroutine

subroutine T_SOLVE(mm1,mm2,mm3,xx,cc)
    real(8), intent(in):: mm1(0:), mm2(0:), mm3(0:), cc(0:)
    real(8), intent(inout):: xx(0:)
    real(8), dimension(0:ith(2)):: m1, m2, m3, c1
    real(8):: temp
    integer:: nn

    m1(0:) = mm1(0:)
    m2(0:) = mm2(0:)
    m3(0:) = mm3(0:)
    c1(0:) = cc(0:)

    do nn = 1, ith(2)
        temp = m1(nn)/m2(nn-1)
        m2(nn) = m2(nn)-temp*m3(nn-1)
        c1(nn) = c1(nn)-temp*c1(nn-1)
    end do

    xx(ith(2)) = c1(ith(2))/m2(ith(2))
    do nn = ith(2)-1, 0, -1
    xx(nn) = (c1(nn)-m3(nn)*xx(nn+1))/m2(nn)
    end do

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine AVG_T_FUEL(tt1,t_f)
    implicit none
    real(8), intent(in):: tt1(0:ith(1))
    real(8), intent(inout):: t_f

    ! heuristic value
    t_f = (3D-1*tt1(0)+7D-1*tt1(ith(1)))*K_B

!    ! average fuel temperature
!    t_f = 0
!    do ii = 1, ith(1)
!        t_f = t_f + (rth(ii)*rth(ii)-rth(ii-1)*rth(ii-1))*(tt1(ii)+tt1(ii-1))
!    end do
!    t_f = t_f / (rth(ith(1))*rth(ith(1))) * 5D-1

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine TEMP_DISTRIBUTE
    use MPI, only: MPI_REAL8, MPI_COMM_WORLD
    implicit none
    integer:: TP, WOR, IERR

    TP = MPI_REAL8
    WOR = MPI_COMM_WORLD

    call MPI_BCAST(t_fuel,1,TP,score,WOR,IERR)
    call MPI_BCAST(t_bulk,1,TP,score,WOR,IERR)
    if ( allocated(t_clad) ) &
    call MPI_BCAST(t_clad,1,TP,score,WOR,IERR)

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine TEMP_CONVERGE
    real(8):: err
    
    if ( icore /= score ) return

    err = norm2(t_fuel-t_save)/norm2(t_fuel)
    t_save = t_fuel

end subroutine


! =============================================================================
! 
! =============================================================================
function T2H(tmp)
    use TH_HEADER, only: h_cool
    implicit none
    real(8):: T2H
    real(8), intent(in):: tmp
    integer:: pt1, pt2, pt3
    real(8):: slope

    ! index search
    pt1 = 1
    pt2 = 11
    do
        if ( pt2 - pt1 == 1 ) exit
        pt3 = (pt2+pt1)/2
        if ( tmp >= h_cool(1,pt3) ) then
            pt1 = pt3
        else
            pt2 = pt3
        end if
    end do
    pt3 = pt1
    slope = (tmp-h_cool(1,pt3))/(h_cool(1,pt3+1)-h_cool(1,pt3))
    T2H = h_cool(2,pt3) + slope * (h_cool(2,pt3+1)-h_cool(2,pt3))
    T2H = T2H

end function

! =============================================================================
! 
! =============================================================================
function T2K(mat,tmp)
    use TH_HEADER, only: k_fuel, k_cool, k_clad
    implicit none
    real(8):: T2K
    integer, intent(in):: mat
    real(8), intent(in):: tmp
    integer:: pt1, pt2, pt3
    real(8):: slope

    ! index search
    pt1 = 1
    select case(mat)
    case(1)
        pt2 = size(k_fuel(1,:))
        do
            if ( pt2 - pt1 == 1 ) exit
            pt3 = (pt2+pt1)/2
            if ( tmp >= k_fuel(1,pt3) ) then
                pt1 = pt3
            else
                pt2 = pt3
            end if
        end do
        pt3 = pt1
        slope = (tmp-k_fuel(1,pt3))/(k_fuel(1,pt3+1)-k_fuel(1,pt3))
        T2K = k_fuel(2,pt3) + slope * (k_fuel(2,pt3+1)-k_fuel(2,pt3))

    case(2)
        pt2 = size(k_clad(1,:))
        do
            if ( pt2 - pt1 == 1 ) exit
            pt3 = (pt2+pt1)/2
            if ( tmp >= k_clad(1,pt3) ) then
                pt1 = pt3
            else
                pt2 = pt3
            end if
        end do
        pt3 = pt1
        slope = (tmp-k_clad(1,pt3))/(k_clad(1,pt3+1)-k_clad(1,pt3))
        T2K = k_clad(2,pt3) + slope * (k_clad(2,pt3+1)-k_clad(2,pt3))

    case(3)
        pt2 = size(k_cool(1,:))
        do
            if ( pt2 - pt1 == 1 ) exit
            pt3 = (pt2+pt1)/2
            if ( tmp >= k_cool(1,pt3) ) then
                pt1 = pt3
            else
                pt2 = pt3
            end if
        end do
        pt3 = pt1
        slope = (tmp-k_cool(1,pt3))/(k_cool(1,pt3+1)-k_cool(1,pt3))
        T2K = k_cool(2,pt3) + slope * (k_cool(2,pt3+1)-k_cool(2,pt3))
    end select

end function

! =============================================================================
! 
! =============================================================================
function H2T(etp)
    real(8):: H2T
    real(8), intent(in):: etp
    integer:: pt1, pt2, pt3
    real(8):: slope

    ! index search
    pt1 = 1
    pt2 = 11
    do
        if ( pt2 - pt1 == 1 ) exit
        pt3 = (pt2+pt1)/2
        if ( etp >= h_cool(2,pt3) ) then
            pt1 = pt3
        else
            pt2 = pt3
        end if
    end do
    pt3 = pt1
    slope = (etp-h_cool(2,pt3))/(h_cool(2,pt3+1)-h_cool(2,pt3))
    H2T = h_cool(1,pt3) + slope * (h_cool(1,pt3+1)-h_cool(1,pt3))
    H2T = H2T

end function

! =============================================================================
! 
! =============================================================================
function T2U(t_b)
    real(8):: T2U
    real(8), intent(in):: t_b
    integer:: pt1, pt2, pt3
    real(8):: slope

    pt1 = 1
    pt2 = 39
    do
        if ( pt2 - pt1 == 1 ) exit
        pt3 = (pt2+pt1)/2
        if ( t_b >= u_cool(1,pt3) ) then
            pt1 = pt3
        else
            pt2 = pt3
        end if
    end do
    pt3 = pt1
    slope = (t_b-u_cool(1,pt3))/(u_cool(1,pt3+1)-u_cool(1,pt3))
    T2U = u_cool(2,pt3) + slope * (u_cool(2,pt3+1)-u_cool(2,pt3))

end function

! =============================================================================
! 
! =============================================================================
function T2C(t_b)
    implicit none
    real(8):: T2C
    real(8), intent(in):: t_b
    integer:: pt1, pt2, pt3
    real(8):: slope

    pt1 = 1
    pt2 = 39
    do
        if ( pt2 - pt1 == 1 ) exit
        pt3 = (pt2+pt1)/2
        if ( t_b >= c_cool(1,pt3) ) then
            pt1 = pt3
        else
            pt2 = pt3
        end if
    end do
    pt3 = pt1
    slope = (t_b-c_cool(1,pt3))/(c_cool(1,pt3+1)-c_cool(1,pt3))
    T2C = c_cool(2,pt3) + slope * (c_cool(2,pt3+1)-c_cool(2,pt3))

end function

subroutine TEMP_UPDATE_BU(istep_burnup)
implicit none
integer, intent(in) :: istep_burnup
if(icore==score) print *, 'CHK;', allocated(t_fuel_bu), temp_grid_on
if( .not. allocated( t_fuel_bu ) ) return
if( .not. temp_grid_on ) return

t_fuel = t_fuel_bu(:,:,:, istep_burnup)
t_clad = t_clad_bu(:,:,:, istep_burnup)
t_bulk = t_bulk_bu(:,:,:, istep_burnup)
rho_bulk = rho_bulk_bu(:,:,:, istep_burnup)
if( icore == score ) print *, 'TFUEL', t_fuel/K_B
if( icore == score ) print *, 'TCOOL', t_bulk/K_B
if( icore == score ) print *, 'RHOCOOL', cool_dens, rho_bulk
if (icore==score) print *, 'Updated temperature grids'
end subroutine

end module
