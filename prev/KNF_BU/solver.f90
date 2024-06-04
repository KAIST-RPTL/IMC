module SOLVERS
    use mpi
    use FMFD_HEADER, only: nfm, ncm, fcr, fcz, zigzagon
    use VARIABLES, only: curr_cyc
    implicit none

    contains

! =============================================================================
!         Hepta diagonal matrix solvers
! =============================================================================
function BiCGSTAB_PRE(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: x, r, rs, v, p, s, t, &
                                               y, z, kk
    real(8), parameter :: e = 1D-10
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: i, j, k, iter

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    ! Jacobi preconditioner
    kk = 0
    where ( M(:,:,:,4) /= 0 ) kk(:,:,:) = 1D0/M(:,:,:,4)

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( ( norm_r .GT. norm_b) .and. (iter < 5D2) )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y        = kk*p

        v(:,:,:) = 0
        ! $omp parallel do default(shared) private(i,j,k)
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*y(i-1,j,k) ! x0
            if ( i /= nfm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,5)*y(i+1,j,k) ! x1
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,2)*y(i,j-1,k) ! y0
            if ( j /= nfm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*y(i,j+1,k) ! y1
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*y(i,j,k-1) ! z0
            if ( k /= nfm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,7)*y(i,j,k+1) ! z1
                               v(i,j,k) = v(i,j,k) + M(i,j,k,4)*y(i,j,k)
        end do
        end do
        end do
        ! $omp end parallel do
         
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z     = kk*s
        t(:,:,:) = 0
        ! $omp parallel do default(shared) private(i,j,k)
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*z(i-1,j,k)
            if ( i /= nfm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,5)*z(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,2)*z(i,j-1,k)
            if ( j /= nfm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*z(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*z(i,j,k-1)
            if ( k /= nfm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,7)*z(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,4)*z(i,j,k)
        end do
        end do
        end do
        ! $omp end parallel do

        omega  = sum(kk*t*s)/sum(kk*t*t)
        x      = x + alpha*y + omega*z
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   
    
end function BiCGSTAB_PRE


! =============================================================================
!         Hepta diagonal matrix solvers
! =============================================================================
function BiCGSTAB_ILU(M,Q) result(x)
    use FMFD_HEADER, only: fn, fx0, fx1, fy0, fy1, fz0, fz1, n_fnodes
    use PRECONDITIONER, only: ILU_FSOLVER, FORWARD_FSUB
    implicit none
    real(8), intent(in) :: M (:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8):: x(1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), dimension(n_fnodes):: xx, r, rs, v, p, s, t, kt, ks
    real(8), dimension(0:n_fnodes+1):: y, z
    !real(8), parameter :: e = 1D-7
    real(8), parameter :: e = 1D-12
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: it = 0
    integer :: mm, nn, oo, iter

    ! neutron source
    do oo = 1, nfm(3)
    do nn = 1, nfm(2)
    do mm = 1, nfm(1)
        if ( fn(mm,nn,oo) == 0 ) cycle
        r(fn(mm,nn,oo)) = Q(mm,nn,oo)
    end do
    end do
    end do

    x     = 0.0
    xx    = 0.0
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0
    y(0)  = 0.0
    z(0)  = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( norm_r .GT. norm_b .and. iter < 3D2 )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y(1:)    = ILU_FSOLVER(p)
        do mm = 1, n_fnodes
            v(mm) = y(fx0(mm))*M(mm,3) &
                  + y(fx1(mm))*M(mm,5) &
                  + y(fy0(mm))*M(mm,2) &
                  + y(fy1(mm))*M(mm,6) &
                  + y(fz0(mm))*M(mm,1) &
                  + y(fz1(mm))*M(mm,7) &
                  + y(mm)*M(mm,4)
        end do

        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z(1:) = ILU_FSOLVER(s)
        do mm = 1, n_fnodes
            t(mm) = z(fx0(mm))*M(mm,3) &
                  + z(fx1(mm))*M(mm,5) &
                  + z(fy0(mm))*M(mm,2) &
                  + z(fy1(mm))*M(mm,6) &
                  + z(fz0(mm))*M(mm,1) &
                  + z(fz1(mm))*M(mm,7) &
                  + z(mm)*M(mm,4)
        end do
        
        kt     = FORWARD_FSUB(t)
        ks     = FORWARD_FSUB(s)
        omega  = dot_product(kt,ks)/dot_product(kt,kt)
        xx     = xx + alpha*y(1:n_fnodes)+ omega*z(1:n_fnodes)
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   

    ! conversion
    do oo = 1, nfm(3)
    do nn = 1, nfm(2)
    do mm = 1, nfm(1)
        if ( fn(mm,nn,oo) == 0 ) cycle
        x(mm,nn,oo) = xx(fn(mm,nn,oo))
    end do
    end do
    end do

    
end function BiCGSTAB_ILU

! =============================================================================
! BICG_G
! =============================================================================
function BiCG_G(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(ncm(1),ncm(2),ncm(3)):: x, r, rs, v, p, s, t, pr
    real(8), dimension(0:ncm(1)+1,0:ncm(2)+1,0:ncm(3)+1):: y, z
    real(8), parameter :: e = 1D-11
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: mm, nn, oo, iter

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0
    pr    = 1D0/M(:,:,:,4)
    y     = 0.0
    z     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( norm_r .GT. norm_b .and. iter < 3D2 )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y(1:ncm(1),1:ncm(2),1:ncm(3)) = pr*p

        do mm = 1, ncm(1)
        do nn = 1, ncm(2)
        do oo = 1, ncm(3)
            v(mm,nn,oo) = y(mm-1,nn,oo)*M(mm,nn,oo,3) &
                        + y(mm+1,nn,oo)*M(mm,nn,oo,5) &
                        + y(mm,nn-1,oo)*M(mm,nn,oo,2) &
                        + y(mm,nn+1,oo)*M(mm,nn,oo,6) &
                        + y(mm,nn,oo-1)*M(mm,nn,oo,1) &
                        + y(mm,nn,oo+1)*M(mm,nn,oo,7) &
                        + y(mm,nn,oo)*M(mm,nn,oo,4)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z(1:ncm(1),1:ncm(2),1:ncm(3)) = pr*s
        do mm = 1, ncm(1)
        do nn = 1, ncm(2)
        do oo = 1, ncm(3)
            t(mm,nn,oo) = z(mm-1,nn,oo)*M(mm,nn,oo,3) &
                        + z(mm+1,nn,oo)*M(mm,nn,oo,5) &
                        + z(mm,nn-1,oo)*M(mm,nn,oo,2) &
                        + z(mm,nn+1,oo)*M(mm,nn,oo,6) &
                        + z(mm,nn,oo-1)*M(mm,nn,oo,1) &
                        + z(mm,nn,oo+1)*M(mm,nn,oo,7) &
                        + z(mm,nn,oo)*M(mm,nn,oo,4)
        end do
        end do
        end do
        
        omega  = sum(pr*t*s)/sum(pr*t*t)
        x      = x + alpha*y(1:ncm(1),1:ncm(2),1:ncm(3)) &
               + omega*z(1:ncm(1),1:ncm(2),1:ncm(3))
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   

    
end function BiCG_G


! =============================================================================
! BICG_L
! =============================================================================
function BICG_L(M,Q,C) result(x1)
    use FMFD_HEADER, only: anode, ax, ay, az, bs0, i_para0, i_para1
    use VARIABLES, only: ncore, icore, score
    use MPI, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_SUM
    use OMP_LIB
    implicit none
    real(8), intent(in) :: M(1:anode,1:fcr,1:fcr,1:fcz,1:7)
    real(8), intent(in) :: Q(1:anode,1:fcr,1:fcr,1:fcz)
    real(8), intent(in) :: C(1:anode,1:fcr,1:fcr,1:fcz)
    real(8) :: x0 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8) :: x1 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), dimension(1:fcr,1:fcr,1:fcz):: xx, r, rs, v, p, s, t, pr
    real(8), dimension(1:fcr,1:fcr,1:fcz,1:7):: MT
    real(8), dimension(0:fcr+1,0:fcr+1,0:fcz+1):: y, z
    real(8), parameter :: e = 1D-15
    !real(8), parameter :: e = 1D-13 use to be for a long time
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: ii, mm, nn, oo
    integer :: id(3), id1(3)
    integer :: iter, ista, iend
    real(8) :: mpie

    x0 = 0.0
    !$omp parallel default(private) shared(Q,M,C,i_para0,i_para1,x0,fcr,fcz,ax,ay,az)
    !$omp do
    do ii = i_para0, i_para1
    xx    = 0.0
    r(:,:,:)    = Q(ii,1:fcr,1:fcr,1:fcz)
    MT(:,:,:,:) = M(ii,1:fcr,1:fcr,1:fcz,1:7)
    pr(:,:,:)   = C(ii,1:fcr,1:fcr,1:fcz)
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0
    y     = 0.0
    z     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e

    iter = 1
    do while ( norm_r > norm_b .and. iter < 3D2 ) 
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y(1:fcr,1:fcr,1:fcz) = pr*p

        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            v(mm,nn,oo) = y(mm-1,nn,oo)*MT(mm,nn,oo,3) &
                        + y(mm+1,nn,oo)*MT(mm,nn,oo,5) &
                        + y(mm,nn-1,oo)*MT(mm,nn,oo,2) &
                        + y(mm,nn+1,oo)*MT(mm,nn,oo,6) &
                        + y(mm,nn,oo-1)*MT(mm,nn,oo,1) &
                        + y(mm,nn,oo+1)*MT(mm,nn,oo,7) &
                        + y(mm,nn,oo)*MT(mm,nn,oo,4)
        end do
        end do
        end do

        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z(1:fcr,1:fcr,1:fcz) = pr*s
        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            t(mm,nn,oo) = z(mm-1,nn,oo)*MT(mm,nn,oo,3) &
                        + z(mm+1,nn,oo)*MT(mm,nn,oo,5) &
                        + z(mm,nn-1,oo)*MT(mm,nn,oo,2) &
                        + z(mm,nn+1,oo)*MT(mm,nn,oo,6) &
                        + z(mm,nn,oo-1)*MT(mm,nn,oo,1) &
                        + z(mm,nn,oo+1)*MT(mm,nn,oo,7) &
                        + z(mm,nn,oo)*MT(mm,nn,oo,4)
        end do
        end do
        end do
        
        omega  = sum(pr*t*s)/sum(pr*t*t)
        if ( isnan(omega) ) omega = 0
        xx     = xx + alpha*y(1:fcr,1:fcr,1:fcz) + omega*z(1:fcr,1:fcr,1:fcz)
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1

    end do
    x0(ax(ii)+1:ax(ii)+fcr,ay(ii)+1:ay(ii)+fcr,az(ii)+1:az(ii)+fcz) = &
        xx(1:fcr,1:fcr,1:fcz)
    end do
    !$omp end do
    !$omp end parallel

    call MPI_REDUCE(x0,x1,nfm(1)*nfm(2)*nfm(3),MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,mpie)

end function BiCG_L


! =============================================================================
! BICG_G_ILU
! =============================================================================
function BiCG_G_ILU(M,Q) result(x)
    use FMFD_HEADER, only: gx0, gx1, gy0, gy1, gz0, gz1, anode, gcn
    use PRECONDITIONER, only: ILU_GSOLVER, FORWARD_GSUB
    use MPI
    implicit none
    real(8), intent(in) :: M (:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8):: x(1:ncm(1),1:ncm(2),1:ncm(3))
    real(8), dimension(anode):: xx, r, rs, v, p, s, t, kt, ks
    real(8), dimension(0:anode+1):: y, z
    !real(8), parameter :: e = 1D-10
    real(8), parameter :: e = 1D-15
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: mm, nn, oo, iter

    
    do oo = 1, ncm(3)
    do nn = 1, ncm(2)
    do mm = 1, ncm(1)
        if ( gcn(mm,nn,oo) == 0 ) cycle
        r(gcn(mm,nn,oo)) = Q(mm,nn,oo)
    end do
    end do
    end do

    x     = 0.0
    xx    = 0.0
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0
    y(0)  = 0.0
    z(0)  = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( norm_r .GT. norm_b .and. iter < 3D2 )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y(1:)    = ILU_GSOLVER(p)
        ! $omp parallel do default(shared) private(mm)
        do mm = 1, anode
            v(mm) = y(gx0(mm))*M(mm,3) &
                  + y(gx1(mm))*M(mm,5) &
                  + y(gy0(mm))*M(mm,2) &
                  + y(gy1(mm))*M(mm,6) &
                  + y(gz0(mm))*M(mm,1) &
                  + y(gz1(mm))*M(mm,7) &
                  + y(mm)*M(mm,4)
        end do
        ! $omp end parallel do

        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z(1:) = ILU_GSOLVER(s)
        ! $omp parallel do default(shared) private(mm)
        do mm = 1, anode
            t(mm) = z(gx0(mm))*M(mm,3) &
                  + z(gx1(mm))*M(mm,5) &
                  + z(gy0(mm))*M(mm,2) &
                  + z(gy1(mm))*M(mm,6) &
                  + z(gz0(mm))*M(mm,1) &
                  + z(gz1(mm))*M(mm,7) &
                  + z(mm)*M(mm,4)
        end do
        ! $omp end parallel do
        
        kt     = FORWARD_GSUB(t)
        ks     = FORWARD_GSUB(s)
        omega  = dot_product(kt,ks)/dot_product(kt,kt)
        xx     = xx + alpha*y(1:anode)+ omega*z(1:anode)
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   


    ! conversion
    do oo = 1, ncm(3)
    do nn = 1, ncm(2)
    do mm = 1, ncm(1)
        if ( gcn(mm,nn,oo) == 0 ) cycle
        x(mm,nn,oo) = xx(gcn(mm,nn,oo))
    end do
    end do
    end do

    
end function BiCG_G_ILU


! =============================================================================
! BICG_L_ILU
! =============================================================================
function BICG_L_ILU(M,Q) result(x1)
    use FMFD_HEADER, only: anode, ax, ay, az, bs0, n_lnodes, bs0, &
                    lx0, lx1, ly0, ly1, lz0, lz1, i_para0, i_para1, n_nodes
    !use VARIABLES, only: ncore, icore, score
    use VARIABLES   ! ***
    use MPI, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_SUM
    use OMP_LIB
    use PRECONDITIONER, only: ILU_SOLVER, FORWARD_SUB
    use BICG_HEADER, only: bm, jb, ib, ln, li, un, ui
    implicit none
    real(8), intent(in) :: M(:,:,:)
    real(8), intent(in) :: Q(:,:)
    real(8) :: x0 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8) :: x1 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8):: MT(n_lnodes,7)
    real(8), dimension(n_lnodes):: xx, r, rs, v, p, s, t, kt, ks
    real(8), dimension(0:n_lnodes):: y, z
    !real(8), parameter :: e = 1D-10
    !real(8), parameter :: e = 1D-15 ! log0608
    real(8), parameter :: e = 1D-15    ! log0508
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: ii, mm, nn, oo
    integer :: id(3), id1(3)
    integer :: iter, ista, iend
    real(8) :: mpie

    x0 = 0

    !$omp parallel default(private) shared(Q,M,i_para0,i_para1,x0,fcr,fcz) &
       shared(ax,ay,az,n_lnodes,lx0,lx1,ly0,ly1,lz0,lz1)
    !$omp do
    do ii = i_para0, i_para1
    xx    = 0.0
    r     = Q(ii,:)
    MT    = M(ii,:,:)
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0
    y(0)  = 0.0
    z(0)  = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e

    iter = 1
    do while ( norm_r > norm_b .and. iter < 3D2 ) 
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)
        y(1:)    = ILU_SOLVER(ii,p)
        do mm = 1, n_lnodes
            v(mm) = y(lx0(mm))*MT(mm,3) &
                  + y(lx1(mm))*MT(mm,5) &
                  + y(ly0(mm))*MT(mm,2) &
                  + y(ly1(mm))*MT(mm,6) &
                  + y(lz0(mm))*MT(mm,1) &
                  + y(lz1(mm))*MT(mm,7) &
                  + y(mm)*MT(mm,4)
        end do


        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        z(1:) = ILU_SOLVER(ii,s)
        do mm = 1, n_lnodes
            t(mm) = z(lx0(mm))*MT(mm,3) &
                  + z(lx1(mm))*MT(mm,5) &
                  + z(ly0(mm))*MT(mm,2) &
                  + z(ly1(mm))*MT(mm,6) &
                  + z(lz0(mm))*MT(mm,1) &
                  + z(lz1(mm))*MT(mm,7) &
                  + z(mm)*MT(mm,4)
        end do


        kt = FORWARD_SUB(ii,t)
        ks = FORWARD_SUB(ii,s)
        omega = dot_product(kt,ks)/dot_product(kt,kt)
        if ( isnan(omega) ) omega = 0
        xx     = xx + alpha*y(1:n_lnodes) + omega*z(1:n_lnodes)
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1

    end do

    x0(ax(ii)+1:ax(ii)+fcr,ay(ii)+1:ay(ii)+fcr,az(ii)+1:az(ii)+fcz) = &
        reshape(xx(:),[fcr,fcr,fcz])
    end do
    !$omp end do
    !$omp end parallel

    call MPI_REDUCE(x0(1:nfm(1),1:nfm(2),1:nfm(3)), &
            x1(1:nfm(1),1:nfm(2),1:nfm(3)),n_nodes,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,mpie)

end function BiCG_L_ILU


function OUT_OF_ZZL(io,jo)
    use FMFD_HEADER, only: zzc0, zz_div, zzc1, zzc2
    implicit none
    logical:: OUT_OF_ZZL
    integer, intent(in):: io, jo
    integer:: mo, no

    if ( .not. zigzagon ) then
        OUT_OF_ZZL = .false.
        return
    end if
    
    do mo = 1, zz_div
    if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
        no = mo
        exit
    end if
    end do

    if ( zzc1(no) < jo .and. jo <= zzc2(no) ) then
        OUT_OF_ZZL = .false.
    else
        OUT_OF_ZZL = .true.
    end if

end function

end module
