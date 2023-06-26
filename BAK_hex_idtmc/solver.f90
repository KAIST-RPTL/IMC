module SOLVERS
    use FMFD_HEADER, only: nfm, ncm, fcr, fcz, zigzagon
    implicit none

    contains

! =============================================================================
!         Hepta diagonal matrix solvers
! =============================================================================
function BiCGStab_hepta(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: x, r, rs, v, p, s, t
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: i, j, k, iter

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( ( norm_r .GT. norm_b) .and. (iter < 3D2) )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*p(i-1,j,k) ! x0
            if ( i /= nfm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,5)*p(i+1,j,k) ! x1
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,2)*p(i,j-1,k) ! y0
            if ( j /= nfm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*p(i,j+1,k) ! y1
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*p(i,j,k-1) ! z0
            if ( k /= nfm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,7)*p(i,j,k+1) ! z1
                               v(i,j,k) = v(i,j,k) + M(i,j,k,4)*p(i,j,k)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do i = 1, nfm(1)
        do j = 1, nfm(2)
        do k = 1, nfm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*s(i-1,j,k)
            if ( i /= nfm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,5)*s(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,2)*s(i,j-1,k)
            if ( j /= nfm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*s(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*s(i,j,k-1)
            if ( k /= nfm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,7)*s(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,4)*s(i,j,k)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        x      = x + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   
    
end function BiCGStab_hepta     


! =============================================================================
! BICG_G
! =============================================================================
function BiCG_G(M,Q) result(x)
    real(8), intent(in) :: M (:,:,:,:)
    real(8), intent(in) :: Q (:,:,:)
    real(8), dimension(ncm(1),ncm(2),ncm(3)):: x, r, rs, v, p, s, t
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    real(8) :: summesion, temp
    integer :: it = 0
    integer :: i, j, k, iter

    x     = 0.0
    r     = Q
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e
    
    iter = 1
    do while ( norm_r .GT. norm_b .and. iter < 3D2 )
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do i = 1, ncm(1)
        do j = 1, ncm(2)
        do k = 1, ncm(3)
            if ( i /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,3)*p(i-1,j,k) ! x0
            if ( i /= ncm(1) ) v(i,j,k) = v(i,j,k) + M(i,j,k,5)*p(i+1,j,k) ! x1
            if ( j /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,2)*p(i,j-1,k) ! y0
            if ( j /= ncm(2) ) v(i,j,k) = v(i,j,k) + M(i,j,k,6)*p(i,j+1,k) ! y1
            if ( k /= 1 )      v(i,j,k) = v(i,j,k) + M(i,j,k,1)*p(i,j,k-1) ! z0
            if ( k /= ncm(3) ) v(i,j,k) = v(i,j,k) + M(i,j,k,7)*p(i,j,k+1) ! z1
                               v(i,j,k) = v(i,j,k) + M(i,j,k,4)*p(i,j,k)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do i = 1, ncm(1)
        do j = 1, ncm(2)
        do k = 1, ncm(3)
            if ( i /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,3)*s(i-1,j,k)
            if ( i /= ncm(1) ) t(i,j,k) = t(i,j,k) + M(i,j,k,5)*s(i+1,j,k)
            if ( j /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,2)*s(i,j-1,k)
            if ( j /= ncm(2) ) t(i,j,k) = t(i,j,k) + M(i,j,k,6)*s(i,j+1,k)
            if ( k /= 1 )      t(i,j,k) = t(i,j,k) + M(i,j,k,1)*s(i,j,k-1)
            if ( k /= ncm(3) ) t(i,j,k) = t(i,j,k) + M(i,j,k,7)*s(i,j,k+1)
                               t(i,j,k) = t(i,j,k) + M(i,j,k,4)*s(i,j,k)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        x      = x + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1
    
    end do   

    
end function BiCG_G

! =============================================================================
! BICG_L
! =============================================================================
function BICG_L(M,Q) result(x)
    real(8), intent(in) :: M (1:nfm(1),1:nfm(2),1:nfm(3),1:7)
    real(8), intent(in) :: Q (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8) :: x (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), dimension(1:fcr,1:fcr,1:fcz):: xx, r, rs, v, p, s, t
    real(8), dimension(1:fcr,1:fcr,1:fcz,1:7):: MT
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: ii, jj, kk, mm, nn, oo
    integer :: id0(3), id1(3)
    integer :: iter

    x = 0
    do jj = 1, ncm(2); id0(2) = (jj-1)*fcr; id1(2) = id0(2) + fcr
    do ii = 1, ncm(1); id0(1) = (ii-1)*fcr; id1(1) = id0(1) + fcr
    do kk = 1, ncm(3); id0(3) = (kk-1)*fcz; id1(3) = id0(3) + fcz

    if ( OUT_OF_ZZL(ii,jj) ) cycle

    xx    = 0.0
    r(:,:,:) = Q(id0(1)+1:id1(1),id0(2)+1:id1(2),id0(3)+1:id1(3))
    MT(:,:,:,:) = M(id0(1)+1:id1(1),id0(2)+1:id1(2),id0(3)+1:id1(3),1:7)
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e

    iter = 1
    do while ( norm_r > norm_b .and. iter < 1D2 ) 
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            if ( mm /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm-1,nn,oo) & ! x0
                                * MT(mm,nn,oo,3)
            if ( mm /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm+1,nn,oo) & ! x1
                                * MT(mm,nn,oo,5)
            if ( nn /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn-1,oo) & ! y0
                                * MT(mm,nn,oo,2)
            if ( nn /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn+1,oo) & ! y1
                                * MT(mm,nn,oo,6)
            if ( oo /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo-1) & ! z0
                                * MT(mm,nn,oo,1)
            if ( oo /= fcz ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo+1) & ! z1
                                * MT(mm,nn,oo,7)
                             v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo) &
                                * MT(mm,nn,oo,4)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            if ( mm /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm-1,nn,oo) & ! x0
                                * MT(mm,nn,oo,3)
            if ( mm /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm+1,nn,oo) & ! x1
                                * MT(mm,nn,oo,5)
            if ( nn /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn-1,oo) & ! y0
                                * MT(mm,nn,oo,2)
            if ( nn /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn+1,oo) & ! y1
                                * MT(mm,nn,oo,6)
            if ( oo /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo-1) & ! z0
                                * MT(mm,nn,oo,1)
            if ( oo /= fcz ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo+1) & ! z1
                                * MT(mm,nn,oo,7)
                             t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo) &
                                * MT(mm,nn,oo,4)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        xx     = xx + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1

    end do   
    x(id0(1)+1:id1(1),id0(2)+1:id1(2),id0(3)+1:id1(3)) = xx(1:fcr,1:fcr,1:fcz)
    end do
    end do
    end do

end function BiCG_L

! =============================================================================
! BICG_L
! =============================================================================
function BICG_LP(M,Q) result(x1)
    use FMFD_HEADER, only: anode, ax, ay, az, bs0
    use VARIABLES, only: ncore, icore, score
    use MPI, only: MPI_COMM_WORLD, MPI_DOUBLE_PRECISION, MPI_SUM
    implicit none
    real(8), intent(in) :: M (1:anode,1:fcr,1:fcr,1:fcz,1:7)
    real(8), intent(in) :: Q (1:anode,1:fcr,1:fcr,1:fcz)
    real(8) :: x0 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8) :: x1 (1:nfm(1),1:nfm(2),1:nfm(3))
    real(8), dimension(1:fcr,1:fcr,1:fcz):: xx, r, rs, v, p, s, t
    real(8), dimension(1:fcr,1:fcr,1:fcz,1:7):: MT
    real(8), parameter :: e = 1D-8
    real(8) :: rho      , rho_prev
    real(8) :: alpha    , omega   , beta
    real(8) :: norm_r   , norm_b
    integer :: ii, jj, kk, mm, nn, oo
    integer :: id0(3), id1(3)
    integer :: iter, ista, iend
    real(8) :: mpie

    x0 = 0
    call MPI_RANGE(anode,ncore,icore,ista,iend)
    do ii = ista, iend
    xx    = 0.0
    r(:,:,:)    = Q(ii,1:fcr,1:fcr,1:fcz)
    MT(:,:,:,:) = M(ii,1:fcr,1:fcr,1:fcz,1:7)
    rs    = r
    rho   = 1.0
    alpha = 1.0
    omega = 1.0
    v     = 0.0
    p     = 0.0

    norm_r = sum(r*r)
    norm_b = norm_r*e

    iter = 1
    do while ( norm_r > norm_b .and. iter < 1D2 ) 
        rho_prev = rho
        rho      = sum(rs*r)
        beta     = (rho/rho_prev) * (alpha/omega)
        p        = r + beta * (p - omega*v)

        v(:,:,:) = 0
        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            if ( mm /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm-1,nn,oo) & ! x0
                                * MT(mm,nn,oo,3)
            if ( mm /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm+1,nn,oo) & ! x1
                                * MT(mm,nn,oo,5)
            if ( nn /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn-1,oo) & ! y0
                                * MT(mm,nn,oo,2)
            if ( nn /= fcr ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn+1,oo) & ! y1
                                * MT(mm,nn,oo,6)
            if ( oo /= 1 )   v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo-1) & ! z0
                                * MT(mm,nn,oo,1)
            if ( oo /= fcz ) v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo+1) & ! z1
                                * MT(mm,nn,oo,7)
                             v(mm,nn,oo) = v(mm,nn,oo) + p(mm,nn,oo) &
                                * MT(mm,nn,oo,4)
        end do
        end do
        end do
        
        alpha = rho/sum(rs*v)
        s     = r - alpha*v
        t(:,:,:) = 0
        do mm = 1, fcr
        do nn = 1, fcr
        do oo = 1, fcz
            if ( mm /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm-1,nn,oo) & ! x0
                                * MT(mm,nn,oo,3)
            if ( mm /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm+1,nn,oo) & ! x1
                                * MT(mm,nn,oo,5)
            if ( nn /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn-1,oo) & ! y0
                                * MT(mm,nn,oo,2)
            if ( nn /= fcr ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn+1,oo) & ! y1
                                * MT(mm,nn,oo,6)
            if ( oo /= 1 )   t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo-1) & ! z0
                                * MT(mm,nn,oo,1)
            if ( oo /= fcz ) t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo+1) & ! z1
                                * MT(mm,nn,oo,7)
                             t(mm,nn,oo) = t(mm,nn,oo) + s(mm,nn,oo) &
                                * MT(mm,nn,oo,4)
        end do
        end do
        end do
        
        omega  = sum(t*s)/sum(t*t)
        xx     = xx + alpha*p + omega*s
        r      = s - omega*t
        norm_r = sum(r*r)
        iter   = iter + 1

    end do
    x0(ax(ii)+1:ax(ii)+fcr,ay(ii)+1:ay(ii)+fcr,az(ii)+1:az(ii)+fcz) = &
        xx(1:fcr,1:fcr,1:fcz)
    end do

    call MPI_REDUCE(x0,x1,nfm(1)*nfm(2)*nfm(3),15,MPI_SUM,score,0,mpie)

end function BiCG_LP

subroutine MPI_RANGE(nn,ncore,icore,ista,iend)
    integer:: iwork1, iwork2
    integer, intent(in):: nn, ncore, icore
    integer, intent(inout):: ista, iend

    iwork1 = nn / ncore
    iwork2 = mod(nn,ncore)
    ista = icore * iwork1 + 1 + min(icore,iwork2)
    iend = ista + iwork1 - 1
    if ( iwork2 > icore ) iend = iend + 1

end subroutine

!! =============================================================================
!! SOR
!! =============================================================================
!subroutine SOR(m0,ss,ff)
!    implicit none
!    real(8), intent(in   ):: m0(:,:,:,:)
!    real(8), intent(in   ):: ss(:,:,:)
!    real(8), intent(inout):: ff(:,:,:)
!    integer:: ii, jj, kk, ee, mm
!    real(8):: temp
!    real(8):: relax = 1.4D0
!    integer:: n_inner = 5
!
!    do mm=1, n_inner
!    do kk=1, nfm(3)
!    do jj=1, nfm(2)
!    do ii=1, nfm(1)
!       temp = ss(ii,jj,kk)
!       if ( ii /= 1 )      temp = temp - m0(ii,jj,kk,3)*ff(ii-1,jj,kk)
!       if ( ii /= nfm(1) ) temp = temp - m0(ii,jj,kk,5)*ff(ii+1,jj,kk)
!       if ( jj /= 1 )      temp = temp - m0(ii,jj,kk,2)*ff(ii,jj-1,kk)
!       if ( jj /= nfm(2) ) temp = temp - m0(ii,jj,kk,6)*ff(ii,jj+1,kk)
!       if ( kk /= 1 )      temp = temp - m0(ii,jj,kk,1)*ff(ii,jj,kk-1)
!       if ( kk /= nfm(3) ) temp = temp - m0(ii,jj,kk,7)*ff(ii,jj,kk+1)
!       ff(ii,jj,kk) = (1D0-relax)*ff(ii,jj,kk)+relax*temp/m0(ii,jj,kk,4)
!    end do
!    end do
!    end do
!    end do
!
!end subroutine
!
!
!! =============================================================================
!! CG is a matrix solver by the SOR
!! =============================================================================
!subroutine SORG(m0,ss,ff)
!    implicit none
!    real(8), intent(in   ):: m0(1:ncm(1),1:ncm(2),1:ncm(3),1:7)
!    real(8), intent(in   ):: ss(1:ncm(1),1:ncm(2),1:ncm(3))
!    real(8), intent(inout):: ff(1:ncm(1),1:ncm(2),1:ncm(3))
!    integer:: ii, jj, kk, ee, mm
!    real(8):: temp
!    integer:: id(3)
!    real(8):: relax = 1.4D0
!    integer:: n_inner = 5
!
!    do mm=1, n_inner
!    do kk=1, ncm(3)
!    do jj=1, ncm(2)
!    do ii=1, ncm(1)
!       if ( OUT_OF_ZZL(ii,jj) ) cycle
!       temp = ss(ii,jj,kk)
!       if ( ii /= 1 )      temp = temp - m0(ii,jj,kk,3)*ff(ii-1,jj,kk)
!       if ( ii /= ncm(1) ) temp = temp - m0(ii,jj,kk,5)*ff(ii+1,jj,kk)
!       if ( jj /= 1 )      temp = temp - m0(ii,jj,kk,2)*ff(ii,jj-1,kk)
!       if ( jj /= ncm(2) ) temp = temp - m0(ii,jj,kk,6)*ff(ii,jj+1,kk)
!       if ( kk /= 1 )      temp = temp - m0(ii,jj,kk,1)*ff(ii,jj,kk-1)
!       if ( kk /= ncm(3) ) temp = temp - m0(ii,jj,kk,7)*ff(ii,jj,kk+1)
!       ff(ii,jj,kk) = (1D0-relax)*ff(ii,jj,kk)+relax*temp/m0(ii,jj,kk,4)
!    end do
!    end do
!    end do
!    end do
!
!end subroutine
!
!subroutine SORL(m0,ss,ff)
!    use FMFD_HEADER, only: zigzagon
!    implicit none
!    real(8), intent(in   ):: m0(1:nfm(1),1:nfm(2),1:nfm(3),1:7)
!    real(8), intent(in   ):: ss(1:nfm(1),1:nfm(2),1:nfm(3))
!    real(8), intent(inout):: ff(1:nfm(1),1:nfm(2),1:nfm(3))
!    integer:: ii, jj, kk, mm, nn, oo, ll
!    real(8):: temp
!    integer:: id0(3), id(3)
!    real(8):: relax = 1.4D0
!    integer:: n_inner = 5
!
!    do kk=1, ncm(3); id0(3) = (kk-1)*fcz
!    do jj=1, ncm(2); id0(2) = (jj-1)*fcr
!    do ii=1, ncm(1); id0(1) = (ii-1)*fcr
!    if ( OUT_OF_ZZL(ii,jj) ) cycle
!    do ll=1, n_inner
!      do mm = 1, fcr; id(1) = id0(1)+mm
!      do nn = 1, fcr; id(2) = id0(2)+nn
!      do oo = 1, fcz; id(3) = id0(3)+oo
!        temp = ss(id(1),id(2),id(3))
!        if ( mm /= 1 )   temp = temp - m0(id(1),id(2),id(3),3)*ff(id(1)-1,id(2),id(3))
!        if ( mm /= fcr ) temp = temp - m0(id(1),id(2),id(3),5)*ff(id(1)+1,id(2),id(3))
!        if ( nn /= 1 )   temp = temp - m0(id(1),id(2),id(3),2)*ff(id(1),id(2)-1,id(3))
!        if ( nn /= fcr ) temp = temp - m0(id(1),id(2),id(3),6)*ff(id(1),id(2)+1,id(3))
!        if ( oo /= 1 )   temp = temp - m0(id(1),id(2),id(3),1)*ff(id(1),id(2),id(3)-1)
!        if ( oo /= fcz ) temp = temp - m0(id(1),id(2),id(3),7)*ff(id(1),id(2),id(3)+1)
!        ff(id(1),id(2),id(3)) = &
!            (1D0-relax)*ff(id(1),id(2),id(3))+relax*temp/m0(id(1),id(2),id(3),4)
!      end do
!      end do
!      end do
!    end do
!    end do
!    end do
!    end do
!
!end subroutine

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
!    !$omp parallel default(shared) private(ii,jj,kk,ll,id0,id,temp)
!    !$omp do
!    !$omp end do
!    !$omp end parallel
