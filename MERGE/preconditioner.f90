module PRECONDITIONER
    use FMFD_HEADER, only: n_lnodes, anode, n_fnodes
    use BICG_HEADER
    implicit none

    contains

! =============================================================================
! 
! =============================================================================
subroutine ILU_INITIAL
    use FMFD_HEADER, only: fcr, fcz, fc1, fc2, ax, ay, az, mvec1, &
                        i_para0, i_para1, anode, ncm, gx0, gx1, gy0, gy1, &
                        gz0, gz1, zigzagon, zz_div, zzc0, zzc1, zzc2, gcn
    use VARIABLES, only: icore, score, ierr
    implicit none
    integer:: n_entry
    integer:: num0, num1   ! matrix number, entry number
    integer:: num2         ! number in a row 
    integer:: io, jo, ko, mo, no

    ! ILU initialization
    n_entry =7*n_lnodes-4*fc1-2*fc2 
    allocate(am(n_entry),ia(n_lnodes+1),ja(n_entry))
    ia(n_lnodes+1) = n_entry+1

    maxfill  = 3
    tol      = 1D-8
    ipar(31) = 1
    dpar(31) = 1D-6
    dpar(32) = 1D-4

    allocate(ln(i_para0:i_para1,n_lnodes),un(i_para0:i_para1,n_lnodes))
    allocate(li(i_para0:i_para1,n_lnodes,maxfill+1))
    allocate(ui(i_para0:i_para1,n_lnodes,maxfill+1))
    li = 0; ui = 0

    allocate(bm(i_para0:i_para1,(2*maxfill+1)*n_lnodes-maxfill*(maxfill+1)+1))
    allocate(jb(i_para0:i_para1,(2*maxfill+1)*n_lnodes-maxfill*(maxfill+1)+1))
    allocate(ib(i_para0:i_para1,n_lnodes+1))


    ! entry setting
    num0 = 0
    num1 = 0
    do ko = 1, fcz
    do jo = 1, fcr
    do io = 1, fcr

        num0 = num0 + 1; num2 = 0
    
        if ( ko /= 1 ) then
            num1 = num1 + 1
            num2 = num2 + 1
            ja(num1) = num0 - fc2
        end if
        if ( jo /= 1 ) then
            num1 = num1 + 1
            num2 = num2 + 1
            ja(num1) = num0 - fcr
        end if
        if ( io /= 1 ) then
            num1 = num1 + 1
            num2 = num2 + 1
            ja(num1) = num0 - 1
        end if
    
        num1 = num1 + 1
        ia(num0) = num1-num2
        ja(num1) = num0
    
        if ( io /= fcr ) then
            num1 = num1 + 1
            ja(num1) = num0 + 1
        end if
        if ( jo /= fcr ) then
            num1 = num1 + 1
            ja(num1) = num0 + fcr
        end if
        if ( ko /= fcz ) then
            num1 = num1 + 1
            ja(num1) = num0 + fc2
        end if

    end do
    end do
    end do

    ! global parameters
    n_entry = anode*7-4*ncm(1)*ncm(3)-2*(anode/ncm(3))
    allocate(gam(n_entry),gia(anode+1),gja(n_entry))
    gia(anode+1) = n_entry+1

    allocate(gln(anode),gun(anode))
    allocate(gli(anode,maxfill+1))
    allocate(gui(anode,maxfill+1))
    gli = 0; gui = 0

    allocate(gbm((2*maxfill+1)*anode-maxfill*(maxfill+1)+1))
    allocate(gjb((2*maxfill+1)*anode-maxfill*(maxfill+1)+1))
    allocate(gib(anode+1))

    allocate(gx0(anode),gx1(anode),gy0(anode))
    allocate(gy1(anode),gz0(anode),gz1(anode))
    gx0 = 0; gx1 = 0; gy0 = 0; gy1 = 0; gz0 = 0; gz1 = 0


    if ( .not. allocated(gcn) ) then
        allocate(gcn(0:ncm(1)+1,0:ncm(2)+1,0:ncm(3)+1))

        gcn = 0; num1 = 0
        do ko = 1, ncm(3)
        do jo = 1, ncm(2)
        do io = 1, ncm(1)
        
            if ( zigzagon ) then
            do mo = 1, zz_div
            if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
                no = mo
                exit
            end if
            end do
            if ( zzc1(no) >= jo .or. jo > zzc2(no) ) cycle
            end if

            num1 = num1 + 1
            gcn(io,jo,ko) = num1

        end do
        end do
        end do
    
        num1 = 0
        do ko = 1, ncm(3)
        do jo = 1, ncm(2)
        do io = 1, ncm(1)
            
            if ( gcn(io,jo,ko) == 0 ) cycle

            num2 = 0
            if ( gcn(io,jo,ko-1) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                gja(num1) = gcn(io,jo,ko-1)
                gz0(gcn(io,jo,ko)) = gcn(io,jo,ko-1)
            end if
            if ( gcn(io,jo-1,ko) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                gja(num1) = gcn(io,jo-1,ko)
                gy0(gcn(io,jo,ko)) = gcn(io,jo-1,ko)
            end if
            if ( gcn(io-1,jo,ko) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                gja(num1) = gcn(io-1,jo,ko)
                gx0(gcn(io,jo,ko)) = gcn(io-1,jo,ko)
            end if
            num1 = num1 + 1
            gia(gcn(io,jo,ko)) = num1-num2
            gja(num1) = gcn(io,jo,ko)
            if ( gcn(io+1,jo,ko) /= 0 ) then
                num1 = num1 + 1
                gja(num1) = gcn(io+1,jo,ko)
                gx1(gcn(io,jo,ko)) = gcn(io+1,jo,ko)
            end if
            if ( gcn(io,jo+1,ko) /= 0 ) then
                num1 = num1 + 1
                gja(num1) = gcn(io,jo+1,ko)
                gy1(gcn(io,jo,ko)) = gcn(io,jo+1,ko)
            end if
            if ( gcn(io,jo,ko+1) /= 0 ) then
                num1 = num1 + 1
                gja(num1) = gcn(io,jo,ko+1)
                gz1(gcn(io,jo,ko)) = gcn(io,jo,ko+1)
            end if
    
        end do
        end do
        end do

    end if


end subroutine

! =============================================================================
! 
! =============================================================================
subroutine ILU_DECOMPOSE(Mfm)
    use FMFD_HEADER, only: fcr, fcz, anode, fc1, fc2, ax, ay, az, mvec1, &
                        i_para0, i_para1
    use VARIABLES, only: icore, score, ierr
    use mpi, only: mpi_wtime
    implicit none
    real(8), intent(inout):: Mfm(:,:,:,:)   ! (nfm(1),nfm(2),nfm(3),7)
    integer:: num0, num1   ! matrix number, entry number
    integer:: num2         ! number in a row 
    integer:: io, jo, ko, mo
    real(8) :: bmtmp((2*maxfill+1)*n_lnodes-maxfill*(maxfill+1)+1)
    integer :: ibtmp(n_lnodes+1)
    integer :: jbtmp((2*maxfill+1)*n_lnodes-maxfill*(maxfill+1)+1)

    ! first element    : ib(ii)+jo-1
    ! diagonal element : ib(ii)+ln(ii)+jo-2

    do mo = i_para0, i_para1

        ! entry setting
        num0 = 0
        num1 = 0
        do ko = az(mo)+1, az(mo)+fcz
        do jo = ay(mo)+1, ay(mo)+fcr
        do io = ax(mo)+1, ax(mo)+fcr

            num0 = num0 + 1
    
            if ( ko /= az(mo)+1 ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,1)
            end if
            if ( jo /= ay(mo)+1 ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,2)
            end if
            if ( io /= ax(mo)+1 ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,3)
            end if
    
            num1 = num1 + 1
            am(num1) = Mfm(io,jo,ko,4)
            mvec1(mo,num0,:) = Mfm(io,jo,ko,:)
    
            if ( io /= ax(mo)+fcr ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,5)
            end if
            if ( jo /= ay(mo)+fcr ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,6)
            end if
            if ( ko /= az(mo)+fcz ) then
                num1 = num1 + 1
                am(num1) = Mfm(io,jo,ko,7)
            end if

        end do
        end do
        end do

        ! ILU decomposition
        call DCSRILUT(n_lnodes,am,ia,ja,bmtmp,ibtmp,jbtmp,tol,maxfill,ipar,dpar,ierr)
        bm(mo,:) = bmtmp
        ib(mo,:) = ibtmp
        jb(mo,:) = jbtmp

        ! indice of ILU decomposed matrix
        do io = 1, n_lnodes
            ln(mo,io) = count(jb(mo,ib(mo,io):ib(mo,io+1)-1)<=io)
            un(mo,io) = count(jb(mo,ib(mo,io):ib(mo,io+1)-1)>=io)
            do jo = 1, ln(mo,io)
                li(mo,io,jo) = ib(mo,io)+jo-1
            end do
            do jo = 1, un(mo,io)
                ui(mo,io,jo) = ib(mo,io)+ln(mo,io)+jo-2
            end do
        end do
    
    end do

end subroutine



! =============================================================================
! 
! =============================================================================
subroutine GLOBAL_ILU
    use VARIABLES, only: n_inact, curr_cyc    ! ***
    use FMFD_HEADER, only: ncm, anode, Mcm, gcn, mvec2
    use FMFD_HEADER, only: deltc0, deltc1
    implicit none
    integer:: num0, num1, ii, jj, kk, ierr

    num0 = 0; num1 = 0
    do kk = 1, ncm(3)
    do jj = 1, ncm(2)
    do ii = 1, ncm(1)
        
        if ( gcn(ii,jj,kk) == 0 ) cycle
        num0 = num0 + 1

        if ( gcn(ii,jj,kk-1) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,1)
        end if
        if ( gcn(ii,jj-1,kk) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,2)
        end if
        if ( gcn(ii-1,jj,kk) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,3)
        end if
    
        num1 = num1 + 1
        gam(num1) = Mcm(ii,jj,kk,4)
        mvec2(num0,:) = Mcm(ii,jj,kk,:)
    
        if ( gcn(ii+1,jj,kk) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,5)
        end if
        if ( gcn(ii,jj+1,kk) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,6)
        end if
        if ( gcn(ii,jj,kk+1) /= 0 ) then
            num1 = num1 + 1
            gam(num1) = Mcm(ii,jj,kk,7)
        end if

    end do
    end do
    end do
    print *, 'GLOBALILU:MVEC', sum(mvec2), sum(Mcm)

    call DCSRILUT(anode,gam,gia,gja,gbm,gib,gjb,tol,maxfill,ipar,dpar,ierr)

    do ii = 1, anode
        gln(ii) = count(gjb(gib(ii):gib(ii+1)-1)<=ii)
        gun(ii) = count(gjb(gib(ii):gib(ii+1)-1)>=ii)

        do jj = 1, gln(ii)
            gli(ii,jj) = gib(ii)+jj-1
        end do
        do jj = 1, gun(ii)
            gui(ii,jj) = gib(ii)+gln(ii)+jj-2
        end do
    end do

end subroutine


! =============================================================================
! LOCAL
! =============================================================================
function FORWARD_SUB(mo,bb) result(xx)
    implicit none
    integer, intent(in):: mo
    real(8), intent(in):: bb(:)
    real(8):: xx(n_lnodes)
    integer:: io

    ! forward substitution
    xx(1) = bb(1)
    do io = 2, n_lnodes
        xx(io) = bb(io) &
            - sum(bm(mo,li(mo,io,1:ln(mo,io)-1))*xx(jb(mo,li(mo,io,1:ln(mo,io)-1))))
    end do

end function

function ILU_SOLVER(mo,bb) result(xx)
    use VARIABLES, only: icore
    implicit none
    integer, intent(in):: mo
    real(8), intent(in):: bb(:)
    real(8):: xx(n_lnodes), yy(n_lnodes)
    integer:: io

    ! forward substitution
    yy(1) = bb(1)
    do io = 2, n_lnodes
        yy(io) = bb(io) &
            - sum(bm(mo,li(mo,io,1:ln(mo,io)-1))*yy(jb(mo,li(mo,io,1:ln(mo,io)-1))))
    end do

    ! backward substitution
    xx(n_lnodes) = yy(n_lnodes)/bm(mo,ib(mo,n_lnodes+1)-1)
    do io = n_lnodes-1, 1, -1
        xx(io) = (yy(io)-sum(bm(mo,ui(mo,io,2:un(mo,io))) &
            *xx(jb(mo,ui(mo,io,2:un(mo,io))))))/bm(mo,ui(mo,io,1))
    end do

end function


! =============================================================================
! GLOBAL
! =============================================================================
function FORWARD_GSUB(bb) result(xx)
    implicit none
    real(8), intent(in):: bb(:)
    real(8):: xx(anode)
    integer:: io

    ! forward substitution
    xx(1) = bb(1)
    do io = 2, anode
        xx(io) = bb(io) &
            - sum(gbm(gli(io,1:gln(io)-1))*xx(gjb(gli(io,1:gln(io)-1))))
    end do

end function


function ILU_GSOLVER(bb) result(xx)
    implicit none
    real(8), intent(in):: bb(:)
    real(8):: xx(anode), yy(anode)
    integer:: io

    ! forward substitution
    yy(1) = bb(1)
    do io = 2, anode
        yy(io) = bb(io) &
            - sum(gbm(gli(io,1:gln(io)-1))*yy(gjb(gli(io,1:gln(io)-1))))
    end do

    ! backward substitution
    xx(anode) = yy(anode)/gbm(gib(anode+1)-1)
    do io = anode-1, 1, -1
        xx(io) = (yy(io)-sum(gbm(gui(io,2:gun(io))) &
            *xx(gjb(gui(io,2:gun(io))))))/gbm(gui(io,1))
    end do

end function

!function BACK_SUB(aa,bb,mm) result(xx)
!    implicit none
!    real(8):: aa(:), bb(:), xx(n_lnodes)
!    integer:: ii, jj, mm
!
!    ! backward substitution
!    xx(n_lnodes) = bb(n_lnodes)/aa(ib(n_lnodes+1)-1)
!    do ii = n_lnodes-1, 1, -1
!        xx(ii) = (bb(ii)-sum(aa(ui(mm,ii,1:un(mm,ii))) &
!            *bb(jb(ui(mm,ii,1:un(mm,ii))))))/aa(ui(mm,ii,1))
!    end do
!
!end function




! =============================================================================
! =============================================================================
! =============================================================================
! 
! =============================================================================
subroutine FINE_ILU_INITIAL
    use FMFD_HEADER, only: bs0, nfm, fx0, fx1, fy0, fy1, fz0, fz1, &
                        zigzagon, fn, zz_div, zzc0, zzc1, zzc2, mvec3, fphi1
    implicit none
    integer:: n_entry
    integer:: num0, num1   ! matrix number, entry number
    integer:: num2         ! number in a row 
    integer:: io, jo, ko, mo, no

    ! ILU initialization
    n_fnodes = bs0
    n_entry = n_fnodes*7-4*nfm(1)*nfm(3)-2*(n_fnodes/nfm(3))
    allocate(fam(n_entry),fia(n_fnodes+1),fja(n_entry))
    fia(n_fnodes+1) = n_entry+1

    maxfill  = 3
    tol      = 1D-8
    ipar(31) = 1
    dpar(31) = 1D-6
    dpar(32) = 1D-4

    allocate(fln(n_fnodes),fun(n_fnodes))
    allocate(fli(n_fnodes,maxfill+1))
    allocate(fui(n_fnodes,maxfill+1))
    fli = 0; fui = 0

    allocate(fbm((2*maxfill+1)*n_fnodes-maxfill*(maxfill+1)+1))
    allocate(fjb((2*maxfill+1)*n_fnodes-maxfill*(maxfill+1)+1))
    allocate(fib(n_fnodes+1))

    allocate(fx0(n_fnodes),fx1(n_fnodes),fy0(n_fnodes))
    allocate(fy1(n_fnodes),fz0(n_fnodes),fz1(n_fnodes))
    fx0 = 0; fx1 = 0; fy0 = 0; fy1 = 0; fz0 = 0; fz1 = 0

    allocate(mvec3(n_fnodes,7))


    ! entry setting
    if ( .not. allocated(fn) ) then
        allocate(fn(0:nfm(1)+1,0:nfm(2)+1,0:nfm(3)+1))

        fn = 0; num1 = 0
        do ko = 1, nfm(3)
        do jo = 1, nfm(2)
        do io = 1, nfm(1)
            if ( fphi1(io,jo,ko) == 0 ) cycle
            num1 = num1 + 1
            fn(io,jo,ko) = num1
        end do
        end do
        end do
    
        num1 = 0
        do ko = 1, nfm(3)
        do jo = 1, nfm(2)
        do io = 1, nfm(1)
            
            if ( fn(io,jo,ko) == 0 ) cycle

            num2 = 0
            if ( fn(io,jo,ko-1) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                fja(num1) = fn(io,jo,ko-1)
                fz0(fn(io,jo,ko)) = fn(io,jo,ko-1)
            end if
            if ( fn(io,jo-1,ko) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                fja(num1) = fn(io,jo-1,ko)
                fy0(fn(io,jo,ko)) = fn(io,jo-1,ko)
            end if
            if ( fn(io-1,jo,ko) /= 0 ) then
                num1 = num1 + 1
                num2 = num2 + 1
                fja(num1) = fn(io-1,jo,ko)
                fx0(fn(io,jo,ko)) = fn(io-1,jo,ko)
            end if
            num1 = num1 + 1
            fia(fn(io,jo,ko)) = num1-num2
            fja(num1) = fn(io,jo,ko)
            if ( fn(io+1,jo,ko) /= 0 ) then
                num1 = num1 + 1
                fja(num1) = fn(io+1,jo,ko)
                fx1(fn(io,jo,ko)) = fn(io+1,jo,ko)
            end if
            if ( fn(io,jo+1,ko) /= 0 ) then
                num1 = num1 + 1
                fja(num1) = fn(io,jo+1,ko)
                fy1(fn(io,jo,ko)) = fn(io,jo+1,ko)
            end if
            if ( fn(io,jo,ko+1) /= 0 ) then
                num1 = num1 + 1
                fja(num1) = fn(io,jo,ko+1)
                fz1(fn(io,jo,ko)) = fn(io,jo,ko+1)
            end if
    
        end do
        end do
        end do

    end if


end subroutine

! =============================================================================
! 
! =============================================================================
subroutine FINE_ILU
    use FMFD_HEADER, only: nfm, Mfm, fn, mvec3
    implicit none
    integer:: num0, num1, ii, jj, kk, ierr

    num0 = 0; num1 = 0
    do kk = 1, nfm(3)
    do jj = 1, nfm(2)
    do ii = 1, nfm(1)
        
        if ( fn(ii,jj,kk) == 0 ) cycle
        num0 = num0 + 1

        if ( fn(ii,jj,kk-1) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,1)
        end if
        if ( fn(ii,jj-1,kk) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,2)
        end if
        if ( fn(ii-1,jj,kk) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,3)
        end if
    
        num1 = num1 + 1
        fam(num1) = Mfm(ii,jj,kk,4)
        mvec3(num0,:) = Mfm(ii,jj,kk,:)
    
        if ( fn(ii+1,jj,kk) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,5)
        end if
        if ( fn(ii,jj+1,kk) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,6)
        end if
        if ( fn(ii,jj,kk+1) /= 0 ) then
            num1 = num1 + 1
            fam(num1) = Mfm(ii,jj,kk,7)
        end if

    end do
    end do
    end do

    call DCSRILUT(n_fnodes,fam,fia,fja,fbm,fib,fjb,tol,maxfill,ipar,dpar,ierr)

    do ii = 1, n_fnodes
        fln(ii) = count(fjb(fib(ii):fib(ii+1)-1)<=ii)
        fun(ii) = count(fjb(fib(ii):fib(ii+1)-1)>=ii)

        do jj = 1, fln(ii)
            fli(ii,jj) = fib(ii)+jj-1
        end do
        do jj = 1, fun(ii)
            fui(ii,jj) = fib(ii)+fln(ii)+jj-2
        end do
    end do

end subroutine


! =============================================================================
! FINE
! =============================================================================
function FORWARD_FSUB(bb) result(xx)
    implicit none
    real(8), intent(in):: bb(:)
    real(8):: xx(n_fnodes)
    integer:: io

    ! forward substitution
    xx(1) = bb(1)
    do io = 2, n_fnodes
        xx(io) = bb(io) &
            - sum(fbm(fli(io,1:fln(io)-1))*xx(fjb(fli(io,1:fln(io)-1))))
    end do

end function


function ILU_FSOLVER(bb) result(xx)
    implicit none
    real(8), intent(in):: bb(:)
    real(8):: xx(n_fnodes), yy(n_fnodes)
    integer:: io

    ! forward substitution
    yy(1) = bb(1)
    do io = 2, n_fnodes
        yy(io) = bb(io) &
            - sum(fbm(fli(io,1:fln(io)-1))*yy(fjb(fli(io,1:fln(io)-1))))
    end do

    ! backward substitution
    xx(n_fnodes) = yy(n_fnodes)/fbm(fib(n_fnodes+1)-1)
    do io = n_fnodes-1, 1, -1
        xx(io) = (yy(io)-sum(fbm(fui(io,2:fun(io))) &
            *xx(fjb(fui(io,2:fun(io))))))/fbm(fui(io,1))
    end do

end function


end module
