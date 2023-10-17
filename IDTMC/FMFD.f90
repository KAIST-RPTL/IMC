module FMFD
    use CONSTANTS,          only: eVtoJoule
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use particle_header,    only : particle
    use FMFD_HEADER
    use CMFD, only: ONE_NODE_CMFD, CMFD_CALCULATION
    use TH_HEADER, only: th_on, pp
    !use STATISTICS, only: STD_MAT
    
    implicit none
    
    contains 


! =============================================================================
! FMFD_initialize_1
! =============================================================================
subroutine FMFD_ALLOCATION()
    !use ENTROPY, only: entrp3
    use PERTURBATION, only: perton
    use COSAMPLING, only: n_pert
    use PCMFD, only: OUT_OF_ZZ, OUT_OF_ZZ0
    use PRECONDITIONER, only: ILU_INITIAL
    implicit none
    integer:: num

    ! parameters allocation
    allocate(k_fmfd(n_batch,n_totcyc))
    if ( icore == score ) allocate(p_fmfd(n_batch,n_act,nfm(1),nfm(2),nfm(3)))
    allocate(fsd_MC(nfm(1),nfm(2),nfm(3)))
    allocate(fsd(nfm(1),nfm(2),nfm(3)))
    if ( dual_fmfd ) then
    allocate(k_fmfd2(n_batch,n_totcyc))
    allocate(fsd2(nfm(1),nfm(2),nfm(3)))
    if ( icore == score ) then
    allocate(p_fmfd2(n_batch,n_act,nfm(1),nfm(2),nfm(3)))
    end if
    end if
    allocate(fm(nfm(1),nfm(2),nfm(3)))
    if(perton) allocate(k_real(n_batch,n_totcyc,n_pert))

    bs0 = nfm(1)*nfm(2)*nfm(3)
    if ( cmfdon ) then
    if ( .not. wholecore ) then
        ! parallel calculation for 1-CMFD
        anode = 0
        if ( zigzagon ) then
            do ii = 1, zz_div
                anode = anode + 2*zzc1(ii)*(zzc0(ii+1)-zzc0(ii))
            end do
        end if
        anode = (ncm(1)*ncm(2) - anode) * ncm(3)
        allocate(mvec(anode,fcr,fcr,fcz,7),svec(anode,fcr,fcr,fcz))
        allocate(mvec1(anode,n_lnodes,7),svec1(anode,n_lnodes))
        allocate(mvec2(anode,7),svec2(anode))
        allocate(pcd(anode,fcr,fcr,fcz))
        allocate(ax(anode),ay(anode),az(anode))

        anode = 0
        do kk = 1, ncm(3)
        do jj = 1, ncm(2)
        do ii = 1, ncm(1)
            if ( OUT_OF_ZZ(ii,jj) ) cycle
            anode = anode + 1
            ax(anode) = (ii-1)*fcr
            ay(anode) = (jj-1)*fcr
            az(anode) = (kk-1)*fcz
        end do
        end do
        end do

        bs0 = anode*fcr*fcr*fcz
        bs1 = bs0*7

    else
        ! parallel calculation for active 1-CMFD
        anode = 0
        if ( zigzagon ) then
            do ii = 1, zz_div
                anode = anode + 2*zzc1(ii)*(zzc0(ii+1)-zzc0(ii))
            end do
        end if
        anode = anode + 4*(ncm(1)-1)
        anode = (ncm(1)*ncm(2) - anode) * (ncm(3)-2)
        allocate(mvec(anode,fcr,fcr,fcz,7),svec(anode,fcr,fcr,fcz))
        allocate(pcd(anode,fcr,fcr,fcz))
        allocate(ax(anode),ay(anode),az(anode))

        anode = 0
        do kk = 1, ncm(3)
        do jj = 1, ncm(2)
        do ii = 1, ncm(1)
            if ( OUT_OF_ZZ0(ii,jj,kk) ) cycle
            anode = anode + 1
            ax(anode) = (ii-1)*fcr
            ay(anode) = (jj-1)*fcr
            az(anode) = (kk-1)*fcz
        end do
        end do
        end do
        bs0 = anode*fcr*fcr*fcz
        bs1 = bs0*7
    end if

    allocate(lx0(n_lnodes),lx1(n_lnodes),ly0(n_lnodes))
    allocate(ly1(n_lnodes),lz0(n_lnodes),lz1(n_lnodes))
    lx0 = 0; lx1 = 0; ly0 = 0; ly1 = 0; lz0 = 0; lz1 = 0;
    num = 0
    do kk = 1, fcz
    do jj = 1, fcr
    do ii = 1, fcr
        num = num + 1
        if ( kk /= 1 )   lz0(num) = num - fc2 
        if ( jj /= 1 )   ly0(num) = num - fcr
        if ( ii /= 1 )   lx0(num) = num - 1
        if ( ii /= fcr ) lx1(num) = num + 1
        if ( jj /= fcr ) ly1(num) = num + fcr
        if ( kk /= fcz ) lz1(num) = num + fc2
    end do
    end do
    end do
    end if

    allocate(acc(n_acc)) 
    do ii = 1, n_acc 
        allocate(acc(ii)%fm(nfm(1),nfm(2),nfm(3)))
    enddo

    allocate(Mfm(nfm(1),nfm(2),nfm(3),7)); Mfm = 0
    allocate(fphi1(nfm(1),nfm(2),nfm(3)),fm_s(nfm(1),nfm(2),nfm(3)))

    if ( cmfdon ) then
        call MPI_RANGE(anode,ncore,icore,i_para0,i_para1)
        call ILU_INITIAL
    end if

    ! depletion
    !if ( DTMCBU ) call DTMC_BU_INITIAL

!    if ( icore /= score ) return

    allocate(fm_avg(nfm(1),nfm(2),nfm(3)))
    allocate(fsd_FM(nfm(1),nfm(2),nfm(3)))
    

    ! FMFD parameters
    allocate(fm_t   (nfm(1),nfm(2),nfm(3)), &
             fmD    (nfm(1),nfm(2),nfm(3)), &
             fm_a   (nfm(1),nfm(2),nfm(3)), &
             fm_nf  (nfm(1),nfm(2),nfm(3)), &
             kappa  (nfm(1),nfm(2),nfm(3)), &
             fphi0  (nfm(1),nfm(2),nfm(3)), &
             fmJ0   (nfm(1),nfm(2),nfm(3),6), &
             fmJ1   (nfm(1),nfm(2),nfm(3),6), &
             fmJn   (nfm(1),nfm(2),nfm(3),6), &
             fmDt   (nfm(1),nfm(2),nfm(3),6), &
             fmDh   (nfm(1),nfm(2),nfm(3),6))
    allocate(ptJn   (500,nfm(1),nfm(2),nfm(3),6))
    allocate(ptJ0   (500,nfm(1),nfm(2),nfm(3),6))
    allocate(ptJ1   (500,nfm(1),nfm(2),nfm(3),6))
    !allocate(ptphi  (500,nfm(1),nfm(2),nfm(3)))
    
    ! area and volume of mesh cell
    a_fm(1:4) = dfm(1)*dfm(3)
    a_fm(5:6) = dfm(1)*dfm(2)
    v_fm    = dfm(1)*dfm(2)*dfm(3)


    ! CMFD parameters
    if ( cmfdon ) then
        allocate(cm_t   (ncm(1),ncm(2),ncm(3)), &
                 cmD    (ncm(1),ncm(2),ncm(3)), &
                 cm_a   (ncm(1),ncm(2),ncm(3)), &
                 cm_nf  (ncm(1),ncm(2),ncm(3)), &
                 cm_s   (ncm(1),ncm(2),ncm(3)), &
                 cphi0  (ncm(1),ncm(2),ncm(3)), &
                 cphi1  (ncm(1),ncm(2),ncm(3)), &
                 deltf0 (nfm(1),nfm(2),nfm(3),6), &
                 deltf1 (nfm(1),nfm(2),nfm(3),6), &
                 deltc0 (ncm(1),ncm(2),ncm(3),6), &
                 deltc1 (ncm(1),ncm(2),ncm(3),6), &
                 jsrc   (nfm(1),nfm(2),nfm(3),6), &
                 fsrc   (nfm(1),nfm(2),nfm(3),6), &
                 cmJ0   (ncm(1),ncm(2),ncm(3),6), &
                 cmJ1   (ncm(1),ncm(2),ncm(3),6), &
                 cmJn   (ncm(1),ncm(2),ncm(3),6), &
                 cmF    (ncm(1),ncm(2),ncm(3),6), &
                 cmDt   (ncm(1),ncm(2),ncm(3),6), &
                 cmDh   (ncm(1),ncm(2),ncm(3),6), &
                 Mcm    (ncm(1),ncm(2),ncm(3),7))
        cm_t    = 0
        cmD     = 0
        cm_a    = 0
        cm_nf   = 0
        cm_s    = 0
        cphi0   = 0
        cphi1   = 0
        deltf0  = 0
        deltf1  = 0
        deltc0  = 0
        deltc1  = 0
        jsrc    = 0
        fsrc    = 0
        cmJ0    = 0
        cmJ1    = 0
        cmJn    = 0
        cmF     = 0
        cmDt    = 0
        cmDh    = 0
        Mcm     = 0

    end if

!    if ( DTMCBU ) then
!        allocate(p_dtmc(nfm(1),nfm(2),nfm(3),n_rings))
!        allocate(f_dtmc(nfm(1),nfm(2),nfm(3),n_rings))
!    end if

end subroutine


! =============================================================================
! FMFD_INITIALIZE_2 initializes parameters used in the FMFD calculaiton
! =============================================================================
subroutine FMFD_initialize()
    implicit none
    integer:: ij

    ! initialization
    fm(:,:,:) % phi     = 0 
    fm(:,:,:) % sig_t   = 0 
    fm(:,:,:) % sig_a   = 0 
    fm(:,:,:) % nusig_f = 0 
    fm(:,:,:) % kappa   = 0 
    do ij = 1, 6
        fm(:,:,:) % J0(ij)   = 0 
        fm(:,:,:) % J1(ij)   = 0 
    enddo 

    fsd_MC = 0

!    if ( DTMCBU ) then
!        !intra_phi = 0
!        !intra_pow = 0
!        do ii = 1, nfm(1)
!        do jj = 1, nfm(2)
!        do kk = 1, nfm(3)
!        cede(ii,jj,kk)%flux(:) = 0
!        cede(ii,jj,kk)%kapa(:) = 0
!        end do
!        end do
!        end do
!    end if

    if ( DTMCBU .and. curr_cyc == n_totcyc ) then
        allocate(tmflux(nfm(1),nfm(2),nfm(3),nrings))
        !allocate(tmkapa(nfm(1),nfm(2),nfm(3),nrings))
        tmflux = 0
        !tmkapa = 0
    end if
        
end subroutine

! =============================================================================
! FMFD_INITIALIZE_THREAD initializes thread-wise parameters
! =============================================================================
subroutine FMFD_initialize_thread() 
    integer:: ij
    integer:: mm, nn, oo
    
    if ( .not. allocated(fm_thread) ) allocate(fm_thread(nfm(1),nfm(2),nfm(3)))
    
    fm_thread(:,:,:) % phi     = 0 
    fm_thread(:,:,:) % sig_t   = 0 
    fm_thread(:,:,:) % sig_a   = 0 
    fm_thread(:,:,:) % nusig_f = 0 
    fm_thread(:,:,:) % kappa   = 0 
    do ij = 1, 6
    fm_thread(:,:,:) % J0(ij)  = 0 
    fm_thread(:,:,:) % J1(ij)  = 0 
    enddo 

    if ( DTMCBU ) then
        ! allocation
        if ( .not. allocated(thflux) ) then
            allocate(thflux(nfm(1),nfm(2),nfm(3),nrings))
            !allocate(thkapa(nfm(1),nfm(2),nfm(3),nrings))
        end if

        ! initialization
        if ( curr_cyc == 1 ) then
        thflux = 0
        !thkapa = 0
        end if

    end if
    
end subroutine
    

!! =============================================================================
!! 
!! =============================================================================
!subroutine DTMC_BU_INITIAL
!    use CONSTANTS, only: pi
!    implicit none
!
!    ! allocation
!    allocate(v_ring(n_rings))
!    allocate(intra_phi(nfm(1),nfm(2),nfm(3),n_rings+1))
!    allocate(intra_pow(nfm(1),nfm(2),nfm(3),n_rings))
!    if ( iscore ) then
!    allocate(intra_kap(acc_skip+1:n_totcyc,nfm(1),nfm(2),nfm(3),n_rings))
!    allocate(f_dep_mc(acc_skip+1:n_totcyc,nfm(1),nfm(2),nfm(3),n_rings))
!    allocate(p_dep_mc(acc_skip+1:n_totcyc,nfm(1),nfm(2),nfm(3),n_rings))
!    allocate(f_dep_dt(n_act,nfm(1),nfm(2),nfm(3),n_rings))
!    allocate(p_dep_dt(n_act,nfm(1),nfm(2),nfm(3),n_rings))
!    end if
!    !$omp parallel
!    allocate(dth_phi(nfm(1),nfm(2),nfm(3),n_rings+1))
!    allocate(dth_pow(nfm(1),nfm(2),nfm(3),n_rings))
!    !$omp end parallel
!
!
!    ! volume of the rings
!    v_ring(1) = pi*rings(1)*rings(1)*dfm(3)
!    do ii = 2, n_rings
!    v_ring(ii) = pi*(rings(ii)*rings(ii)-rings(ii-1)*rings(ii-1))*dfm(3)
!    end do
!
!end subroutine


! =============================================================================
! ISINCOMING determines if the particle is coming to the FMFD grid
! =============================================================================
subroutine INCOMING(xyz,income)
    integer, intent(in):: xyz(:)
    integer, intent(out):: income

    ! x0
    if ( xyz(1) == 0 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 2
        end if
        end if
    end if
    ! x1
    if ( xyz(1) == nfm(1)+1 ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
            income = 1
        end if
        end if
    end if
    ! y0
    if ( xyz(2) == 0 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 4
        end if
        end if
    end if
    ! y1
    if ( xyz(2) == nfm(2)+1 ) then
        if ( 1 <= xyz(3) .and. xyz(3) <= nfm(3) ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
            income = 3
        end if
        end if
    end if
    ! z0
    if ( xyz(3) == 0 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 6
        end if
        end if
    end if
    ! z1
    if ( xyz(3) == nfm(3)+1 ) then
        if ( 1 <= xyz(1) .and. xyz(1) <= nfm(1) ) then
        if ( 1 <= xyz(2) .and. xyz(2) <= nfm(2) ) then
            income = 5
        end if
        end if
    end if

end subroutine

! =============================================================================
! FMFD_TRK calculates FMFD parameters such as flux, group contstans by 
! track-length estiamtor
! =============================================================================
subroutine FMFD_TRK(wgt,distance,macro_xs,id)
    implicit none
    type(Particle):: p
    real(8), intent(in) :: wgt
    real(8), intent(in) :: distance
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(1:3)
    real(8) :: flux
    
    flux = wgt * distance

    fm_thread(id(1),id(2),id(3)) % phi = & 
    fm_thread(id(1),id(2),id(3)) % phi + flux
    fm_thread(id(1),id(2),id(3)) % sig_t = &
    fm_thread(id(1),id(2),id(3)) % sig_t + flux*macro_xs(1)
    fm_thread(id(1),id(2),id(3)) % sig_a = &
    fm_thread(id(1),id(2),id(3)) % sig_a + flux*macro_xs(2)
    fm_thread(id(1),id(2),id(3)) % nusig_f = &
    fm_thread(id(1),id(2),id(3)) % nusig_f + flux*macro_xs(4)
    fm_thread(id(1),id(2),id(3)) % kappa = &
    fm_thread(id(1),id(2),id(3)) % kappa + flux*macro_xs(5)
    
end subroutine

! =============================================================================
! FMFD_COL calculates FMFD parameters such as flux, group contstans by 
! collision estimator
! =============================================================================
subroutine FMFD_COL(wgt, macro_xs,id)
    real(8), intent(in) :: wgt
    real(8), intent(in) :: macro_xs(5)
    integer, intent(in) :: id(3)
    real(8) :: flux
    
    flux = wgt / macro_xs(1)
    
    fm_thread(id(1),id(2),id(3)) % phi = & 
    fm_thread(id(1),id(2),id(3)) % phi + flux
    fm_thread(id(1),id(2),id(3)) % sig_t = &
    fm_thread(id(1),id(2),id(3)) % sig_t + wgt
    fm_thread(id(1),id(2),id(3)) % sig_a = &
    fm_thread(id(1),id(2),id(3)) % sig_a + flux*macro_xs(2)
    fm_thread(id(1),id(2),id(3)) % nusig_f = &
    fm_thread(id(1),id(2),id(3)) % nusig_f + flux*macro_xs(4)
    fm_thread(id(1),id(2),id(3)) % kappa = &
    fm_thread(id(1),id(2),id(3)) % kappa + flux*macro_xs(5)
    
end subroutine


! =============================================================================
! FMFD_SURF calculates FMFD surface parameters like net and particle current
! =============================================================================
subroutine FMFD_SURF (inside,income, is, id, uvw, wgt, bc)
    logical, intent(in) :: inside
    integer, intent(in) :: income
    integer, intent(in) :: is, id(3)
    real(8), intent(in) :: uvw(3)
    real(8), intent(in) :: wgt
    integer, intent(in) :: bc
    integer:: dir

    if ( .not. fmfdon ) return
    
    ! print '(I1, I3, I3, I2, L2, I2, I2, F8.3, I2)', icore, id(1:3), inside, income, is, wgt, bc
    ! inner nodes
    if ( inside ) then 
        ! surface partial current
        select case(is)
        case(1,3,5)
            fm_thread(id(1),id(2),id(3))%J0(is) = &
            fm_thread(id(1),id(2),id(3))%J0(is) + wgt
        case(2,4,6)
            fm_thread(id(1),id(2),id(3))%J1(is) = &
            fm_thread(id(1),id(2),id(3))%J1(is) + wgt
        end select

        ! boundary condition
        if ( bc == 2 ) then
        select case(is)
        case(1,3,5)
            fm_thread(id(1),id(2),id(3))%J1(is) = &
            fm_thread(id(1),id(2),id(3))%J1(is) + wgt
        case(2,4,6)
            fm_thread(id(1),id(2),id(3))%J0(is) = &
            fm_thread(id(1),id(2),id(3))%J0(is) + wgt
        end select
        end if
        return
    end if

    ! boundary nodes
    select case(income)
    case(1)
        fm_thread(1,id(2),id(3))%J1(1) = &
        fm_thread(1,id(2),id(3))%J1(1) + wgt
    case(2)
        fm_thread(nfm(1),id(2),id(3))%J0(2) = &
        fm_thread(nfm(1),id(2),id(3))%J0(2) + wgt
    case(3)
        fm_thread(id(1),1,id(3))%J1(3) = &
        fm_thread(id(1),1,id(3))%J1(3) + wgt
    case(4)
        fm_thread(id(1),nfm(2),id(3))%J0(4) = &
        fm_thread(id(1),nfm(2),id(3))%J0(4) + wgt
    case(5)
        fm_thread(id(1),id(2),1)%J1(5) = &
        fm_thread(id(1),id(2),1)%J1(5) + wgt
    case(6)
        fm_thread(id(1),id(2),nfm(3))%J0(6) = &
        fm_thread(id(1),id(2),nfm(3))%J0(6) + wgt
    end select
            
end subroutine

!! =============================================================================
!! 
!! =============================================================================
!subroutine WHICH_RING(xyz,ring)
!    implicit none
!    real(8), intent(in) :: xyz(3)
!    integer, intent(out):: ring
!    integer:: rr
!
!    ring = n_rings+1
!    do rr = 1, n_rings
!        if ( xyz(1)*xyz(1)+xyz(2)*xyz(2)-rings(rr)*rings(rr) < 0 ) then
!            ring = rr
!            exit
!        end if
!    end do
!
!end subroutine

! =============================================================================
! 
! =============================================================================
subroutine DTMC_BU_TRK(id,id_ce,flux,macro_xs)
    use GEOMETRY_HEADER, only: cells
    implicit none
    integer, intent(in):: id(3), id_ce
    real(8), intent(in):: flux
    real(8), intent(in):: macro_xs

    if ( id_ce == 0 ) return

    ! flux
    thflux(id(1),id(2),id(3),id_ce) = &
    thflux(id(1),id(2),id(3),id_ce) + flux
!    ! kappa x sigma_f
!    thkapa(id(1),id(2),id(3),id_ce) = &
!    thkapa(id(1),id(2),id(3),id_ce) + flux * macro_xs

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine DTMC_BU_COL(id,id_ce,wgt,macro_xs)
    use GEOMETRY_HEADER, only: cells
    implicit none
    integer, intent(in):: id(3), id_ce
    real(8), intent(in):: wgt
    real(8), intent(in):: macro_xs

    if ( id_ce == 0 ) return

    ! flux
    thflux(id(1),id(2),id(3),id_ce) = &
    thflux(id(1),id(2),id(3),id_ce) + wgt
!    ! kappa x sigma_f
!    thkapa(id(1),id(2),id(3),id_ce) = &
!    thkapa(id(1),id(2),id(3),id_ce) + wgt * macro_xs

end subroutine


! =============================================================================
! NORM_FMFD normalizes cycle-wise FMFD parameters
! =============================================================================
subroutine NORM_FMFD(cyc)
    implicit none
    integer, intent(in):: cyc
    integer:: i, j, k, ij

    !> gather thread FMFD parameters
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
    fm(i,j,k)%phi     = fm(i,j,k)%phi     + fm_thread(i,j,k)%phi
    fm(i,j,k)%sig_t   = fm(i,j,k)%sig_t   + fm_thread(i,j,k)%sig_t 
    fm(i,j,k)%sig_a   = fm(i,j,k)%sig_a   + fm_thread(i,j,k)%sig_a 
    fm(i,j,k)%nusig_f = fm(i,j,k)%nusig_f + fm_thread(i,j,k)%nusig_f 
    fm(i,j,k)%kappa   = fm(i,j,k)%kappa   + fm_thread(i,j,k)%kappa 
    do mm = 1, 6
    fm(i,j,k)%J0(mm)  = fm(i,j,k)%J0(mm)  + fm_thread(i,j,k)%J0(mm)
    fm(i,j,k)%J1(mm)  = fm(i,j,k)%J1(mm)  + fm_thread(i,j,k)%J1(mm)
    end do
    end do
    end do
    end do
    

    !> gather thread intra-pin FMFD parameters for depletion
    if ( DTMCBU .and. cyc == n_totcyc ) then
        tmflux = tmflux + thflux
        !tmkapa = tmkapa + thkapa
    end if

end subroutine


! =============================================================================
! PROCESS_FMFD deals with MPI process and average quantities
! =============================================================================
subroutine PROCESS_FMFD(bat,cyc)
    use GEOMETRY_HEADER, only: universes
    use TALLY, only: tallyon
    use CMFD, only: L_DTILDA
    use PCMFD, only: L_PDHAT
    implicit none
    !> MPI derived type reduce parameters 
    integer, intent(in):: bat
    integer, intent(in):: cyc
    integer :: dsize    ! data size
    real(8) :: aa, bb
    integer :: ij, lc
    real(8), allocatable:: fsd_MC0(:,:,:)
    real(8), allocatable:: intra_flux(:,:,:,:), intra_kapa(:,:,:,:)
    type(FMFD_accumulation), pointer:: ac
    real(8):: tt0, tt1
    real(8):: Jn(nfm(1),nfm(2),nfm(3),6)
    integer:: length, which_idtmc
    integer:: reactor_type
    integer:: ncell


    ! -------------------------------------------------------------------------
    ! data transmission I
    allocate(fsd_MC0(nfm(1),nfm(2),nfm(3)))

    dsize = nfm(1)*nfm(2)*nfm(3)
    lc = mod(cyc-1,n_acc)+1

    ac => acc(lc)
    tt0 = MPI_WTIME()
!    call MPI_REDUCE(fm(:,:,:)%sig_a,ac%fm(:,:,:)%sig_a,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fm(:,:,:)%sig_t,ac%fm(:,:,:)%sig_t,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fm(:,:,:)%nusig_f,ac%fm(:,:,:)%nusig_f,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fm(:,:,:)%kappa,ac%fm(:,:,:)%kappa,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fm(:,:,:)%phi,ac%fm(:,:,:)%phi,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fsd_MC,fsd_MC0,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fm(:,:,:)%sig_a,ac%fm(:,:,:)%sig_a,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fm(:,:,:)%sig_t,ac%fm(:,:,:)%sig_t,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fm(:,:,:)%nusig_f,ac%fm(:,:,:)%nusig_f,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fm(:,:,:)%kappa,ac%fm(:,:,:)%kappa,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fm(:,:,:)%phi,ac%fm(:,:,:)%phi,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!!    call MPI_REDUCE(fsd_MC,fsd_MC0,dsize,15,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    do ij = 1, 6
!    call MPI_REDUCE(fm(:,:,:)%J0(ij),ac%fm(:,:,:)%J0(ij),dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    call MPI_REDUCE(fm(:,:,:)%J1(ij),ac%fm(:,:,:)%J1(ij),dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
!    end do

    call MPI_REDUCE(fm(:,:,:), ac%fm(:,:,:), dsize*17, MPI_REAL8, MPI_SUM, score, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(fsd_MC,fsd_MC0,dsize,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)

    tt1 = MPI_WTIME()

    ! Depletion parameter allocation
    if ( DTMCBU .and. cyc ==  n_totcyc ) then
        if ( .not. allocated(buflux) ) allocate(buflux(nfm(1),nfm(2),nfm(3),nrings))
        buflux = 0
        call MPI_REDUCE(tmflux,buflux,dsize*nrings,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        deallocate(tmflux)
    end if

    if ( .not. iscore ) then
        deallocate(fsd_MC0)
        nullify(ac)
        !if ( DTMCBU .and. cyc == n_totcyc ) deallocate(tmflux)
        !if ( DTMCBU .and. cyc == n_totcyc ) deallocate(tmflux,tmkapa)
        return
    end if

    ! -------------------------------------------------------------------------
    ! data transmission II
    fsd_MC = fsd_MC0
    deallocate(fsd_MC0)

    if ( tallyon ) call FMFD_TO_MC(bat,cyc,fm)

    ! current sweeping
    ! 1 (+) 0 (-)
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 )       ac%fm(ii,jj,kk)%J1(1) = ac%fm(ii-1,jj,kk)%J1(2)
        if ( ii /= nfm(1) )  ac%fm(ii,jj,kk)%J0(2) = ac%fm(ii+1,jj,kk)%J0(1)
        if ( jj /= 1 )       ac%fm(ii,jj,kk)%J1(3) = ac%fm(ii,jj-1,kk)%J1(4)
        if ( jj /= nfm(2) )  ac%fm(ii,jj,kk)%J0(4) = ac%fm(ii,jj+1,kk)%J0(3)
        if ( kk /= 1 )       ac%fm(ii,jj,kk)%J1(5) = ac%fm(ii,jj,kk-1)%J1(6)
        if ( kk /= nfm(3) )  ac%fm(ii,jj,kk)%J0(6) = ac%fm(ii,jj,kk+1)%J0(5)
    end do
    end do
    end do

    ! initilization
    fm_avg(:,:,:)%phi     = 0 
    fm_avg(:,:,:)%sig_t   = 0 
    fm_avg(:,:,:)%sig_a   = 0 
    fm_avg(:,:,:)%nusig_f = 0 
    fm_avg(:,:,:)%kappa   = 0 
    do ii = 1, 6 
      fm_avg(:,:,:)%J0(ii) = 0 
      fm_avg(:,:,:)%J1(ii) = 0 
    enddo 


    ! cycle length
    which_idtmc = 1
    !if ( inactive_CMFD ) then
    select case(which_idtmc)
    case(1)
    ! iDTMC1
        if ( cyc > acc_skip ) then
        length = cyc-acc_skip
        else
        length = 1
        end if
        !if ( cyc >= n_inact+1 ) print*, "iDTMC1 : ", length
    case(2)
    ! iDTMC2
    if ( cyc <= n_inact ) then
        if ( cyc > acc_skip ) then
        length = cyc-acc_skip
        else
        length = 1
        end if
    else
        length = n_inact-acc_skip+1
        if ( length < cyc-n_inact ) then
        length = cyc-n_inact
        end if
        if ( cyc >= n_inact+1 ) print*, "iDTMC2 : ", length
    end if
    end select
!    else
!        length = n_acc-acc_skip
!    end if

    ! accumulation
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    do mm = cyc-length+1, cyc
       if ( cyc > n_acc ) then
       nn = mod(mm,n_acc)+1
       else
       nn = mm
       end if
       fm_avg(ii,jj,kk)%phi     = fm_avg(ii,jj,kk)%phi     &
                                + acc(nn)%fm(ii,jj,kk)%phi
       fm_avg(ii,jj,kk)%sig_t   = fm_avg(ii,jj,kk)%sig_t   &
                                + acc(nn)%fm(ii,jj,kk)%sig_t
       fm_avg(ii,jj,kk)%sig_a   = fm_avg(ii,jj,kk)%sig_a   &
                                + acc(nn)%fm(ii,jj,kk)%sig_a
       fm_avg(ii,jj,kk)%nusig_f = fm_avg(ii,jj,kk)%nusig_f &
                                + acc(nn)%fm(ii,jj,kk)%nusig_f
       fm_avg(ii,jj,kk)%kappa   = fm_avg(ii,jj,kk)%kappa &
                                + acc(nn)%fm(ii,jj,kk)%kappa

       fm_avg(ii,jj,kk)%J0(:)   = fm_avg(ii,jj,kk)%J0(:)   &
                                + acc(nn)%fm(ii,jj,kk)%J0(:)
       fm_avg(ii,jj,kk)%J1(:)   = fm_avg(ii,jj,kk)%J1(:)   &
                                + acc(nn)%fm(ii,jj,kk)%J1(:)
    end do
    end do
    end do
    end do

!    ! accumulation
!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!    if ( inactive_CMFD ) then
!    !if ( cyc <= n_inact ) then
!    do mm = acc_skip+1, cyc
!       fm_avg(ii,jj,kk)%phi     = fm_avg(ii,jj,kk)%phi     &
!                                + acc(mm)%fm(ii,jj,kk)%phi
!       fm_avg(ii,jj,kk)%sig_t   = fm_avg(ii,jj,kk)%sig_t   &
!                                + acc(mm)%fm(ii,jj,kk)%sig_t
!       fm_avg(ii,jj,kk)%sig_a   = fm_avg(ii,jj,kk)%sig_a   &
!                                + acc(mm)%fm(ii,jj,kk)%sig_a
!       fm_avg(ii,jj,kk)%nusig_f = fm_avg(ii,jj,kk)%nusig_f &
!                                + acc(mm)%fm(ii,jj,kk)%nusig_f
!
!       fm_avg(ii,jj,kk)%J0(:)   = fm_avg(ii,jj,kk)%J0(:)   &
!                                + acc(mm)%fm(ii,jj,kk)%J0(:)
!       fm_avg(ii,jj,kk)%J1(:)   = fm_avg(ii,jj,kk)%J1(:)   &
!                                + acc(mm)%fm(ii,jj,kk)%J1(:)
!    end do
!    length = cyc-acc_skip
!!    else
!!    do mm = cyc-n_inact+acc_skip, cyc
!!       fm_avg(ii,jj,kk)%phi     = fm_avg(ii,jj,kk)%phi     &
!!                                + acc(mm)%fm(ii,jj,kk)%phi
!!       fm_avg(ii,jj,kk)%sig_t   = fm_avg(ii,jj,kk)%sig_t   &
!!                                + acc(mm)%fm(ii,jj,kk)%sig_t
!!       fm_avg(ii,jj,kk)%sig_a   = fm_avg(ii,jj,kk)%sig_a   &
!!                                + acc(mm)%fm(ii,jj,kk)%sig_a
!!       fm_avg(ii,jj,kk)%nusig_f = fm_avg(ii,jj,kk)%nusig_f &
!!                                + acc(mm)%fm(ii,jj,kk)%nusig_f
!!
!!       fm_avg(ii,jj,kk)%J0(:)   = fm_avg(ii,jj,kk)%J0(:)   &
!!                                + acc(mm)%fm(ii,jj,kk)%J0(:)
!!       fm_avg(ii,jj,kk)%J1(:)   = fm_avg(ii,jj,kk)%J1(:)   &
!!                                + acc(mm)%fm(ii,jj,kk)%J1(:)
!!    end do
!!    length = n_inact-acc_skip+1
!!    end if
!    else
!    do mm = acc_skip+1, n_acc
!       fm_avg(ii,jj,kk)%phi     = fm_avg(ii,jj,kk)%phi     &
!                                + acc(mm)%fm(ii,jj,kk)%phi
!       fm_avg(ii,jj,kk)%sig_t   = fm_avg(ii,jj,kk)%sig_t   &
!                                + acc(mm)%fm(ii,jj,kk)%sig_t
!       fm_avg(ii,jj,kk)%sig_a   = fm_avg(ii,jj,kk)%sig_a   &
!                                + acc(mm)%fm(ii,jj,kk)%sig_a
!       fm_avg(ii,jj,kk)%nusig_f = fm_avg(ii,jj,kk)%nusig_f &
!                                + acc(mm)%fm(ii,jj,kk)%nusig_f
!
!       fm_avg(ii,jj,kk)%J0(:)   = fm_avg(ii,jj,kk)%J0(:)   &
!                                + acc(mm)%fm(ii,jj,kk)%J0(:)
!       fm_avg(ii,jj,kk)%J1(:)   = fm_avg(ii,jj,kk)%J1(:)   &
!                                + acc(mm)%fm(ii,jj,kk)%J1(:)
!    end do
!    length = n_acc-acc_skip
!    end if
!    end do
!    end do
!    end do

    where ( fm_avg(:,:,:)%phi /= 0 )
    fm_avg(:,:,:) % sig_t   = fm_avg(:,:,:) % sig_t   / fm_avg(:,:,:) % phi
    fm_avg(:,:,:) % sig_a   = fm_avg(:,:,:) % sig_a   / fm_avg(:,:,:) % phi
    fm_avg(:,:,:) % nusig_f = fm_avg(:,:,:) % nusig_f / fm_avg(:,:,:) % phi
    fm_avg(:,:,:) % kappa   = fm_avg(:,:,:) % kappa   / fm_avg(:,:,:) % phi
    fm_avg(:,:,:) % phi     = fm_avg(:,:,:) % phi     / (dble(ngen)*v_fm*2D0*length)
    end where

    ! surface quantity normalization
    do ii = 1, 6
    bb = dble(ngen)*a_fm(ii)*length
    fm_avg(:,:,:) % J0(ii)   = fm_avg(:,:,:) % J0(ii) / bb
    fm_avg(:,:,:) % J1(ii)   = fm_avg(:,:,:) % J1(ii) / bb
    end do


    ! zigzag at the corner
    if ( zigzagon ) call SET_ZERO_FLUX(fm_avg(:,:,:)%phi)

!    ! intra pin power
!    if ( DTMCBU .and. cyc == n_totcyc ) then
!        call INTRA_PIN_MC
!        deallocate(tmflux)
!        !deallocate(tmflux,tmkapa)
!    end if

end subroutine 

!! =============================================================================
!! INTRA_PIN_POWER
!! =============================================================================
!subroutine INTRA_PIN_MC
!    implicit none
!    integer:: mm, nn, oo
!    integer:: nring
!
!    ! normalization for the sum to be unity
!    do mm = 1, nfm(1)
!    do nn = 1, nfm(2)
!    do oo = 1, nfm(3)
!        nring = buring(mm,nn,oo)
!        buflux(mm,nn,oo,1:nring) = &
!        buflux(mm,nn,oo,1:nring) / sum(buflux(mm,nn,oo,1:nring))
!!        bukapa(mm,nn,oo,1:nring) = &
!!        bukapa(mm,nn,oo,1:nring) / sum(bukapa(mm,nn,oo,1:nring))
!    end do
!    end do
!    end do
!
!end subroutine

!
!! =============================================================================
!! 
!! =============================================================================
!subroutine BU_TALLY_GROUPING(flux,kappa)
!    implicit none
!    integer:: n_group
!    real(8):: n_type
!    real(8):: flux(:,:,:,:), kappa(:,:,:,:)
!    real(8), allocatable:: group0(:,:), group1(:,:)
!
!    n_group = maxval(butype)
!
!    allocate(group0(n_group,n_rings))
!    allocate(group1(n_group,n_rings))
!    group0 = 0; group1 = 0
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!        group0(butype(ii,jj,kk),:) = group0(butype(ii,jj,kk),:) + &
!            flux(ii,jj,kk,:)
!        group1(butype(ii,jj,kk),:) = group1(butype(ii,jj,kk),:) + &
!            kappa(ii,jj,kk,:)
!    end do
!    end do
!    end do
!
!    do ii = 1, n_group
!    n_type = count(butype==ii)
!    group0(ii,:) = group0(ii,:) / n_type
!    group1(ii,:) = group1(ii,:) / n_type
!    end do
!
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!        flux(ii,jj,kk,:)  = group0(butype(ii,jj,kk),:)
!        kappa(ii,jj,kk,:) = group1(butype(ii,jj,kk),:)
!    end do
!    end do
!    end do
!
!end subroutine


! =============================================================================
! SET_ZERO_FLUX
! =============================================================================
subroutine SET_ZERO_FLUX(mat)
    implicit none
    real(8), intent(out):: mat(:,:,:)   ! for given matrix

    do ii = 1, zz_div
        mat(zzf0(ii)+1:zzf0(ii+1),1:zzf1(ii),:)        = 0D0
        mat(zzf0(ii)+1:zzf0(ii+1),zzf2(ii)+1:nfm(2),:) = 0D0
    end do

end subroutine

! =============================================================================
! FMFD_TO_MC transfers the tally parameters from FMFD to MC
! =============================================================================
subroutine FMFD_TO_MC(bat,cyc,fm)
    use TALLY, only: ttally, MC_tally, n_type
    use ENTROPY, only: fetnusigf
    implicit none
    integer, intent(in):: bat, cyc
    type(FMFD_parameters), intent(in):: fm(:,:,:)
    integer:: acyc

    if ( cyc <= n_inact ) return
    if ( bat == 0 ) return
    acyc = cyc - n_inact

    do ii = 1, n_type
    select case(ttally(ii))
    case(2)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%sig_t
    case(3)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%sig_a
    case(4)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%nusig_f
    case(12)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%sig_t/fm(:,:,:)%phi
    case(13)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%sig_a/fm(:,:,:)%phi
    case(14)
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%nusig_f/fm(:,:,:)%phi
    case default
        MC_tally(bat,acyc,ii,1,:,:,:) = fm(:,:,:)%phi
    end select
    end do

    !fetnusigf(:,:,:) = fm(:,:,:)%nusig_f

end subroutine

! =============================================================================
! FMFD_SOLVE solves FMFD eigenvalue problem
! =============================================================================
subroutine FMFD_SOLVE(bat,cyc)
    use FMFD_HEADER, only: acc
    implicit none
    integer, intent(in)    :: bat           ! batch number
    integer, intent(in)    :: cyc           ! cycle number
    real(8):: phi1(nfm(1),nfm(2),nfm(3))

    call BASE_FMFD_CALCULATION(bat,cyc,phi1)
    if ( dual_fmfd ) call DUAL_FMFD_CALCULATION(bat,cyc,phi1)

end subroutine

! =============================================================================
! BASE_FMFD_CALCULATION
! =============================================================================
subroutine BASE_FMFD_CALCULATION(bat,cyc,phi1)
    use ENTROPY, only: mprupon
    use PERTURBATION, only: perton, PERTURBATION_KEFF
    use CMFD, only: ONE_NODE_ADJOINT
    implicit none
    type(FMFD_accumulation), pointer:: ac
    integer, intent(in):: bat, cyc
    real(8), intent(inout):: phi1(:,:,:)
    real(8):: k_eff
    integer:: acyc
    integer:: lc
    real(8):: aa
    real(8):: tt0, tt1, tt2


    ! parameter initialization
    !tt1 = MPI_WTIME()
    if ( icore == score ) then
    k_eff = keff
    ! copy parameters
    if ( inactive_cmfd .and. cyc <= n_inact ) then
        lc = mod(cyc-1,n_acc)+1
        ac => acc(lc)

        fphi1(:,:,:) = ac%fm(:,:,:)%phi
        if ( zigzagon ) call SET_ZERO_FLUX(fphi1)

        where ( fphi1 /= 0 )
        fm_t(:,:,:)    = ac%fm(:,:,:)%sig_t
        fm_a(:,:,:)    = ac%fm(:,:,:)%sig_a
        fm_nf(:,:,:)   = ac%fm(:,:,:)%nusig_f 
        end where
        do ii = 1, 6
        aa = dble(ngen)*a_fm(ii)
        where ( fphi1 /= 0 ) 
        fmJ0(:,:,:,ii) = ac%fm(:,:,:)%J0(ii)/aa
        fmJ1(:,:,:,ii) = ac%fm(:,:,:)%J1(ii)/aa
        end where
        end do
        nullify(ac)
    
        where ( fphi1 /= 0 ) 
        fm_t   = fm_t  / fphi1
        fm_a   = fm_a  / fphi1
        fm_nf  = fm_nf / fphi1
        fphi1  = fphi1 / (dble(ngen)*v_fm*2D0)
        end where


    else
        fphi1(:,:,:)    = fm_avg(:,:,:)%phi
        where ( fphi1 /= 0 )
        fm_t(:,:,:)     = fm_avg(:,:,:)%sig_t
        fm_a(:,:,:)     = fm_avg(:,:,:)%sig_a
        fm_nf(:,:,:)    = fm_avg(:,:,:)%nusig_f 
        end where
        do ii = 1, 6
        where ( fphi1 /= 0 )
        fmJ0(:,:,:,ii)  = fm_avg(:,:,:)%J0(ii)
        fmJ1(:,:,:,ii)  = fm_avg(:,:,:)%J1(ii)
        end where 
        end do
    end if
    fmJn = fmJ1-fmJ0
    fmD  = 1D0 / (3D0 * fm_t)
    where( fphi1 == 0 ) fmD = 0
    where( fphi1 == 0 .or. fphi1 < 1E-13 ) fm_nf = 0


!    if ( cyc == n_inact+1 ) then
!        call RECALL(1,k_eff)
!        pause
!        stop
!    end if

    end if
    !tt2 = MPI_WTIME()
    !if ( iscore ) print*, " - FMFD parameter reading : ", tt2-tt1

    if ( inactive_cmfd ) then
        if ( curr_cyc <= n_inact ) then
        call CMFD_CALCULATION(k_eff)
        else
        call ONE_NODE_CMFD(k_eff,1D-9)
        end if
    elseif ( cmfdon ) then
        if ( icore == score ) then
        if ( .not. allocated(fmF) ) allocate(fmF(nfm(1),nfm(2),nfm(3),6))
        fmF = 2D0*(fmJ1+fmJ0)
        end if
        call ONE_NODE_CMFD(k_eff,1D-9)
    else
        if ( iscore ) then
        call D_TILDA_CALCULATION
        if ( .not. pfmfdon ) then
            call D_HAT_CALCULATION
            call FMFD_MATRIX
        else
            call D_PHAT_CALCULATION
            call PFMFD_MATRIX
            write(*,*) 'PFMFD', sum(fmDh)
        endif 
        if ( zigzagon ) call SET_ZERO_M
        call POWER(k_eff)
        end if
    end if

    ! weight update
    call WEIGHT_UPDATE(bat,cyc,k_eff)
    if(iscore) print *, 'keff_nopert', k_eff
   
    ! error quantification by 1st order perturbation
    tt1 = MPI_WTIME()
    if ( perton )  call PERTURBATION_KEFF(bat,k_eff,cyc)
    

    tt2 = MPI_WTIME()
    if ( iscore ) print*, " - perturbation total : ", tt2-tt1

    if ( icore /= score ) return
    print*, "keff ", k_eff
    !if (perton) write(*,*) "COSAMPLING", AVG(k_real(bat,cyc,:))

    !> CMFD feedback (modulation)
    acyc = cyc - n_inact
    !if ( DTMCBU .and. acyc > 0 ) call INTRA_PIN_DTMC(acyc)
    if ( DTMCBU .and. acyc > 0 ) &
        p_dep_dt(acyc,:,:,:) = fm_avg(:,:,:)%kappa*fphi1(:,:,:)
    if ( dual_fmfd ) phi1 = fphi1
    if ( mprupon ) k_mprup = k_eff
    if ( bat > 0 ) k_fmfd(bat,cyc) = k_eff
    if ( bat /= 0 .and. acyc > 0 ) &
        p_fmfd(bat,acyc,:,:,:) = fm_nf(:,:,:)*fphi1(:,:,:)

end subroutine

! =============================================================================
! 
! =============================================================================
subroutine INTRA_PIN_DTMC
    use constants
    use MATERIAL_HEADER, only: n_materials, Material_CE, materials
    implicit none
    real(8):: factor
    ! -----
    integer:: nz, nring
    integer:: id0(4), id1(4)
    real(8):: smflux
    real(8), allocatable:: tmring(:,:,:,:)
    type(Material_CE), pointer:: mat


    ! -------------------------------------------------------------------------
    ! Normalized flux calculation
    ! normalization factor
    ! kappa (MeV)
    ! power (MW)
    factor = Nominal_Power / &
        sum(fm_avg(:,:,:)%kappa*fphi1(:,:,:)*eVtoJoule)

    ! DTMC volume-weighted flux for normalized for real power
    fphi1(:,:,:) = factor*fphi1(:,:,:)
    print *, 'FPHI', factor, Nominal_Power


    ! -------------------------------------------------------------------------
    ! Real flux in materials
    ! intra pin flux distribution for reconstruction
    allocate(tmring(nfm(1),nfm(2),nfm(3),nrings))
    tmring = 0
    do ii = 1, n_materials
        if ( .not. materials(ii)%depletable ) cycle
        mat => materials(ii)
        nz = size(mat%dtmc) / 4

        id0(1:4) = mat%dtmc(1:4)
        do jj = 1, nz
            id1(1:4) = mat%dtmc(4*(jj-1)+1:4*jj)
            nring = buring(id1(1),id1(2),id1(3))
            do kk = 1, nring
            tmring(id0(1),id0(2),id0(3),kk) = tmring(id0(1),id0(2),id0(3),kk) &
                + buflux(id1(1),id1(2),id1(3),kk)
            end do
        end do
    end do

    ! normalization of intra pin flux distribution
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        nring = buring(ii,jj,kk)
        if ( tmring(ii,jj,kk,1) == 0 ) cycle
        tmring(ii,jj,kk,1:nring) = &
        tmring(ii,jj,kk,1:nring) / sum(tmring(ii,jj,kk,1:nring))
    end do
    end do
    end do

    ! material-wise reconstructed flux
    do ii = 1, n_materials
        if ( .not. materials(ii)%depletable ) cycle
        mat => materials(ii)
        nz = size(mat%dtmc) / 4
        smflux = 0

        id0(1:4) = mat%dtmc(1:4)
        do jj = 1, nz
            id1(1:4) = mat%dtmc(4*(jj-1)+1:4*jj)
            smflux = smflux + fphi1(id1(1),id1(2),id1(3))
        end do
        mat%flux = smflux/dble(mat%vol)*tmring(id0(1),id0(2),id0(3),id0(4))
        print *, 'MAT', mat%mat_name, mat%flux
    end do
    if ( allocated(tmring) ) deallocate(tmring)


!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!        buflux(ii,jj,kk,:) = buflux(ii,jj,kk,:) * fphi1(ii,jj,kk)
!!        bukapa(ii,jj,kk,:) = bukapa(ii,jj,kk,:) &
!!            * fphi1(ii,jj,kk) * fm_avg(ii,jj,kk)%kappa * eVtoJoule
!    end do
!    end do
!    end do
!    end if

!    allocate(pavg(nfm(1),nfm(2),nfm(3),n_rings))
!
!    ! DTMC flux
!    f_dep_dt(acyc,:,:,:,1) = factor*fphi1(:,:,:)
!
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!    do rr = 1, n_rings
!    pavg(ii,jj,kk,rr) = sum(f_dep_mc(acc_skip+1:cyc,ii,jj,kk,rr))
!    end do
!    f_dep_dt(acyc,ii,jj,kk,:) = pavg(ii,jj,kk,:) &
!        /sum(pavg(ii,jj,kk,:))*n_rings*f_dep_dt(acyc,ii,jj,kk,1)
!    end do
!    end do
!    end do
!
!    ! DTMC power
!    p_dep_dt(acyc,:,:,:,1) = factor*fm_avg(:,:,:)%kappa*fphi1(:,:,:)
!
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!    do rr = 1, n_rings
!    pavg(ii,jj,kk,rr) = sum(p_dep_mc(acc_skip+1:cyc,ii,jj,kk,rr))
!    end do
!    p_dep_dt(acyc,ii,jj,kk,:) = pavg(ii,jj,kk,:)*v_fm &
!        /sum(pavg(ii,jj,kk,:))*p_dep_dt(acyc,ii,jj,kk,1)
!    end do
!    end do
!    end do





!    ! DTMC flux
!    f_dep_dt(acyc,:,:,:,1) = factor/1.602E-19*fphi1(:,:,:)
!
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!    do rr = 1, n_rings
!    pavg(ii,jj,kk,rr) = sum(f_dep_mc(n_inact+1:cyc,ii,jj,kk,rr))
!    end do
!    f_dep_dt(acyc,ii,jj,kk,:) = pavg(ii,jj,kk,:)*v_ring(:) &
!        /sum(pavg(ii,jj,kk,:))*n_rings*f_dep_dt(acyc,ii,jj,kk,1)
!    end do
!    end do
!    end do
!
!    ! DTMC power
!    p_dep_dt(acyc,:,:,:,1) = factor*fm_avg(:,:,:)%kappa*fphi1(:,:,:)
!
!    do kk = 1, nfm(3)
!    do jj = 1, nfm(2)
!    do ii = 1, nfm(1)
!    do rr = 1, n_rings
!    pavg(ii,jj,kk,rr) = sum(p_dep_mc(n_inact+1:cyc,ii,jj,kk,rr))
!    end do
!    p_dep_dt(acyc,ii,jj,kk,:) = pavg(ii,jj,kk,:)*v_fm &
!        /sum(pavg(ii,jj,kk,:))*p_dep_dt(acyc,ii,jj,kk,1)
!    end do
!    end do
!    end do

end subroutine


! =============================================================================
! 
! =============================================================================
subroutine RECALL(opt,k_eff)
    implicit none
    integer, intent(in):: opt
    real(8):: k_eff

    
    open(21,file='FMFD_PARAMETERS.out')

    select case(opt)
    case(1) ! ===== saving =====
    write(21,1), fphi1
    write(21,1), fm_t
    write(21,1), fm_a
    write(21,1), fm_nf
    write(21,1), fmJ0
    write(21,1), fmJ1
    write(21,*), k_eff

    case(2) ! ===== loading =====
    read(21,1), fphi1
    read(21,1), fm_t
    read(21,1), fm_a
    read(21,1), fm_nf
    read(21,1), fmJ0
    read(21,1), fmJ1
    read(21,*), k_eff
    fmJn = fmJ1-fmJ0
    fmD = 1D0 / (3D0*fm_t)
    where ( fphi1 == 0 ) fmD = 0
    where ( fphi1 == 0 .or. fphi1 < 1E-12 ) fm_nf = 0

    end select

    1 format(<nfm(1)>ES17.9)
    close(21)

end subroutine

! =============================================================================
! DUAL_FMFD_CALCULATION
! =============================================================================
subroutine DUAL_FMFD_CALCULATION(bat,cyc,phi2)
    use ENTROPY, only: mprupon
    implicit none
    type(FMFD_accumulation), pointer:: ac
    integer, intent(in):: bat, cyc
    real(8), intent(inout):: phi2(:,:,:)
    real(8):: k_eff
    integer:: acyc
    integer:: lc
    real(8):: aa

    ! parameter initialization
    if ( icore == score ) then
    ! copy parameters
    k_eff = keff
    lc = mod(cyc,n_acc)+1
    ac => acc(lc)
    fm_t(:,:,:)    = ac%fm(:,:,:)%sig_t
    fm_a(:,:,:)    = ac%fm(:,:,:)%sig_a
    fm_nf(:,:,:)   = ac%fm(:,:,:)%nusig_f 
    fphi1(:,:,:)   = ac%fm(:,:,:)%phi
    do ii = 1, 6
    fmJ0(:,:,:,ii) = ac%fm(:,:,:)%J0(ii)
    fmJ1(:,:,:,ii) = ac%fm(:,:,:)%J1(ii)
    end do
    nullify(ac)

    where ( fphi1 /= 0 ) 
    fm_t   = fm_t / fphi1
    fm_a   = fm_a / fphi1
    fm_nf  = fm_nf / fphi1
    fphi1  = fphi1 / (dble(ngen)*v_fm*2D0)
    end where

    do ii = 1, 6
    aa = dble(ngen)*a_fm(ii)
    fmJ0(:,:,:,ii) = fmJ0(:,:,:,ii) / aa
    fmJ1(:,:,:,ii) = fmJ1(:,:,:,ii) / aa
    end do

    if ( zigzagon ) call SET_ZERO_FLUX(fphi1)
    fmD(:,:,:) = 1D0 / (3D0 * fm_t(:,:,:))
    where( fphi1 == 0 ) fmD = 0
    where( fphi1 == 0 .or. fphi1 < 1E-13 ) fm_nf = 0
    fmJn = fmJ1-fmJ0
    end if

    if ( inactive_cmfd ) then
        if ( cyc > n_inact ) call ONE_NODE_CMFD(k_eff,1D-9,fphi1)
    elseif ( cmfdon ) then
        if ( icore == score ) fmF = 2D0*(fmJ1+fmJ0)
        call ONE_NODE_CMFD(k_eff,1D-9,fphi1)
    else
        if ( icore == score ) then
        call D_TILDA_CALCULATION
        if ( .not. pfmfdon ) then
            call D_HAT_CALCULATION
            call FMFD_MATRIX
        else
            call D_PHAT_CALCULATION
            call PFMFD_MATRIX
        endif 

        if ( zigzagon ) call SET_ZERO_M
        call POWER(k_eff,fphi1)
        end if
    end if

    ! weight update
    call WEIGHT_UPDATE(bat,cyc,k_eff)

    if ( icore /= score ) return

    !> CMFD feedback (modulation)
    print*, "keff2", k_eff
    acyc = cyc - n_inact
    if ( mprupon ) k_mprup = k_eff
    if ( bat > 0 ) k_fmfd2(bat,cyc) = k_eff
    if ( bat /= 0 .and. acyc > 0 ) &
        p_fmfd2(bat,acyc,:,:,:) = fm_nf(:,:,:)*fphi1(:,:,:)
    
end subroutine

! =============================================================================
! D_TILDA_CALCULATION
! =============================================================================
!subroutine D_TILDA_CALCULATION
!    implicit none
!    
!    fmDt = 0
!
!    ! inner region
!    do ii = 1, nfm(1)
!    do jj = 1, nfm(2)
!    do kk = 1, nfm(3)
!        if ( ii /= 1 )      fmDt(ii,jj,kk,1) = 2D0*fmD(ii,jj,kk) &
!            *fmD(ii-1,jj,kk)/((fmD(ii,jj,kk)+fmD(ii-1,jj,kk))*dfm(1))
!        if ( ii /= nfm(1) ) fmDt(ii,jj,kk,2) = 2D0*fmD(ii+1,jj,kk) &
!            *fmD(ii,jj,kk)/((fmD(ii+1,jj,kk)+fmD(ii,jj,kk))*dfm(1))
!        if ( jj /= 1 )      fmDt(ii,jj,kk,3) = 2D0*fmD(ii,jj,kk) &
!            *fmD(ii,jj-1,kk)/((fmD(ii,jj,kk)+fmD(ii,jj-1,kk))*dfm(2))
!        if ( jj /= nfm(2) ) fmDt(ii,jj,kk,4) = 2D0*fmD(ii,jj+1,kk) &
!            *fmD(ii,jj,kk)/((fmD(ii,jj+1,kk)+fmD(ii,jj,kk))*dfm(2))
!        if ( kk /= 1 )      fmDt(ii,jj,kk,5) = 2D0*fmD(ii,jj,kk) &
!            *fmD(ii,jj,kk-1)/((fmD(ii,jj,kk)+fmD(ii,jj,kk-1))*dfm(3))
!        if ( kk /= nfm(3) ) fmDt(ii,jj,kk,6) = 2D0*fmD(ii,jj,kk+1) &
!            *fmD(ii,jj,kk)/((fmD(ii,jj,kk+1)+fmD(ii,jj,kk))*dfm(3))
!    end do
!    end do
!    end do
!
!    where ( isnan(fmDt) ) fmDt = 0
!
!end subroutine
!
!
!! =============================================================================
!! D_HAT_CALCULATION calculates correction factors
!! =============================================================================
!subroutine D_HAT_CALCULATION
!    integer :: i, j, k
!
!    ! inner region
!    do i = 1, nfm(1)
!    do j = 1, nfm(2)
!    do k = 1, nfm(3)
!        if ( i /= 1 ) &      ! x0
!        fmDh(i,j,k,1) = (fmJn(i,j,k,1)+fmDt(i,j,k,1) &
!            *(fphi1(i,j,k)-fphi1(i-1,j,k)))/(fphi1(i,j,k)+fphi1(i-1,j,k))
!        if ( i /= nfm(1) ) & ! x1
!        fmDh(i,j,k,2) = (fmJn(i,j,k,2)+fmDt(i,j,k,2) &
!            *(fphi1(i+1,j,k)-fphi1(i,j,k)))/(fphi1(i+1,j,k)+fphi1(i,j,k))
!        if ( j /= 1 ) &      ! y0
!        fmDh(i,j,k,3) = (fmJn(i,j,k,3)+fmDt(i,j,k,3) &
!            *(fphi1(i,j,k)-fphi1(i,j-1,k)))/(fphi1(i,j,k)+fphi1(i,j-1,k))
!        if ( j /= nfm(2) ) & ! y1
!        fmDh(i,j,k,4) = (fmJn(i,j,k,4)+fmDt(i,j,k,4) &
!            *(fphi1(i,j+1,k)-fphi1(i,j,k)))/(fphi1(i,j+1,k)+fphi1(i,j,k))
!        if ( k /= 1 ) &      ! y0
!        fmDh(i,j,k,5) = (fmJn(i,j,k,5)+fmDt(i,j,k,5) &
!            *(fphi1(i,j,k)-fphi1(i,j,k-1)))/(fphi1(i,j,k)+fphi1(i,j,k-1))
!        if ( k /= nfm(3) ) & ! y1
!        fmDh(i,j,k,6) = (fmJn(i,j,k,6)+fmDt(i,j,k,6) &
!            *(fphi1(i,j,k+1)-fphi1(i,j,k)))/(fphi1(i,j,k+1)+fphi1(i,j,k))
!    end do
!    end do
!    end do
!
!    ! boundary
!    i = 1;      fmDh(i,:,:,1) = fmJn(i,:,:,1)/fphi1(i,:,:)
!    i = nfm(1); fmDh(i,:,:,2) = fmJn(i,:,:,2)/fphi1(i,:,:)
!    j = 1;      fmDh(:,j,:,3) = fmJn(:,j,:,3)/fphi1(:,j,:)
!    j = nfm(2); fmDh(:,j,:,4) = fmJn(:,j,:,4)/fphi1(:,j,:)
!    k = 1;      fmDh(:,:,k,5) = fmJn(:,:,k,5)/fphi1(:,:,k)
!    k = nfm(3); fmDh(:,:,k,6) = fmJn(:,:,k,6)/fphi1(:,:,k)
!
!end subroutine
!
!
!! =============================================================================
!! D_HAT_CALCULATION calculates correction factors
!! =============================================================================
!subroutine D_PHAT_CALCULATION
!    integer :: i, j, k
!
!    ! inner region (outgoing direction)
!    do i = 1, nfm(1)
!    do j = 1, nfm(2)
!    do k = 1, nfm(3)
!        if ( i /= 1 ) &      ! x0
!        fmDh(i,j,k,1) = (fmJ0(i,j,k,1)-5D-1*fmDt(i,j,k,1) &
!                    *(fphi1(i,j,k)-fphi1(i-1,j,k)))/fphi1(i,j,k)
!        if ( i /= nfm(1) ) & ! x1
!        fmDh(i,j,k,2) = (fmJ1(i,j,k,2)+5D-1*fmDt(i,j,k,2) &
!                    *(fphi1(i+1,j,k)-fphi1(i,j,k)))/fphi1(i,j,k)
!        if ( j /= 1 ) &      ! y0
!        fmDh(i,j,k,3) = (fmJ0(i,j,k,3)-5D-1*fmDt(i,j,k,3) &
!                    *(fphi1(i,j,k)-fphi1(i,j-1,k)))/fphi1(i,j,k)
!        if ( j /= nfm(2) ) & ! y1
!        fmDh(i,j,k,4) = (fmJ1(i,j,k,4)+5D-1*fmDt(i,j,k,4) &
!                    *(fphi1(i,j+1,k)-fphi1(i,j,k)))/fphi1(i,j,k)
!        if ( k /= 1 ) &      ! y0
!        fmDh(i,j,k,5) = (fmJ0(i,j,k,5)-5D-1*fmDt(i,j,k,5) &
!                    *(fphi1(i,j,k)-fphi1(i,j,k-1)))/fphi1(i,j,k)
!        if ( k /= nfm(3) ) & ! y1
!        fmDh(i,j,k,6) = (fmJ1(i,j,k,6)+5D-1*fmDt(i,j,k,6) &
!                    *(fphi1(i,j,k+1)-fphi1(i,j,k)))/fphi1(i,j,k)
!    end do
!    end do
!    end do
!
!    ! boundary
!    if ( .not. zigzagon ) then
!    i = 1;      fmDh(i,:,:,1) = -fmJn(i,:,:,1)/fphi1(i,:,:)
!    i = nfm(1); fmDh(i,:,:,2) = +fmJn(i,:,:,2)/fphi1(i,:,:)
!    j = 1;      fmDh(:,j,:,3) = -fmJn(:,j,:,3)/fphi1(:,j,:)
!    j = nfm(2); fmDh(:,j,:,4) = +fmJn(:,j,:,4)/fphi1(:,j,:)
!    else
!    do i = 1, zz_div
!        fmDh(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) = &
!            -fmJn(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) &
!            /fphi1(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:)
!        fmDh(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) = &
!            +fmJn(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) &
!            /fphi1(zzf2(i),zzf0(i)+1:zzf0(i+1),:)
!        fmDh(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) = &
!            -fmJn(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) &
!            /fphi1(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:)
!        fmDh(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) = &
!            +fmJn(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) &
!            /fphi1(zzf0(i)+1:zzf0(i+1),zzf2(i),:)
!    end do
!    end if
!    k = 1;      fmDh(:,:,k,5) = -fmJn(:,:,k,5)/fphi1(:,:,k)
!    k = nfm(3); fmDh(:,:,k,6) = +fmJn(:,:,k,6)/fphi1(:,:,k)
!
!end subroutine


! ========================================================= !
!  Miscel. functions for CMFD_solve
! ========================================================= !

! function which produces M matrix from xs & D_hat
subroutine FMFD_MATRIX
    integer :: i, j, k
    
    ! initialization
    Mfm = 0
    
    ! Mfm matrix set 
    do k = 1, nfm(3)
    do j = 1, nfm(2)
    do i = 1, nfm(1)
        if ( i /= 1 ) &         ! x0
            Mfm(i,j,k,3) = -(fmDt(i,j,k,1)+fmDh(i,j,k,1))/dfm(1)
        if ( i /= nfm(1) ) &    ! x1
            Mfm(i,j,k,5) = -(fmDt(i,j,k,2)-fmDh(i,j,k,2))/dfm(1)
        if ( j /= 1 ) &         ! y0
            Mfm(i,j,k,2) = -(fmDt(i,j,k,3)+fmDh(i,j,k,3))/dfm(2)
        if ( j /= nfm(2) ) &    ! y1
            Mfm(i,j,k,6) = -(fmDt(i,j,k,4)-fmDh(i,j,k,4))/dfm(2)
        if ( k /= 1 ) &         ! z0
            Mfm(i,j,k,1) = -(fmDt(i,j,k,5)+fmDh(i,j,k,5))/dfm(3)
        if ( k /= nfm(3) ) &    ! z1
            Mfm(i,j,k,7) = -(fmDt(i,j,k,6)-fmDh(i,j,k,6))/dfm(3)
        
        Mfm(i,j,k,4) = &
            +(fmDt(i,j,k,1)-fmDh(i,j,k,1)+fmDt(i,j,k,2)+fmDh(i,j,k,2))/dfm(1) &
            +(fmDt(i,j,k,3)-fmDh(i,j,k,3)+fmDt(i,j,k,4)+fmDh(i,j,k,4))/dfm(2) &
            +(fmDt(i,j,k,5)-fmDh(i,j,k,5)+fmDt(i,j,k,6)+fmDh(i,j,k,6))/dfm(3) &
            +fm_a(i,j,k)

    enddo
    enddo
    enddo

end subroutine FMFD_MATRIX


! function which produces M matrix from xs & D_hat
!subroutine PFMFD_MATRIX
!    integer :: i, j, k
!    
!    ! initialization
!    Mfm = 0
!    
!    ! Mfm matrix set 
!    do k = 1, nfm(3)
!    do j = 1, nfm(2)
!    do i = 1, nfm(1)
!        if ( i /= 1 ) &         ! x0
!            Mfm(i,j,k,3) = -(fmDt(i,j,k,1)+fmDh(i-1,j,k,2))/dfm(1)
!        if ( i /= nfm(1) ) &    ! x1
!            Mfm(i,j,k,5) = -(fmDt(i,j,k,2)+fmDh(i+1,j,k,1))/dfm(1)
!        if ( j /= 1 ) &         ! y0
!            Mfm(i,j,k,2) = -(fmDt(i,j,k,3)+fmDh(i,j-1,k,4))/dfm(2)
!        if ( j /= nfm(2) ) &    ! y1
!            Mfm(i,j,k,6) = -(fmDt(i,j,k,4)+fmDh(i,j+1,k,3))/dfm(2)
!        if ( k /= 1 ) &         ! z0
!            Mfm(i,j,k,1) = -(fmDt(i,j,k,5)+fmDh(i,j,k-1,6))/dfm(3)
!        if ( k /= nfm(3) ) &    ! z1
!            Mfm(i,j,k,7) = -(fmDt(i,j,k,6)+fmDh(i,j,k+1,5))/dfm(3)
!        
!        Mfm(i,j,k,4) = &
!            +(fmDt(i,j,k,1)+fmDh(i,j,k,1)+fmDt(i,j,k,2)+fmDh(i,j,k,2))/dfm(1) &
!            +(fmDt(i,j,k,3)+fmDh(i,j,k,3)+fmDt(i,j,k,4)+fmDh(i,j,k,4))/dfm(2) &
!            +(fmDt(i,j,k,5)+fmDh(i,j,k,5)+fmDt(i,j,k,6)+fmDh(i,j,k,6))/dfm(3) &
!            +fm_a(i,j,k)
!
!    enddo
!    enddo
!    enddo
!
!end subroutine PFMFD_MATRIX
!
!subroutine SET_ZERO_M
!    implicit none
!    integer:: ij
!
!    where( fphi1 == 0 ) 
!    Mfm(:,:,:,1) = 0
!    Mfm(:,:,:,2) = 0
!    Mfm(:,:,:,3) = 0
!    Mfm(:,:,:,4) = 0
!    Mfm(:,:,:,5) = 0
!    Mfm(:,:,:,6) = 0
!    Mfm(:,:,:,7) = 0
!    end where
!
!    do ij = 1, zz_div
!        Mfm(zzf1(ij)+1,zzf0(ij)+1:zzf0(ij+1),:,3) = 0
!        Mfm(zzf2(ij),zzf0(ij)+1:zzf0(ij+1),:,5) = 0
!        Mfm(zzf0(ij)+1:zzf0(ij+1),zzf1(ij)+1,:,2) = 0
!        Mfm(zzf0(ij)+1:zzf0(ij+1),zzf2(ij),:,6) = 0
!    end do
!
!end subroutine
!
!
subroutine POWER (k_eff,phi2)
    use SOLVERS, only: BiCGSTAB_PRE, BiCGSTAB_ILU
    use PRECONDITIONER, only: FINE_ILU_INITIAL, FINE_ILU
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in), optional:: phi2(:,:,:)
    real(8), parameter:: ONE = 1D0
    real(8):: phi0(nfm(1),nfm(2),nfm(3))
    real(8):: sorphi(nfm(1),nfm(2),nfm(3))
    integer:: iter, iter_max = 1D5
    real(8):: err
    real(8):: tt0, tt1, idx
    real(8):: kpre
    real(8):: ww = 1.7

    if ( .not. allocated(fn) ) call FINE_ILU_INITIAL
    call FINE_ILU

    err = ONE
    iter = 1
!    tt0 = MPI_WTIME()
!    idx = 1
    if ( dual_fmfd .and. present(phi2) ) fphi1 = phi2
    do while ( ( err > 1D-9 ) .and. ( iter < iter_max ) )
        !call CPU_TIME(tt0)
        iter = iter + 1
        !sorphi = ww*fphi1+(1D0-ww)*phi0
        phi0 = fphi1
        kpre = k_eff
        fm_s = fm_nf(:,:,:)*fphi1(:,:,:)/k_eff
        !fm_s = fm_nf(:,:,:)*sorphi(:,:,:)/k_eff
        !fphi1 = BiCGSTAB_PRE(Mfm(:,:,:,:),fm_s(:,:,:))
        fphi1 = BiCGSTAB_ILU(mvec3(:,:),fm_s(:,:,:))
        k_eff = k_eff*sum(fm_nf*fphi1*fm_nf*fphi1) &
               / sum(fm_nf*fphi1*fm_nf*phi0)
        !err = abs(k_eff-kpre)/k_eff
!        if ( err < 1D-3 ) then
!            tt1 = MPI_WTIME()
!            print*, tt1-tt0
!            exit
!        end if
        err = sum(abs(fphi1-phi0))/sum(fphi1)
        !print*, iter, k_eff, err
        !call CPU_TIME(tt1)
    enddo
!    tt1 = MPI_WTIME()
!    print*, "time", tt1-tt0

    
end subroutine 


! =============================================================================
! 
! =============================================================================
subroutine WEIGHT_UPDATE(bat,cyc,k_eff,phi2)
    use BANK_HEADER, only: fission_bank, bank_size
    use TALLY, only: FM_ID, CM_ID
    implicit none
    integer, intent(in):: bat, cyc
    real(8), intent(in):: k_eff
    real(8), intent(in), optional:: phi2(:,:,:)
    real(8), allocatable:: fsd3(:,:,:), fsd_MC3(:,:,:), &
                           fsd_FM3(:,:,:) ! for inactive CMFD
    integer:: ix0, ix1, iy0, iy1, iz0, iz1
    logical:: update

    if ( inactive_cmfd .and. .not. allocated(fsd3) ) then
        allocate(fsd3(ncm(1),ncm(2),ncm(3)))
        allocate(fsd_MC3(ncm(1),ncm(2),ncm(3)))
        allocate(fsd_FM3(ncm(1),ncm(2),ncm(3)))
    end if


    if ( icore == score ) then
    update = .true.
    if ( isnan(k_eff) .or. ( k_eff < 1D-2 .or. k_eff > 2D0 ) ) update = .false.
    if ( .not. fmfd2mc .and. n_batch == 1 .and. cyc > n_inact ) update = .false.
    if ( .not. fmfd2mc .and. n_batch > 1 .and. bat >= 1 ) update = .false.

    if ( update ) then
    if ( inactive_CMFD ) then
        do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1
            fsd_MC3(ii,jj,kk) = sum(fsd_MC(ix0:ix1,iy0:iy1,iz0:iz1))
        end do
        end do
        end do

        fsd_FM3 = cm_nf*cphi1
        fsd_MC3 = fsd_MC3 / sum(fsd_MC3)
        fsd_FM3 = fsd_FM3 / sum(fsd_FM3)
        fsd3 = fsd_FM3 / fsd_MC3
    else
        fsd_FM = fm_nf*fphi1
        fsd_MC = fsd_MC / sum(fsd_MC)
        fsd_FM = fsd_FM / sum(fsd_FM)
        fsd = fsd_FM / fsd_MC
    if ( dual_fmfd .and. present(phi2) ) then
        fsd_FM = fm_nf*phi2
        fsd_FM = fsd_FM / sum(fsd_FM)
        fsd2 = fsd_FM / fsd_MC
    end if
    end if
    end if
    end if
    call MPI_BCAST(update,1,MPI_LOGICAL,score,MPI_COMM_WORLD,ierr)

    if ( update ) then
    if ( inactive_cmfd .and. cyc <= n_inact ) then
    call MPI_BCAST(fsd3,ncm(1)*ncm(2)*ncm(3),MPI_REAL8,score,MPI_COMM_WORLD,ierr)
    do ii = 1, bank_size
        id = CM_ID(fission_bank(ii)%xyz)
        if ( id(1) < 1 .or. id(1) > ncm(1) ) cycle
        if ( id(2) < 1 .or. id(2) > ncm(2) ) cycle
        if ( id(3) < 1 .or. id(3) > ncm(3) ) cycle
        fission_bank(ii)%wgt = fission_bank(ii)%wgt * fsd3(id(1),id(2),id(3))
    enddo
    else
    call MPI_BCAST(fsd,n_nodes,MPI_REAL8,score,MPI_COMM_WORLD,ierr)
    if ( dual_fmfd ) then
    call MPI_BCAST(fsd2,n_nodes,MPI_REAL8,score,MPI_COMM_WORLD,ierr)
    allocate(fwgt(bank_size))
    fwgt(:) = fission_bank(:)%wgt
    do ii = 1, bank_size
        id = FM_ID(fission_bank(ii)%xyz)
        if ( id(1) < 1 .or. id(1) > nfm(1) ) cycle
        if ( id(2) < 1 .or. id(2) > nfm(2) ) cycle
        if ( id(3) < 1 .or. id(3) > nfm(3) ) cycle
        fission_bank(ii)%wgt = fission_bank(ii)%wgt * fsd2(id(1),id(2),id(3))
        if ( dual_fmfd ) &
        fwgt(ii) = fwgt(ii) * fsd(id(1),id(2),id(3))
    enddo
    return
    end if

    do ii = 1, bank_size
        id = FM_ID(fission_bank(ii)%xyz)
        if ( id(1) < 1 .or. id(1) > nfm(1) ) cycle
        if ( id(2) < 1 .or. id(2) > nfm(2) ) cycle
        if ( id(3) < 1 .or. id(3) > nfm(3) ) cycle
        fission_bank(ii)%wgt = fission_bank(ii)%wgt * fsd(id(1),id(2),id(3))
    enddo
    end if
    end if

end subroutine

!subroutine FMFD_ENTRP(fsd,cyc)
!    use ENTROPY, only: entrp3
!    implicit none
!    real(8), intent(in) :: fsd(:,:,:)
!    integer, intent(in) :: cyc
!    real(8):: ee0(nfm(1),nfm(2),nfm(3))
!
!
!    ee0 = fsd/sum(fsd)
!    where ( ee0 /= 0 ) ee0 = -ee0*log(ee0)/log(2D0)
!    entrp3(cyc) = sum(ee0)
!
!end subroutine


! =============================================================================
! DET_POWER
! =============================================================================
subroutine DET_POWER(pd)
    use TH_HEADER, only: nth, pp
    implicit none
    real(8), intent(in):: pd(:,:,:) ! power distribution
    integer:: ds(1:3)

    ds(:) = nfm(:)/nth(:)
    
    do ii = 1, nth(1)
    do jj = 1, nth(2)
    do kk = 1, nth(3)
        pp(ii,jj,kk) = sum(pd((ii-1)*ds(1)+1:ii*ds(1), &
                              (jj-1)*ds(2)+1:jj*ds(2), &
                              (kk-1)*ds(3)+1:kk*ds(3)))
    end do
    end do
    end do

end subroutine


! =============================================================================
! 
! =============================================================================
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

end module
