module FMFD
    use omp_lib 
    use mpi 
    use geometry_header,    only : lattices, find_lat_idx, lattice_coord
    use variables
    use constants,          only : INFINITY, TiNY_BIT, prt_keff
    use particle_header,    only : particle
    use FMFD_HEADER
    
    implicit none
    
    contains 


! =============================================================================
! FMFD_initialize_1
! =============================================================================
subroutine FMFD_ALLOCATION()
    use ENTROPY, only: entrp3
    !use CMFD, only: OUT_OF_ZZ
    use PCMFD, only: OUT_OF_ZZ
    implicit none

    ! parameters allocation
    if ( n_batch == 1 ) then
    allocate(k_fmfd(n_batch,n_totcyc))
    allocate(p_fmfd(n_batch,n_act,nfm(1),nfm(2),nfm(3)))
    else
    allocate(k_fmfd(0:n_batch,n_totcyc))
    allocate(p_fmfd(0:n_batch,n_act,nfm(1),nfm(2),nfm(3)))
    end if
    allocate(e_fmfd(nfm(1),nfm(2),nfm(3)))
    allocate(fm(nfm(1),nfm(2),nfm(3)))
    allocate(fm_avg(nfm(1),nfm(2),nfm(3)))
    allocate(fsd_MC(nfm(1),nfm(2),nfm(3)))
    allocate(fsd_FM(nfm(1),nfm(2),nfm(3)))
    allocate(fsd(nfm(1),nfm(2),nfm(3)))
    allocate(acc(n_acc)) 
    do ii = 1, n_acc 
        allocate(acc(ii)%fm(nfm(1),nfm(2),nfm(3)))
        acc(ii)%fm(:,:,:)%phi     = 0
        acc(ii)%fm(:,:,:)%sig_t   = 0
        acc(ii)%fm(:,:,:)%sig_a   = 0
        acc(ii)%fm(:,:,:)%nusig_f = 0
        do jj = 1, 6
        acc(ii)%fm(:,:,:)%Jn(jj)   = 0
        acc(ii)%fm(:,:,:)%J0(jj)   = 0
        acc(ii)%fm(:,:,:)%J1(jj)   = 0
        end do
    enddo
    allocate(entrp3(n_totcyc))
    entrp3 = 0
    
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
                 cm_phi0(ncm(1),ncm(2),ncm(3)), &
                 cm_phi1(ncm(1),ncm(2),ncm(3)), &
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
                 Mcm    (ncm(1),ncm(2),ncm(3),7), &
                 Mfm    (nfm(1),nfm(2),nfm(3),7))
        cm_t    = 0
        cmD     = 0
        cm_a    = 0
        cm_nf   = 0
        cm_s    = 0
        cm_phi0 = 0
        cm_phi1 = 0
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
        Mfm     = 0


        ! parallel calculation for 1-CMFD
        anode = 0
        if ( zigzagon ) then
            do ii = 1, zz_div
                anode = anode + 2*zzc1(ii)*(zzc0(ii+1)-zzc0(ii))
            end do
        end if
        anode = (ncm(1)*ncm(2) - anode) * ncm(3)
        allocate(mvec(anode,fcr,fcr,fcz,7),svec(anode,fcr,fcr,fcz))
        allocate(ax(anode),ay(anode),az(anode))

        anode = 0
        do ii = 1, ncm(1)
        do jj = 1, ncm(2)
        do kk = 1, ncm(3)
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

    end if

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
    do ij = 1, 6
        fm(:,:,:) % Jn(ij)   = 0 
        fm(:,:,:) % J0(ij)   = 0 
        fm(:,:,:) % J1(ij)   = 0 
    enddo 

    fsd_MC = 0
        
end subroutine

! =============================================================================
! FMFD_INITIALIZE_THREAD initializes thread-wise parameters
! =============================================================================
subroutine FMFD_initialize_thread() 
    integer:: ij
    
    if ( .not. allocated(fm_thread) ) allocate(fm_thread(nfm(1),nfm(2),nfm(3)))
    
    fm_thread(:,:,:) % phi     = 0 
    fm_thread(:,:,:) % sig_t   = 0 
    fm_thread(:,:,:) % sig_a   = 0 
    fm_thread(:,:,:) % nusig_f = 0 
    do ij = 1, 6
    fm_thread(:,:,:) % Jn(ij)   = 0 
    fm_thread(:,:,:) % J0(ij)   = 0 
    fm_thread(:,:,:) % J1(ij)   = 0 
    enddo 
    
end subroutine
    


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
    fm_thread(id(1),id(2),id(3)) % sig_t + flux*macro_xs(1)
    fm_thread(id(1),id(2),id(3)) % sig_a = &
    fm_thread(id(1),id(2),id(3)) % sig_a + flux*macro_xs(2)
    fm_thread(id(1),id(2),id(3)) % nusig_f = &
    fm_thread(id(1),id(2),id(3)) % nusig_f + flux*macro_xs(4)
    
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

! =============================================================================
! NORM_FMFD normalizes cycle-wise FMFD parameters
! =============================================================================
subroutine NORM_FMFD()
    implicit none
    integer:: i, j, k

    !> gather thread FMFD parameters
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
    fm(i,j,k)%phi     = fm(i,j,k)%phi     + fm_thread(i,j,k)%phi
    fm(i,j,k)%sig_t   = fm(i,j,k)%sig_t   + fm_thread(i,j,k)%sig_t 
    fm(i,j,k)%sig_a   = fm(i,j,k)%sig_a   + fm_thread(i,j,k)%sig_a 
    fm(i,j,k)%nusig_f = fm(i,j,k)%nusig_f + fm_thread(i,j,k)%nusig_f 
    do mm = 1, 6
    fm(i,j,k)%J0(mm)  = fm(i,j,k)%J0(mm)  + fm_thread(i,j,k)%J0(mm)
    fm(i,j,k)%J1(mm)  = fm(i,j,k)%J1(mm)  + fm_thread(i,j,k)%J1(mm)
    end do
    end do
    end do
    end do

end subroutine


! =============================================================================
! PROCESS_FMFD deals with MPI process and average quantities
! =============================================================================
subroutine PROCESS_FMFD(bat,cyc)
    use TALLY, only: tallyon
    implicit none
    !> MPI derived type reduce parameters 
    integer, intent(in):: bat
    integer, intent(in):: cyc
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: &
        sig_t0, sig_t1, &
        sig_a0, sig_a1, &
        sig_f0, sig_f1, &
        flux0,  flux1, &
        fsd_MC0
    real(8), dimension(nfm(1),nfm(2),nfm(3),6):: &
        Jp0, Jp1, &
        Jn0, Jn1
    integer :: dsize
    real(8) :: aa, bb
    integer :: ij

    ! data transmission I
    sig_t0(:,:,:) = fm(:,:,:)%sig_t
    sig_a0(:,:,:) = fm(:,:,:)%sig_a
    sig_f0(:,:,:) = fm(:,:,:)%nusig_f
    flux0(:,:,:)  = fm(:,:,:)%phi
    do ij = 1, 6
    Jp0(:,:,:,ij) = fm(:,:,:)%J1(ij)
    Jn0(:,:,:,ij) = fm(:,:,:)%J0(ij)
    end do

    ! data gathering
    dsize = nfm(1)*nfm(2)*nfm(3)
    call MPI_REDUCE(fsd_MC,fsd_MC0,dsize,15,MPI_SUM,score,0,ierr)
    call MPI_REDUCE(sig_t0,sig_t1,dsize,15,MPI_SUM,score,0,ierr)
    call MPI_REDUCE(sig_a0,sig_a1,dsize,15,MPI_SUM,score,0,ierr)
    call MPI_REDUCE(sig_f0,sig_f1,dsize,15,MPI_SUM,score,0,ierr)
    call MPI_REDUCE(flux0,flux1,dsize,15,MPI_SUM,score,0,ierr)
    do ij = 1, 6
    call MPI_REDUCE(Jp0(:,:,:,ij),Jp1(:,:,:,ij),dsize,15,MPI_SUM,score,0,ierr)
    call MPI_REDUCE(Jn0(:,:,:,ij),Jn1(:,:,:,ij),dsize,15,MPI_SUM,score,0,ierr)
    end do

    if ( icore /= score ) return 

    ! data transmission II
    fsd_MC = fsd_MC0
    fm(:,:,:)%sig_t   = sig_t1(:,:,:)
    fm(:,:,:)%sig_a   = sig_a1(:,:,:)
    fm(:,:,:)%nusig_f = sig_f1(:,:,:)
    fm(:,:,:)%phi     = flux1(:,:,:)
    do ij = 1, 6
    fm(:,:,:)%J1(ij) = Jp1(:,:,:,ij)
    fm(:,:,:)%J0(ij) = Jn1(:,:,:,ij)
    end do

    if ( tallyon ) call FMFD_TO_MC(bat,cyc,fm)

    ! current swapping
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 )       fm(ii,jj,kk)%J1(1) = fm(ii-1,jj,kk)%J1(2)
        if ( ii /= nfm(1) )  fm(ii,jj,kk)%J0(2) = fm(ii+1,jj,kk)%J0(1)
        if ( jj /= 1 )       fm(ii,jj,kk)%J1(3) = fm(ii,jj-1,kk)%J1(4)
        if ( jj /= nfm(2) )  fm(ii,jj,kk)%J0(4) = fm(ii,jj+1,kk)%J0(3)
        if ( kk /= 1 )       fm(ii,jj,kk)%J1(5) = fm(ii,jj,kk-1)%J1(6)
        if ( kk /= nfm(3) )  fm(ii,jj,kk)%J0(6) = fm(ii,jj,kk+1)%J0(5)
    end do
    end do
    end do

    ! zigzag at the corner
    if ( zigzagon ) call SET_ZERO_FLUX

    ! group constant
    where ( fm(:,:,:)%phi /= 0 ) 
    fm(:,:,:) % sig_t   = fm(:,:,:) % sig_t   / fm(:,:,:) % phi
    fm(:,:,:) % sig_a   = fm(:,:,:) % sig_a   / fm(:,:,:) % phi
    fm(:,:,:) % nusig_f = fm(:,:,:) % nusig_f / fm(:,:,:) % phi
    fm(:,:,:) % phi     = fm(:,:,:) % phi     / (dble(ngen)*v_fm)
    end where
!    where ( fm(:,:,:)%phi == 0 ) 
!    fm(:,:,:) % sig_t   = fm_avg(:,:,:) % sig_t
!    fm(:,:,:) % sig_a   = fm_avg(:,:,:) % sig_a
!    fm(:,:,:) % nusig_f = fm_avg(:,:,:) % nusig_f
!    end where

    ! surface quantity normalization
    do ii = 1, 6
    bb = dble(ngen)*a_fm(ii)
    fm(:,:,:) % J0(ii)   = fm(:,:,:) % J0(ii) / bb
    fm(:,:,:) % J1(ii)   = fm(:,:,:) % J1(ii) / bb
    end do

    ! net current
    do ii = 1, 6
    fm(:,:,:) % Jn(ii) = fm(:,:,:) % J1(ii) - fm(:,:,:) % J0(ii)
    end do

    ! save the next accumulation
    do ii = n_acc, 2, -1
        acc(ii)%fm(:,:,:) = acc(ii-1)%fm(:,:,:)
    enddo 
    acc(1)%fm(:,:,:) = fm(:,:,:)
    
    !> average with the accumulated parameters
    fm_avg(:,:,:)%phi     = 0 
    fm_avg(:,:,:)%sig_t   = 0 
    fm_avg(:,:,:)%sig_a   = 0 
    fm_avg(:,:,:)%nusig_f = 0 
    do ii = 1, 6 
      fm_avg(:,:,:)%Jn(ii) = 0 
      fm_avg(:,:,:)%J0(ii) = 0 
      fm_avg(:,:,:)%J1(ii) = 0 
    enddo 

    ! accumulation
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
    do mm = 1, n_acc
       fm_avg(ii,jj,kk)%phi     = fm_avg(ii,jj,kk)%phi     &
                                + acc(mm)%fm(ii,jj,kk)%phi
       fm_avg(ii,jj,kk)%sig_t   = fm_avg(ii,jj,kk)%sig_t   &
                                + acc(mm)%fm(ii,jj,kk)%sig_t
       fm_avg(ii,jj,kk)%sig_a   = fm_avg(ii,jj,kk)%sig_a   &
                                + acc(mm)%fm(ii,jj,kk)%sig_a
       fm_avg(ii,jj,kk)%nusig_f = fm_avg(ii,jj,kk)%nusig_f &
                                + acc(mm)%fm(ii,jj,kk)%nusig_f

       fm_avg(ii,jj,kk)%Jn(:)   = fm_avg(ii,jj,kk)%Jn(:)   &
                                + acc(mm)%fm(ii,jj,kk)%Jn(:)
       fm_avg(ii,jj,kk)%J0(:)   = fm_avg(ii,jj,kk)%J0(:)   &
                                + acc(mm)%fm(ii,jj,kk)%J0(:)
       fm_avg(ii,jj,kk)%J1(:)   = fm_avg(ii,jj,kk)%J1(:)   &
                                + acc(mm)%fm(ii,jj,kk)%J1(:)
    end do
    end do
    end do
    end do
    
    ! average
    fm_avg(:,:,:)%phi     = fm_avg(:,:,:)%phi     / dble(n_acc)
    fm_avg(:,:,:)%sig_t   = fm_avg(:,:,:)%sig_t   / dble(n_acc)
    fm_avg(:,:,:)%sig_a   = fm_avg(:,:,:)%sig_a   / dble(n_acc)
    fm_avg(:,:,:)%nusig_f = fm_avg(:,:,:)%nusig_f / dble(n_acc)
    do ii = 1, 6
    fm_avg(:,:,:)%Jn(ii)  = fm_avg(:,:,:)%Jn(ii)  / dble(n_acc)
    fm_avg(:,:,:)%J0(ii)  = fm_avg(:,:,:)%J0(ii)  / dble(n_acc)
    fm_avg(:,:,:)%J1(ii)  = fm_avg(:,:,:)%J1(ii)  / dble(n_acc)
    enddo

end subroutine 

! =============================================================================
! SET_ZERO_FLUX
! =============================================================================
subroutine SET_ZERO_FLUX
    implicit none
    
    do ii = 1, zz_div
        fm(zzf0(ii)+1:zzf0(ii+1),1:zzf1(ii),:)%phi        = 0D0
        fm(zzf0(ii)+1:zzf0(ii+1),zzf2(ii)+1:nfm(2),:)%phi = 0D0
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

    fetnusigf(:,:,:) = fm(:,:,:)%nusig_f

end subroutine

! =============================================================================
! FMFD_SOLVE solves FMFD eigenvalue problem
! =============================================================================
subroutine FMFD_SOLVE(keff,fsd,bat,cyc)
    use CMFD, only: ONE_NODE_CMFD
    use TH_HEADER, only: th_on, pp
    implicit none
    real(8), intent(inout) :: keff          ! multiplication factor
    real(8), intent(inout) :: fsd(:,:,:)    ! fission source distribution
    integer, intent(in)    :: bat           ! batch number
    integer, intent(in)    :: cyc           ! cycle number
    real(8) :: M(nfm(1),nfm(2),nfm(3),7)    ! FMFD matrix
    real(8), dimension(nfm(1),nfm(2),nfm(3)):: &
        sig_t, &      ! total
        sig_a, &      ! absorption
        nusig_f, &    ! nu X fission
        D, &          ! diffusion coefficient
        phi0, &       ! neutron flux (previous)
        phi1          ! neutron flux (current)
    real(8), dimension(nfm(1),nfm(2),nfm(3),6):: &
        D_tilda, &    ! D tilda
        D_hat, &      ! correction factor
        Jn, &         ! net current
        J0, &         ! partial current -
        J1, &         ! partial current +
        sphi          ! surface flux
    integer:: acyc

    if ( icore == score ) then
    ! copy parameters
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        sig_t(ii,jj,kk)   = fm_avg(ii,jj,kk)%sig_t
        sig_a(ii,jj,kk)   = fm_avg(ii,jj,kk)%sig_a
        nusig_f(ii,jj,kk) = fm_avg(ii,jj,kk)%nusig_f 
        D(ii,jj,kk)       = 1D0 / (3D0 * fm_avg(ii,jj,kk)%sig_t) 
        Jn(ii,jj,kk,:)    = fm_avg(ii,jj,kk)%Jn(:) 
        J0(ii,jj,kk,:)    = fm_avg(ii,jj,kk)%J0(:) 
        J1(ii,jj,kk,:)    = fm_avg(ii,jj,kk)%J1(:) 
        phi1(ii,jj,kk)    = fm_avg(ii,jj,kk)%phi
    enddo 
    enddo
    enddo
    where( phi1 == 0 ) D = 0
    end if

    if ( cmfdon ) then
        if ( icore == score ) then
        !print*, "CMFD"
        do ii = 1, 6
            sphi(:,:,:,ii) = 2D0*J1(:,:,:,ii)+2D0*J0(:,:,:,ii)
        end do
        end if

        call ONE_NODE_CMFD(keff,sig_t,sig_a,nusig_f,D,phi1,J0,J1,Jn,sphi)

    else
        if ( icore == score ) then
        call D_TILDA_CALCULATION(D,D_tilda)
        if ( .not. pfmfdon ) then
            call D_HAT_CALCULATION(D_tilda,Jn,phi1,D_hat)
            M = getM(D_tilda,D_hat,sig_a)
        else
            call D_PHAT_CALCULATION(D_tilda,J0,J1,Jn,phi1,D_hat)
            M = getPM(D_tilda,D_hat,sig_a)
        endif 

        if ( zigzagon ) call SET_ZERO_M(phi1,M)
        call POWER (keff, M, phi0, phi1, nusig_f)
        end if
    end if

    if ( icore /= score ) return

    !> CMFD feedback (modulation)
    print*, "keff", keff
    !pause

    if ( ( .not. fmfd2mc .and. cyc > n_inact ) &
        .or. isnan(keff) .or. ( keff < 1D-2 .or. keff > 2D0 ) ) then
        fsd = 1D0
        if ( th_on ) call DET_POWER(nusig_f)
    else
        where ( phi1 == 0 ) fsd_MC = 0
        fsd_MC(:,:,:) = fsd_MC(:,:,:) / sum(fsd_MC)
        fsd_FM(:,:,:) = nusig_f(:,:,:)*phi1(:,:,:)
        fsd_FM(:,:,:) = fsd_FM(:,:,:) / sum(fsd_FM)
        if ( th_on ) call DET_POWER(fsd_FM)
        where ( fsd_MC(:,:,:) /= 0 ) &
        fsd(:,:,:)    = fsd_FM(:,:,:) / fsd_MC(:,:,:)
    end if

    call FMFD_ENTRP(fsd_FM,cyc)
    acyc = cyc - n_inact
    if ( acyc > 0 ) p_fmfd(bat,acyc,:,:,:) = nusig_f(:,:,:)*phi1(:,:,:)
    
end subroutine

! =============================================================================
! D_TILDA_CALCULATION
! =============================================================================
subroutine D_TILDA_CALCULATION(D,Dt)
    implicit none
    real(8), intent(in) :: D(:,:,:)
    real(8), intent(out):: Dt(:,:,:,:)
    
    Dt = 0

    ! inner region
    do ii = 1, nfm(1)
    do jj = 1, nfm(2)
    do kk = 1, nfm(3)
        if ( ii /= 1 ) &      ! x0
        Dt(ii,jj,kk,1) = &
            2D0*D(ii,jj,kk)*D(ii-1,jj,kk)/((D(ii,jj,kk)+D(ii-1,jj,kk))*dfm(1))
        if ( ii /= nfm(1) ) & ! x1
        Dt(ii,jj,kk,2) = &
            2D0*D(ii+1,jj,kk)*D(ii,jj,kk)/((D(ii+1,jj,kk)+D(ii,jj,kk))*dfm(1))
        if ( jj /= 1 ) &      ! y0
        Dt(ii,jj,kk,3) = &
            2D0*D(ii,jj,kk)*D(ii,jj-1,kk)/((D(ii,jj,kk)+D(ii,jj-1,kk))*dfm(2))
        if ( jj /= nfm(2) ) & ! y1
        Dt(ii,jj,kk,4) = &
            2D0*D(ii,jj+1,kk)*D(ii,jj,kk)/((D(ii,jj+1,kk)+D(ii,jj,kk))*dfm(2))
        if ( kk /= 1 ) &      ! y0
        Dt(ii,jj,kk,5) = &
            2D0*D(ii,jj,kk)*D(ii,jj,kk-1)/((D(ii,jj,kk)+D(ii,jj,kk-1))*dfm(3))
        if ( kk /= nfm(3) ) & ! y1
        Dt(ii,jj,kk,6) = &
            2D0*D(ii,jj,kk+1)*D(ii,jj,kk)/((D(ii,jj,kk+1)+D(ii,jj,kk))*dfm(3))
    end do
    end do
    end do

    where ( isnan(Dt) ) Dt = 0

end subroutine


! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_HAT_CALCULATION(Dt,Jn,phi,Dh)
    real(8), intent(in) :: Dt(:,:,:,:), Jn(:,:,:,:), phi(:,:,:)
    real(8), intent(out):: Dh(:,:,:,:)
    integer :: i, j, k

    ! inner region
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        Dh(i,j,k,1) = (Jn(i,j,k,1)+Dt(i,j,k,1)*(phi(i,j,k)-phi(i-1,j,k)))/ &
                    (phi(i,j,k)+phi(i-1,j,k))
        if ( i /= nfm(1) ) & ! x1
        Dh(i,j,k,2) = (Jn(i,j,k,2)+Dt(i,j,k,2)*(phi(i+1,j,k)-phi(i,j,k)))/ &
                    (phi(i+1,j,k)+phi(i,j,k))
        if ( j /= 1 ) &      ! y0
        Dh(i,j,k,3) = (Jn(i,j,k,3)+Dt(i,j,k,3)*(phi(i,j,k)-phi(i,j-1,k)))/ &
                    (phi(i,j,k)+phi(i,j-1,k))
        if ( j /= nfm(2) ) & ! y1
        Dh(i,j,k,4) = (Jn(i,j,k,4)+Dt(i,j,k,4)*(phi(i,j+1,k)-phi(i,j,k)))/ &
                    (phi(i,j+1,k)+phi(i,j,k))
        if ( k /= 1 ) &      ! y0
        Dh(i,j,k,5) = (Jn(i,j,k,5)+Dt(i,j,k,5)*(phi(i,j,k)-phi(i,j,k-1)))/ &
                    (phi(i,j,k)+phi(i,j,k-1))
        if ( k /= nfm(3) ) & ! y1
        Dh(i,j,k,6) = (Jn(i,j,k,6)+Dt(i,j,k,6)*(phi(i,j,k+1)-phi(i,j,k)))/ &
                    (phi(i,j,k+1)+phi(i,j,k))
    end do
    end do
    end do

    ! boundary
    i = 1;      Dh(i,:,:,1) = Jn(i,:,:,1)/phi(i,:,:)
    i = nfm(1); Dh(i,:,:,2) = Jn(i,:,:,2)/phi(i,:,:)
    j = 1;      Dh(:,j,:,3) = Jn(:,j,:,3)/phi(:,j,:)
    j = nfm(2); Dh(:,j,:,4) = Jn(:,j,:,4)/phi(:,j,:)
    k = 1;      Dh(:,:,k,5) = Jn(:,:,k,5)/phi(:,:,k)
    k = nfm(3); Dh(:,:,k,6) = Jn(:,:,k,6)/phi(:,:,k)

end subroutine


! =============================================================================
! D_HAT_CALCULATION calculates correction factors
! =============================================================================
subroutine D_PHAT_CALCULATION(Dt,J0,J1,Jn,phi,Dh)
    real(8), intent(in) :: Dt(:,:,:,:), J0(:,:,:,:), J1(:,:,:,:), &
                           Jn(:,:,:,:), phi(:,:,:)
    real(8), intent(out):: Dh(:,:,:,:)
    integer :: i, j, k

    ! inner region (outgoing direction)
    do i = 1, nfm(1)
    do j = 1, nfm(2)
    do k = 1, nfm(3)
        if ( i /= 1 ) &      ! x0
        Dh(i,j,k,1) = (J0(i,j,k,1)-5D-1*Dt(i,j,k,1) &
                    *(phi(i,j,k)-phi(i-1,j,k)))/phi(i,j,k)
        if ( i /= nfm(1) ) & ! x1
        Dh(i,j,k,2) = (J1(i,j,k,2)+5D-1*Dt(i,j,k,2) &
                    *(phi(i+1,j,k)-phi(i,j,k)))/phi(i,j,k)
        if ( j /= 1 ) &      ! y0
        Dh(i,j,k,3) = (J0(i,j,k,3)-5D-1*Dt(i,j,k,3) &
                    *(phi(i,j,k)-phi(i,j-1,k)))/phi(i,j,k)
        if ( j /= nfm(2) ) & ! y1
        Dh(i,j,k,4) = (J1(i,j,k,4)+5D-1*Dt(i,j,k,4) &
                    *(phi(i,j+1,k)-phi(i,j,k)))/phi(i,j,k)
        if ( k /= 1 ) &      ! y0
        Dh(i,j,k,5) = (J0(i,j,k,5)-5D-1*Dt(i,j,k,5) &
                    *(phi(i,j,k)-phi(i,j,k-1)))/phi(i,j,k)
        if ( k /= nfm(3) ) & ! y1
        Dh(i,j,k,6) = (J1(i,j,k,6)+5D-1*Dt(i,j,k,6) &
                    *(phi(i,j,k+1)-phi(i,j,k)))/phi(i,j,k)
    end do
    end do
    end do

    ! boundary
    if ( .not. zigzagon ) then
    i = 1;      Dh(i,:,:,1) = -Jn(i,:,:,1)/phi(i,:,:)
    i = nfm(1); Dh(i,:,:,2) = +Jn(i,:,:,2)/phi(i,:,:)
    j = 1;      Dh(:,j,:,3) = -Jn(:,j,:,3)/phi(:,j,:)
    j = nfm(2); Dh(:,j,:,4) = +Jn(:,j,:,4)/phi(:,j,:)
    else
    do i = 1, zz_div
        Dh(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) = &
            -Jn(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:,1) &
            /phi(zzf1(i)+1,zzf0(i)+1:zzf0(i+1),:)
        Dh(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) = &
            +Jn(zzf2(i),zzf0(i)+1:zzf0(i+1),:,2) &
            /phi(zzf2(i),zzf0(i)+1:zzf0(i+1),:)
        Dh(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) = &
            -Jn(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:,3) &
            /phi(zzf0(i)+1:zzf0(i+1),zzf1(i)+1,:)
        Dh(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) = &
            +Jn(zzf0(i)+1:zzf0(i+1),zzf2(i),:,4) &
            /phi(zzf0(i)+1:zzf0(i+1),zzf2(i),:)
    end do
    end if
    k = 1;      Dh(:,:,k,5) = -Jn(:,:,k,5)/phi(:,:,k)
    k = nfm(3); Dh(:,:,k,6) = +Jn(:,:,k,6)/phi(:,:,k)

end subroutine


! ========================================================= !
!  Miscel. functions for CMFD_solve
! ========================================================= !

! function which produces M matrix from xs & D_hat
function getM (Dt,Dh,sig_a) result(M)
    real(8), intent(in) :: Dt(:,:,:,:)
    real(8), intent(in) :: Dh(:,:,:,:)
    real(8), intent(in) :: sig_a(:,:,:)
    real(8):: M(nfm(1),nfm(2),nfm(3),7)
    integer :: i, j, k
    
    ! initialization
    M = 0
    
    ! M matrix set 
    do k = 1, nfm(3)
    do j = 1, nfm(2)
    do i = 1, nfm(1)
        if ( i /= 1 )      M(i,j,k,3) = -(Dt(i,j,k,1)+Dh(i,j,k,1))/dfm(1) ! x0
        if ( i /= nfm(1) ) M(i,j,k,5) = -(Dt(i,j,k,2)-Dh(i,j,k,2))/dfm(1) ! x1
        if ( j /= 1 )      M(i,j,k,2) = -(Dt(i,j,k,3)+Dh(i,j,k,3))/dfm(2) ! y0
        if ( j /= nfm(2) ) M(i,j,k,6) = -(Dt(i,j,k,4)-Dh(i,j,k,4))/dfm(2) ! y1
        if ( k /= 1 )      M(i,j,k,1) = -(Dt(i,j,k,5)+Dh(i,j,k,5))/dfm(3) ! z0
        if ( k /= nfm(3) ) M(i,j,k,7) = -(Dt(i,j,k,6)-Dh(i,j,k,6))/dfm(3) ! z1
        
        M(i,j,k,4) = &
            +(Dt(i,j,k,1)-Dh(i,j,k,1)+Dt(i,j,k,2)+Dh(i,j,k,2))/dfm(1) &
            +(Dt(i,j,k,3)-Dh(i,j,k,3)+Dt(i,j,k,4)+Dh(i,j,k,4))/dfm(2) &
            +(Dt(i,j,k,5)-Dh(i,j,k,5)+Dt(i,j,k,6)+Dh(i,j,k,6))/dfm(3) &
            +sig_a(i,j,k)

    enddo
    enddo
    enddo

end function getM


! function which produces M matrix from xs & D_hat
function getPM (Dt,Dh,sig_a) result(M)
    real(8), intent(in) :: Dt(:,:,:,:)
    real(8), intent(in) :: Dh(:,:,:,:)
    real(8), intent(in) :: sig_a(:,:,:)
    real(8):: M(nfm(1),nfm(2),nfm(3),7)
    integer :: i, j, k
    
    ! initialization
    M = 0
    
    ! M matrix set 
    do k = 1, nfm(3)
    do j = 1, nfm(2)
    do i = 1, nfm(1)
        if ( i /= 1 )      M(i,j,k,3) = -(Dt(i,j,k,1)+Dh(i-1,j,k,2))/dfm(1) ! x0
        if ( i /= nfm(1) ) M(i,j,k,5) = -(Dt(i,j,k,2)+Dh(i+1,j,k,1))/dfm(1) ! x1
        if ( j /= 1 )      M(i,j,k,2) = -(Dt(i,j,k,3)+Dh(i,j-1,k,4))/dfm(2) ! y0
        if ( j /= nfm(2) ) M(i,j,k,6) = -(Dt(i,j,k,4)+Dh(i,j+1,k,3))/dfm(2) ! y1
        if ( k /= 1 )      M(i,j,k,1) = -(Dt(i,j,k,5)+Dh(i,j,k-1,6))/dfm(3) ! z0
        if ( k /= nfm(3) ) M(i,j,k,7) = -(Dt(i,j,k,6)+Dh(i,j,k+1,5))/dfm(3) ! z1
        
        M(i,j,k,4) = &
            +(Dt(i,j,k,1)+Dh(i,j,k,1)+Dt(i,j,k,2)+Dh(i,j,k,2))/dfm(1) &
            +(Dt(i,j,k,3)+Dh(i,j,k,3)+Dt(i,j,k,4)+Dh(i,j,k,4))/dfm(2) &
            +(Dt(i,j,k,5)+Dh(i,j,k,5)+Dt(i,j,k,6)+Dh(i,j,k,6))/dfm(3) &
            +sig_a(i,j,k)

    enddo
    enddo
    enddo

end function getPM

subroutine SET_ZERO_M(phi1,M)
    implicit none
    real(8), intent(in)::  phi1(:,:,:)
    real(8), intent(out):: M(:,:,:,:)
    integer:: ij

    where( phi1 == 0 ) 
    M(:,:,:,1) = 0
    M(:,:,:,2) = 0
    M(:,:,:,3) = 0
    M(:,:,:,4) = 1
    M(:,:,:,5) = 0
    M(:,:,:,6) = 0
    M(:,:,:,7) = 0
    end where

    do ij = 1, zz_div
        M(zzf1(ij)+1,zzf0(ij)+1:zzf0(ij+1),:,3) = 0
        M(zzf2(ij),zzf0(ij)+1:zzf0(ij+1),:,5) = 0
        M(zzf0(ij)+1:zzf0(ij+1),zzf1(ij)+1,:,2) = 0
        M(zzf0(ij)+1:zzf0(ij+1),zzf2(ij),:,6) = 0
    end do

end subroutine


subroutine POWER (k_eff, M, phi0, phi1, nusig_f)
    use SOLVERS, only: BICGStab_hepta !, SOR
    implicit none
    real(8), intent(inout):: k_eff
    real(8), intent(in)   :: M(:,:,:,:), nusig_f(:,:,:)
    real(8), intent(inout):: phi0(:,:,:), phi1(:,:,:)
    real(8), parameter:: ONE = 1D0
    real(8) :: F(1:nfm(1),1:nfm(2),1:nfm(3))
    integer :: iter, iter_max = 2D2 
    real(8) :: err
    real(8):: tt0, tt1
    real(8):: kpre
    
    err = ONE
    iter = 1
    do while ( ( err > 1D-10 ) .and. ( iter < iter_max ) )
        !call CPU_TIME(tt0)
        print*, k_eff
        iter = iter + 1
        phi0 = phi1
        kpre = k_eff
        F = nusig_f(:,:,:)*phi1(:,:,:)/k_eff
        !call SOR(M(:,:,:,:),F(:,:,:),phi1(:,:,:))
        phi1 = BiCGStab_hepta(M(:,:,:,:),F(:,:,:))
        k_eff = k_eff*sum(nusig_f*phi1*nusig_f*phi1) &
               / sum(nusig_f*phi1*nusig_f*phi0)
        !err = maxval(abs(phi1(:,:,:)-phi0(:,:,:))/phi1(:,:,:))
        err = abs(k_eff-kpre)/k_eff
        !call CPU_TIME(tt1)
        !write(7,*), err, tt1-tt0
    enddo
    
end subroutine 

subroutine FMFD_ENTRP(fsd,cyc)
    use ENTROPY, only: entrp3
    implicit none
    real(8), intent(in) :: fsd(:,:,:)
    integer, intent(in) :: cyc
    real(8):: ee0(nfm(1),nfm(2),nfm(3))


    ee0 = fsd/sum(fsd)
    where ( ee0 /= 0 ) ee0 = -ee0*log(ee0)/log(2D0)
    entrp3(cyc) = sum(ee0)

end subroutine


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

end module
