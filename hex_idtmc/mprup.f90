module MPRUP
    use VARIABLES, only: icore, score, ngen, n_totcyc, n_inact, n_act
    use ENTROPY
    use FMFD_HEADER, only: fmfdon, k_fmfd
    implicit none

    ! * Convergence criteria
    !   1 convergence : crt1
    !    (1) relative error
    !    (2) length of cycles
    !   2 convergence : crt2
    !    (1) relative error
    !    (2) final generation size

    contains

! =============================================================================
! 
! =============================================================================
subroutine GENSIZE(bat,cyc)
    use FMFD_HEADER, only: acc, n_acc, k_mprup
    use ENTROPY, only: mprupon
    implicit none
    integer, intent(in) :: bat,cyc

    ! FMFD divergence
    if ( cyc > 1 .and. fmfdon ) then
    if ( isnan(k_mprup) .or. k_mprup < 1D-1 .or. k_mprup > 2D0 ) then
        call GENSIZEUP
        return
    end if
    end if

    if ( genup ) then
    ! 1 convergence test
    select case(ccrt)
    case(1); if ( dshannon > crt1 ) return
    case(2); if ( cyc == 1 .or. mod(cyc-1,int(crt1)) /= 0 ) return
    end select

    ! 2 stationary point
    select case(scrt)
    case(1)
        if ( dual_fmfd ) then
        entrp1 = sum(shannon(1:elength))/dble(elength)
        entrp2 = sum(shannon_(1:elength))/dble(elength)
        dentrp = abs(entrp2-entrp1)/entrp2
        !print*, "dentrp", dentrp, crt2
        if ( dentrp < crt2 ) then
            genup = .false.
            return
        end if
        else
        entrp2 = sum(shannon(1:elength))/dble(elength)
        dentrp = abs(entrp2-entrp1)/entrp2
        if ( dentrp < crt2 ) then
            genup = .false.
            return
        end if
        entrp1 = entrp2
        end if

        call GENSIZEUP
        return

    case(2)
        if ( ngen >= crt2 ) then
            genup = .false.
            return
        end if
        call GENSIZEUP
        return

    end select

    ! 3 stopping test
    else
        !print*, dshannon, crt3
        if ( dshannon < crt3 .or. curr_cyc == crt4c ) then
            up_sign = .true.
            mprupon = .false.
            genup = .true.
            return
        end if
    end if


end subroutine


! =============================================================================
! GENSIZEUP increases the generation size and reallocate the fission neutron
! =============================================================================
subroutine GENSIZEUP
    use BANK_HEADER,    only: source_bank
    use FMFD_HEADER,    only: fmfdon
    implicit none

    source_bank(:)%wgt = (ngen+rampup)/dble(ngen)
    ngen = ngen + rampup
    up_sign = .true.
    if ( scrt == 2 .and. ngen > crt2c ) ngen = crt2c

    ! update criteria
    if ( fmfdon ) then
        if ( ccrt == 1 ) then
            crt1 = crt1c/sqrt(dble(ngen))
        end if
        if ( scrt == 1 ) then
            crt2 = crt2c/sqrt(dble(ngen))
        end if
        crt3 = crt3c/sqrt(dble(ngen))
    else
        if ( ccrt == 1 ) then
            crt1 = crt1c/sqrt(dble(ngen))
        end if
        if ( scrt == 1 ) then
            crt2 = crt2c/sqrt(dble(ngen))
        end if
        crt3 = crt3c/sqrt(dble(ngen))
    end if

end subroutine

! =============================================================================
! CYCLECHANGE goes to active cycle
! =============================================================================
subroutine CYCLECHANGE(cyc)
    use VARIABLES, only: n_batch
    use FMFD_HEADER, only: fake_MC, n_fake, inactive_CMFD
    implicit none
    integer, intent(in):: cyc

    if ( .not. fake_MC ) then
    n_inact = cyc + 3
    else
    if ( inactive_CMFD ) then
    if ( .not. do_burn ) then
    n_inact = cyc + 8
    n_fake  = cyc + 8
    else
    n_inact = cyc + 3
    n_fake  = cyc + 3
    end if
    else
    n_inact = cyc + 7
    n_fake  = cyc
    end if
    end if
    genup = .false.
    if ( n_batch == 1 ) then
        n_totcyc = n_inact + n_act
    else
        n_totcyc = n_inact
    end if

end subroutine

! =============================================================================
! MPRUP_DIST distributes the necessary parameters to the nodes
! =============================================================================
subroutine MPRUP_DIST(sz,source_bank)
    use MPI,         only: MPI_REAL8, MPI_COMM_WORLD
    use BANK_HEADER, only: Bank
    implicit none
    integer, intent(in):: sz
    type(Bank):: source_bank(:)
    integer:: TP
    integer:: WOR
    integer:: ii

    TP = MPI_REAL8
    WOR = MPI_COMM_WORLD

    do ii = 1, sz
    call MPI_BCAST(source_bank(ii)%wgt,1,TP,score,WOR,ierr)
    end do
    call MPI_BCAST(ngen,1,TP,score,WOR,ierr)
    call MPI_BCAST(genup,1,TP,score,WOR,ierr)
    call MPI_BCAST(mprupon,1,TP,score,WOR,ierr)
   
end subroutine

end module
