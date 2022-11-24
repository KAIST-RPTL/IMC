module COSAMPLING
    use FMFD_HEADER
    use mpi
    use VARIABLES, only: n_totcyc, iscore
    implicit none
    real(8), allocatable, dimension(:,:,:,:):: acc_sigt, acc_siga
    real(8), allocatable, dimension(:,:,:,:):: acc_nufi, acc_phi, acc_kap
    real(8), allocatable, dimension(:,:,:,:,:):: acc_J0, acc_J1, acc_Jn
    real(8), allocatable:: fphi2(:,:,:)
    real(8):: sqrt2, sqrt2pi
    real(8):: cyclen
    real(8), allocatable:: avg_(:), std_(:), corr(:,:), cor(:,:), chol(:,:), invchol(:,:), tmpchol(:,:)
    real(8), allocatable:: lhs(:,:), unlhs(:,:)
    real(8), allocatable:: courn(:,:,:,:,:), coxs(:,:,:,:,:)

    integer, parameter:: n_pert = 100
    real(8):: k_pert(n_pert)    ! by perturbation with forward flux
    real(8):: k_pert2(n_pert)   ! with adjoint flux
    real(8):: k_pert3(n_pert)   ! direct calcultion
    real(8), allocatable:: p_pert(:,:,:,:)
    real(8), allocatable:: s_pert(:)
    integer:: types, n_data
    real(8) :: corcrit = 1d-2

    ! correlation between nodes
    real(8), allocatable:: avgn(:), stdn(:)


    contains

    ! =========================================================================
    ! 
    ! =========================================================================
    subroutine COSAMPLING_MAIN(cyc)
        use STATISTICS, only: AVG
        implicit none 
        integer, intent(in):: cyc
        real(8):: avg_cor(3), oavg(3), ostd(3), incor(3,3), tmp, ucor(3,3), diffcor(3,3), prev(3,3)
        integer:: nodes1, nodes2, i, it
        integer :: aa
        ! parameter allocation
        if ( .not. allocated(avg_)) then
            sqrt2 = sqrt(2.0)
            allocate(avg_(3),std_(3),corr(3,3),cor(3,3),chol(3,3),invchol(3,3),tmpchol(3,3))
            allocate(lhs(n_pert,3),unlhs(n_pert,3))
            allocate(courn(nfm(1),nfm(2),nfm(3),n_pert,3))
            allocate(coxs(nfm(1),nfm(2),nfm(3),n_pert,3))
            allocate(p_pert(n_pert,nfm(1),nfm(2),nfm(3)))
        end if
        
        !call DATA_INO(cyc)
    
        cyclen = dble(cyc-acc_skip)
!        nodes1 = 0
!        nodes2 = 0
!        avg_cor(1:3) = 0
        do kk = 1, nfm(3)
        do jj = 1, nfm(2)
        do ii = 1, nfm(1)
            if ( fphi1(ii,jj,kk) == 0 ) cycle
!            ! correlation matrix
            call AVG_N_STD(acc_sigt(acc_skip+1:cyc,ii,jj,kk), &
                           acc_siga(acc_skip+1:cyc,ii,jj,kk), &
                           acc_nufi(acc_skip+1:cyc,ii,jj,kk),cyclen)
            if ( std_(3) == 0 ) then
                types = 2
            else
                types = 3
            end if
            call CORRELATION_MATRIX(acc_sigt(acc_skip+1:cyc,ii,jj,kk), & 
                                    acc_siga(acc_skip+1:cyc,ii,jj,kk), &
                                    acc_nufi(acc_skip+1:cyc,ii,jj,kk),cyclen)

            oavg(1:types) = avg_(1:types)
            ostd(1:types) = std_(1:types)
            incor = cor

            ! random number generation by LHS
            call LHS_SAMPLING
            
            ! inverse normal CDF
            do mm = 1, n_pert
                lhs(mm,:) = ICDF_NORMAL(lhs(mm,:))
            end do

            ! Cholesky decomposition
            !corcvg = .false.
            !do while( .not. corcvg)
            call DECHOLESKY(cor, chol)
            tmpchol = chol
            
            call AVG_N_STD(lhs(:,1), &
                           lhs(:,2), &
                           lhs(:,3),dble(n_pert))
            call CORRELATION_MATRIX(lhs(:,1), &
                           lhs(:,2), &
                           lhs(:,3),dble(n_pert))
            
            call DECHOLESKY(cor, chol)
            
            select case(types)
            case(2)
                invchol(1:2,1:2) = matinv2(chol)
            case(3)
                invchol(1:3,1:3) = matinv3(chol)
            end select

            do mm = 1,n_pert
            do nn = 1,types
                unlhs(mm,nn) = sum(lhs(mm,1:types)*invchol(1:types,nn))
            enddo
            enddo
            
            do mm = 1, n_pert
            do nn = 1, types
                courn(ii,jj,kk,mm,nn) = sum(unlhs(mm,1:types)*tmpchol(1:types,nn))
            end do
            end do

            call AVG_N_STD(courn(ii,jj,kk,:,1), & 
                           courn(ii,jj,kk,:,2), &
                           courn(ii,jj,kk,:,3),dble(n_pert))
            call CORRELATION_MATRIX(courn(ii,jj,kk,:,1), & 
                                    courn(ii,jj,kk,:,2), &
                                    courn(ii,jj,kk,:,3),dble(n_pert))

            ucor = cor
            ! CFD for correlated urn generation
            do mm = 1, n_pert
            courn(ii,jj,kk,mm,:) = CDF_COURN(courn(ii,jj,kk,mm,:))
            end do
            call AVG_N_STD(courn(ii,jj,kk,:,1), & 
                           courn(ii,jj,kk,:,2), &
                           courn(ii,jj,kk,:,3),dble(n_pert))
            call CORRELATION_MATRIX(courn(ii,jj,kk,:,1), & 
                                    courn(ii,jj,kk,:,2), &
                                    courn(ii,jj,kk,:,3),dble(n_pert))
            
            ! inverse CDF for XS
            coxs(ii,jj,kk,:,:) = ICDF_XS(acc_sigt(acc_skip+1:cyc,ii,jj,kk), & 
                                         acc_siga(acc_skip+1:cyc,ii,jj,kk), &
                                         acc_nufi(acc_skip+1:cyc,ii,jj,kk), &
                                         courn(ii,jj,kk,:,:))
            ! inverse CDF for Jn & phi for BC
            if(iscore) ptphi(:,ii,jj,kk)  = ICDF_PHI(acc_phi(acc_skip+1:cyc,ii,jj,kk))

            call AVG_N_STD(coxs(ii,jj,kk,:,1), & 
                           coxs(ii,jj,kk,:,2), &
                           coxs(ii,jj,kk,:,3),dble(n_pert))

            do i = 1,types
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) - avg_(i)
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) * ostd(i) / std_(i)
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) + oavg(i)
            enddo
            
            call AVG_N_STD(coxs(ii,jj,kk,:,1), & 
                           coxs(ii,jj,kk,:,2), &
                           coxs(ii,jj,kk,:,3),dble(n_pert))

            call CORRELATION_MATRIX(coxs(ii,jj,kk,:,1), & 
                                    coxs(ii,jj,kk,:,2), &
                                    coxs(ii,jj,kk,:,3), dble(n_pert))

            ! ADJUST COR: COR = 2COR-INCOR
            ! IF DET(COR) <= 0: Just return
            ! IF DET(COR) >  0: Perform other CORSAMPLING
            if(cyc==131 .and. iscore) then
            print *, 'MESH', ii, jj, kk
            print *, 'INCOR', incor
            print *, 'RESCOR', cor
            print *, 'PERTCOR', incor + (incor-cor)*1d0
            endif

            prev = incor

            do it = 1, 5
            cor = prev + (incor - cor) * 2d-1
            prev= cor

            if(cor(1,1)*(cor(2,2)*cor(3,3)-cor(2,3)*cor(3,2)) &
              -cor(1,2)*(cor(2,1)*cor(3,3)-cor(2,3)*cor(3,1)) &
              +cor(1,3)*(cor(2,1)*cor(3,2)-cor(2,2)*cor(3,1)) <= 0) exit

            ! random number generation by LHS
            call LHS_SAMPLING
            
            ! inverse normal CDF
            do mm = 1, n_pert
                lhs(mm,:) = ICDF_NORMAL(lhs(mm,:))
            end do

            ! Cholesky decomposition
            !corcvg = .false.
            !do while( .not. corcvg)
            call DECHOLESKY(cor, chol)
            tmpchol = chol
            
            call AVG_N_STD(lhs(:,1), &
                           lhs(:,2), &
                           lhs(:,3),dble(n_pert))
            call CORRELATION_MATRIX(lhs(:,1), &
                           lhs(:,2), &
                           lhs(:,3),dble(n_pert))
            
            call DECHOLESKY(cor, chol)
            
            select case(types)
            case(2)
                invchol(1:2,1:2) = matinv2(chol)
            case(3)
                invchol(1:3,1:3) = matinv3(chol)
            end select

            do mm = 1,n_pert
            do nn = 1,types
                unlhs(mm,nn) = sum(lhs(mm,1:types)*invchol(1:types,nn))
            enddo
            enddo
            
            do mm = 1, n_pert
            do nn = 1, types
                courn(ii,jj,kk,mm,nn) = sum(unlhs(mm,1:types)*tmpchol(1:types,nn))
            end do
            end do

            call AVG_N_STD(courn(ii,jj,kk,:,1), & 
                           courn(ii,jj,kk,:,2), &
                           courn(ii,jj,kk,:,3),dble(n_pert))
            call CORRELATION_MATRIX(courn(ii,jj,kk,:,1), & 
                                    courn(ii,jj,kk,:,2), &
                                    courn(ii,jj,kk,:,3),dble(n_pert))

            ! CFD for correlated urn generation
            do mm = 1, n_pert
            courn(ii,jj,kk,mm,:) = CDF_COURN(courn(ii,jj,kk,mm,:))
            end do
            call AVG_N_STD(courn(ii,jj,kk,:,1), & 
                           courn(ii,jj,kk,:,2), &
                           courn(ii,jj,kk,:,3),dble(n_pert))
            call CORRELATION_MATRIX(courn(ii,jj,kk,:,1), & 
                                    courn(ii,jj,kk,:,2), &
                                    courn(ii,jj,kk,:,3),dble(n_pert))
            
            ! inverse CDF for XS
            coxs(ii,jj,kk,:,:) = ICDF_XS(acc_sigt(acc_skip+1:cyc,ii,jj,kk), & 
                                         acc_siga(acc_skip+1:cyc,ii,jj,kk), &
                                         acc_nufi(acc_skip+1:cyc,ii,jj,kk), &
                                         courn(ii,jj,kk,:,:))

            call AVG_N_STD(coxs(ii,jj,kk,:,1), & 
                           coxs(ii,jj,kk,:,2), &
                           coxs(ii,jj,kk,:,3),dble(n_pert))

            do i = 1,types
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) - avg_(i)
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) * ostd(i) / std_(i)
                coxs(ii,jj,kk,:,i) = coxs(ii,jj,kk,:,i) + oavg(i)
            enddo
            
            call AVG_N_STD(coxs(ii,jj,kk,:,1), & 
                           coxs(ii,jj,kk,:,2), &
                           coxs(ii,jj,kk,:,3),dble(n_pert))

            call CORRELATION_MATRIX(coxs(ii,jj,kk,:,1), & 
                                    coxs(ii,jj,kk,:,2), &
                                    coxs(ii,jj,kk,:,3), dble(n_pert))
            if(cyc==131 .and. iscore) print *, 'PRVCOS', it, prev
            if(cyc==131 .and. iscore) print *, 'AFPCOS', it, cor
            enddo
    
        end do
        end do
        end do
    end subroutine


    subroutine DATA_INO(cyc)
        implicit none
        integer:: cyc

!        open(1,file='accxs.out')
!        write(1,1), acc_sigt(acc_skip+1:cyc,:,:,:)
!        write(1,1), acc_siga(acc_skip+1:cyc,:,:,:)
!        write(1,1), acc_nufi(acc_skip+1:cyc,:,:,:)
!        write(1,1), acc_phi(acc_skip+1:cyc,:,:,:)
!        close(1)
        open(1,file='accxs.out')
        read(1,1), acc_sigt(acc_skip+1:cyc,:,:,:)
        read(1,1), acc_siga(acc_skip+1:cyc,:,:,:)
        read(1,1), acc_nufi(acc_skip+1:cyc,:,:,:)
        read(1,1), acc_phi(acc_skip+1:cyc,:,:,:)
        close(1)

        1 format(<nfm(1)>ES16.8)

    end subroutine


    ! =========================================================================
    ! 
    ! =========================================================================
    subroutine PEARSON_CORRELATION(xst,xsa,xsf)
        implicit none
        real(8), intent(in):: xst(:), xsa(:), xsf(:)
        real(8), allocatable:: vma(:,:) ! variable minus average

        cyclen = size(xst)
        avg_(1) = sum(xst(:))/cyclen
        avg_(2) = sum(xsa(:))/cyclen
        avg_(3) = sum(xsf(:))/cyclen

        allocate(vma(int(cyclen),3))
        vma(:,1) = xst-avg_(1)
        vma(:,2) = xsa-avg_(2)
        vma(:,3) = xsf-avg_(3)

        do mm = 1, 3
        do nn = mm+1, 3, 1
        cor(mm,nn) = sum(vma(:,mm)*vma(:,nn)) &
            /(sqrt(sum(vma(:,mm)*vma(:,mm)))*sqrt(sum(vma(:,nn)*vma(:,nn))))
        end do
        end do

        cor(1,1) = 1
        cor(2,2) = 1
        cor(3,3) = 1
        cor(2,1) = cor(1,2)
        cor(3,1) = cor(1,3)
        cor(3,2) = cor(2,3)
        !print*, cor

        deallocate(vma)

    end subroutine
    
    subroutine RANDOM_CORR
        use RANDOMS
        implicit none

        cor = 1
        cor(1,2) = 2*rang()-1
        cor(1,3) = 2*rang()-1
        cor(2,3) = 2*rang()-1
        cor(2,1) = cor(1,2)
        cor(3,1) = cor(1,3)
        cor(3,2) = cor(2,3)

    end subroutine

    ! =========================================================================
    ! AVG_N_STD
    ! =========================================================================
    subroutine AVG_N_STD(xst,xsa,xsf, cyclen)
        implicit none
        real(8), intent(in):: xst(:), xsa(:), xsf(:), cyclen

        avg_(1) = sum(xst(:))/cyclen
        avg_(2) = sum(xsa(:))/cyclen
        avg_(3) = sum(xsf(:))/cyclen

        std_(1) = sqrt(sum((xst(:)-avg_(1))*(xst(:)-avg_(1)))/(cyclen))
        std_(2) = sqrt(sum((xsa(:)-avg_(2))*(xsa(:)-avg_(2)))/(cyclen))
        std_(3) = sqrt(sum((xsf(:)-avg_(3))*(xsf(:)-avg_(3)))/(cyclen))

!        if ( std_(3) == 0 ) then
!            types = 2
!        else
!            types = 3
!        end if

    end subroutine


    ! =========================================================================
    ! CORRELATION_MATRIX
    ! =========================================================================
    subroutine CORRELATION_MATRIX(xst,xsa,xsf,cyclen)
        implicit none
        real(8), intent(in):: xst(:), xsa(:), xsf(:), cyclen


        corr(1,2) = sum(xst*xsa)/cyclen
        corr(1,3) = sum(xst*xsf)/cyclen
        corr(2,3) = sum(xsa*xsf)/cyclen

        ! correlation matrix
        do mm = 1, types
        do nn = mm+1, types
        cor(mm,nn) = (corr(mm,nn)-avg_(mm)*avg_(nn))/(std_(mm)*std_(nn))
        cor(nn,mm) = cor(mm,nn)
        end do
        end do
        cor(1,1) = 1
        cor(2,2) = 1
        cor(3,3) = 1

    end subroutine

    function CORR_VAL(x1,x2)
        implicit none
        real(8), intent(in):: x1(:), x2(:)
        real(8) :: CORR_VAL
        real(8) :: c, a1, a2, s1, s2
        integer :: n
        
        n = size(x1)
        if(n/=size(x2)) then
            write(*,*) '?', n, size(x1), size(x2)
            CORR_VAL = 0D0
            return
        endif
        a1= sum(x1)/dble(n); a2 = sum(x2)/dble(n)
        s1= sqrt(sum((x1-a1)**2)/dble(n))
        s2= sqrt(sum((x2-a2)**2)/dble(n))

        c = sum(x1*x2)/dble(n)
        CORR_VAL = (c-a1*a2)/(s1*s2)

    end function

    ! =========================================================================
    ! LHS_SAPMLING samples random number from the Latin hypercube sampling
    ! =========================================================================
    subroutine LHS_SAMPLING
        use RANDOMS, only: rang
        implicit none
        integer:: pos
        real(8):: temp

!        do mm = 1, n_pert
!            do nn = 1, types
!                lhs(mm,nn) = nn-rang()
!            end do
!            do nn = types, 2, -1
!                pos = int(nn*rang())+1
!                temp = lhs(mm,pos)
!                lhs(mm,pos) = lhs(mm,nn)
!                lhs(mm,nn) = temp
!            end do
!        end do
!
!        lhs = lhs / dble(types)
!!        do mm = 1, n_pert
!!            print*, lhs(mm,:)
!!        end do
!!        stop

        do mm = 1, n_pert
            do nn = 1, types
                lhs(mm,nn) = mm/dble(n_pert) -rang()/dble(n_pert)
            enddo
        enddo

        do nn = 1, types
            do mm = 1, n_pert
                pos = int(mm*rang()) + 1
                temp= lhs(pos,nn)
                lhs(pos,nn) = lhs(mm,nn)
                lhs(mm,nn)  = temp
            enddo
        enddo

    end subroutine


    ! =========================================================================
    ! ICDF returns the value from the inverse CDF
    ! =========================================================================
    function ICDF_NORMAL(urn) result(lhs_)
        implicit none
        real(8), intent(in):: urn(:)
        real(8):: lhs_(3)
        real(8):: err(3), inverr(3)

!        do nn = 1, types
!            if ( urn(nn) < 0.5 ) then
!                lhs_(nn) = log(2D0*urn(nn))
!            else
!                lhs_(nn) = -log(2D0-2D0*urn(nn))
!            end if
!        end do

        err(1:types) = 2D0*urn(1:types)-1D0
        call VDERFINV(1,err(1),inverr(1))
        call VDERFINV(1,err(2),inverr(2))
        call VDERFINV(1,err(3),inverr(3))
        lhs_(1:types) = sqrt2*inverr(1:types)

    end function


    ! =========================================================================
    ! DECHOLESKY
    ! =========================================================================
    subroutine DECHOLESKY(cor, chol)
        implicit none
        real(8), intent(inout)  :: cor(3,3)
        real(8), intent(out) :: chol(3,3)
        real(8):: detm  ! determinant

        do
        chol = 0
        do mm = 1, types
            chol(mm,mm) = sqrt(cor(mm,mm)-sum(chol(1:mm-1,mm)*chol(1:mm-1,mm)))
            do nn = mm+1, types
            chol(mm,nn) = (cor(nn,mm)-sum(chol(1:mm-1,mm)*chol(1:mm-1,nn)))/chol(mm,mm)
            end do
        end do

        detm = cor(1,1)*(cor(2,2)*cor(3,3)-cor(2,3)*cor(3,2)) &
              -cor(1,2)*(cor(2,1)*cor(3,3)-cor(2,3)*cor(3,1)) &
              +cor(1,3)*(cor(2,1)*cor(3,2)-cor(2,2)*cor(3,1))
        
        if ( detm < 0 ) then
            do mm = 1, 3
                cor(mm,mm) = cor(mm,mm) + 1D-4
            end do
        else
            exit
        end if
        end do

    end subroutine


    ! =========================================================================
    ! 
    ! =========================================================================
    function CDF_COURN(urn) result(cdf)
        implicit none
        real(8), intent(in):: urn(:)
        real(8):: cdf(3)

        cdf(1:types) = 5D-1*(1D0+erf(urn(1:types)/sqrt2))

    end function

        
    ! =========================================================================
    ! 
    ! =========================================================================
    function ICDF_XS(xst,xsa,xsf,urnxs) result(cxs)
        use RANDOMS, only: rang
        implicit none
        real(8), intent(in):: xst(:), xsa(:), xsf(:)
        real(8), intent(in):: urnxs(:,:)
        real(8):: cxs(n_pert,3)
        real(8):: v_min, v_max, v_diff
        real(8):: val(0:10)
        real(8):: cdf(0:10)
        integer:: steps = 1D1
        integer:: xx, yy
        real(8):: urn


        ! ===== total cross section =====
        ! minimum & maximum
        v_min = minval(xst(:))
        v_max = maxval(xst(:))
        v_diff = (v_max-v_min)/real(steps-1)
        val(0) = v_min-v_diff*5D-1
        do xx = 1, steps
            val(xx) = val(xx-1)+v_diff
        end do
    
        ! cumulative density function (CDF)
        cdf = 0
        do xx = 1, cyclen
        do yy = 1, steps
        if ( xst(xx) < val(yy) ) then
            cdf(yy:) = cdf(yy:) + 1
            exit
        end if
        end do
        end do
        cdf = cdf / cdf(steps)
    
        ! sampling
        do yy = 1, n_pert
!        urn = rang()
        do xx = 1, steps
        if ( urnxs(yy,1) < cdf(xx) ) then
            cxs(yy,1) = val(xx-1)+v_diff*(urnxs(yy,1)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
            exit
        end if
!        if ( urn < cdf(xx) ) then
!            cxs(yy,1) = val(xx-1)+v_diff*(urn-cdf(xx-1)) &
!                        /(cdf(xx)-cdf(xx-1))
!            exit
!        end if
        end do
        end do
    

        ! ===== absorption cross section =====
        ! minimum & maximum
        v_min = minval(xsa(:))
        v_max = maxval(xsa(:))
        v_diff = (v_max-v_min)/real(steps-1)
        val(0) = v_min-v_diff*5D-1
        do xx = 1, steps
            val(xx) = val(xx-1)+v_diff
        end do
    
        ! cumulative density function (CDF)
        cdf = 0
        do xx = 1, cyclen
        do yy = 1, steps
        if ( xsa(xx) < val(yy) ) then
            cdf(yy:) = cdf(yy:) + 1
            exit
        end if
        end do
        end do
        cdf = cdf / cdf(steps)
    
        ! sampling
        do yy = 1, n_pert
!        urn = rang()
        do xx = 1, steps
        if ( urnxs(yy,2) < cdf(xx) ) then
            cxs(yy,2) = val(xx-1)+v_diff*(urnxs(yy,2)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
            exit
        end if
!        if ( urn < cdf(xx) ) then
!            cxs(yy,2) = val(xx-1)+v_diff*(urn-cdf(xx-1)) &
!                        /(cdf(xx)-cdf(xx-1))
!            exit
!        end if
        end do
        end do
    
        !if ( xsf(1) == 0 ) then
        if ( types == 2 ) then
            cxs(:,3) = 0
            return
        end if

        ! ===== nu fission cross section =====
        ! minimum & maximum
        v_min = minval(xsf(:))
        v_max = maxval(xsf(:))
        v_diff = (v_max-v_min)/real(steps-1)
        val(0) = v_min-v_diff*5D-1
        do xx = 1, steps
            val(xx) = val(xx-1)+v_diff
        end do
    
        ! cumulative density function (CDF)
        cdf = 0
        do xx = 1, cyclen
        do yy = 1, steps
        if ( xsf(xx) < val(yy) ) then
            cdf(yy:) = cdf(yy:) + 1
            exit
        end if
        end do
        end do
        cdf = cdf / cdf(steps)

        ! sampling
        do yy = 1, n_pert
!        urn = rang()
        do xx = 1, steps
        if ( urnxs(yy,3) < cdf(xx) ) then
            cxs(yy,3) = val(xx-1)+v_diff*(urnxs(yy,3)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
            exit
        end if
!        if ( urn < cdf(xx) ) then
!            cxs(yy,3) = val(xx-1)+v_diff*(urn-cdf(xx-1)) &
!                        /(cdf(xx)-cdf(xx-1))
!            exit
!        end if
        end do
        end do

    end function

        
    ! =========================================================================
    ! 
    ! =========================================================================
    function ICDF_J(ptJ0) result(ptJ1)
        use RANDOMS, only: rang
        implicit none
        real(8), intent(in):: ptJ0(:,:)
        real(8):: ptJ1(n_pert,6)
        real(8):: v_min, v_max, v_diff
        real(8):: val(0:10)
        real(8):: cdf(0:10)
        integer:: steps = 1D1
        integer:: xx, yy, ss
        real(8):: urn


        do ss = 1, 6
        ! minimum & maximum
        v_min = minval(ptJ0(:,ss))
        v_max = maxval(ptJ0(:,ss))
        v_diff = (v_max-v_min)/real(steps-1)
        val(0) = v_min-v_diff*5D-1
        do xx = 1, steps
            val(xx) = val(xx-1)+v_diff
        end do
    
        ! cumulative density function (CDF)
        cdf = 0
        do xx = 1, cyclen
        do yy = 1, steps
        if ( ptJ0(xx,ss) < val(yy) ) then
            cdf(yy:) = cdf(yy:) + 1
            exit
        end if
        end do
        end do
        cdf = cdf / cdf(steps)
    
        ! sampling
        do yy = 1, n_pert
        urn = rang()
        do xx = 1, steps
        if ( urn < cdf(xx) ) then
            ptJ1(yy,ss) = val(xx-1)+v_diff*(urn-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
            exit
        end if
        end do
        end do
        end do

    end function

        
    ! =========================================================================
    ! 
    ! =========================================================================
    function ICDF_PHI(ptphi0) result(ptphi1)
        use RANDOMS, only: rang
        implicit none
        real(8), intent(in):: ptphi0(:)
        real(8):: ptphi1(n_pert)
        real(8):: v_min, v_max, v_diff
        real(8):: val(0:10)
        real(8):: cdf(0:10)
        integer:: steps = 1D1
        integer:: xx, yy
        real(8):: urn


        ! minimum & maximum
        v_min = minval(ptphi0(:))
        v_max = maxval(ptphi0(:))
        v_diff = (v_max-v_min)/real(steps-1)
        val(0) = v_min-v_diff*5D-1
        do xx = 1, steps
            val(xx) = val(xx-1)+v_diff
        end do
    
        ! cumulative density function (CDF)
        cdf = 0
        do xx = 1, cyclen
        do yy = 1, steps
        if ( ptphi0(xx) < val(yy) ) then
            cdf(yy:) = cdf(yy:) + 1
            exit
        end if
        end do
        end do
        cdf = cdf / cdf(steps)
    
        ! sampling
        do yy = 1, n_pert
        urn = rang()
        do xx = 1, steps
        if ( urn < cdf(xx) ) then
            ptphi1(yy) = val(xx-1)+v_diff*(urn-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
            exit
        end if
        end do
        end do

    end function


    ! #########################################################################

    ! =========================================================================
    ! 
    ! =========================================================================
    subroutine COSAMPLING_NODE(cyc)
        use VARIABLES, only: n_inact
        implicit none 
        integer, intent(in):: cyc
        integer:: ix0, ix1, iy0, iy1, iz0, iz1
        integer:: idm, idn, ido, id0, id1, id2

        cyclen = dble(cyc-acc_skip)
        if ( cyc <= n_inact ) then
            cyclen = dble(cyc-acc_skip)
        else
            cyclen = dble(n_inact-acc_skip)
        end if
        ! parameter allocation
        if ( .not. allocated(avg_) ) then
            n_data = n_lnodes*3
            allocate(avg_(n_data),std_(n_data),corr(n_data,n_data))
            allocate(cor(n_data,n_data),chol(n_data,n_data))
            allocate(lhs(n_pert,n_data))
            allocate(courn(nfm(1),nfm(2),nfm(3),n_pert,3))
            allocate(coxs(nfm(1),nfm(2),nfm(3),n_pert,3))
        end if

        do ii = 1, ncm(1); ix1 = ii*fcr; ix0 = ix1-fcr+1
        do jj = 1, ncm(2); iy1 = jj*fcr; iy0 = iy1-fcr+1
        do kk = 1, ncm(3); iz1 = kk*fcz; iz0 = iz1-fcz+1

            if ( OUT_ZZ(ii,jj) ) cycle

!            print*, ix0, ix1, iy0
!            print*, n_data

            ! correlation matrix
            call NODE_AVG_STD( &
                acc_sigt(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1), &
                acc_siga(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1), &
                acc_nufi(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1))
!            write(1,*), "average & SD"
!            write(1,1), avg_
!            write(1,1), std_
!            write(1,*)
            call NODE_COR_MATRIX( &
                acc_sigt(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1), &
                acc_siga(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1), &
                acc_nufi(cyc-cyclen+1:cyc,ix0:ix1,iy0:iy1,iz0:iz1))

!            write(1,*), "correlation"
!            write(1,2), cor
!            write(1,*)

            ! random number generation by LHS
            call LHS_SAMPLING_NODE

!            write(1,*), "LHS"
!            write(1,2), lhs

            ! inverse normal CDF
            lhs(:,:) = ICDF_NORMAL_NODE(lhs(:,:))

!            write(1,*), "inverse LHS"
!            write(1,2), lhs

            ! Cholesky decomposition
            call DECHOLESKY_NODE

!            write(1,*), "cholesky"
!            write(1,2), chol

            ! correlated uniform random number
!            do mm = 1, fcr; id0 = ix0+mm-1
!            do nn = 1, fcr; id1 = iy0+nn-1
!            do oo = 1, fcz; id2 = iz0+oo-1
!
!            end do
!            end do
!            end do
            lhs = matmul(lhs,chol)
            ! --- conversion
            do mm = 1, fcr; id0 = ix0+mm-1
            do nn = 1, fcr; id1 = iy0+nn-1
            do oo = 1, fcz; id2 = iz0+oo-1
                idm = mm+(nn-1)*fcr+(oo-1)*fcr*fcr
                idn = n_lnodes+idm
                ido = n_lnodes+idn
                courn(id0,id1,id2,:,1) = lhs(:,idm)
                courn(id0,id1,id2,:,2) = lhs(:,idn)
                courn(id0,id1,id2,:,3) = lhs(:,ido)
            end do
            end do
            end do

!            write(1,*), "correlated inverse urn"
!            write(1,1), courn(1:8,1:8,1,1:3,1)
!            write(1,1), courn(1:8,1:8,1,1:3,2)
!            write(1,1), courn(1:8,1:8,1,1:3,3)

            ! CFD for correlated urn generation
!            courn(ix0:ix1,iy0:iy1,iz0:iz1,:,:) = &
!                CDF_COURN_NODE(courn(ix0:ix1,iy0:iy1,iz0:iz1,:,:))
            courn(ix0:ix1,iy0:iy1,iz0:iz1,:,:) = &
                5D-1*(1D0+erf(courn(ix0:ix1,iy0:iy1,iz0:iz1,:,:)/sqrt2))

!            write(1,*), "correlated urn"
!            write(1,1), courn(1:8,1:8,1,1:3,1)
!            write(1,1), courn(1:8,1:8,1,1:3,2)
!            write(1,1), courn(1:8,1:8,1,1:3,3)

            ! inverse CDF for XS
!            print*, "1"
            do mm = 1, fcr; id0 = ix0+mm-1
            do nn = 1, fcr; id1 = iy0+nn-1
            do oo = 1, fcz; id2 = iz0+oo-1
                coxs(id0,id1,id2,:,:) = ICDF_XS( &
                    acc_sigt(cyc-cyclen+1:cyc,id0,id1,id2), & 
                    acc_siga(cyc-cyclen+1:cyc,id0,id1,id2), & 
                    acc_nufi(cyc-cyclen+1:cyc,id0,id1,id2), & 
                    courn(id0,id1,id2,:,:))
!                    print*, coxs(id0,id1,id2,1,:)
            end do
            end do
            end do

!            write(1,*), "cross section"
!            write(1,1), coxs(1:8,1:8,1,1:3,1)
!            write(1,1), coxs(1:8,1:8,1,1:3,2)
!            write(1,1), coxs(1:8,1:8,1,1:3,3)
!            pause
!            stop

        end do
        end do
        end do
        1 format(1000es15.7)
        2 format(<n_data>es15.7)

    end subroutine


    ! =============================================================================
    ! OUT_OF_ZZ determines if a region is in or out of the zigzag boundary region
    ! =============================================================================
    function OUT_ZZ(io,jo)
        logical:: OUT_ZZ
        integer, intent(in):: io, jo
        integer:: mo, no
        
        if ( .not. zigzagon ) then
            OUT_ZZ = .false.
            return
        end if
    
        do mo = 1, zz_div
        if ( zzc0(mo) < io .and. io <= zzc0(mo+1) ) then
            no = mo
            exit
        end if
        end do
    
        if ( zzc1(no) < jo .and. jo <= zzc2(no) ) then
            OUT_ZZ = .false.
        else
            OUT_ZZ = .true.
        end if
    
    end function


    ! =========================================================================
    ! 
    ! =========================================================================
    subroutine NODE_AVG_STD(xst,xsa,xsf)
        implicit none
        real(8), intent(in):: xst(:,:,:,:), xsa(:,:,:,:), xsf(:,:,:,:)
        integer:: id0, id1, id2

        do oo = 1, fcz
        do nn = 1, fcr
        do mm = 1, fcr

            id0 = mm+fcr*(nn-1)+fcr*fcr*(oo-1)
            id1 = n_lnodes+id0
            id2 = n_lnodes+id1

            avg_(id0) = sum(xst(:,mm,nn,oo))/dble(cyclen)
            avg_(id1) = sum(xsa(:,mm,nn,oo))/dble(cyclen)
            avg_(id2) = sum(xsf(:,mm,nn,oo))/dble(cyclen)

            std_(id0) = sqrt(sum((xst(:,mm,nn,oo)-avg_(id0)) &
                      * (xst(:,mm,nn,oo)-avg_(id0)))/dble(cyclen))
            std_(id1) = sqrt(sum((xsa(:,mm,nn,oo)-avg_(id1)) &
                      * (xsa(:,mm,nn,oo)-avg_(id1)))/dble(cyclen))
            std_(id2) = sqrt(sum((xsf(:,mm,nn,oo)-avg_(id2)) &
                      * (xsf(:,mm,nn,oo)-avg_(id2)))/dble(cyclen))

        end do
        end do
        end do

    end subroutine



    ! =========================================================================
    ! CORRELATION_MATRIX
    ! =========================================================================
    subroutine NODE_COR_MATRIX(xst,xsa,xsf)
        implicit none
        real(8), intent(in):: xst(:,:,:,:), xsa(:,:,:,:), xsf(:,:,:,:)
        integer:: xx, yy, zz
        integer:: idm, idn, ido
        integer:: idx, idy, idz

        ! sum product
        do oo = 1, fcz
        do nn = 1, fcr
        do mm = 1, fcr

            idm = mm+fcr*(nn-1)+fcr*fcr*(oo-1)
            idn = n_lnodes+idm
            ido = n_lnodes+idn

            do zz = 1, fcz
            do yy = 1, fcr
            do xx = 1, fcr
    
                idx = xx+fcr*(yy-1)+fcr*fcr*(zz-1)
                idy = n_lnodes+idx
                idz = n_lnodes+idy

                ! ---
                if ( idm <= idx ) corr(idm,idx) = &
                    sum(xst(:,mm,nn,oo)*xst(:,xx,yy,zz))/dble(cyclen)
                if ( idm < idy ) corr(idm,idy) = &
                    sum(xst(:,mm,nn,oo)*xsa(:,xx,yy,zz))/dble(cyclen)
                if ( idm < idz ) corr(idm,idz) = &
                    sum(xst(:,mm,nn,oo)*xsf(:,xx,yy,zz))/dble(cyclen)
                ! ---
                if ( idn <= idy ) corr(idn,idy) = &
                    sum(xsa(:,mm,nn,oo)*xsa(:,xx,yy,zz))/dble(cyclen)
                if ( idn < idz ) corr(idn,idz) = &
                    sum(xsa(:,mm,nn,oo)*xsf(:,xx,yy,zz))/dble(cyclen)
                ! ---
                if ( ido <= idz ) corr(ido,idz) = &
                    sum(xsf(:,mm,nn,oo)*xsf(:,xx,yy,zz))/dble(cyclen)
    
            end do
            end do
            end do

        end do
        end do
        end do


        ! correlation matrix
        cor = 0
        do mm = 1, n_data
            !cor(mm,mm) = 1
        !do nn = mm+1, n_data
        do nn = mm, n_data
            if ( std_(mm)*std_(nn) == 0 ) then
            if ( mm == nn ) then
            cor(mm,nn) = 1
            else
            cor(mm,nn) = 0
            end if
            else
            cor(mm,nn) = (corr(mm,nn)-avg_(mm)*avg_(nn))/(std_(mm)*std_(nn))
            end if
            cor(nn,mm) = cor(mm,nn)
        end do
        end do


!        ! correlation matrix
!        do oo = 1, fcz
!        do nn = 1, fcr
!        do mm = 1, fcr
!
!            idm = mm+fcr*(nn-1)+fcr*fcr*(oo-1)
!            idn = n_lnodes+idm
!            ido = n_lnodes+idn
!
!            do zz = 1, fcz
!            do yy = 1, fcr
!            do xx = 1, fcr
!    
!                idx = xx+fcr*(yy-1)+fcr*fcr*(zz-1)
!                idy = n_lnodes+idx
!                idz = n_lnodes+idy
!
!                ! ---
!                if ( idm <= idx ) then; cor(idm,idx) = (sum(xst(:,mm,nn,oo) &
!                    *xst(:,xx,yy,zz))-cyclen*avg_(idm)*avg_(idx)) &
!                    /(sqrt(sum(xst(:,mm,nn,oo)*xst(:,mm,nn,oo))-cyclen &
!                    *avg_(idm)*avg_(idm))*sqrt(sum(xst(:,xx,yy,zz) &
!                    *xst(:,xx,yy,zz))-cyclen*avg_(idx)*avg_(idx)))
!                    cor(idx,idm) = cor(idm,idx)
!                end if
!                if ( idm <  idy ) then; cor(idm,idy) = (sum(xst(:,mm,nn,oo) &
!                    *xsa(:,xx,yy,zz))-cyclen*avg_(idm)*avg_(idy)) &
!                    /(sqrt(sum(xst(:,mm,nn,oo)*xst(:,mm,nn,oo))-cyclen &
!                    *avg_(idm)*avg_(idm))*sqrt(sum(xsa(:,xx,yy,zz) &
!                    *xsa(:,xx,yy,zz))-cyclen*avg_(idy)*avg_(idy)))
!                    cor(idy,idm) = cor(idm,idy)
!                end if
!                if ( idm <  idz ) then; cor(idm,idz) = (sum(xst(:,mm,nn,oo) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(idm)*avg_(idz)) &
!                    /(sqrt(sum(xst(:,mm,nn,oo)*xst(:,mm,nn,oo))-cyclen &
!                    *avg_(idm)*avg_(idm))*sqrt(sum(xsf(:,xx,yy,zz) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(idz)*avg_(idz)))
!                    cor(idz,idm) = cor(idm,idz)
!                end if
!                ! ---
!                if ( idn <= idy ) then; cor(idn,idy) = (sum(xsa(:,mm,nn,oo) &
!                    *xsa(:,xx,yy,zz))-cyclen*avg_(idn)*avg_(idy)) &
!                    /(sqrt(sum(xsa(:,mm,nn,oo)*xsa(:,mm,nn,oo))-cyclen &
!                    *avg_(idn)*avg_(idn))*sqrt(sum(xsa(:,xx,yy,zz) &
!                    *xsa(:,xx,yy,zz))-cyclen*avg_(idy)*avg_(idy)))
!                    cor(idy,idn) = cor(idn,idy)
!                end if
!                if ( idn <  idz ) then; cor(idn,idz) = (sum(xsa(:,mm,nn,oo) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(idn)*avg_(idz)) &
!                    /(sqrt(sum(xsa(:,mm,nn,oo)*xsa(:,mm,nn,oo))-cyclen &
!                    *avg_(idn)*avg_(idn))*sqrt(sum(xsf(:,xx,yy,zz) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(idz)*avg_(idz)))
!                    cor(idz,idn) = cor(idn,idz)
!                end if
!                ! ---
!                if ( ido <= idz ) then; cor(ido,idz) = (sum(xsf(:,mm,nn,oo) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(ido)*avg_(idz)) &
!                    /(sqrt(sum(xsf(:,mm,nn,oo)*xsf(:,mm,nn,oo))-cyclen &
!                    *avg_(ido)*avg_(ido))*sqrt(sum(xsf(:,xx,yy,zz) &
!                    *xsf(:,xx,yy,zz))-cyclen*avg_(idz)*avg_(idz)))
!                    cor(idz,ido) = cor(ido,idz)
!                end if
!    
!            end do
!            end do
!            end do
!
!        end do
!        end do
!        end do
!
!        where ( isnan(cor) ) cor = 0
!
!        write(2,2), cor
!        write(2,*)
!        pause
!        stop

    end subroutine



    ! =========================================================================
    ! LHS_SAPMLING samples random number from the Latin hypercube sampling
    ! =========================================================================
    subroutine LHS_SAMPLING_NODE
        use RANDOMS, only: rang
        implicit none
        integer:: pos
        real(8):: temp

        do mm = 1, n_pert
            do nn = 1, n_data
                lhs(mm,nn) = nn-rang()
            end do
            do nn = n_data, 2, -1
                pos = int(nn*rang())+1
                temp = lhs(mm,pos)
                lhs(mm,pos) = lhs(mm,nn)
                lhs(mm,nn) = temp
            end do
        end do

        lhs = lhs / dble(n_data)

    end subroutine


    ! =========================================================================
    ! ICDF returns the value from the inverse CDF
    ! =========================================================================
    function ICDF_NORMAL_NODE(urn) result(lhs_)
        implicit none
        real(8), intent(in):: urn(:,:)
        real(8):: lhs_(n_pert,n_data)
        real(8):: err(1), inverr(1)

        do mm = 1, n_pert
        do nn = 1, n_data
            err(1) = 2D0*urn(mm,nn)-1D0
            call VDERFINV(1,err(1),inverr(1))
            lhs_(mm,nn) = sqrt2*inverr(1)
        end do
        end do

    end function


    ! =========================================================================
    ! DECHOLESKY
    ! =========================================================================
    subroutine DECHOLESKY_NODE
        implicit none
        logical:: cor_nan

        do

        chol = 0; cor_nan = .false.
        do mm = 1, n_data
            chol(mm,mm) = sqrt(cor(mm,mm)-sum(chol(mm,1:mm-1)*chol(mm,1:mm-1)))
            if ( isnan(chol(mm,mm)) ) then
                cor_nan = .true.
                exit
            end if
        end do

        if ( cor_nan ) then
            print*, "NAN CHOLESKY"
            print*, cor(1,1)
            do mm = 1, n_data
                cor(mm,mm) = cor(mm,mm) + 1D-5
            end do
            cycle
        end if

        do mm = 1, n_data
        do nn = mm+1, n_data
        chol(nn,mm) = (cor(nn,mm)-sum(chol(mm,1:mm-1)*chol(nn,1:mm-1)))/chol(mm,mm)
        end do
        end do
        exit

        end do

    end subroutine


    ! =========================================================================
    ! 
    ! =========================================================================
    function CDF_COURN_NODE(urn) result(cdf)
        implicit none
        real(8), intent(in):: urn(:,:,:,:,:)
        real(8):: cdf(fcr,fcr,fcz,n_pert,3)

        cdf = 5D-1*(1D0+erf(urn/sqrt2))

    end function
	  pure function matinv2(A) result(B)
		!! Performs a direct calculation of the inverse of a 2×2 matrix.
		real(8), intent(in) :: A(2,2)   !! Matrix
		real(8)             :: B(2,2)   !! Inverse matrix
		real(8)             :: detinv

		! Calculate the inverse determinant of the matrix
		detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

		! Calculate the inverse of the matrix
		B(1,1) = +detinv * A(2,2)
		B(2,1) = -detinv * A(2,1)
		B(1,2) = -detinv * A(1,2)
		B(2,2) = +detinv * A(1,1)
	  end function        
	  pure function matinv3(A) result(B)
		!! Performs a direct calculation of the inverse of a 3×3 matrix.
		real(8), intent(in) :: A(3,3)   !! Matrix
		real(8)             :: B(3,3)   !! Inverse matrix
		real(8)             :: detinv

		! Calculate the inverse determinant of the matrix
		detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
				  - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
				  + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

		! Calculate the inverse of the matrix
		B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
		B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
		B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
		B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
		B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
		B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
		B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
		B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
		B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
	  end function

!    ! =========================================================================
!    ! 
!    ! =========================================================================
!    function ICDF_XS_NODE1(xst,urnxs) result(cxs)
!        implicit none
!        real(8), intent(in):: xst(:)
!        real(8), intent(in):: urnxs(:)
!        real(8):: cxs(n_pert)
!        real(8):: v_min, v_max, v_diff
!        real(8):: val(0:10)
!        real(8):: cdf(0:10)
!        integer:: steps = 1D1
!        integer:: xx, yy
!
!
!        ! ===== total cross section =====
!        ! minimum & maximum
!        v_min = minval(xst(:))
!        v_max = maxval(xst(:))
!        v_diff = (v_max-v_min)/real(steps-1)
!        val(0) = v_min-v_diff*5D-1
!        do xx = 1, steps
!            val(xx) = val(xx-1)+v_diff
!        end do
!    
!        ! cumulative density function (CDF)
!        cdf = 0
!        do xx = 1, cyclen
!        do nn = 1, steps
!        if ( xst(xx) < val(nn) ) then
!            cdf(nn:) = cdf(nn:) + 1
!            exit
!        end if
!        end do
!        end do
!        cdf = cdf / cdf(steps)
!    
!        ! sampling
!        do yy = 1, n_pert
!        do xx = 1, steps
!        if ( urnxs(yy) < cdf(xx) ) then
!            cxs(yy) = val(xx-1)+v_diff*(urnxs(yy)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
!            exit
!        end if
!        end do
!        end do
!
!    end function
!
!
!    ! =========================================================================
!    ! 
!    ! =========================================================================
!    function ICDF_XS_NODE2(xsa,urnxs) result(cxs)
!        implicit none
!        real(8), intent(in):: xsa(:)
!        real(8), intent(in):: urnxs(:)
!        real(8):: cxs(n_pert)
!        real(8):: v_min, v_max, v_diff
!        real(8):: val(0:10)
!        real(8):: cdf(0:10)
!        integer:: steps = 1D1
!        integer:: xx, yy
!
!
!        ! ===== absorption cross section =====
!        ! minimum & maximum
!        v_min = minval(xsa(:))
!        v_max = maxval(xsa(:))
!        v_diff = (v_max-v_min)/real(steps-1)
!        val(0) = v_min-v_diff*5D-1
!        do xx = 1, steps
!            val(xx) = val(xx-1)+v_diff
!        end do
!    
!        ! cumulative density function (CDF)
!        cdf = 0
!        do xx = 1, cyclen
!        do nn = 1, steps
!        if ( xsa(xx) < val(nn) ) then
!            cdf(nn:) = cdf(nn:) + 1
!            exit
!        end if
!        end do
!        end do
!        cdf = cdf / cdf(steps)
!    
!        ! sampling
!        do yy = 1, n_pert
!        do xx = 1, steps
!        if ( urnxs(yy) < cdf(xx) ) then
!            cxs(yy) = val(xx-1)+v_diff*(urnxs(yy)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
!            exit
!        end if
!        end do
!        end do
!
!    end function
!
!
!    ! =========================================================================
!    ! 
!    ! =========================================================================
!    function ICDF_XS_NODE3(xsf,urnxs) result(cxs)
!        implicit none
!        real(8), intent(in):: xsf(:)
!        real(8), intent(in):: urnxs(:)
!        real(8):: cxs(n_pert)
!        real(8):: v_min, v_max, v_diff
!        real(8):: val(0:10)
!        real(8):: cdf(0:10)
!        integer:: steps = 1D1
!        integer:: xx, yy
!
!
!        ! ===== nu X fission cross section =====
!        ! minimum & maximum
!        v_min = minval(xsf(:))
!        v_max = maxval(xsf(:))
!        v_diff = (v_max-v_min)/real(steps-1)
!        val(0) = v_min-v_diff*5D-1
!        do xx = 1, steps
!            val(xx) = val(xx-1)+v_diff
!        end do
!    
!        ! cumulative density function (CDF)
!        cdf = 0
!        do xx = 1, cyclen
!        do nn = 1, steps
!        if ( xsf(xx) < val(nn) ) then
!            cdf(nn:) = cdf(nn:) + 1
!            exit
!        end if
!        end do
!        end do
!        cdf = cdf / cdf(steps)
!    
!        ! sampling
!        do yy = 1, n_pert
!        do xx = 1, steps
!        if ( urnxs(yy) < cdf(xx) ) then
!            cxs(yy) = val(xx-1)+v_diff*(urnxs(yy)-cdf(xx-1))/(cdf(xx)-cdf(xx-1))
!            exit
!        end if
!        end do
!        end do
!
!    end function


end module
