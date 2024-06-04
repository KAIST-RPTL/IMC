module STATISTICS
    use VARIABLES, only: n_batch, n_totcyc, n_act, t_totcyc
    implicit none
    
    contains

    function AVG(val)
        real(8):: avg
        real(8), intent(in):: val(:)
    
        avg = sum(val)/size(val)
    
    end function
    
    function PCM(val)
        real(8):: pcm
        real(8), intent(in):: val
    
        pcm = val*1E5
    
    end function
    
    function STD_M(val) ! STD of the sample mean
        real(8):: std_m
        real(8), intent(in):: val(:)
        integer:: length
        real(8):: avg
    
        length = size(val)
        avg = sum(val)/length
        std_m = sqrt(dot_product((val-avg),(val-avg))/(length*(length-1)))
        if ( isnan(std_m) ) std_m = 0
    
    end function
    
    function STD_S(val) ! sample STD
        real(8):: std_s
        real(8), intent(in):: val(:)
        integer:: length
        real(8):: avg
    
        length = size(val)
        avg = sum(val)/length
        std_s = sqrt(dot_product((val-avg),(val-avg))/(length-1))
        if ( isnan(std_s) ) std_s = 0
    
    end function
    
    function STD_P(val) ! STD for the power distribution
        real(8):: std_p
        real(8), intent(in):: val(:,:,:,:)
        real(8), allocatable:: std_(:,:,:)
        integer:: mm, nn, oo
        integer:: s1, s2, s3
       
        s1 = size(val,2)
        s2 = size(val,3)
        s3 = size(val,4)
    
        allocate(std_(s1,s2,s3))
    
        do mm = 1, s1
        do nn = 1, s2
        do oo = 1, s3
            std_(mm,nn,oo) = std_m(val(:,mm,nn,oo))
        end do
        end do
        end do
    
        std_p = sum(std_)/(s1*s2*s3)
        if ( isnan(std_p) ) std_p = 0
    
        deallocate(std_)
    
    end function
    
    function STD_PS(val) ! sample STD for the power distribution
        real(8):: std_ps
        real(8), intent(in):: val(:,:,:,:)
        real(8), allocatable:: std_(:,:,:)
        integer:: mm, nn, oo
        integer:: s1, s2, s3
    
        s1 = size(val,2)
        s2 = size(val,3)
        s3 = size(val,4)

        allocate(std_(s1,s2,s3))
    
        do mm = 1, s1
        do nn = 1, s2
        do oo = 1, s3
            std_(mm,nn,oo) = std_s(val(:,mm,nn,oo))
        end do
        end do
        end do
    
        std_ps = sum(std_)/(s1*s2*s3)
        if ( isnan(std_ps) ) std_ps = 0
    
        deallocate(std_)
    
    end function


    function STD_MAT(val) ! sample STD for the power distribution
        real(8):: std_mat
        real(8), intent(in):: val(:,:,:)
        real(8), allocatable:: arr(:)
        integer:: mm, nn, oo
        integer:: s1, s2, s3
        integer:: num

        s1 = size(val,1)
        s2 = size(val,2)
        s3 = size(val,3)


        num = 0 
        do mm = 1, s1
        do nn = 1, s2
        do oo = 1, s3
            if ( val(mm,nn,oo) == 0 ) cycle
            num = num + 1
        end do
        end do
        end do

        allocate(arr(num))

        num = 0
        do mm = 1, s1
        do nn = 1, s2
        do oo = 1, s3
            if ( val(mm,nn,oo) == 0 ) cycle
            num = num + 1
            arr(num) = val(mm,nn,oo)
        end do
        end do
        end do

        std_mat = std_s(arr(:))

        deallocate(arr)
    
    end function

    function COV_M(val1,val2) ! covariance of the sample means
        real(8):: cov_m
        real(8), intent(in):: val1(:), val2(:)
        integer:: length
        real(8):: avg1, avg2
    
        length = size(val1)
        avg1 = sum(val1)/length
        avg2 = sum(val2)/length
        cov_m = dot_product((val1-avg1),(val2-avg2))/(length*(length-1))
        if ( isnan(cov_m) ) cov_m = 0
    
    end function
!
!    subroutine TIME_MEASURE
!        use SIMULATION_HEADER, only: t_MC, t_det, t_tot
!        implicit none
!        if ( n_batch > 1 ) then
!        allocate(t_MC(n_batch,t_totcyc))
!        allocate(t_det(n_batch,t_totcyc))
!        allocate(t_tot(n_batch,t_totcyc))
!        else
!        allocate(t_MC(n_batch,n_totcyc))
!        allocate(t_det(n_batch,n_totcyc))
!        allocate(t_tot(n_batch,n_totcyc))
!        end if
!        t_MC  = 0
!        t_det = 0
!        t_tot = 0
!    
!    end subroutine
    
    subroutine NORM_DIST(dist)
        real(8), intent(inout):: dist(1:,1:,:,:,:) ! (bat,cyc,x,y,z)
        integer:: ii, jj
    
        do ii = 1, n_batch
        do jj = 1, n_act
            dist(ii,jj,:,:,:) = dist(ii,jj,:,:,:)/AVG_P(dist(ii,jj,:,:,:))
        end do
        end do
    
    end subroutine
    
    function AVG_P(val)
        real(8):: avg_p
        real(8), intent(in):: val(:,:,:)
        integer:: sz
    
        sz = size(val)
        avg_p = sum(val)/dble(sz)
    
    end function

end module
