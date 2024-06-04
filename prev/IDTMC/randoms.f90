module randoms 
    use omp_lib
    use constants, only: PI
    use, intrinsic :: ISO_C_BINDING
    
    implicit none
    !integer(C_INT) :: seed
    !! !$omp threadprivate(seed)
    !integer :: seed = 123456789
    
    interface
    function rang() result(r) bind(C)
        use ISO_C_BINDING
        implicit none
        real(C_DOUBLE) :: r
    end function rang
    end interface
    
    contains 
    
    !function rang() result (r) 
    !    real(8) :: r 
    !    call random_number ( r )
    !end function

    !function rang() result (r)    
    !    real(8) :: r
    !
    !    seed1 = mod ( seed1, 65536 )
    !    seed1 = mod ( ( 3125 * seed1 ), 65536 )
    !    r = real ( seed1, kind = 8 ) / 65536.0D+00
    !    return
    !end function 
    
    !function rang() result (r) 
    !    real(8) :: r 
    !    call random_value ( r )
    !end function
    
    
    function rand_vec() result(uvw) 
        real(8) :: uvw(3) 
        real(8) :: z, theta 
        
        z       = 2*rang() - 1
        theta = 2*pi*rang()
        
        uvw(1) = sqrt(1-z**2)*cos(theta) 
        uvw(2) = sqrt(1-z**2)*sin(theta) 
        uvw(3) = z
    end function 
    
    
    function rand_vec_2d() result(uvw) 
        real(8) :: uvw(3) 
        real(8) :: z, theta 
        
        theta = 2*pi*rang()
        
        uvw(1) = cos(theta) 
        uvw(2) = sin(theta) 
        uvw(3) = 0
        
    end function 
    
    !subroutine random_value (r)
    !
    !implicit none
    !
    !real ( kind = 8 ) r
    !!integer ( kind = 4 ) seed
    !integer :: tid 
    !
    !tid = omp_get_thread_num()
    !seed = mod ( seed, 65536 )
    !seed = mod ( ( 3125 * seed ), 65536 )
    !r = real ( seed, kind = 8 ) / 65536.0D+00
    !return
    !end subroutine


end module 
