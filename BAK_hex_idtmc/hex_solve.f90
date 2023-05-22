module hex_solve
    
    !! NOTE !!
    ! neither of the BiCG algorithms work for larger matrices
    ! SOR algorithm appended at the end as a temporary patch
    
    contains
    
	! multiply square multidiagonal array array(:,:), diagonal offsets diags(:)
    ! with vector vec(:); return resulting vector
	function hex_mult(array, diags, vec)
	    real(8), intent(in) :: array(:,:)
		integer, intent(in) :: diags(:)
		real(8), intent(in) :: vec(:)
		real(8), allocatable :: hex_mult(:)
		integer :: vec_dim, num_diags
		integer :: i, j
		
		! allocate running variables
		vec_dim = size(vec)
		num_diags = size(diags)
		allocate(hex_mult(vec_dim))
		hex_mult(:) = 0d0
		
		! calculate vector multiplication
		do i = 1, vec_dim
		    do j = 1, num_diags
			    if (i + diags(j) < 1) then
				    cycle
				else if (i + diags(j) > vec_dim) then
				    cycle
				end if
				hex_mult(i) = hex_mult(i) + array(i, j) * vec(i + diags(j))
			end do
		end do
	end function hex_mult
    
    ! solves the linear system of equations Ax=b by successive overrelaxation
    ! matrix A array(:,:) with diagonal offsets diags(:), vector b vec(:) and x to sol(:)
    subroutine hex_sor(array, diags, vec, sol)
        real(8), intent(in) :: array(:,:)
        integer, intent(in) :: diags(:)
        real(8), intent(in) :: vec(:)
        real(8), intent(inout) :: sol(:)
        
        real(8) :: omega
        real(8) :: rel_err
        
        real(8) :: sigma
        integer :: h, w, main
        integer :: i, j, iter, idx
        
        omega = 5d-1
        
        h = size(sol)
        w = size(diags)
        
        main = -1
        do j = 1, w
            if (diags(j) .eq. 0) then
                main = j
            end if
        end do
        if (main .eq. -1) then
            print *, "SOR error: no primary diagonal!"
        end if
        
        do iter = 1, 500
            do i = 1, h
                if (array(i, main) == 0d0) then
                    cycle
                end if
                sigma = 0d0
                do j = 1, w
                    if (j .ne. main) then
                        idx = i + diags(j)
                        if ((idx .ge. 1) .and. (idx .le. h)) then
                            sigma = sigma + array(i, j) * sol(idx)
                        end if
                    end if
                end do
                sol(i) = (1d0 - omega) * sol(i) + omega * (vec(i) - sigma) / array(i, main)
            end do
            rel_err = sum(abs(hex_mult(array, diags, sol) - vec)) / sum(abs(vec))
            if (rel_err < 1e-8) then
                exit
            end if
        end do
		print *, "SOR REL ERROR ", rel_err
    end subroutine hex_sor
end module hex_solve