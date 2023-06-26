module hex_cmfd
    
	use hex_variables
	use hex_solve

    contains
	
	subroutine hc_allocate
        
        allocate(hc_t  (ncm(1), ncm(2), ncm(3)))
        allocate(hcD   (ncm(1), ncm(2), ncm(3)))
        allocate(hc_a  (ncm(1), ncm(2), ncm(3)))
        allocate(hc_nf (ncm(1), ncm(2), ncm(3)))
        allocate(hcphi0(ncm(1), ncm(2), ncm(3)))
        allocate(hcphi1(ncm(1), ncm(2), ncm(3)))
        allocate(hcJ0  (ncm(1), ncm(2), ncm(3), 8))
        allocate(hcJ1  (ncm(1), ncm(2), ncm(3), 8))
        allocate(hcJn  (ncm(1), ncm(2), ncm(3), 8))
        allocate(hcDt  (ncm(1), ncm(2), ncm(3), 8))
        allocate(hcDh  (ncm(1), ncm(2), ncm(3), 8))
		
	end subroutine hc_allocate
    
    subroutine hc_homogenisation
        integer :: i, j, k, ii, jj, kk, x, y, z, b
        integer :: boundaries(8)
        
        hc_sr = dcm(1) / sqrt(3.0)
        hc_sa = (dcm(1) ** 2.0) * (sqrt(3.0) / 2.0)
        hc_v = hc_sa * dcm(3) ! volume of cylinder is base area times height
        
        hc_t = 0.0
        hcD = 0.0
        hc_a = 0.0
        hc_nf = 0.0
        hcphi0 = 0.0
        hcphi1 = 0.0
        hcJ0 = 0.0
        hcJ1 = 0.0
        hcJn = 0.0
        
	    do i = 1, ncm(1)
		    do j = 1, ncm(2)
                if (.not. hc_in_zz(i,j)) cycle ! only consider assemblies that exist
			    do k = 1, ncm(3)
					do ii = 1, (2*fcr-1)
					    do jj = 1, (2*fcr-1)
						    ! check if the designated pin actually exists, and if it's on a boundary
                            ! 0 for assembly edges, 1 for internal boundaries
							boundaries(:) = 1
						    if (ii + jj < fcr + 1) then
							    cycle
							else if (ii + jj == fcr + 1) then
							    boundaries(1) = 0
							else if (ii + jj > 3*fcr - 1) then
							    cycle
							else if (ii + jj == 3*fcr - 1) then
							    boundaries(8) = 0
							end if
							if (ii == 1) then
							    boundaries(2) = 0
							else if (ii == (2*fcr - 1)) then
							    boundaries(7) = 0
							end if
							if (jj == 1) then
							    boundaries(3) = 0
							else if (jj == (2*fcr - 1)) then
							    boundaries(6) = 0
							end if
							! compute global fine mesh coordinates
							x = (fcr - 1) * (ncm(2) - 1) + ii + (i - 1) * fcr + (j - 1) * (1 - fcr)
						    y = jj + (i - 1) * (fcr - 1) + (j - 1) * (2 * fcr - 1)
						    do kk = 1, fcz
                                boundaries(4:5) = 1
                                if (kk == 1) boundaries(4) = 0
                                if (kk == z_max) boundaries(5) = 0
							    z = kk + (k - 1) * fcz
                                hc_t(i,j,k) = hc_t(i,j,k) + v_hf(x,y,z) * fphi1(x,y,z) * fm_t(x,y,z)
                                hc_a(i,j,k) = hc_a(i,j,k) + v_hf(x,y,z) * fphi1(x,y,z) * fm_a(x,y,z)
                                hc_nf(i,j,k) = hc_nf(i,j,k) + v_hf(x,y,z) * fphi1(x,y,z) * fm_nf(x,y,z)
                                hcphi1(i,j,k) = hcphi1(i,j,k) + v_hf(x,y,z) * fphi1(x,y,z)
                                do b = 1, 8
                                    if (boundaries(b) == 1) cycle
                                    hcJ0(i,j,k,b) = hcJ0(i,j,k,b) + fmJ0(x,y,z,b) * s_hf(x,y,z,b)
                                    hcJ1(i,j,k,b) = hcJ1(i,j,k,b) + fmJ1(x,y,z,b) * s_hf(x,y,z,b)
                                end do
                            end do
                        end do
                    end do
                    ! normalisation
                    hc_t(i,j,k) = hc_t(i,j,k) / hcphi1(i,j,k)
                    hc_a(i,j,k) = hc_a(i,j,k) / hcphi1(i,j,k)
                    hc_nf(i,j,k) = hc_nf(i,j,k) / hcphi1(i,j,k)
                    hcphi1(i,j,k) = hcphi1(i,j,k) / hc_v
                    hcD(i,j,k) = 1.0 / (3.0 * hc_t(i,j,k)) ! calculate diffusion length
                    hcJ0(i,j,k,1:3) = hcJ0(i,j,k,1:3) / hc_sr
                    hcJ0(i,j,k,4) = hcJ0(i,j,k,4) / hc_sa
                    hcJ1(i,j,k,5) = hcJ1(i,j,k,5) / hc_sa
                    hcJ1(i,j,k,6:8) = hcJ1(i,j,k,6:8) / hc_sr
                end do
            end do
        end do
        hcJn = hcJ1 - hcJ0
    end subroutine hc_homogenisation
    
    subroutine hc_diffusion
        integer :: i, j, k, ii, jj, kk, b
        logical :: boundaries(8)
        
        hcDt = 0.0
        hcDh = 0.0
        do i = 1, ncm(1)
            do j = 1, ncm(2)
                if (.not. hc_in_zz(i,j)) cycle ! consider only nodes that exist
                do k = 1, ncm(3)
                    ! see if neighbouring cells exist
                    boundaries(:) = .false.
                    if (.not. hc_in_zz(i-1,j))   boundaries(1) = .true.
                    if (.not. hc_in_zz(i-1,j+1)) boundaries(2) = .true.
                    if (.not. hc_in_zz(i,j-1))   boundaries(3) = .true.
                    if (k == 1)                  boundaries(4) = .true.
                    if (k == ncm(3))             boundaries(5) = .true.
                    if (.not. hc_in_zz(i,j+1))   boundaries(6) = .true.
                    if (.not. hc_in_zz(i+1,j-1)) boundaries(7) = .true.
                    if (.not. hc_in_zz(i+1,j))   boundaries(8) = .true.
                    do b = 1, 8
                        if (.not. boundaries(b)) then ! process internal boundaries
                            ! get index of neighbour node
                            ii = i
                            jj = j
                            kk = k
                            if (b == 1) ii = i - 1
                            if (b == 2) then
                                ii = i - 1
                                jj = j + 1
                            end if
                            if (b == 3) jj = j - 1
                            if (b == 4) kk = k - 1
                            if (b == 5) kk = k + 1
                            if (b == 6) jj = j + 1
                            if (b == 7) then
                                ii = i + 1
                                jj = j - 1
                            end if
                            if (b == 8) ii = i + 1
                            hcDt(i,j,k,b) = 2.0 * hcD(i,j,k) * hcD(ii,jj,kk) / (hcD(i,j,k) + hcD(ii,jj,kk)) ! diffusion coefficient at boundary
                            if (b == 4 .or. b == 5) then
                                hcDt(i,j,k,b) = hcDt(i,j,k,b) / ncm(3)
                            else
                                hcDt(i,j,k,b) = hcDt(i,j,k,b) / ncm(1)
                            end if
                        else ! process reactor boundaries
                            if (b == 4 .or. b == 5) then
                                hcDt(i,j,k,b) = hcD(i,j,k) / (ncm(3) + 2.0 * hcD(i,j,k))
                            else
                                hcDt(i,j,k,b) = hcD(i,j,k) / (ncm(1) + 2.0 * hcD(i,j,k))
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end subroutine hc_diffusion
    
    subroutine hc_correction
        integer :: i, j, k, b
        
        hcDh = 0.0
        do i = 1, ncm(1)
            do j = 1, ncm(2)
                if (.not. hc_in_zz(i,j)) cycle ! only consider cells that exist
                do k = 1, ncm(3)
                    do b = 1, 8
                        if (hcDt(i,j,k,b) == 0.0) cycle
                        if (b .ge. 5) then
                            hcDh(i,j,k,b) = hcJ1(i,j,k,b) / hcphi1(i,j,k) - hcDt(i,j,k,b)
                        else
                            hcDh(i,j,k,b) = hcJ0(i,j,k,b) / hcphi1(i,j,k) - hcDt(i,j,k,b)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine hc_correction
    
    subroutine hc_matrix
        integer :: i, j, k, ii, jj, kk, b, bb
        logical :: boundaries(8)
        
        Mhc = 0.0
        do i = 1, ncm(1)
            do j = 1, ncm(2)
                if (.not. hc_in_zz(i,j)) cycle
                do k = 1, ncm(3)
                    Mhc(i,j,k,5) = hc_a(i,j,k) ! ABSORPTION
                    ! see if neighbouring cells exist
                    boundaries(:) = .false.
                    if (.not. hc_in_zz(i-1,j))   boundaries(1) = .true.
                    if (.not. hc_in_zz(i-1,j+1)) boundaries(2) = .true.
                    if (.not. hc_in_zz(i,j-1))   boundaries(3) = .true.
                    if (k == 1)                  boundaries(4) = .true.
                    if (k == ncm(3))             boundaries(5) = .true.
                    if (.not. hc_in_zz(i,j+1))   boundaries(6) = .true.
                    if (.not. hc_in_zz(i+1,j-1)) boundaries(7) = .true.
                    if (.not. hc_in_zz(i+1,j))   boundaries(8) = .true.
                    do b = 1, 8
                        ! get index of neighbour node
                        ii = i
                        jj = j
                        kk = k
                        if (b == 1) ii = i - 1
                        if (b == 2) then
                            ii = i - 1
                            jj = j + 1
                        end if
                        if (b == 3) jj = j - 1
                        if (b == 4) kk = k - 1
                        if (b == 5) kk = k + 1
                        if (b == 6) jj = j + 1
                        if (b == 7) then
                            ii = i + 1
                            jj = j - 1
                        end if
                        if (b == 8) ii = i + 1
                        ! get diagonal index
                        bb = b
                        if (b .ge. 5) bb = b + 1
                        ! build matrix
                        if (b == 4 .or. b == 5) then
                            if (.not. boundaries(b)) Mhc(i,j,k,bb) = Mhc(i,j,k,bb) - hc_sa * (hcDt(ii,jj,kk,9-b) + hcDh(ii,jj,kk,9-b)) / hc_v ! INCOMING PARTIAL CURRENT
                            Mhc(i,j,k,5) = Mhc(i,j,k,5) + hc_sa * (hcDt(i,j,k,b) + hcDh(i,j,k,b)) / hc_v ! OUTGOING PARTIAL CURRENT
                        else
                            if (.not. boundaries(b)) Mhc(i,j,k,bb) = Mhc(i,j,k,bb) - hc_sr * (hcDt(ii,jj,kk,9-b) + hcDh(ii,jj,kk,9-b)) / hc_v ! INCOMING PARTIAL CURRENT
                            Mhc(i,j,k,5) = Mhc(i,j,k,5) + hc_sr * (hcDt(i,j,k,b) + hcDh(i,j,k,b)) / hc_v ! OUTGOING PARTIAL CURRENT
                        end if
                    end do
                end do
            end do
        end do
    end subroutine hc_matrix
    
    subroutine hc_init
        integer :: n_cells
        
		n_cells = ncm(1) * ncm(2) * ncm(3)
		if (.not. allocated(A_hc)) then
            allocate(A_hc(n_cells, 9))
		    allocate(b_hc(n_cells))
		    allocate(b_hc_old(n_cells))
		    allocate(x_hc(n_cells))
        end if
        if (.not. allocated(Mhc)) allocate(Mhc(ncm(1), ncm(2), ncm(3), 9))
		diags_hc(5) =   0
		diags_hc(4) = - 1
		diags_hc(6) =   1
		diags_hc(3) = - ncm(3)
		diags_hc(7) =   ncm(3)
		diags_hc(2) = - (ncm(2) - 1) * ncm(3)
		diags_hc(8) =   (ncm(2) - 1) * ncm(3)
		diags_hc(1) = - ncm(2) * ncm(3)
		diags_hc(9) =   ncm(2) * ncm(3)
        Mhc = 0.0
    end subroutine hc_init
    
    subroutine hc_expansion
	    integer :: x, y, z, idx
		
        A_hc = 0d0
        x_hc = 0d0
		do x = 1, ncm(1)
		    do y = 1, ncm(2)
			    do z = 1, ncm(3)
				    idx = (x - 1) * ncm(2) * ncm(3) + (y - 1) * ncm(3) + z
				    A_hc(idx, :) = Mhc(x, y, z, :)
					x_hc(idx) = hcphi1(x, y, z)
			    end do
			end do
        end do
    end subroutine hc_expansion
    
    subroutine hc_interpretation
	    integer :: x, y, z, idx
        
		do x = 1, ncm(1)
		    do y = 1, ncm(2)
			    do z = 1, ncm(3)
				    idx = (x - 1) * ncm(2) * ncm(3) + (y - 1) * ncm(3) + z
					Mhc(x, y, z, :) = A_hc(idx, :)
					hcphi1(x, y, z) = x_hc(idx)
				end do
			end do
		end do
    end subroutine hc_interpretation
    
    subroutine hc_source
        integer :: x, y, z, idx
        
		do x = 1, ncm(1)
		    do y = 1, ncm(2)
			    do z = 1, ncm(3)
				    idx = (x - 1) * ncm(2) * ncm(3) + (y - 1) * ncm(3) + z
		    		b_hc(idx) = hcphi1(x, y, z) * hc_nf(x, y, z)
		    	end do
		    end do
		end do
    end subroutine hc_source
    
    subroutine hc_cmfd(k_eff, corr)
        logical, intent(in) :: corr
	    real(8), intent(inout) :: k_eff
	    real(8) :: k_eff_old
		real(8) :: rel_err
        
        integer :: i, j, k, ii, jj, kk, x, y, z
        
		if ( .not. allocated(hc_t) ) call hc_allocate
        call hc_init
		
		i = 0
		k_eff = 1d0
		rel_err = 1d0
        call hc_homogenisation
		hcphi0 = hcphi1
		call hc_source ! set initial source
		b_hc_old = b_hc
	    call hc_diffusion
		if (corr) call hc_correction
		call hc_matrix
		
		! TEST REMOVE
		do i = 1, ncm(1)
		    do j = 1, ncm(2)
			    if (.not. hc_in_zz(i,j)) cycle
			    do k = 1, ncm(3)
				    print *, i, j, k
					print "(A, 8E12.4)", "CURR", hcJ0(i,j,k,1:4), hcJ1(i,j,k,5:8)
					print "(A, 4E12.4)", "CMXS", hcphi1(i,j,k), hc_t(i,j,k), hc_a(i,j,k), hc_nf(i,j,k)
				end do
			end do
		end do
		
		do while (rel_err > 1d-9 .and. i < 500)
		    !if (corr) call hc_correction
			!if (corr) call hc_matrix
		    i = i + 1
			call hc_expansion
            b_hc_old = b_hc
			b_hc = b_hc / k_eff
			call hex_sor(A_hc, diags_hc, b_hc, x_hc) ! SOLVE
			call hc_interpretation
            call hc_source
            k_eff_old = k_eff
			k_eff = k_eff * sum(b_hc * b_hc_old) / sum(b_hc_old * b_hc_old) ! K UPDATE
            rel_err = abs(k_eff - k_eff_old) ! CONVERGENCE TEST
        end do
    end subroutine hc_cmfd
	
end module hex_cmfd