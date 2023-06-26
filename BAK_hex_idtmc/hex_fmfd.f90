module hex_fmfd

    use variables
	use MPI
    
	use FMFD_header
    use hex_variables
	use hex_geometry, only: hc_in_zz, hc_cmfd_coords
	use hex_cmfd

    ! use ncm(3) (input)
	! use fcr, fcz (input)
	! use dfm (input) --> edge-to-edge width/height of fine mesh cells
	! use fmD (input)
	
	! use fmDt (running variable)
	! use deltf0 (output?)
	! use Dh (running variable)
	
	! use fmJ0, fmJ1, fmJn (inputs) --> fmJ1 partial current through boundary in positive direction, fmJ0 the same in negative direction, fmJn net current
	! use Mfm(:,:,:,:) (running variable)
	
	use hex_solve
		
    real(8) :: s_s, s_l, s_e, s_v! radial surface areas (internal, long, edge external, vertex external boundaries)
    real(8) :: a_i, a_e, a_v ! axial surface areas
    real(8) :: v_i, v_e, v_v ! volumes (internal, edge, vertex cells)
	
	real(8), allocatable :: A_hf(:,:)
	real(8), allocatable :: b_hf(:), x_hf(:), b_hf_old(:)
	integer :: diags_hf(9)
	real(8), allocatable :: J0_tmp(:,:,:,:), J1_tmp(:,:,:,:), Jn_tmp(:,:,:,:)
	! integer :: x_max, y_max, z_max ! put in a header file somewhere!!

    contains
	
	subroutine hf_allocate ! equivalent to FMFD_ALLOCATION()
        use PERTURBATION, only: perton
        use COSAMPLING, only: n_pert
        use PRECONDITIONER, only: ILU_INITIAL
		
        implicit none
	    
	    integer :: i, j, k, ii, jj, kk, x, y, z ! coarse-mesh, fine-mesh, global mesh indices for iteration
		integer :: s ! surface iteration
        integer :: num
		integer :: boundaries(8)
		
	    x_max = (fcr - 1) * (ncm(2) - 1) + (2 * fcr - 1) + (ncm(1) - 1) * fcr ! calculate necessary size of the array
	    y_max = (2 * fcr - 1) + (ncm(1) - 1) * (fcr - 1) + (ncm(2) - 1) * (2 * fcr - 1)
	    z_max = fcz * ncm(3)

        nfm(1) = x_max
        nfm(2) = y_max
        nfm(3) = z_max ! TODO: Properly adjust
		
        ! parameters allocation
        allocate(k_fmfd(n_batch,n_totcyc))
        if ( icore == score ) allocate(p_fmfd(n_batch,n_act,nfm(1),nfm(2),nfm(3)))
        allocate(fsd_MC(x_max, y_max, z_max))
        allocate(fsd(x_max, y_max, z_max))
        !if ( dual_fmfd ) then
        !    allocate(k_fmfd2(n_batch,n_totcyc))
        !    allocate(fsd2(x_max, y_max, z_max))
        !    if ( icore == score ) then
        !        allocate(p_fmfd2(n_batch,n_act,x_max, y_max, z_max))
        !    end if
        !end if
        allocate(fm(x_max, y_max, z_max))
        if(perton) allocate(k_real(n_batch,n_totcyc,n_pert))
		
		! ...I have no idea what this part does; copied over from original
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
            allocate(mvec(anode,fcr,fcr,fcz,9),svec(anode,fcr,fcr,fcz))
            allocate(mvec1(anode,n_lnodes,9),svec1(anode,n_lnodes))
            allocate(mvec2(anode,9),svec2(anode))
            allocate(pcd(anode,fcr,fcr,fcz))
            allocate(ax(anode),ay(anode),az(anode))
	    
            anode = 0
            do kk = 1, ncm(3)
            do jj = 1, ncm(2)
            do ii = 1, ncm(1)
                if ( .not. hc_in_zz(ii,jj) ) cycle
                anode = anode + 1
                ax(anode) = (ii-1)*fcr
                ay(anode) = (jj-1)*fcr
                az(anode) = (kk-1)*fcz
            end do
            end do
            end do
	    
            bs0 = anode*fcr*fcr*fcz
            bs1 = bs0*9
	    
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
            allocate(mvec(anode,fcr,fcr,fcz,9),svec(anode,fcr,fcr,fcz))
            allocate(pcd(anode,fcr,fcr,fcz))
            allocate(ax(anode),ay(anode),az(anode))
	    
            anode = 0
            do kk = 1, ncm(3)
            do jj = 1, ncm(2)
            do ii = 1, ncm(1)
                if ( .not. hc_in_zz0(ii,jj,kk) ) cycle
                anode = anode + 1
                ax(anode) = (ii-1)*fcr
                ay(anode) = (jj-1)*fcr
                az(anode) = (kk-1)*fcz
            end do
            end do
            end do
            bs0 = anode*fcr*fcr*fcz
            bs1 = bs0*9
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
            allocate(acc(ii)%fm(x_max, y_max, z_max))
        enddo
	    
        allocate(Mfm(nfm(1),nfm(2),nfm(3),9)); Mfm = 0
        allocate(fphi1(x_max, y_max, z_max),fm_s(x_max, y_max, z_max))
	    
		! again, what does this do?!?
        !if ( cmfdon ) then
        !    call MPI_RANGE(anode,ncore,icore,i_para0,i_para1)
        !    call ILU_INITIAL
        !end if
	    
        allocate(fm_avg(x_max, y_max, z_max))
        allocate(fsd_FM(x_max, y_max, z_max))
	    
	    allocate(fm_t (x_max, y_max, z_max))
	    allocate(fmD  (x_max, y_max, z_max))
	    allocate(fm_a (x_max, y_max, z_max))
	    allocate(fm_nf(x_max, y_max, z_max))
	    allocate(kappa(x_max, y_max, z_max))
	    allocate(fphi0(x_max, y_max, z_max))
	    allocate(fmJ0 (x_max, y_max, z_max,8))
	    allocate(fmJ1 (x_max, y_max, z_max,8))
	    allocate(fmJn (x_max, y_max, z_max,8))
	    allocate(fmDt (x_max, y_max, z_max,8))
	    allocate(fmDh (x_max, y_max, z_max,8))
		fmDh = 0.0
	    
        allocate(ptJn (500,x_max, y_max, z_max,8))
        allocate(ptJ0 (500,x_max, y_max, z_max,8))
        allocate(ptJ1 (500,x_max, y_max, z_max,8))
	    
	    allocate(s_hf (x_max, y_max, z_max, 8))
	    allocate(v_hf (x_max, y_max, z_max))
	    
		! create arrays of cell volumes and surface areas
		s_s = dfm(1) / sqrt(3.0) ! radial surfaces, base length
		s_l = (dfm(1) / sqrt(12.0)) + dduct
		s_e = dfm(1)
		s_v = (dfm(1) / 2.0) + (dduct / sqrt(3.0))
		s_s = s_s * dfm(3) ! radial surfaces, multiply by height
		s_l = s_l * dfm(3)
		s_e = s_e * dfm(3)
		s_v = s_v * dfm(3)
		a_i = (dfm(1) ** 2.0) * (sqrt(3.0) / 2.0) ! axial surfaces
		a_e = (dfm(1) ** 2.0) * (sqrt(3.0) / 4.0) + dfm(1) * dduct
		a_v = (dfm(1) ** 2.0) * (sqrt(3.0) / 6.0) + dfm(1) * dduct + (dduct ** 2.0) / sqrt(3.0)
		v_i = a_i * dfm(3) ! volume of a cylinder is area of base times height
		v_e = a_e * dfm(3)
		v_v = a_v * dfm(3)
	    do i = 1, ncm(1)
		    do j = 1, ncm(2)
                if (.not. hc_in_zz(i,j)) cycle
			    do k = 1, ncm(3)
					do ii = 1, (2*fcr-1)
					    do jj = 1, (2*fcr-1)
						    ! check if the fine mesh cell actually exists, and if it's on an assembly edge
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
							    z = kk + (k - 1) * fcz
								! radial surface areas
								if (boundaries(1) .eq. 0 .and. boundaries(2) .ne. 0) then
								    if (boundaries(3) .eq. 1) then
									    s_hf(x, y, z, 1) = s_e
										s_hf(x, y, z, 2) = s_l
									    s_hf(x, y, z, 6) = s_s
										s_hf(x, y, z, 7) = s_l
										s_hf(x, y, z, 8) = s_s
									else
									    s_hf(x, y, z, 1) = s_v
										s_hf(x, y, z, 2) = s_l
										s_hf(x, y, z, 3) = s_v
										s_hf(x, y, z, 6) = s_s
										s_hf(x, y, z, 8) = s_l
									end if
								else if (boundaries(3) .eq. 0 .and. boundaries(1) .ne. 0) then
								    if (boundaries(7) .eq. 1) then
									    s_hf(x, y, z, 1) = s_l
										s_hf(x, y, z, 2) = s_s
										s_hf(x, y, z, 3) = s_e
										s_hf(x, y, z, 6) = s_s
										s_hf(x, y, z, 8) = s_l
									else
									    s_hf(x, y, z, 1) = s_l
										s_hf(x, y, z, 2) = s_s
										s_hf(x, y, z, 3) = s_v
										s_hf(x, y, z, 6) = s_l
										s_hf(x, y, z, 7) = s_v
									end if
								else if (boundaries(7) .eq. 0 .and. boundaries(3) .ne. 0) then
								    if (boundaries(8) .eq. 1) then
									    s_hf(x, y, z, 1) = s_s
										s_hf(x, y, z, 2) = s_s
										s_hf(x, y, z, 3) = s_l
										s_hf(x, y, z, 6) = s_l
										s_hf(x, y, z, 7) = s_e
									else
									    s_hf(x, y, z, 1) = s_s
										s_hf(x, y, z, 2) = s_l
										s_hf(x, y, z, 3) = s_l
										s_hf(x, y, z, 7) = s_v
										s_hf(x, y, z, 8) = s_v
									end if
								else if (boundaries(8) .eq. 0 .and. boundaries(7) .ne. 0) then
								    if (boundaries(6) .eq. 1) then
									    s_hf(x, y, z, 1) = s_s
										s_hf(x, y, z, 2) = s_l
										s_hf(x, y, z, 3) = s_s
										s_hf(x, y, z, 7) = s_l
										s_hf(x, y, z, 8) = s_e
									else
									    s_hf(x, y, z, 1) = s_l
										s_hf(x, y, z, 3) = s_s
										s_hf(x, y, z, 6) = s_v
										s_hf(x, y, z, 7) = s_l
										s_hf(x, y, z, 8) = s_v
									end if
								else if (boundaries(6) .eq. 0 .and. boundaries(8) .ne. 0) then
								    if (boundaries(2) .eq. 1) then
									    s_hf(x, y, z, 1) = s_l
										s_hf(x, y, z, 3) = s_s
										s_hf(x, y, z, 6) = s_e
										s_hf(x, y, z, 7) = s_s
										s_hf(x, y, z, 8) = s_l
									else
									    s_hf(x, y, z, 2) = s_v
										s_hf(x, y, z, 3) = s_l
										s_hf(x, y, z, 6) = s_v
										s_hf(x, y, z, 7) = s_s
										s_hf(x, y, z, 8) = s_l
									end if
								else if (boundaries(2) .eq. 0 .and. boundaries(6) .ne. 0) then
								    if (boundaries(1) .eq. 1) then
									    s_hf(x, y, z, 2) = s_e
										s_hf(x, y, z, 3) = s_l
										s_hf(x, y, z, 6) = s_l
										s_hf(x, y, z, 7) = s_s
										s_hf(x, y, z, 8) = s_s
									else
									    s_hf(x, y, z, 1) = s_v
										s_hf(x, y, z, 2) = s_v
										s_hf(x, y, z, 6) = s_l
										s_hf(x, y, z, 7) = s_l
										s_hf(x, y, z, 8) = s_s
									end if
								else if ((sum(boundaries(1:3)) + sum(boundaries(6:8))) .eq. 6) then
								    s_hf(x, y, z, 1:3) = s_s
									s_hf(x, y, z, 6:8) = s_s
								else
								    print *, "something has gone badly wrong"
                                    print *, [x, y, z]
                                    print *, [ii, jj, kk]
                                    print *, boundaries
									pause
								end if
								! volume and axial surface areas
								if (sum(boundaries(:)) == 6) then ! vertex case
								    v_hf(x, y, z) = v_v
									s_hf(x, y, z, 4:5) = a_v
								else if (sum(boundaries(:)) == 7) then ! edge case
								    v_hf(x, y, z) = v_e
									s_hf(x, y, z, 4:5) = a_e
								else if (sum(boundaries(:)) == 8) then ! internal case
								    v_hf(x, y, z) = v_i
									s_hf(x, y, z, 4:5) = a_i
								else ! error case
								    print *, "ERROR (HF_ALLOCATE): INVALID BOUNDARY CONDITIONS"
									stop
								end if
							end do
						end do
					end do
				end do
			end do
		end do
	end subroutine hf_allocate
	
    subroutine hf_FMFD_TO_MC(bat,cyc,fm) ! TODO: verify this!
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

	subroutine hf_process(bat,cyc) ! equivalent to PROCESS_FMFD
        use GEOMETRY_HEADER, only: universes
        use TALLY, only: tallyon
        use CMFD, only: L_DTILDA
        use PCMFD, only: L_PDHAT
        implicit none
        !> MPI derived type reduce parameters 
        integer, intent(in):: bat
        integer, intent(in):: cyc
        real(8) :: aa, bb
        integer :: ij, lc
        real(8), allocatable:: fsd_MC0(:,:,:)
        real(8), allocatable:: intra_flux(:,:,:,:), intra_kapa(:,:,:,:)
        type(FMFD_accumulation), pointer:: ac
        real(8):: tt0, tt1
        real(8):: Jn(x_max, y_max, z_max,8)
        integer:: length, which_idtmc
        integer:: reactor_type
        integer:: ncell
		
		integer :: i, j, k, ii, jj, kk, x, y, z
	    
        ! -------------------------------------------------------------------------
        ! data transmission I
        allocate(fsd_MC0(x_max, y_max, z_max))
	    
        n_cells_hf = x_max * y_max * z_max
        lc = mod(cyc-1,n_acc)+1
	    
        ac => acc(lc)
        tt0 = MPI_WTIME()
        call MPI_REDUCE(fm(:,:,:)%sig_a,ac%fm(:,:,:)%sig_a,n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(fm(:,:,:)%sig_t,ac%fm(:,:,:)%sig_t,n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(fm(:,:,:)%nusig_f,ac%fm(:,:,:)%nusig_f,n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(fm(:,:,:)%kappa,ac%fm(:,:,:)%kappa,n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        call MPI_REDUCE(fm(:,:,:)%phi,ac%fm(:,:,:)%phi,n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        do ij = 1, 8
            call MPI_REDUCE(fm(:,:,:)%J0(ij),ac%fm(:,:,:)%J0(ij),n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
            call MPI_REDUCE(fm(:,:,:)%J1(ij),ac%fm(:,:,:)%J1(ij),n_cells_hf,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
        end do
        tt1 = MPI_WTIME()

	    
        if ( .not. iscore ) then
            deallocate(fsd_MC0)
            nullify(ac)
            return
        end if
	    
        ! -------------------------------------------------------------------------
        ! data transmission II
        !fsd_MC = fsd_MC0
        !deallocate(fsd_MC0)
	    
        if ( tallyon ) call hf_FMFD_TO_MC(bat,cyc,fm)
	    
        ! current sweeping
        ! 1 (+) 0 (-)
        do i = 1, x_max
            do j = 1, y_max
                do k = 1, z_max
                    ! set incoming partial fluxes
					if (i .ne. 1)                    ac%fm(i,j,k)%J1(1) = ac%fm(i-1,j,k)%J1(8)
					if (i .ne. 1 .and. j .ne. y_max) ac%fm(i,j,k)%J1(2) = ac%fm(i-1,j+1,k)%J1(7)
					if (j .ne. 1)                    ac%fm(i,j,k)%J1(3) = ac%fm(i,j-1,k)%J1(6)
					if (k .ne. 1)                    ac%fm(i,j,k)%J1(4) = ac%fm(i,j,k-1)%J1(5)
					if (k .ne. z_max)                ac%fm(i,j,k)%J0(5) = ac%fm(i,j,k+1)%J0(4)
					if (j .ne. y_max)                ac%fm(i,j,k)%J0(6) = ac%fm(i,j+1,k)%J0(3)
					if (i .ne. x_max .and. j .ne. 1) ac%fm(i,j,k)%J0(7) = ac%fm(i+1,j-1,k)%J0(2)
					if (i .ne. x_max)                ac%fm(i,j,k)%J0(8) = ac%fm(i+1,j,k)%J0(1)
                    if (i .ne. 1)                    fmJ1(i, j, k, 1) = fmJ1(i-1, j, k, 8)
                    if (i .ne. 1 .and. j .ne. y_max) fmJ1(i, j, k, 2) = fmJ1(i-1, j+1, k, 7)
                    if (j .ne. 1)                    fmJ1(i, j, k, 3) = fmJ1(i, j-1, k, 6)
                    if (k .ne. 1)                    fmJ1(i, j, k, 4) = fmJ1(i, j, k-1, 5)
                    if (k .ne. z_max)                fmJ0(i, j, k, 5) = fmJ0(i, j, k+1, 4)
                    if (j .ne. y_max)                fmJ0(i, j, k, 6) = fmJ0(i, j+1, k, 3)
                    if (i .ne. x_max .and. j .ne. 1) fmJ0(i, j, k, 7) = fmJ0(i+1, j-1, k, 2)
                    if (i .ne. x_max)                fmJ0(i, j, k, 8) = fmJ0(i+1, j, k, 1)
                    ! calculate net fluxes ==> done in hf_fmfd_calc() !
                end do
            end do
        end do
	    
        ! initialisation
        fm_avg(:,:,:)%phi     = 0 
        fm_avg(:,:,:)%sig_t   = 0 
        fm_avg(:,:,:)%sig_a   = 0 
        fm_avg(:,:,:)%nusig_f = 0 
        fm_avg(:,:,:)%kappa   = 0 
        do ii = 1, 8 
          fm_avg(:,:,:)%J0(ii) = 0 
          fm_avg(:,:,:)%J1(ii) = 0 
        enddo 
	    
	    
        ! cycle length
        which_idtmc = 1
        select case(which_idtmc)
        case(1)
        ! iDTMC1
            if ( cyc > acc_skip ) then
            length = cyc-acc_skip
            else
            length = 1
            end if
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
	    
        where ( fm_avg(:,:,:)%phi /= 0 )
        fm_avg(:,:,:) % sig_t   = fm_avg(:,:,:) % sig_t   / fm_avg(:,:,:) % phi
        fm_avg(:,:,:) % sig_a   = fm_avg(:,:,:) % sig_a   / fm_avg(:,:,:) % phi
        fm_avg(:,:,:) % nusig_f = fm_avg(:,:,:) % nusig_f / fm_avg(:,:,:) % phi
        fm_avg(:,:,:) % kappa   = fm_avg(:,:,:) % kappa   / fm_avg(:,:,:) % phi
        end where
		do i = 1, x_max
		    do j = 1, y_max
			    do k = 1, z_max
				    if (fm_avg(i, j, k) % phi .eq. 0) cycle
					fm_avg(i, j, k) % phi = fm_avg(i, j, k) % phi / (dble(ngen) * v_hf(i, j, k) * 2d0 * length)
				end do
			end do
		end do
	    
        ! surface quantity normalization
		do i = 1, x_max
		    do j = 1, y_max
			    do k = 1, z_max
				    do ii = 1, 8
					    bb = dble(ngen) * s_hf(i, j, k, ii) * length
						if (bb == 0.0) then
						    fm_avg(i,j,k) % J0(ii) = 0.0
							fm_avg(i,j,k) % J1(ii) = 0.0
							cycle
						end if
						fm_avg(i,j,k) % J0(ii)   = fm_avg(i,j,k) % J0(ii) / bb
						fm_avg(i,j,k) % J1(ii)   = fm_avg(i,j,k) % J1(ii) / bb
					end do
				end do
			end do
		end do
		
		! SET CROSS SECTIONS AND FLUX OF NONEXISTENT CELLS TO ZERO
		do i = 1, x_max
		    do j = 1, y_max
			    do k = 1, z_max
				    if (v_hf(i,j,k) /= 0.0) then
						cycle
					else
					fm_avg(i,j,k)%sig_t = 0.0
					fm_avg(i,j,k)%sig_a = 0.0
					fm_avg(i,j,k)%nusig_f = 0.0
					fm_avg(i,j,k)%kappa = 0.0
					fm_avg(i,j,k)%phi = 0.0
					end if
				end do
			end do
		end do
	    
        ! zigzag at the corner
        if ( zigzagon ) call hf_SET_ZERO_FLUX(fm_avg(:,:,:)%phi)
    end subroutine
	
    subroutine hf_diffusion
        integer :: i, j, k, b, ii, jj, kk, bb ! cell (i, j, k) with boundary b facing (ii, jj, kk), rep. at diag. bb in matrix
        integer :: boundary ! 1 for internal, 2 for assembly boundary, 3 for reactor boundary
        real(8) :: d_temp
        
        fmDt = 0.0
        do i = 1, x_max
            do j = 1, y_max
                do k = 1, z_max
                    if (v_hf(i,j,k) == 0.0) cycle ! only consider cells that exist...
                    do b = 1, 8
                        if (s_hf(i,j,k,b) == 0.0) then
						    fmDt(i,j,k,b) = 0.0
							cycle
                        end if ! ...and boundaries that exist
                        ! what type of surface is this?
                        ! code assumes that internal boundaries and assembly-edge radial boundaries have different areas; throw an error if not
                        if (s_e == s_s .or. s_e == s_l) print *, "WARNING AREAS OF DIFFERENT BOUNDARY TYPES MATCH"                        
                        if (s_v == s_s .or. s_v == s_l) print *, "WARNING AREAS OF DIFFERENT BOUNDARY TYPES MATCH"
                        if (b == 4) then ! code assumes that coarse-mesh and fine-mesh axial boundaries are identical in area
                            boundary = 1
                            if (k == 1) boundary = 3
                        else if (b == 5) then
                            boundary = 1
                            if (k == z_max) boundary = 3
                        else if (s_hf(i,j,k,b) == s_e .or. s_hf(i,j,k,b) == s_v) then
                            boundary = 2
                        else if (s_hf(i,j,k,b) == s_s .or. s_hf(i,j,k,b) == s_l) then
                            boundary = 1
                        else
                            print *, "ERROR SURFACE AREA WEIRD ", i, j, k, b
                            stop
                        end if
                        ! find the index of the cell on the other side of the boundary
                        ii = i
                        jj = j
                        kk = k
                        if (b == 1) then
                            ii = i - 1
                        else if (b == 2) then
                            ii = i - 1
                            jj = j + 1
                        else if (b == 3) then
                            jj = j - 1
                        else if (b == 4) then
                            kk = k - 1
                        else if (b == 5) then
                            kk = k + 1
                        else if (b == 6) then
                            jj = j + 1
                        else if (b == 7) then
                            ii = i + 1
                            jj = j - 1
                        else if (b == 8) then
                            ii = i + 1
                        else
                            print *, "ERROR INVALID SURFACE INDEX", i, j, k, b
                            stop
                        end if
                        ! if radial assembly edge, check if this is also the reactor boundary
                        if (boundary == 2) then
                            ! check if this cell actually exists; if not, then it's a reactor boundary
                            if (ii < 1 .or. ii > x_max) then
                                boundary = 3
                            else if (jj < 1 .or. jj > y_max) then
                                boundary = 3
                            else if (kk < 1 .or. kk > z_max) then
                                boundary = 3
                            else if (v_hf(ii, jj, kk) == 0.0) then
                                boundary = 3
                            end if
                        end if
                        ! diffusion coefficient
                        if (boundary == 1 .or. boundary == 2) then
                            d_temp = 2.0 * fmD(i,j,k) * fmD(ii,jj,kk) / (fmD(i,j,k) + fmD(ii,jj,kk))
                        else if (boundary == 3) then
                            d_temp = 0.0 ! fmD(i,j,k)
                        else
                            print *, "ERROR BOUNDARY NOT SET", i, j, k, b, boundary
                            stop
                        end if
                        ! pin pitch adjusted diffusion coefficient
                        if (boundary == 1) then
                            if (b == 4 .or. b == 5) then
                                d_temp = d_temp / dfm(3)
                            else
                                d_temp = d_temp / dfm(1)
                            end if
                        else if (boundary == 2) then
                            if (b == 4 .or. b == 5) then
                                print *, "AXIAL BOUNDARY ODDNESS", i, j, k, b
                            else
                                d_temp = d_temp / (2.0 * dduct)
                            end if
                        else
                            if (b == 4 .or. b == 5) then
                                d_temp = d_temp / (0.5 * dfm(3) + 2.0 * d_temp)
                            else
                                d_temp = d_temp / (dduct + 2.0 * d_temp)
                            end if
                        end if
                        ! save results
                        fmDt(i,j,k,b) = d_temp
                    end do
                end do
            end do
        end do
    end subroutine hf_diffusion
	
    subroutine hf_correction
        integer :: i, j, k, b
        
        do i = 1, x_max
            do j = 1, y_max
                do k = 1, z_max
                    if (v_hf(i, j, k) == 0.0) cycle ! skip nonexistent cells
                    do b = 1, 8
                        if (s_hf(i,j,k,b) == 0) then
						    fmDh(i,j,k,b) = 0.0
						    cycle
                        end if						! skip surfaces with no diffusion
                        if (b .ge. 5) then
                            fmDh(i,j,k,b) = fmJ1(i,j,k,b) / fphi1(i,j,k) - fmDt(i,j,k,b)
                        else
                            fmDh(i,j,k,b) = fmJ0(i,j,k,b) / fphi1(i,j,k) - fmDt(i,j,k,b)
                        end if
                    end do
                end do
            end do
        end do
    end subroutine hf_correction
	
    subroutine hf_matrix
        integer :: i, j, k, b, ii, jj, kk, bb
        integer :: boundary
        
        Mfm = 0.0
        do i = 1, x_max
            do j = 1, y_max
                do k = 1, z_max
                    if (v_hf(i,j,k) == 0.0) cycle ! if the cell doesn't exist, cycle
                    Mfm(i,j,k,5) = fm_a(i,j,k) ! ABSORPTION
                    do b = 1, 8
                        if (s_hf(i,j,k,b) == 0.0) cycle ! if the boundary doesn't exist, cycle
                        ! find the index of the cell on the other side of the boundary
                        ii = i
                        jj = j
                        kk = k
                        if (b == 1) then
                            ii = i - 1
                        else if (b == 2) then
                            ii = i - 1
                            jj = j + 1
                        else if (b == 3) then
                            jj = j - 1
                        else if (b == 4) then
                            kk = k - 1
                        else if (b == 5) then
                            kk = k + 1
                        else if (b == 6) then
                            jj = j + 1
                        else if (b == 7) then
                            ii = i + 1
                            jj = j - 1
                        else if (b == 8) then
                            ii = i + 1
                        else
                            print *, "ERROR INVALID SURFACE INDEX", i, j, k, b
                            stop
                        end if
                        ! check if the cell on the other side exists
                        boundary = 1
                        if (ii < 1 .or. ii > x_max) boundary = 0
                        if (jj < 1 .or. jj > y_max) boundary = 0
                        if (kk < 1 .or. kk > z_max) boundary = 0
                        if (boundary == 1) then
                            if (v_hf(ii, jj, kk) == 0.0) boundary = 0
                        end if
                        bb = b
                        if (b .ge. 5) bb = b + 1 ! save space for the middle diagonal
                        if (boundary == 1) Mfm(i,j,k,bb) = Mfm(i,j,k,bb) - s_hf(i,j,k,b) * (fmDt(ii,jj,kk,9-b) + fmDh(ii,jj,kk,9-b)) / v_hf(i,j,k) ! INCOMING PARTIAL FLUX
                        Mfm(i,j,k,5) = Mfm(i,j,k,5) + s_hf(i,j,k,b) * (fmDt(i,j,k,b) + fmDh(i,j,k,b)) / v_hf(i,j,k) ! OUTGOING PARTIAL FLUX
                    end do
                end do
            end do
        end do
    end subroutine hf_matrix
	
    subroutine hf_init
        integer :: n_cells
        
		n_cells = x_max * y_max * z_max
		if (.not. allocated(A_hf)) then
            allocate(A_hf(n_cells, 9))
		    allocate(b_hf(n_cells))
		    allocate(b_hf_old(n_cells))
		    allocate(x_hf(n_cells))
        end if
		if (.not. allocated(Mfm)) allocate(Mfm(x_max, y_max, z_max, 9))
		diags_hf(5) =   0
		diags_hf(4) = - 1
		diags_hf(6) =   1
		diags_hf(3) = - z_max
		diags_hf(7) =   z_max
		diags_hf(2) = - (y_max - 1) * z_max
		diags_hf(8) =   (y_max - 1) * z_max
		diags_hf(1) = - y_max * z_max
		diags_hf(9) =   y_max * z_max
        Mfm = 0.0
    end subroutine hf_init
	
	! expand 4-D matrix with fine mesh coordinates into a 2-D matrix with mesh cell indices
	subroutine hf_expansion
	    integer :: x, y, z, idx
		
        A_hf = 0d0
        x_hf = 0d0
		do x = 1, x_max
		    do y = 1, y_max
			    do z = 1, z_max
				    idx = (x - 1) * y_max * z_max + (y - 1) * z_max + z
				    A_hf(idx, :) = Mfm(x, y, z, :)
					x_hf(idx) = fphi1(x, y, z)
			    end do
			end do
	    end do
	end subroutine hf_expansion
	
	! reverse of hf_expansion
	subroutine hf_interpretation
	    integer :: x, y, z, idx
		do x = 1, x_max
		    do y = 1, y_max
			    do z = 1, z_max
				    idx = (x - 1) * y_max * z_max + (y - 1) * z_max + z
					Mfm(x, y, z, :) = A_hf(idx, :)
					fphi1(x, y, z) = x_hf(idx)
				end do
			end do
		end do
    end subroutine hf_interpretation
	
    ! source vector update
    subroutine hf_source
        integer :: x, y, z, idx
		do x = 1, x_max
		    do y = 1, y_max
		    	do z = 1, z_max
		    		idx = (x - 1) * y_max * z_max + (y - 1) * z_max + z
		    		b_hf(idx) = fphi1(x, y, z) * fm_nf(x, y, z)
		    	end do
		    end do
		end do
    end subroutine hf_source
	
	! DEBUGGING FUNCTION
	subroutine hf_printvals
	    integer :: i, j, k, ii, jj, kk, x, y, z
		
	    do i = 1, ncm(1)
		    do j = 1, ncm(2)
			    if (.not. hc_in_zz(i, j)) cycle
			    do k = 1, ncm(3)
					do ii = 1, (2*fcr-1)
					    do jj = 1, (2*fcr-1)
						    ! check if the fine mesh cell actually exists
						    if (ii + jj < fcr + 1) then
							    cycle
							else if (ii + jj > 3*fcr - 1) then
							    cycle
							end if
							! compute global fine mesh coordinates
							x = (fcr - 1) * (ncm(2) - 1) + ii + (i - 1) * fcr + (j - 1) * (1 - fcr)
						    y = jj + (i - 1) * (fcr - 1) + (j - 1) * (2 * fcr - 1)
						    do kk = 1, fcz
							    z = kk + (k - 1) * fcz
								print *, "CELL ID"
								print "(9I3)", i, j, k, ii, jj, kk, x, y, z
								print *, "PARTIAL CURRENTS"
								print "(8E12.5)", fmJ0(x,y,z,1:4), fmJ1(x,y,z,5:8)
								print *, "FLUX"
								print "(E12.5)", fphi1(x,y,z)
								print *, "CROSS-SECTIONS"
								print "(5E12.5)", fm_t(x,y,z), fmD(x,y,z), fm_a(x,y,z), fm_nf(x,y,z), kappa(x,y,z)
								!print *, "CELL VOLUME"
								!print "(E12.5)", v_hf(x,y,z)
								!print *, "MATRIX COEFFICIENTS"
								!print "(9E12.4)", A_hf((x - 1) * y_max * z_max + (y - 1) * z_max + z, :)
							end do
						end do
					end do
				end do
			end do
		end do
	end subroutine hf_printvals
	
    ! MAIN FMFD SUBOUTINE
    subroutine hf_fmfd(k_eff, corr)
        logical, intent(in) :: corr
	    real(8), intent(inout) :: k_eff
	    real(8) :: k_eff_old
		real(8) :: rel_err
        
        integer :: i, j, k, ii, jj, kk, x, y, z
        
        call hf_init
		
		i = 0
		k_eff = 1d0
		rel_err = 1d0
		call hf_source ! set initial source
		b_hf_old = b_hf
	    call hf_diffusion
		if (corr) call hf_correction
		call hf_matrix
		fphi0 = fphi1 ! save initial MC flux
		do while (rel_err > 1d-8 .and. i < 500)
		    i = i + 1
			call hf_expansion
            b_hf_old = b_hf
			b_hf = b_hf / k_eff
			call hex_sor(A_hf, diags_hf, b_hf, x_hf) ! SOLVE
			call hf_interpretation
            call hf_source
            k_eff_old = k_eff
			k_eff = k_eff * sum(b_hf * b_hf_old) / sum(b_hf_old * b_hf_old) ! K UPDATE
            rel_err = sum(abs(b_hf - b_hf_old)) / sum(b_hf_old) ! CONVERGENCE TEST
        end do
		!do i = 1, x_max ! TEST REMOVE
		!    do j = 1, y_max
		!	    do k = 1, z_max
		!		    if (v_hf(i,j,k) == 0.0) cycle
		!			print "(2E12.4)", fphi0(i,j,k) / (sum(fphi0 * v_hf) / sum(v_hf)), fphi1(i,j,k) / (sum(fphi1 * v_hf) / sum(v_hf))
		!		end do
		!	end do
		!end do
	end subroutine hf_fmfd
	
    subroutine hf_WEIGHT_UPDATE(bat,cyc,k_eff,phi2) ! TODO: verify this subroutine!
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
		integer :: i, j, k, ii, jj, kk, x, y, z
    
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
        if ( .not. inactive_CMFD .and. cyc <= n_inact) update = .false.
        !        ADDITIONAL: no update for no inactive CMFD

        if ( update ) then
            if ( inactive_CMFD ) then
			
			    fsd_MC3 = 0.0 ! LINKPOINT - sum up the Monte Carlo fission source distribution!
			    do i = 1, ncm(1)
				    do j = 1, ncm(2)
					    if (.not. hc_in_zz(i, j)) cycle ! skip nonexistent fuel assemblies
					    do k = 1, ncm(3)
						    do ii = 1, (2*fcr-1)
							    do jj = 1, (2*fcr-1)
						            if (ii + jj < fcr + 1) then ! skip nonexistent fuel pins
							            cycle
							        else if (ii + jj > 3*fcr - 1) then
							            cycle
							        end if
								    do kk = 1, fcz
							            x = (fcr - 1) * (ncm(2) - 1) + ii + (i - 1) * fcr + (j - 1) * (1 - fcr)
						                y = jj + (i - 1) * (fcr - 1) + (j - 1) * (2 * fcr - 1)
							            z = kk + (k - 1) * fcz
										fsd_MC3(i,j,k) = fsd_MC3(i,j,k) + fsd_MC(x,y,z)
									end do
								end do
							end do
						end do
					end do
				end do
        
                fsd_FM3 = hc_nf*hcphi1
				fsd_MC3 = hc_nf*hcphi0
                fsd_FM3 = fsd_FM3 / sum(fsd_FM3)
                fsd_MC3 = fsd_MC3 / sum(fsd_MC3)
                fsd3 = fsd_FM3 / fsd_MC3
				
				do i = 1, ncm(1)
				    do j = 1, ncm(2)
					    if (.not. hc_in_zz(i,j)) fsd3(i,j,:) = 1.0
					end do
				end do
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
                id = hc_cmfd_coords(fission_bank(ii)%xyz, fm0, dcm, fcr) ! CM_ID(fission_bank(ii)%xyz)
                if ( id(1) < 1 .or. id(1) > ncm(1) ) cycle
                if ( id(2) < 1 .or. id(2) > ncm(2) ) cycle
                if ( id(3) < 1 .or. id(3) > ncm(3) ) cycle
				if ( .not. hc_in_zz(id(1), id(2))  ) cycle
                fission_bank(ii)%wgt = fission_bank(ii)%wgt * fsd3(id(1),id(2),id(3))
            enddo
        else
            call MPI_BCAST(fsd,n_cells_hf,MPI_REAL8,score,MPI_COMM_WORLD,ierr)
            if ( dual_fmfd ) then
                call MPI_BCAST(fsd2,n_cells_hf,MPI_REAL8,score,MPI_COMM_WORLD,ierr)
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
    
    end subroutine hf_WEIGHT_UPDATE

	subroutine hf_fmfd_calc(bat,cyc,phi1) ! equivalent to BASE_FMFD_CALCULATION
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
		integer :: i, j, k, b, ii, jj, kk

		! parameter initialization
		!tt1 = MPI_WTIME()
		if ( icore == score ) then
		k_eff = keff
		! copy parameters
		if ( inactive_cmfd .and. cyc <= n_inact ) then
			lc = mod(cyc-1,n_acc)+1
			ac => acc(lc)

			fphi1(:,:,:) = ac%fm(:,:,:)%phi
			if ( zigzagon ) call hf_SET_ZERO_FLUX(fphi1)

			where ( fphi1 /= 0 )
			fm_t(:,:,:)    = ac%fm(:,:,:)%sig_t
			fm_a(:,:,:)    = ac%fm(:,:,:)%sig_a
			fm_nf(:,:,:)   = ac%fm(:,:,:)%nusig_f 
			end where
			do ii = 1, 8
			    where ( fphi1 /= 0 ) 
			        fmJ0(:,:,:,ii) = ac%fm(:,:,:)%J0(ii)/(dble(ngen)*s_hf(:,:,:,ii))
			        fmJ1(:,:,:,ii) = ac%fm(:,:,:)%J1(ii)/(dble(ngen)*s_hf(:,:,:,ii))
			    end where
			end do
			nullify(ac)
		
			where ( fphi1 /= 0 ) 
			fm_t   = fm_t  / fphi1
			fm_a   = fm_a  / fphi1
			fm_nf  = fm_nf / fphi1
			fphi1  = fphi1 / (dble(ngen)*v_hf*2D0)
			end where

		else
			fphi1(:,:,:)    = fm_avg(:,:,:)%phi
			where ( fphi1 /= 0 )
			fm_t(:,:,:)     = fm_avg(:,:,:)%sig_t
			fm_a(:,:,:)     = fm_avg(:,:,:)%sig_a
			fm_nf(:,:,:)    = fm_avg(:,:,:)%nusig_f 
			end where
			do ii = 1, 8
			where ( fphi1 /= 0 )
			fmJ0(:,:,:,ii)  = fm_avg(:,:,:)%J0(ii)
			fmJ1(:,:,:,ii)  = fm_avg(:,:,:)%J1(ii)
			end where 
			end do

!            do ii = 1,8
!            do j = 1, y_max
!                print '(A,I2,<x_max>E15.5)', 'Jn', ii, (fmJ1(i,j,1,ii)-fmJ0(i,j,1,ii), i=1,x_max)
!            enddo
!            print *, ' '
!            enddo
            
!            do ii = 1,8
!            do j = 1, y_max
!                print '(A,I2,<x_max>E15.5)', 'J1', ii, (fmJ1(i,j,1,ii), i=1,x_max)
!            enddo
!            print *, ' '
!            enddo
!            print *, 'PHI', fphi1(5,6,1)
!            print '(A,8E15.5)', 'J0', (fmJ0(5,6,1,i), i=1,8)
!            print '(A,8E15.5)', 'J1', (fmJ1(5,6,1,i), i=1,8)
		end if
		! CALCULATE INCOMING PARTIAL FLUXES
		! (TODO): check if redundant
		do i = 1, x_max
		    do j = 1, y_max
			    do k = 1, z_max
				    if (v_hf(i,j,k) == 0.0) then
					    fmJ0(i,j,k,:) = 0.0
						fmJ1(i,j,k,:) = 0.0
					end if
					do b = 1, 8
				        ii = i
					    jj = j
						kk = k
				        if (k == 1 .or. k == 2) ii = i - 1
					    if (k == 7 .or. k == 8) ii = i + 1
						if (k == 4)             kk = k - 1
						if (k == 5)             kk = k + 1
					    if (k == 3 .or. k == 7) jj = j - 1
					    if (k == 2 .or. k == 6) jj = j + 1
						if (ii < 1 .or. ii > x_max .or. jj < 1 .or. jj > y_max .or. kk < 1 .or. kk > z_max) cycle
						if (v_hf(ii,jj,kk) == 0.0) cycle
						if (b .le. 4) fmJ1(i,j,k,b) = fmJ1(ii,jj,kk,9-b)
						if (b .ge. 5) fmJ0(i,j,k,b) = fmJ0(ii,jj,kk,9-b)
					end do
				end do
			end do
		end do
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

		if ( inactive_cmfd .and. curr_cyc <= n_inact ) then
		    if ( iscore ) then ! LINKPOINT
			    call hc_cmfd(k_eff, .true.)
			end if
		else if ( curr_cyc > n_inact) then
			if ( iscore ) then
			    call hf_diffusion ! D_TILDA_CALCULATION
			    if ( pfmfdon ) then
			    	call hf_correction ! D_PHAT_CALCULATION
			    	call hf_matrix ! PFMFD_MATRIX
			    	write(*,*) 'PFMFD', sum(fmDh)
			    endif 
			    if ( zigzagon ) call hf_SET_ZERO_FLUX(fphi1) ! SET_ZERO_M
			    call hf_fmfd(k_eff, .true.) ! POWER
			end if
		end if
		
		! weight update
		call hf_WEIGHT_UPDATE(bat,cyc,k_eff)
		!if(iscore) print *, 'keff_nopert', k_eff ! REMOVE: reimplement these print statements!
	   
		! error quantification by 1st order perturbation
		tt1 = MPI_WTIME()
		if ( perton )  call PERTURBATION_KEFF(bat,k_eff,cyc)
		

		tt2 = MPI_WTIME()
		!if ( iscore ) print*, " - perturbation total : ", tt2-tt1

		if ( icore /= score ) return
		!print*, "keff ", k_eff
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
	end subroutine hf_fmfd_calc
end module hex_fmfd
