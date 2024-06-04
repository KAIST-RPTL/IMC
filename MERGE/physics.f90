module physics
    use omp_lib
    use variables
    use constants
    use particle_header 
    use XS_header 
    use bank_header
    use randoms
	USE tally, ONLY:INSIDE, FM_ID, tallyon, AM_ID, OUT_OF_ZZ_ADJ
	USE FMFD_header, ONLY:fmfdon, fsd_MC
	
    implicit none 
    
    contains
    
    subroutine collision_MG(p)
		implicit none 
        type(Particle), intent(inout) :: p
        real(8) :: sig_tot, rnum, wgt_s, uvw_temp(3)
        integer :: i, i_group, idx_group, n_group, n, bsize
		logical :: delayed
		integer :: pg, ng, nsplit
		real(8) :: temp, beta, lambda_b, speedn, fd
		real(8) :: rn, lambda_d, beta_d, val
		INTEGER :: iMT = 0
		INTEGER :: code_xyz
		INTEGER :: node_xyz(3)
		INTEGER :: node_x, node_y, node_z
		LOGICAL :: is_inside 
		INTEGER :: AHRI
		LOGICAL :: SONA
		INTEGER, ALLOCATABLE :: tmp_Code(:)  !> COPY & PASTE p%ptc_Code (size: p%n_prog)
		REAL(8), ALLOCATABLE :: tmp_PwTL(:)  !> COPY & PASTE p%ptc_PwTL (size: p%n_prog)
		INTEGER, ALLOCATABLE :: idx_STORE(:) !> SELECT INDEX FOR SELECTING [1:nainfo] outof tmp_Code & tmp_PwTL
		REAL(8) :: tmp_FACTOR                !> FACTOR for scaling the weight according to fainfo
		integer :: id(3)
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
		delayed = .false. 
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
        
        !> Collision estimator 
        !$omp atomic
        k_col = k_col + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)
		
		!> RELATED TO ADJOINT FLUX EVALUATION 
		IF(curr_cyc > n_inact .AND. tally_adj_flux) THEN
			! --- VECTORIZE THE POSITION
			node_xyz = AM_ID(p % coord(1) % xyz)
			node_x   = node_xyz(1)
			node_y   = node_xyz(2)
			node_z   = node_xyz(3)
			IF(zigzag_adjflux) THEN
				IF(OUT_OF_ZZ_ADJ(node_x,node_y)) THEN
					GO TO 27
				END IF
			END IF
			! --- INTEGER VALUE FOR THE POSITION (GROUP CONSIDERED)
			IF(num_adj_group == 1) THEN
				code_xyz = Code_node_XYZ(node_xyz)
			ELSE
				code_xyz = Code_node_XYZ(node_xyz,p%G)
			END IF
			IF(ANY(code_xyz .EQ. p%ptc_Code(1:p%n_prog)) .AND. p%n_prog > 0) THEN
				DO i = 1,p%n_prog
					IF(code_xyz .EQ. p%ptc_Code(i)) THEN
						p%ptc_PwTL(i) = p%ptc_PwTL(i) + p%trvlength * p%ptc_wgt0
						EXIT
					END IF
				END DO
			ELSE
				p % n_prog = p % n_prog +1
				p%ptc_Code(p%n_prog) = code_xyz
				p%ptc_PwTL(p%n_prog) = p%trvlength * p%ptc_wgt0
			END IF
		END IF
		
        !> Fission bank add
27      n = int(p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
        
		! fission site for FMFD calculation
		if ( fmfdon ) then
			if ( INSIDE(p%coord(1)%xyz) ) then
				id(:) = FM_ID(p%coord(1)%xyz)
				fsd_MC(id(1),id(2),id(3)) = fsd_MC(id(1),id(2),id(3)) + n
			end if
		end if		

        if (n > 0) then
			!> DELAYED ?
			if (allocated(MGD)) then 
				if (rang() <= sum(MGD(p%material)%beta(:))) delayed = .true. 
			endif 
			!> INCREMENT IN BANK INDEX
			bank_idx = bank_idx + 1
			!> UPDATE BANK WEIGHT & POSITION 
			thread_bank(bank_idx)%wgt = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:))/sig_tot
			thread_bank(bank_idx)%xyz = p%coord(1)%xyz
			thread_bank(bank_idx)%uvw = rand_vec()
			
			!> ADJOINT : pass particle's IFP related info. to bank
			IF(latent>1) thread_bank(bank_idx)%delayedarr(1:latent-1) = p%delayedarr(2:latent)
			IF(latent>1) thread_bank(bank_idx)%delayedlam(1:latent-1) = p%delayedlam(2:latent)
			IF(latent>1) thread_bank(bank_idx)%nlifearr  (1:latent-1) = p%nlifearr  (2:latent)
			thread_bank(bank_idx)%nlifearr(latent) = p%trvltime
			
			if (delayed) then 
				ng = size(MGD(p%material)%beta)
				thread_bank(bank_idx)%delayed 		= .true.
				thread_bank(bank_idx)%time 			= p%time
				thread_bank(bank_idx)%G 			= fission_G(p%material,.true.,iMT)
				thread_bank(bank_idx)%G_delayed     = iMT
				thread_bank(bank_idx)%E 			= real(p%material,8)
				thread_bank(bank_idx)%beta(1:ng)	= MGD(p%material)%beta(:)
				thread_bank(bank_idx)%lambda(1:ng)	= MGD(p%material)%lambda(:)
				! --- IFP ADJOINT RELATED
				thread_bank(bank_idx)%delayedarr(latent) = iMT
				thread_bank(bank_idx)%delayedlam(latent) = MGD(p%material)%lambda(iMT)
			else 
				thread_bank(bank_idx)%G 		= fission_G(p%material,.false.)
				thread_bank(bank_idx)%delayed 	= .false.
				! --- IFP ADJOINT RELATED
				thread_bank(bank_idx)%delayedarr(latent) = 0
				thread_bank(bank_idx)%delayedlam(latent) = 0
			endif 
			
			! +++ APPEND TO BANK_LONG
			IF(do_IFP_long .AND. latent > 1 .AND. curr_cyc > n_inact) THEN
				! --- APPEND THE thread_bank information to thread_bank_ENERGY
				thread_bank_LONG(bank_idx)%wgt        = thread_bank(bank_idx)%wgt
				thread_bank_LONG(bank_idx)%xyz        = thread_bank(bank_idx)%xyz
				thread_bank_LONG(bank_idx)%uvw        = thread_bank(bank_idx)%uvw
				thread_bank_LONG(bank_idx)%E          = thread_bank(bank_idx)%E
				thread_bank_LONG(bank_idx)%G          = thread_bank(bank_idx)%G
				thread_bank_LONG(bank_idx)%G_delayed  = thread_bank(bank_idx)%G_delayed
				thread_bank_LONG(bank_idx)%delayed    = thread_bank(bank_idx)%delayed
				thread_bank_LONG(bank_idx)%time       = thread_bank(bank_idx)%time
				thread_bank_LONG(bank_idx)%beta       = thread_bank(bank_idx)%beta
				thread_bank_LONG(bank_idx)%lambda     = thread_bank(bank_idx)%lambda
				thread_bank_LONG(bank_idx)%delayedarr = thread_bank(bank_idx)%delayedarr
				thread_bank_LONG(bank_idx)%delayedlam = thread_bank(bank_idx)%delayedlam
				thread_bank_LONG(bank_idx)%nlifearr   = thread_bank(bank_idx)%nlifearr
				! --- RELATED TO IFP ADJOINT TALLY
				thread_bank_LONG(bank_idx)%i_Parent = p%i_Parent
				thread_bank_LONG(bank_idx) % arr_code(1:nainfo_src) = 0
				thread_bank_LONG(bank_idx) % arr_PWTL(1:nainfo_src) = 0.d0
				IF(p%n_prog <= nainfo_src) THEN
					thread_bank_LONG(bank_idx) % arr_code(1:p%n_prog) = p%ptc_Code(1:p%n_prog)
					thread_bank_LONG(bank_idx) % arr_PWTL(1:p%n_prog) = p%ptc_PwTL(1:p%n_prog)
				ELSE
					IF(ALLOCATED(idx_STORE)) DEALLOCATE(idx_STORE)
					ALLOCATE(idx_STORE(nainfo_src))
					CALL pick_random_indexes(idx_STORE,p%n_prog,nainfo_src)
					thread_bank_LONG(bank_idx) % arr_code(1:nainfo_src) = p%ptc_Code(idx_STORE)
					thread_bank_LONG(bank_idx) % arr_PWTL(1:nainfo_src) = p%ptc_PwTL(idx_STORE) * REAL(p%n_prog,8) / REAL(nainfo_src,8)
				END IF
			END IF
        endif
        
		!> For Dynamic Neutron Source Initialization 
		if (do_DMC .and. allocated(MGD) .and. curr_cyc > n_inact) then 
			ng = size(MGD(p%material)%beta)
			beta = sum(MGD(p%material)%beta(:))
			temp = 0 
			do i = 1, ng 
				temp = temp + MGD(p%material)%beta(i) / MGD(p%material)%lambda(i)
			enddo 
			lambda_b = beta / temp
			!> Neutron Source Sample for Transient Calculation (not fission source) 
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 			= p%wgt / (MGD(p%material)%vel(p%g)*sig_tot)
			thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 			= p%coord(1)%uvw
			thread_bank_init(init_idx)%delayed 		= .false.
			thread_bank_init(init_idx)%time 		= 0
			thread_bank_init(init_idx)%G 			= p%G
			
			!> Precursor bank add
			temp = p%wgt*(beta/lambda_b)*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)
			nsplit = int(temp/1.0) + 1 
			do i = 1, nsplit
				prec_idx = prec_idx + 1
				prec_thread(prec_idx)%wgt 			= temp/real(nsplit,8)
				prec_thread(prec_idx)%xyz 			= p%coord(1)%xyz
				prec_thread(prec_idx)%G 			= p%G
				prec_thread(prec_idx)%idx 			= p%material
				prec_thread(prec_idx)%time 			= 0
				prec_thread(prec_idx)%beta(1:ng)	= MGD(p%material)%beta(:)
				prec_thread(prec_idx)%lambda(1:ng)	= MGD(p%material)%lambda(:)
			enddo
		endif
		
		!> CHANGE IN THE GROUP AFTER COLLISION & PTC WEIGHT
        rnum = rang()
        do i_group = 1, size(XS_MG(p%material)%sig_scat(p%g,:))
            if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
                idx_group = i_group
                exit
            endif 
        enddo 
        p % wgt = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:))/sig_tot
		
		!> UPDATE PTC INFO
        p % g   = idx_group
        p % last_uvw(:) = p % coord(1)% uvw(:)
        p % coord(1)% uvw(:) = rand_vec()
        
        if (p%wgt < wgt_min) THEN !call Russian_Roulette(p)
            wgt_s = 2*wgt_min
            if ((p%wgt/wgt_s).ge.rang()) then
                p%wgt = wgt_s
            else
                p%alive = .false.
            endif 
        endif 
    end subroutine collision_MG
    
    subroutine collision_MG_DT(p, macro_major)
        type(Particle), intent(inout) :: p
        real(8), intent(in) :: macro_major
        real(8) :: sig_tot, temp, rnum, wgt_s, uvw_temp(3)
        integer :: i, i_group, idx_group, n, bsize, i_source
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
        
        !> Collision estimator 
        !$omp atomic
        k_col = k_col + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)

        
        rnum = rang()
        if (rnum < (XS_MG(p%material)%sig_abs(p%g) - XS_MG(p%material)%sig_fis(p%g)) / sig_tot) then 
            p%wgt   = 0
            p%alive = .false.
            return
            
        elseif (rnum < (XS_MG(p%material)%sig_abs(p%g)) / sig_tot) then 
            !> Fission bank add
            !n = int(p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
            n = int(p%wgt*XS_MG(p%material)%nu(p%g)*(1./keff) + rang())
            do i_source = 1, n
                bank_idx = bank_idx + 1
                thread_bank(bank_idx)%xyz = p%coord(1)%xyz
                thread_bank(bank_idx)%uvw = rand_vec()
                rnum = rang()
                do i_group = 1, size(XS_MG(p%material)%chi(:))
                    if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
                        thread_bank(bank_idx)%G = i_group
                        exit
                    endif 
                enddo
            enddo
            p%wgt   = 0
            p%alive = .false.
            return
             
        
        else
            rnum = rang()
            do i_group = 1, size(XS_MG(p%material)%sig_scat(p%g,:))
                if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
                    idx_group = i_group
                    exit
                endif 
            enddo 
            !p % wgt  = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:))/sig_tot
            p % g    = idx_group
            p % coord(1) % uvw(:) = rand_vec()
        endif

    end subroutine

	function fission_G (i_mat, delayed, G_d) result(G)
		integer, intent(in) :: i_mat
		logical, intent(in) :: delayed
		integer, optional, intent(inout) :: G_d ! --- PRECURSOR GROUP INDEX
		integer :: G
		integer :: i_group, n, ng, G_prec, i
		real(8) :: rn, rnum, temp, beta
		real(8) :: val0, val1
		
		G = 1;
		rn = rang()
		n = n_group 
		G_prec = -1
		if (delayed) then 
			! sample G_prec
			temp = 0
			beta = sum(MGD(i_mat)%beta(:))
			rnum = rang() 
			ng = size(MGD(i_mat)%beta)
			do i = 1, ng
				temp = temp + MGD(i_mat)%beta(i)
				if (rnum < temp/beta) then 
					G_prec = i
					exit 
				endif
			enddo 
			! Score the precursor group index
			IF(present(G_d)) G_d = G_prec
			! Sample the outgoing group
			G = -1; val1 = 0
			val0 = sum(MGD(i_mat)%spectra(G_prec,:))
			do i_group = 1, n
				val1 = val1 + MGD(i_mat)%spectra(G_prec,i_group)
				if (rn < val1/val0) then 
					G = i_group
					exit
				endif 
			enddo
		else 
			val1 = 0
			val0 = sum(XS_MG(i_mat)%chi(:))
			do i_group = 1, n
				val1 = val1 + XS_MG(i_mat)%chi(i_group)
				if (rn < val1/val0) then 
					G = i_group
					exit
				endif 
			enddo
		endif 
		if (G < 0) then 
			print *, "fission_G :: group not selected", delayed
			stop
		endif 
		return 
	end function fission_G
	
	!===============================================================================
	! FUNCITON that converts integer array (x,y,z,g) into single integer
	!===============================================================================
	FUNCTION CODE_node_XYZ(node_xyz,node_g) RESULT (code_xyz)
		USE variables, ONLY:num_adj_group,nfm_adj
		IMPLICIT NONE
		INTEGER, INTENT(IN) 		  :: node_xyz(:)
		INTEGER, INTENT(IN), OPTIONAL :: node_g
		INTEGER             		  :: code_xyz
		INTEGER             		  :: node_size
		INTEGER :: my_g
		node_size = SIZE(node_xyz)
		IF(node_size /= 3) THEN
			WRITE(*,'(A)') 'INVALID CODING OF FINE-MESH NODE POSITION (X,Y,Z) WHILST CALLING CODE_node_XYZ.f90'
			STOP
		END IF
		code_xyz = 0
		IF(PRESENT(node_g)) THEN
			my_g = node_g
		ELSE 
			my_g = 1
		END IF
		IF(my_g > num_adj_group) THEN
			WRITE(*,'(A)') 'INVALID PROVISION OF NODE_G WHILST CALLING Code_node_XYZ.f90'
			STOP
		END IF
		code_xyz = (my_g-1)*nfm_adj(1)*nfm_adj(2)*nfm_adj(3) + (node_xyz(3)-1)*nfm_adj(1)*nfm_adj(2) + (node_xyz(2)-1)*nfm_adj(1) + node_xyz(1)
	END FUNCTION CODE_node_XYZ
	
end module 
