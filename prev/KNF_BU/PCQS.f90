module PCQS

	use transient 
	use bank_header
	use particle_header, only : particle 
	use ace_header, 	only : ace, CrossSectionDataForm, EnergyDist
	use variables
	use tally, only : k_eff 
	
	implicit none 

	real(8) :: PCQS_keff
	integer :: n_pcqs_act
	integer :: n_pcqs_inact
	integer :: n_pcqs_totcyc
	
	
	!> PKE parameters 
	real(8), allocatable :: PKE_beta_d(:), PKE_beta_tally1(:), PKE_beta_tally2(:)
	real(8), allocatable :: PKE_prec(:), PKE_prec_tally1(:), PKE_prec_tally2(:), PKE_prec0(:)
	real(8), allocatable :: PKE_lambda(:), PKE_lambda_tally1(:), PKE_lambda_tally2(:)
	real(8) :: PKE_gen, PKE_gen_tally1, PKE_gen_tally2
	
	!> PKE parameter tally 
	real(8), allocatable :: PKE_tally(:) 
	integer :: N_PKE_tally 
	
	
	
	
	
	real(8) :: PKE_beta 
	real(8) :: PKE_alpha, PKE_alpha0 
	real(8) :: PKE_rho, PKE_rho0 
	real(8) :: PKE_keff, PKE_keff_tally, PKE_keff0 
	real(8) :: PKE_n, PKE_n0 
	real(8) :: PKE_Z, PKE_Z_tally1, PKE_Z_tally2
	real(8) :: PKE_f
	real(8) :: PKE_amp 
	real(8) :: PKE_gamma
	real(8) :: PCQS_abs, PCQS_prod, PCQS_leak 
	
	logical :: corrector = .false. 
	
	real(8), allocatable :: PCQS_keff_cyc(:), PCQS_power_cyc(:)
	real(8), allocatable :: PCQS_power(:,:) 
	
	real(8), allocatable :: PCQS_beta_cyc(:,:)
	real(8), allocatable :: PCQS_lambda_cyc(:,:)
	real(8), allocatable :: PCQS_gen_cyc(:)
	
	contains 
	
	subroutine PKE_init() 
	
	
		N_PKE_tally = npg & ! beta_d*F*flux 
					+ npg & ! beta_d*F*flux/lambda 
					+ 1   & ! F*flux 
					+ 1     ! flux / vel 
	
		allocate(PKE_tally(N_PKE_tally))
		PKE_tally(:) = 0 
		
		allocate(PKE_beta_d(npg),PKE_beta_tally1(npg),PKE_beta_tally2(npg)) 
		allocate(PKE_prec(npg),PKE_prec_tally1(npg),PKE_prec_tally2(npg), PKE_prec0(npg)) 
		allocate(PKE_lambda(npg),PKE_lambda_tally1(npg),PKE_lambda_tally2(npg)) 
		
		n_pcqs_totcyc = N_PCQS_INACT + N_PCQS_ACT
		
		!> Initialize tally bank
		PKE_beta_tally1 = 0;  PKE_beta_tally2 = 0 
		PKE_lambda_tally1 = 0; PKE_lambda_tally2 = 0;
		PKE_prec_tally1 = 0;  PKE_prec_tally2 = 0;
		PKE_gen_tally1 =0;  PKE_gen_tally2 = 0
		PKE_Z_tally1 =0; PKE_Z_tally2 = 0 
		
		
		PKE_beta_d(:) = 0
		PKE_lambda(:) = 0
		PKE_gen = 0
		PKE_Z_tally1 = 0
		
		
		PKE_keff0 = 0
		PKE_keff_tally = 0 
		PKE_f = 1
		PKE_n = 1
		PKE_n0 = 1
		PKE_amp = 1 
		allocate(PCQS_keff_cyc(n_pcqs_act)) 
		allocate(PCQS_power_cyc(n_pcqs_act)) 
		allocate(PCQS_power(0:n_timestep,2)) 
		
		
		allocate(PCQS_beta_cyc(n_pcqs_act,npg))
		allocate(PCQS_lambda_cyc(n_pcqs_act,npg))
		allocate(PCQS_gen_cyc(n_pcqs_act))
				 
				 
		PCQS_power(:,:) = 0 
		PKE_gamma = 0 
		
		n_col_avg = int(n_col_avg / real(n_act,8)) 
		n_cross_avg = int(n_cross_avg / real(n_act,8)) 
		
		
		
		!SSP = ngen*10.0/ real(n_col_avg,8) 
		SSP = 1.0
		
		call MPI_BCAST(SSP, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		
	end subroutine 
	
	
	
	
	subroutine solve_PKE() 
		integer :: i, j, i_micro
		real(8), allocatable :: A(:,:), B(:), x(:) 
		real(8) :: rcv_buf
		real(8) :: alpha, Z, amp
		real(8) :: n0, n1 
		
		
		if (icore==score) then
		
		!> Calculate PKE parameters
			PKE_beta_d(:) = PKE_beta_d(:) / real(n_pcqs_act,8)
			PKE_lambda(:) = PKE_lambda(:) / real(n_pcqs_act,8)
			PKE_gen = PKE_gen / real(n_pcqs_act,8)
			PKE_keff = PKE_keff_tally / real(n_pcqs_act,8)
			PKE_beta = sum(PKE_beta_d) 
			PKE_prec0 = PKE_prec
			
			PKE_rho = 1.0d0 - PKE_keff0 / PKE_keff
			
			PKE_alpha0 = PKE_alpha
			PKE_alpha = (PKE_rho - PKE_beta) / PKE_gen
			
			
			
			
			if (curr_timestep == 0) then 
				PKE_Z_tally2 = PKE_Z_tally1
				PKE_keff0 = 1.0d0 
				PKE_keff = PKE_keff0
				
				PKE_n = 1.0d0 
				PKE_prec0(:) = PKE_beta_d(:) * PKE_n / (PKE_gen * PKE_lambda) 
				PKE_prec = PKE_prec0
				
				PKE_Z = 1.0d0
				PKE_f = 1.0d0 
			else 
			
				PKE_Z = PKE_Z_tally1 / PKE_Z_tally2
				PKE_n0 = PKE_n 
				
				!> Solve PKE
				allocate(A(1+npg, 1+npg), B(1+npg), x(1+npg))
				
				!print '(I,2F14.5)', 0, PKE_n, PKE_Z
				do i_micro = 1, n_microtimestep
					A(:,:) = 0 
					B(:) = 0 
					
					!> linear interpolation of alpha
					alpha = PKE_alpha0 + (PKE_alpha-PKE_alpha0)*i_micro/real(n_microtimestep,8) 
					
					do i = 1, npg
						A(1,i+1)	= -dt_micro * PKE_lambda(i)
						A(i+1,1)	= -dt_micro * PKE_beta_d(i) / PKE_gen 
						A(i+1,i+1)  = 1.0d0 + dt_micro*PKE_lambda(i)
						B(i+1) 		= PKE_prec(i)
					enddo 
					A(1,1) = 1.0d0 - dt_micro*alpha
					B(1)   = PKE_n
					
					!> Solve Ax=B 
					call GEM(A,B,x,1+npg)
					
					
					!> Update PKE_n, PKE_prec 
					PKE_n = x(1) 
					PKE_prec(1:npg) = x(2:npg+1) 
					amp = PKE_n
					
					!write(prt_dynamic, '(I,5E14.5)') i_micro, PKE_n / PKE_Z, amp, PKE_Z, PKE_rho, PKE_rho/PKE_beta  
				enddo 
				deallocate(A,B,x) 
			
			
				
			
				PKE_f = PKE_n / PKE_Z
				PKE_amp = amp 
				
				
			endif 
			
			
			
		endif 
		
		
		
		
		!> Broadcast normalization factors
		call MPI_BCAST(PKE_f, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		call MPI_BCAST(PKE_n, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		call MPI_BCAST(PKE_n0, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		
		
		
		!> Adjust weight
		if (curr_timestep > 0) then 
			prompt_bank(:)%wgt = prompt_bank(:)%wgt * PKE_f 
			delayed_bank(:)%wgt = delayed_bank(:)%wgt * PKE_f 
		endif 
		
		
		PKE_beta_d(:) = 0
		PKE_lambda(:) = 0
		PKE_gen = 0
		PKE_Z_tally1 = 0
		
		
		!> Zero tally bank 
		PKE_beta_tally1 = 0;  PKE_beta_tally2 = 0;
		PKE_lambda_tally1 = 0; PKE_lambda_tally2 = 0;
		PKE_prec_tally1 = 0;  PKE_prec_tally2 = 0;
		PKE_gen_tally1 = 0;  PKE_gen_tally2 = 0;
		PKE_keff_tally = 0 
		PKE_Z_tally1 = 0
		
	end subroutine 
	
	subroutine collision_pcqs_MG(p) 
	
        type(Particle), intent(inout) :: p
        real(8) :: sig_tot, rnum, wgt_s, uvw_temp(3), val 
        integer :: i, i_group, idx_group, n_group, n, imat
		integer :: pg, ng, nsplit
		real(8) :: temp, beta, beta_d, lambda_d
		real(8) :: speedn, wgt_prev , fd, rn
		logical :: delayed
		real(8) :: wgt, sigtot_pcqs 
		
		
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
		
		imat = p%material
        wgt_prev = p%wgt
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
		speedn = MGD(p%material)%vel(p%g)
		temp = 1.0d0 / (speedn*del_t) 
		
		sigtot_pcqs = sig_tot + temp + (PKE_gamma / speedn) 
		wgt = p%wgt * sig_tot / sigtot_pcqs
		
        !> keff estimator 
        !$omp atomic
        PCQS_keff = PCQS_keff + p%wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/(sig_tot))
        
		
		!beta = sum(MGD(p%material)%beta(:))
		!if (curr_cyc > n_pcqs_inact .and. beta > 0) then 
		!	!> PKE parameter tally 
		!	!$omp critical 
		!	do i = 1, npg 
		!		PKE_tally(i) = PKE_tally(i)+ wgt*MGD(p%material)%beta(i)*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g) / sig_tot
		!		PKE_tally(i+npg) = PKE_tally(i+npg) + wgt*MGD(p%material)%beta(i)*XS_MG(p%material)%nu(p%g) &
		!											*XS_MG(p%material)%sig_fis(p%g) / MGD(p%material)%lambda(i) / sig_tot
		!	enddo 
		!	PKE_tally(npg*2+1)     = PKE_tally(npg*2+1) + wgt*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g) / sig_tot
		!	PKE_tally(npg*2+2)     = PKE_tally(npg*2+2) + wgt / (speedn * sig_tot)
		!	!$omp end critical 
		!endif
		
		
		! ===================================================================================================
		!  오래 걸리는 부분 
		!  1. 조건문 내부 작업 최대한 줄이기 
		!  2. 조건문의 조건 단순화
		! ===================================================================================================

		if (rang() < SSP) then 
		
		!> Source tally for the next ITERATION 
		val = (XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/PKE_keff0)
		if (XS_MG(p%material)%sig_fis(p%g) > 0) then 
			ng = size(MGD(p%material)%beta)
			beta = sum(MGD(p%material)%beta(:))
		
			if (rang() < beta) then 
				delayed = .true. 
			else 
				delayed = .false. 
			endif 
			bank_idx = bank_idx + 1 
			thread_bank(bank_idx)%xyz 	= p%coord(1)%xyz
			thread_bank(bank_idx)%uvw 	= rand_vec()
			
			if (delayed) then 
				temp = 0
				rn = rang()
				do i = 1, ng 
					temp = temp + MGD(p%material)%beta(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				lambda_d = MGD(p%material)%lambda(i_group)
				fd = exp(-lambda_d*del_t) * (1-exp(lambda_d*del_t)+lambda_d*del_t*exp(lambda_d*del_t))/(lambda_d*del_t)  ! f_3,d
				thread_bank(bank_idx)%G 	= fission_G(p%material,.true., i_group)
				thread_bank(bank_idx)%wgt 	= wgt * fd * val * (1.0/SSP)
				thread_bank(bank_idx)%delayed 	= delayed
				
			else 
				thread_bank(bank_idx)%G 	= fission_G(p%material,.false.)
				thread_bank(bank_idx)%wgt 	= wgt * val * (1.0/SSP)
				thread_bank(bank_idx)%delayed 	= delayed
			endif 
			
			
		endif 
        
        endif 
        
		
		
		!> Source tally for the next TIME-STEP
		if (curr_cyc == n_pcqs_totcyc ) then 
		
		if (rang() < SSP) then 
			!> For PCQS Neutron Source Initialization 
			! Source 1 
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 			= (1.0/SSP) * p%wgt * exp(PKE_gamma * del_t) * (sigtot_pcqs-sig_tot-(PKE_gamma/speedn)) / sigtot_pcqs 
			thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 			= p%coord(1)%uvw
			thread_bank_init(init_idx)%G 			= p%G  ! incident E group 
			
			
			ng = size(MGD(p%material)%beta)
			beta = sum(MGD(p%material)%beta(:))
			if (beta > 0) then 
				! Source 2
				! select delayed group 
				temp = 0 
				rn = rang()
				do i = 1, ng 
					temp = temp + MGD(p%material)%beta(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				beta_d = MGD(p%material)%beta(i_group)
				lambda_d = MGD(p%material)%lambda(i_group)
				
				val = (beta*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/PKE_keff0)
				
				
				fd = (1-exp(-lambda_d*del_t)*(1+lambda_d*del_t))/(lambda_d*del_t)  ! f_2,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 	= wgt * fd *val * (1.0/SSP)
				thread_bank_init(init_idx)%xyz 	= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 	= rand_vec()
				thread_bank_init(init_idx)%G 	= fission_G(p%material,.true.,i_group)
			
				! Source 3 
				! select delayed group 
				rn = rang()
				temp = 0 
				do i = 1, ng 
					temp = temp + MGD(p%material)%beta(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				
				beta_d = MGD(p%material)%beta(i_group)
				lambda_d = MGD(p%material)%lambda(i_group)
				val = (beta*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/(sig_tot*lambda_d))
				fd = exp(-lambda_d*del_t)   ! f_1,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 			= wgt * fd * lambda_d * val* (1.0/SSP)
				thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 			= rand_vec() 
				thread_bank_init(init_idx)%G 			= fission_G(p%material,.true.,i_group)
				
			endif 
		endif 
		
		
		endif 
		! ===================================================================================================
		! ===================================================================================================
		
		
		
        rnum = rang()
		n = size(XS_MG(p%material)%sig_scat(p%g,:))
        do i_group = 1, n
            if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
                idx_group = i_group
                exit
            endif 
        enddo 
		
		
		
		
		! PCQS weight reduction 
        p % wgt = p % wgt * (sum(XS_MG(p%material)%sig_scat(p%g,:))) / (sig_tot+1.0d0 / (speedn*del_t) )
        p % g   = idx_group
        p % last_uvw(:) = p % coord(1)% uvw(:)
        p % coord(1)% uvw(:) = rand_vec()
		
        if (p%wgt < wgt_min_dyn) THEN !call Russian_Roulette(p)
            wgt_s = 2*wgt_min_dyn
            if ((p%wgt/wgt_s).ge.rang()) then
                p%wgt = wgt_s
            else
                p%alive = .false.
            endif 
        endif 	
	


	end subroutine 

	
	
	
	
	subroutine collision_pcqs_MG_init(p) 
	
        type(Particle), intent(inout) :: p
        real(8) :: sig_tot, rnum, wgt_s, uvw_temp(3)
        integer :: i, i_group, idx_group, n_group, n, imat
		integer :: pg, ng, nsplit
		real(8) :: rn, temp, beta, lambda_d, beta_d
		real(8) :: speedn, wgt_prev , fd, val 
		logical :: delayed
		real(8) :: wgt
		
		
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
		imat = p%material
		speedn = MGD(p%material)%vel(p%g)
        sig_tot = sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g)
		wgt_min_dyn = wgt_min 
		
		wgt = p%wgt !* sig_tot / (sig_tot + 1.0/(speedn*del_t))
		
        !> keff estimator 
        !$omp atomic
        PCQS_keff = PCQS_keff + wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g))/(sig_tot)
		
		
        !> Fission bank add
        n = int(wgt*(XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1/keff) + rang())
        if (n > 0) then
			if (allocated(MGD)) then 
				if (rang() <= sum(MGD(p%material)%beta(:))) delayed = .true. 
			endif 
			
			bank_idx = bank_idx + 1
			thread_bank(bank_idx)%xyz = p%coord(1)%xyz
			thread_bank(bank_idx)%uvw = rand_vec()
			
			if (delayed) then 
				thread_bank(bank_idx)%delayed 		= .true.
				thread_bank(bank_idx)%G 			= fission_G(p%material,.true.)
			else 
				thread_bank(bank_idx)%G 		= fission_G(p%material,.false.)
				thread_bank(bank_idx)%delayed 	= .false.
			endif 
			
        endif
		
		
		if (rang() < SSP) then 
		
		wgt = wgt / SSP
		
		
		!> For PCQS Neutron Source Initialization 
		! Source 1 
		init_idx = init_idx + 1
		thread_bank_init(init_idx)%wgt 			= wgt/(speedn*del_t*sig_tot)
		thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
		thread_bank_init(init_idx)%uvw 			= p%coord(1)%uvw
		thread_bank_init(init_idx)%G 			= p%G  ! incident E group 
		
		
		ng = size(MGD(p%material)%beta)
		beta = sum(MGD(p%material)%beta(:))
		if (beta > 0) then 
			val = (beta*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)!*(1/keff)
			! Source 2
			! select delayed group 
			temp = 0 
			rn = rang()
			do i = 1, ng 
				temp = temp + MGD(p%material)%beta(i)
				if (rn < temp/beta) then 
					i_group = i  
					exit 
				endif 
			enddo 
			lambda_d = MGD(p%material)%lambda(i_group)
			
			fd = (1-exp(-lambda_d*del_t)*(1+lambda_d*del_t))/(lambda_d*del_t)  ! f_2,d
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 	= wgt * fd * val
			thread_bank_init(init_idx)%xyz 	= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 	= rand_vec()
			thread_bank_init(init_idx)%G 	= fission_G(p%material,.true.,i_group)
		
		
			! Source 3 
			! select delayed group 
			rn = rang()
			temp = 0 
			do i = 1, ng 
				temp = temp + MGD(p%material)%beta(i)
				if (rn < temp/beta) then 
					i_group = i  
					exit 
				endif 
			enddo 
			
			beta_d = MGD(p%material)%beta(i_group)
			lambda_d = MGD(p%material)%lambda(i_group)
			val = beta*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/(sig_tot*lambda_d)
			fd = exp(-lambda_d*del_t) ! f_1,d 
        
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 			= wgt * fd * lambda_d * val
			thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 			= rand_vec() 
			thread_bank_init(init_idx)%G 			= fission_G(p%material,.true.,i_group)
			
			
			
		endif 
		
		
		
		
			
		val = (XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)/sig_tot)*(1.0/keff)
		if (XS_MG(p%material)%sig_fis(p%g) > 0) then 
			ng = size(MGD(p%material)%beta)
			beta = sum(MGD(p%material)%beta(:))
		
			if (rang() < beta) then 
				delayed = .true. 
			else 
				delayed = .false. 
			endif 
			split_idx = split_idx + 1 
			split_thread(split_idx)%xyz 	= p%coord(1)%xyz
			split_thread(split_idx)%uvw 	= rand_vec()
			
			if (delayed) then 
				temp = 0
				rn = rang()
				do i = 1, ng 
					temp = temp + MGD(p%material)%beta(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				lambda_d = MGD(p%material)%lambda(i_group)
				fd = exp(-lambda_d*del_t) * (1-exp(lambda_d*del_t)+lambda_d*del_t*exp(lambda_d*del_t))/(lambda_d*del_t)  ! f_3,d
				split_thread(split_idx)%G 	= fission_G(p%material,delayed, i_group)
				split_thread(split_idx)%wgt = wgt * fd * val
				split_thread(split_idx)%delayed = delayed
			else 
				split_thread(split_idx)%G 	= fission_G(p%material,delayed)
				split_thread(split_idx)%wgt = wgt * val
				split_thread(split_idx)%delayed = delayed
			endif 
			
		endif 
		
		endif 
		
        rnum = rang()
		n = size(XS_MG(p%material)%sig_scat(p%g,:))
        do i_group = 1, n
            if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
                idx_group = i_group
                exit
            endif 
        enddo 
		
		
        p % wgt = p % wgt * sum(XS_MG(p%material)%sig_scat(p%g,:)) / sig_tot
        p % g   = idx_group
        p % last_uvw(:) = p % coord(1)% uvw(:)
        p % coord(1)% uvw(:) = rand_vec()
		
		
        if (p%wgt < wgt_min_dyn) THEN !call Russian_Roulette(p)
            wgt_s = 2*wgt_min_dyn
            if ((p%wgt/wgt_s).ge.rang()) then
                p%wgt = wgt_s
            else
                p%alive = .false.
            endif 
        endif 	
	


	end subroutine 	
	
	
	
	
	
	
! ================================================== !
!	collision_PCQS_CE : 
! ================================================== !
subroutine collision_PCQS_CE (p)
    use constants, only: k_b
	use ace_reactions
    implicit none 
    type(particle), intent(inout) :: p
    integer :: iso, i, i_iso, xn, isab
    real(8) :: rn, el, noel, r, sigt_sum, temp, sum1, sum2
    real(8) :: micro_xs(6), macro_xs(5)
    ! * microscopic cross section
    ! 1 : total
    ! 2 : elastic
    ! 3 : absorption
    ! 4 : fission
    ! 5 : nufission
    ! 6 : thermal elastic
    real(8) :: ipfac
    integer :: ierg
    integer :: n_iso
    real(8) :: dtemp
    integer :: ii, jj, kk 
    real(8) :: xs_t(5)
	
	integer :: imat, ng, pt1, pt2, pt3, NE, i_group 
	real(8) :: speedn, sigtot_pcqs, wgt, beta, pdf 
	real(8) :: nu_del, nu 
	real(8) :: beta_g(8), lambda(8)
	real(8) :: val, fd 
	logical :: delayed 

    p%n_collision = p%n_collision + 1
    p % n_coord = 1
    xn = 1
    !===============================================
    ! Sample a target isotope in the mixture
    call WHAT_TEMPERATURE(p)
    macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
    rn = rang(); temp = 0
    do i = 1, materials(p%material)%n_iso
        dtemp = abs(p%kT-ace(materials(p%material)%ace_idx(i))%temp)
        if ( materials(p%material)%db .and. dtemp > K_B .and. p%E < 1D0 ) then
        ! On-the-fly Doppler broadening
        call GET_OTF_DB_MIC(p%kT,materials(p%material)%ace_idx(i),p%E,micro_xs)
        else
        ! point-wise data at the given temperature
        micro_xs = getMicroXS( materials(p%material)%ace_idx(i), p%E)
        end if
        ! S(a,b)
        call GET_SAB_MIC(materials(p%material),i,p%E,micro_xs)
        temp = temp + micro_xs(1)*materials(p%material)%numden(i)*barn
        if ( rn < temp/macro_xs(1) ) then
            iso = materials(p%material)%ace_idx(i)
            isab = materials(p%material) % sablist(i)
            i_iso = i
            if ( materials(p%material)%sab .and. ace(iso)%sab_iso /= 0 &
                .and. p%E < 4D-6 ) then
                p%yes_sab = .true.
            else
                p%yes_sab = .false.
            endif
            exit
        endif
    enddo

		! ===================================================================================================
		!  PCQS PART
		! ===================================================================================================
		speedn = sqrt(2.0d0*p%E*mevj/(m_u*m_n))*1.0d2   ! cm/s
		temp = 1.0d0 / (speedn*del_t) 
		
		sigtot_pcqs = macro_xs(1) + temp + (PKE_gamma / speedn) 
		wgt = p%wgt * macro_xs(1) / sigtot_pcqs  
		
		
		nu_del = getnudel(iso,p%E)
		nu = getnu(iso, p%E)
		ng = ace(iso)%NXS(8)
		if (ng > 8) then 
			print *, "delayed precursor group is larger than 8 :: ", ng, ace(iso)%library 
			stop
		endif 
		
		beta_g(:) = 0 
		lambda(:) = 0 
		! sample precursor group
		do i = 1, ng
			NE = ace(iso) % prcr( i ) % NE
			pt1 = 1; pt2 =  NE
			BS: do 
				if(pt2-pt1 == 1) exit BS
				pt3 = (pt2+pt1)/2 
				if (p%E >= ace(iso)%prcr(i)%E(pt3)) then 
					pt1 = pt3 
				else 
					pt2 = pt3
				endif 
			enddo BS
			
			ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%prcr(i)%E(pt1))/(ace(iso)%prcr(i)%E(pt1+1)-ace(iso)%prcr(i)%E(pt1))))
			pdf = ace(iso)%prcr(i)%F(pt1) + ipfac*(ace(iso)%prcr(i)%F(pt1+1)-ace(iso)%prcr(i)%F(pt1))
			beta_g(i) = pdf * nu_del/nu
			lambda(i) = ace(iso) % prcr( i ) % decay_const
		enddo
		beta = sum(beta_g)
		
		
		if (curr_cyc > n_pcqs_inact) then 
		
			!> PKE parameter tally 
			if (beta > 0) then 
				temp = 0
				rn = rang()
				i_group = ng
				do i = 1, ng
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i
						exit 
					endif 
				enddo 
				
				!여기서 sig_tot은 micro로
				!$omp atomic
				PKE_beta_tally1(i_group) = PKE_beta_tally1(i_group) + wgt*beta_g(i_group)*micro_xs(5) / (micro_xs(1) * w_tot)
				!$omp atomic
				PKE_beta_tally2(i_group) = PKE_beta_tally2(i_group) + wgt*micro_xs(5) / (micro_xs(1) * w_tot)
				!$omp atomic
				PKE_lambda_tally1(i_group) = PKE_lambda_tally1(i_group) + wgt*beta_g(i_group) & 
													* micro_xs(5) / (micro_xs(1)*w_tot)
				!$omp atomic
				PKE_lambda_tally2(i_group) = PKE_lambda_tally2(i_group) + wgt* beta_g(i_group)*micro_xs(5) / (lambda(i_group) * micro_xs(1)*w_tot)
			endif
			!$omp atomic
			PKE_gen_tally1 = PKE_gen_tally1 + wgt / (speedn * micro_xs(1) * w_tot)
			!$omp atomic
			PKE_gen_tally2 = PKE_gen_tally2 + wgt*micro_xs(5) / (micro_xs(1)*PKE_keff0 * w_tot)
			!$omp atomic
			PKE_Z_tally1 = PKE_Z_tally1 + wgt / (speedn * micro_xs(1))
		endif
		
		

		if (rang() < SSP) then 
		
		!> Source tally for the next ITERATION 
		val = (micro_xs(5)/micro_xs(1))*(1/PKE_keff0)
		if (micro_xs(4)>0) then 
		
			if (rang() < beta) then 
				delayed = .true. 
			else 
				delayed = .false. 
			endif 
			bank_idx = bank_idx + 1 
			thread_bank(bank_idx)%xyz 	= p%coord(1)%xyz
			thread_bank(bank_idx)%uvw 	= rand_vec()
			
			if (delayed) then 
				temp = 0
				rn = rang()
				i_group = ng
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				fd = exp(-lambda(i_group)*del_t)  &
					* (1-exp(lambda(i_group)*del_t)+lambda(i_group)*del_t*exp(lambda(i_group)*del_t))/(lambda(i_group)*del_t)  ! f_3,d
				thread_bank(bank_idx)%E 	= fission_E (p%E, iso,.true.,i_group)
				thread_bank(bank_idx)%wgt 	= wgt * fd * val * (1.0/SSP)
				thread_bank(bank_idx)%delayed 	= delayed
				
			else 
				thread_bank(bank_idx)%E 	= fission_E (p%E, iso,.false.)
				thread_bank(bank_idx)%wgt 	= wgt * val * (1.0/SSP)
				thread_bank(bank_idx)%delayed 	= delayed
			endif 
			
			
		endif 
        
        endif 
        
		
		
		!> Source tally for the next TIME-STEP
		if (curr_cyc == n_pcqs_totcyc ) then 
		
		if (rang() < SSP) then 
			!> For PCQS Neutron Source Initialization 
			! Source 1 
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 			= (1.0/SSP) * p%wgt * exp(PKE_gamma * del_t) * (sigtot_pcqs-macro_xs(1)-(PKE_gamma/speedn)) / sigtot_pcqs 
			thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 			= p%coord(1)%uvw
			thread_bank_init(init_idx)%E 			= p%E  ! incident E group 
			

			
			if (beta > 0) then 
				! Source 2
				! select delayed group 
				temp = 0 
				rn = rang()
				i_group = ng
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				
				val = (beta*micro_xs(5)/micro_xs(1))*(1/PKE_keff0)
				
				
				fd = (1-exp(-lambda(i_group)*del_t)*(1+lambda(i_group)*del_t))/(lambda(i_group)*del_t)  ! f_2,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 	= wgt * fd *val * (1.0/SSP)
				thread_bank_init(init_idx)%xyz 	= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 	= rand_vec()
				thread_bank_init(init_idx)%E 	= fission_E (p%E, iso,.true.,i_group)
			

			
			
				! Source 3 
				! select delayed group 
				rn = rang()
				temp = 0 
				i_group = ng
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				
				val = (beta*micro_xs(5)/(micro_xs(1)*lambda(i_group)))
				fd = exp(-lambda(i_group)*del_t)   ! f_1,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 			= wgt * fd * lambda(i_group) * val* (1.0/SSP)
				thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 			= rand_vec() 
				thread_bank_init(init_idx)%E 			= fission_E (p%E, iso,.true.,i_group)
				

				
			endif 
		endif 
		
		
		endif 
		! ===================================================================================================
		! ===================================================================================================	
	
	
	
    !!===============================================
    !Sampling reaction: elastic vs.non-elastic
    el   = micro_xs(2)
    noel = 0 
    call getierg(iso,ierg,p%E)
    ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
    do i = 1, ace(iso)%NXS(5) !> through the reaction types...
        if (abs(ace(iso)%TY(i)) == 19) cycle 
        noel = noel + ace(iso)%sig_MT(i)%cx(ierg) & 
                    + ipfac*(ace(iso)%sig_MT(i)%cx(ierg+1) - ace(iso)%sig_MT(i)%cx(ierg))
    enddo 

    r = rang()*(noel+el)-el
    if( ace(iso)%nxs(5) == 0 .or. r <= 0.0d0 ) then 
        if ( p%yes_sab .and. isab > 0 ) then
            call SAB_CE(p,iso,isab,micro_xs(2),micro_xs(6))
        elseif( p % yes_sab .and. isab < 0) then
            call SAB_THERM_CE(p, iso, abs(isab), micro_xs(2), micro_xs(6))
        else
            call elastic_CE (p, iso)
        end if
    else
        call notElastic_CE (p, iso, xn)
    end if
    
    p%wgt = p%wgt * ((el+noel)/micro_xs(1)) 
    
    !> (n, xn) reaction
    p%wgt = p%wgt * dble(xn)
    
    call absorption_CE(p)
    
end subroutine collision_PCQS_CE
	
	
	
	
	
	
subroutine collision_PCQS_CE_init (p)
    use constants, only: k_b
    implicit none 
    type(particle), intent(inout) :: p
    integer :: iso, i, i_iso, xn, isab
    real(8) :: rn, el, noel, r, sigt_sum, temp, sum1, sum2
    real(8) :: micro_xs(6), macro_xs(5)
    ! * microscopic cross section
    ! 1 : total
    ! 2 : elastic
    ! 3 : absorption
    ! 4 : fission
    ! 5 : nufission
    ! 6 : thermal elastic
    real(8) :: ipfac
    integer :: ierg
    integer :: n_iso
    real(8) :: dtemp
    integer :: ii, jj, kk 
    real(8) :: xs_t(5)
	
	integer :: imat, ng, pt1, pt2, pt3, NE, i_group 
	real(8) :: speedn, sigtot_pcqs, wgt, beta, pdf 
	real(8) :: nu_del, nu 
	real(8) :: beta_g(8), lambda(8)
	real(8) :: val, fd 
	logical :: delayed 
	
    p%n_collision = p%n_collision + 1
    p % n_coord = 1
    xn = 1
	
    !===============================================
    ! Sample a target isotope in the mixture
    call WHAT_TEMPERATURE(p)
    macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
    rn = rang(); temp = 0
    do i = 1, materials(p%material)%n_iso
        dtemp = abs(p%kT-ace(materials(p%material)%ace_idx(i))%temp)
        if ( materials(p%material)%db .and. dtemp > K_B .and. p%E < 1D0 ) then
        ! On-the-fly Doppler broadening
        call GET_OTF_DB_MIC(p%kT,materials(p%material)%ace_idx(i),p%E,micro_xs)
        else
        ! point-wise data at the given temperature
        micro_xs = getMicroXS( materials(p%material)%ace_idx(i), p%E)
        end if
        ! S(a,b)
        call GET_SAB_MIC(materials(p%material),i,p%E,micro_xs)
        temp = temp + micro_xs(1)*materials(p%material)%numden(i)*barn
        if ( rn < temp/macro_xs(1) ) then
            iso = materials(p%material)%ace_idx(i)
            isab = materials(p%material) % sablist(i)
            i_iso = i
            if ( materials(p%material)%sab .and. ace(iso)%sab_iso /= 0 &
                .and. p%E < 4D-6 ) then
                p%yes_sab = .true.
            else
                p%yes_sab = .false.
            endif
            exit
        endif
    enddo
    
    call fissionSite_CE(p, iso, micro_xs)
	
	
	
	
		! ===================================================================================================
		!  PCQS PART
		! ===================================================================================================
		wgt_min_dyn = wgt_min 
		wgt = p%wgt 
		speedn = sqrt(2.0d0*p%E*mevj/(m_u*m_n))*1.0d2   ! cm/s
	
	
		
		
		nu_del = getnudel(iso,p%E)
		nu = getnu(iso, p%E)
		ng = ace(iso)%NXS(8)
		if (ng > 8) then 
			print *, "delayed precursor group is larger than 8 :: ", ng, ace(iso)%library 
			stop
		endif 
		
		beta_g(:) = 0 
		lambda(:) = 0 
		! sample precursor group
		do i = 1, ng
			NE = ace(iso) % prcr( i ) % NE
			pt1 = 1; pt2 =  NE
			BS: do 
				if(pt2-pt1 == 1) exit BS
				pt3 = (pt2+pt1)/2 
				if (p%E >= ace(iso)%prcr(i)%E(pt3)) then 
					pt1 = pt3 
				else 
					pt2 = pt3
				endif 
			enddo BS
			
			ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%prcr(i)%E(pt1))/(ace(iso)%prcr(i)%E(pt1+1)-ace(iso)%prcr(i)%E(pt1))))
			pdf = ace(iso)%prcr(i)%F(pt1) + ipfac*(ace(iso)%prcr(i)%F(pt1+1)-ace(iso)%prcr(i)%F(pt1))
			beta_g(i) = pdf * nu_del/nu
			lambda(i) = ace(iso) % prcr( i ) % decay_const
		enddo
		beta = nu_del/nu
		
				
		
		

		if (rang() < SSP) then 
		
		!> Source tally for the next ITERATION 
		val = (micro_xs(5)/micro_xs(1))
		if (micro_xs(4)>0) then 
		
			if (rang() < beta) then 
				delayed = .true. 
			else 
				delayed = .false. 
			endif 
			split_idx = split_idx + 1 
			split_thread(split_idx)%xyz 	= p%coord(1)%xyz
			split_thread(split_idx)%uvw 	= rand_vec()
			
			if (delayed) then 
				temp = 0
				rn = rang()
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				fd = exp(-lambda(i_group)*del_t)  &
					* (1-exp(lambda(i_group)*del_t)+lambda(i_group)*del_t*exp(lambda(i_group)*del_t))/(lambda(i_group)*del_t)  ! f_3,d
				split_thread(split_idx)%E 		= fission_E (p%E, iso,.true.,i_group)
				split_thread(split_idx)%wgt 	= wgt * fd * val * (1.0/SSP)
				split_thread(split_idx)%delayed = delayed
				
			else 
				split_thread(split_idx)%E 		= fission_E (p%E, iso,.false.)
				split_thread(split_idx)%wgt 	= wgt * val * (1.0/SSP)
				split_thread(split_idx)%delayed = delayed
			endif 
			
		endif 
        
        endif 
        
		
		if (rang() < SSP) then 
			wgt = wgt / SSP 
			!> For PCQS Neutron Source Initialization 
			! Source 1 
			init_idx = init_idx + 1
			thread_bank_init(init_idx)%wgt 			= wgt/(speedn*del_t*macro_xs(1))
			thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
			thread_bank_init(init_idx)%uvw 			= p%coord(1)%uvw
			thread_bank_init(init_idx)%E 			= p%E  ! incident E group 

			if (beta > 0) then 
				! Source 2
				! select delayed group 
				temp = 0 
				rn = rang()
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				
				val = (beta*micro_xs(5)/micro_xs(1))
				
				
				fd = (1-exp(-lambda(i_group)*del_t)*(1+lambda(i_group)*del_t))/(lambda(i_group)*del_t)  ! f_2,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 	= wgt * fd *val 
				thread_bank_init(init_idx)%xyz 	= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 	= rand_vec()
				thread_bank_init(init_idx)%E 	= fission_E (p%E, iso,.true.,i_group)

				! Source 3 
				! select delayed group 
				rn = rang()
				temp = 0 
				do i = 1, ng 
					temp = temp + beta_g(i)
					if (rn < temp/beta) then 
						i_group = i  
						exit 
					endif 
				enddo 
				
				val = (beta*micro_xs(5)/(micro_xs(1)*lambda(i_group)))
				fd = exp(-lambda(i_group)*del_t)   ! f_1,d
				init_idx = init_idx + 1
				thread_bank_init(init_idx)%wgt 			= wgt * fd * lambda(i_group) * val
				thread_bank_init(init_idx)%xyz 			= p%coord(1)%xyz
				thread_bank_init(init_idx)%uvw 			= rand_vec() 
				thread_bank_init(init_idx)%E 			= fission_E (p%E, iso,.true.,i_group)

				endif 
		endif 
		
		
		! ===================================================================================================
		! ===================================================================================================	
	
    !!===============================================
    !Sampling reaction: elastic vs.non-elastic
    el   = micro_xs(2)
    noel = 0 
    call getierg(iso,ierg,p%E)
    ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
    do i = 1, ace(iso)%NXS(5) !> through the reaction types...
        if (abs(ace(iso)%TY(i)) == 19) cycle 
        noel = noel + ace(iso)%sig_MT(i)%cx(ierg) & 
                    + ipfac*(ace(iso)%sig_MT(i)%cx(ierg+1) - ace(iso)%sig_MT(i)%cx(ierg))
    enddo 

    r = rang()*(noel+el)-el
    if( ace(iso)%nxs(5) == 0 .or. r <= 0.0d0 ) then 
        if ( p%yes_sab .and. isab > 0 ) then
            call SAB_CE(p,iso,isab,micro_xs(2),micro_xs(6))
        elseif( p % yes_sab .and. isab < 0) then
            call SAB_THERM_CE(p, iso, abs(isab), micro_xs(2), micro_xs(6))
        else
            call elastic_CE (p, iso)
        end if
    else
        call notElastic_CE (p, iso, xn)
    end if
    
    p%wgt = p%wgt * ((el+noel)/micro_xs(1)) 
    
    !> (n, xn) reaction
    p%wgt = p%wgt * dble(xn)
    
    call absorption_CE(p)
    

end subroutine collision_PCQS_CE_init
	
	
	
	
	
	subroutine set_PCQS_source() 
		integer :: i, isize, jsize, rcv_int, size_delay
		real(8) :: wgt_g(7), wgt_prompt, wgt_delay, rcv_buf, size_prompt, wgt_1, wgt_2
		
		!wgt_g(:) = 0
		!wgt_1 = 0 
		!wgt_2 = 0 
		!
		!
		!wgt_prompt = sum(prompt_bank(:)%wgt)
		!wgt_delay = sum(delayed_bank(:)%wgt)
		!size_prompt = size(prompt_bank)
		!size_delay = size(delayed_bank)
		!
		!
		!do i = 1, size(delayed_bank) 
		!	if (delayed_bank(i)%delayed) then 
		!		wgt_1 = wgt_1 + delayed_bank(i)%wgt
		!	else 
		!		wgt_2 = wgt_2 + delayed_bank(i)%wgt
		!	endif 
		!enddo 
		!call MPI_REDUCE(wgt_1,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
		!wgt_1 = rcv_buf
		!call MPI_REDUCE(wgt_2,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
		!wgt_2 = rcv_buf
		!
		!
		!
		!call MPI_REDUCE(wgt_prompt,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
		!wgt_prompt = rcv_buf
		!call MPI_REDUCE(wgt_delay,rcv_buf,1,MPI_REAL8,MPI_SUM,score,MPI_COMM_WORLD,ierr)
		!wgt_delay = rcv_buf
		!call MPI_REDUCE(size_delay,rcv_int,1,MPI_INTEGER,MPI_SUM,score,MPI_COMM_WORLD,ierr)
		!size_delay = rcv_int
		
		
		
		! prompt_bank는 고정 
		! delayed_bank를 tally 
		! fission_bank에 합침
		
		if (allocated(temp_bank)) deallocate(temp_bank) 
		if (allocated(fission_bank)) deallocate(fission_bank) 
		



		
		if (allocated(delayed_bank)) then 
			!if (icore==score) write(prt_prec,'(I,4E15.5,I,2E15.5)') curr_cyc, wgt_prompt, wgt_delay, wgt_1, wgt_2, size_delay, cyc_power, wgt_min_dyn
			isize = size(delayed_bank) 
			jsize = size(prompt_bank) 
			allocate(temp_bank(isize+jsize)) 
						
			temp_bank(1:isize) = delayed_bank(:)
			temp_bank(isize+1:isize+jsize) = prompt_bank(:)
			
			deallocate(delayed_bank) 
		else 
			!if (icore==score)  write(prt_prec,'(I,2E15.5,2I)') curr_cyc, sum(prompt_bank(:)%wgt), 0, size(prompt_bank), 0
			jsize = size(prompt_bank) 
			allocate(temp_bank(jsize))
			temp_bank = prompt_bank		
		endif 
		
		
		
		!> Combing the source 
		isize = int(ngen/ncore)
		call combing_source(temp_bank, fission_bank, isize)
		deallocate(temp_bank) 
		wgt_min_dyn = 0.25 * sum(fission_bank(:)%wgt) / size(fission_bank)
		
		
		!do i = 1, isize 
		!	if (E_mode == 1 .and. fission_bank(i)%E == 0 ) then 
		!		print *, icore, i, fission_bank(i)%E,fission_bank(i)%wgt,fission_bank(i)%xyz, 'fission_bank'
		!		stop 
		!	endif 
		!enddo 
		
		
	end subroutine 
	
	
	
	subroutine combing_source(src_in, src_out, M) 
		integer, intent(in) :: M 
		type(bank), intent(in) :: src_in(:) 
		type(bank), allocatable :: src_out(:)
		integer :: i, j, isize, idx
		real(8) :: w_av, tooth, val 
		
		
		if (allocated(src_out)) deallocate(src_out) 
		allocate(src_out(M)) 
		isize = size(src_in) 
		w_av = sum(src_in(:)%wgt) / real(M,8) 
		
		
		val = src_in(1)%wgt 
		idx = 1
		tooth = w_av * rang()
		
		do i = 1, M
			do j = idx, isize
				if (tooth < val) then 
					src_out(i) = src_in(j) 
					idx = j 
					exit
				endif 
				
				val = val + src_in(j+1)%wgt
			enddo 
			tooth = tooth + w_av
		enddo 
		src_out(:)%wgt = w_av
		
		
		
		
		
		!val = 0; idx = 1
		!tooth = w_av * rang() 
		!do i = 1, isize 
		!	val = val + src_in(i)%wgt
		!	if (tooth < val) then 
		!		src_out(idx) = src_in(i) 
		!		tooth = tooth + w_av
		!		idx = idx + 1 
		!		if (idx > M) return 
		!	endif 
		!enddo 
		
		
	end subroutine 
	
	
	
	subroutine combing_source_par (src_in, src_out, M) 
		use omp_lib 
		implicit none 
		
		integer, intent(in) :: M 
		type(bank), intent(in) :: src_in(:) 
		type(bank), allocatable :: src_out(:)
		real(8), allocatable :: src_c(:)
		integer :: i, j, isize, idx
		real(8) :: w_av, tooth, val, rn 
		
		integer :: tid, ista, iend, nthreads
		integer :: pt1, pt2, pt3 
		real(8) :: time1, time2, time3

		
		if (allocated(src_out)) deallocate(src_out) 
		allocate(src_out(M)) 
		isize = size(src_in) 
		w_av = sum(src_in(:)%wgt) / real(M,8) 
		
		allocate(src_c(0:isize)) 
		src_c(0) = 0
		do i = 1, isize 
			src_c(i) = src_c(i-1)+src_in(i)%wgt
		enddo 
		
		call random_number(rn)
		
		!$OMP PARALLEL private(tid, ista, iend, pt1, pt2, pt3, val, tooth) 
		tid = omp_get_thread_num() 
		nthreads = omp_get_num_threads() 
		ista = tid * M / nthreads + 1
		iend = (tid+1)*M / nthreads  
		
		tooth = w_av * rn + w_av * (ista-1) 
		!if (icore==score) print *, tid, tooth 
		! binary search로 각 thread마다 시작지점 찾음 (idx) 
        pt1 = 0; pt2 = isize
        do 
            if(pt2-pt1 == 1) exit 
            pt3 = (pt2+pt1)/2 
            if (tooth >= src_c(pt3) ) then 
                pt1 = pt3 
            else 
                pt2 = pt3
            endif
        enddo
		
		val = src_c(pt2) 
		do i = ista, iend
			inner: do j = pt2, isize
				if (tooth < val) then 
					src_out(i) = src_in(j) 
					pt2 = j 
					exit inner
				endif 
				val = val + src_in(j+1)%wgt
			enddo inner
			tooth = tooth + w_av
		enddo 
		
		
		
		!$OMP END PARALLEL
		
		deallocate(src_c)
		src_out(:)%wgt = w_av
		
		
		
	end subroutine	
	
	
	
	
	

	
	
subroutine GEM(A,B,x,n)
    integer :: i, k 
	integer, intent(in) :: n 
	real(8) :: A(n,n), B(n), x(n) 
	
	x = B

    do k = 1, n-1
        a(k+1: n, k) = a(k+1: n, k) / a(k, k)

        a(k+1: n, k+1: n) = a(k+1: n, k+1: n) - &
                matmul(a(k+1: n, k: k), a(k: k, k+1: n))
    end do

    do i = 1, n
        x(i) = x(i) - dot_product(a(i, 1: i-1), x(1: i-1))
    end do

    do i = n, 1, -1
        x(i) = x(i) - dot_product(a(i, i+1: n), x(i+1: n))
        x(i) = x(i) / a(i, i)
    end do

end subroutine GEM
	
	
	
	


end module
