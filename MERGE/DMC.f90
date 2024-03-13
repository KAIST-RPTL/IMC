module DMC 

	use transient 
	use bank_header
	use particle_header, only : particle 
	use ace_header, 	only : ace, CrossSectionDataForm, EnergyDist
	use variables!, 		only : E_mode, ngen, icore,score, n_act
	
	
implicit none 

	

contains 

	
!===============================================================================
! SET_DYNAMC_BANK - Set dynamic_bank for transient simulation source.
!===============================================================================
	subroutine set_dynamic_bank() 
		integer :: i, j, isize, size_bank, ierr
		real(8) :: t, t0, t1, temp, beta, w_nav
		real(8) :: lambda_b
		real(8) :: av_wgtp, av_wgtd
		real(8) :: rn, temp1(npg), val 
		integer :: delayed_group 
		integer :: np, nd
		real(8), allocatable :: fd(:), Tk(:), w_timed(:), time(:) 
		real(8) :: w, w_av, comb
		type(Bank), allocatable :: comb_bank(:)
		type(PrecBank), allocatable :: comb_prec(:)
		integer :: M, idx
		integer, allocatable :: comb_idx(:) 
		
		
		M = int(real(ngen,8) / 2.0)
		size_bank = 2*M 
		
		! 4. ¥Ÿ¿Ω dynamic_bank∏¶ ∏∏µÎ (prompt + forced decay precursor bank)
		
		allocate(dynamic_bank(1:size_bank))
		
		
		if (icore == score) then 
			isize = size(prec_bank)
			allocate(fd(1:npg))
			allocate(Tk(1:isize))
			allocate(w_timed(1:isize))
			allocate(time(1:isize))
			allocate(comb_idx(1:M))
			
			Tk(:) = 0 
			! 4.1 sample src from prec_bank 
			do i = 1, isize 
				t0 = prec_bank(i)%time
				t1 = curr_time !(curr_timestep-1) *del_t
				t  = t1 + rang() * del_t
				time(i) = t
				
				beta = sum(prec_bank(i)%beta(:))
				temp = 0
				do j = 1, npg
					temp = temp + prec_bank(i)%beta(j) / prec_bank(i)%lambda(j)
				enddo 
				lambda_b = beta / temp
				
				! Define fd 
				if (t0 == 0) then
					do j = 1, npg
						fd(j) = (lambda_b/prec_bank(i)%lambda(j))*(prec_bank(i)%beta(j) / beta)
					enddo 
				else 
					do j = 1, npg
						fd(j) = prec_bank(i)%beta(j) / beta
					enddo 
				endif 
				
				do j = 1, npg 
					Tk(i) = Tk(i) + fd(j) * exp(-prec_bank(i)%lambda(j)*(t-t0))
				enddo 
				
				w_timed(i) = Tk(i) * prec_bank(i)%wgt
				
			enddo
			
			w_av = sum(w_timed)/ real(M,8)
			
			allocate(comb_prec(1:M))
			comb = rang()*w_av
			idx = 1 
			do i = 1, M 
				PRC : do j = idx, isize
					w = sum(w_timed(1:j)) 
					if (comb < w) then 
						idx = j 
						exit PRC
					endif 
				enddo PRC
				comb_prec(i) = prec_bank(idx)
				comb_idx(i) = idx
				comb = comb + w_av
				comb_prec(i)%wgt = w_av / Tk(idx)
			enddo 
			
			deallocate(prec_bank) 
			call move_alloc(comb_prec, prec_bank) 
			
			!print *, 'prec comb complete', size(prec_bank), prec_bank(1)%wgt
			
			do i = 1, M
				t0 = prec_bank(i)%time
				t1 = curr_time !(curr_timestep-1) *del_t
				t  = time(comb_idx(i))
				
				
				beta = sum(prec_bank(i)%beta(:))
				temp = 0
				do j = 1, npg
					temp = temp + prec_bank(i)%beta(j) / prec_bank(i)%lambda(j)
				enddo 
				lambda_b = beta / temp
				
				! Define fd 
				if (t0 == 0) then
					do j = 1, npg
						fd(j) = (lambda_b/prec_bank(i)%lambda(j))*(prec_bank(i)%beta(j) / beta)
					enddo 
				else 
					do j = 1, npg
						fd(j) = prec_bank(i)%beta(j) / beta
					enddo 
				endif 
				
				temp = 0; 
				do j = 1, npg
					val = fd(j) * prec_bank(i)%lambda(j)*exp(-prec_bank(i)%lambda(j)*(t-t0))
					temp = temp + val
					temp1(j) = val
				enddo
				
				rn = rang() 
				delayed_group = 0
				do j = 1, npg
					if (rn < sum(temp1(1:j))/temp) then 
						delayed_group = j
						exit
					endif 
				enddo 
				
				dynamic_bank(i)%wgt 	= prec_bank(i)%wgt*del_t*temp
				dynamic_bank(i)%xyz(:) 	= prec_bank(i)%xyz(:)
				dynamic_bank(i)%uvw(:) 	= rand_vec()
				dynamic_bank(i)%time 	= t
				dynamic_bank(i)%delayed	= .true.
				
				!write (prt_dynamic, '(I, 4e14.6)') i, dynamic_bank(i)%wgt, prec_bank(i)%wgt, del_t, temp				
				
				!DMC_prod = DMC_prod + dynamic_bank(i)%wgt
				
				if (E_mode == 1) then
					dynamic_bank(i)%E = fission_E (prec_bank(i)%E,prec_bank(i)%idx,.true.,delayed_group)
				else 
					dynamic_bank(i)%G = fission_G (prec_bank(i)%idx,.true., delayed_group) 
				endif 
				
			enddo
			deallocate(fd, Tk, time, w_timed, comb_idx) 
			!print *, 'delayed neutron sample complete', dynamic_bank(1)%wgt
			
			
			
			
			print '(a,2I,2E13.5)', 'asdfsfas', size(prompt_bank), isize , sum(prompt_bank(:)%wgt)/real(size_bank-M,8), sum(dynamic_bank(1:M)%wgt)/real(M,8)
			
			
			! Combing technique for time source
			psize = size(prompt_bank)
			
			if (psize < 1) then 
				print *, "prompt_bank size is zero"
				stop 
			endif 
			
			w_av = sum(prompt_bank(:)%wgt) / real(M,8)
			allocate(comb_bank(1:M))
			comb = rang()*w_av
			idx = 1 
			do i = 1, M 
				w = sum(prompt_bank(1:idx-1)%wgt)
				PRT : do j = idx, psize
					w = w + prompt_bank(j)%wgt
					if (comb < w) then 
						idx = j 
						exit PRT
					endif 
				enddo PRT
				comb_bank(i) = prompt_bank(idx)
				comb = comb + w_av
			enddo 
			comb_bank(:)%wgt = w_av
			deallocate(prompt_bank) 
			call move_alloc(comb_bank, prompt_bank) 
			!print *, 'time source comb complete', size(prompt_bank), prompt_bank(1)%wgt, size_bank
			
			
			! 4.3 Rest of the dynamic_bank is from prompt_bank
			do i = M+1, size_bank 
				dynamic_bank(i)%wgt 	= prompt_bank(i-M)%wgt
				dynamic_bank(i)%xyz(:) 	= prompt_bank(i-M)%xyz(:)
				dynamic_bank(i)%uvw(:) 	= prompt_bank(i-M)%uvw(:)
				dynamic_bank(i)%E 		= prompt_bank(i-M)%E
				dynamic_bank(i)%G 		= prompt_bank(i-M)%G
				dynamic_bank(i)%time	= prompt_bank(i-M)%time
				dynamic_bank(i)%delayed	= .false.
			enddo 
			deallocate(prompt_bank)
			
			av_wgtp = sum(dynamic_bank(M+1:size_bank)%wgt)/real(size_bank-M,8)
			av_wgtd = sum(dynamic_bank(1:M)%wgt)/real(M,8)
			
			
			
			!print *, av_wgtp, av_wgtd
			
			!> Set weight windows
			!wgt_min_dyn   = 0.25*av_wgtp
			!wgt_split_dyn = av_wgtd
			wgt_min_dyn   = 0.05*av_wgtp
			wgt_split_dyn = 2.0*av_wgtp
			
			! Shuffle dynamic_bank
			allocate(temp_bank(1:size_bank))
			allocate(comb_idx(1:size_bank))
			
			call scramble( comb_idx )
			
			do i = 1, size_bank 
				temp_bank(i) = dynamic_bank(comb_idx(i))
			enddo 
			deallocate(dynamic_bank)
			call move_alloc(temp_bank, dynamic_bank)
			
		else
			deallocate(prec_bank, prompt_bank)
		endif 
		
		
		if (icore/=score) allocate(prec_bank(1:M))
		call MPI_BCAST(dynamic_bank, size_bank, MPI_bank, score, MPI_COMM_WORLD, ierr) 
		call MPI_BCAST(prec_bank, M, MPI_precbank, score, MPI_COMM_WORLD, ierr) 
		
		call MPI_BCAST(wgt_min_dyn, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		call MPI_BCAST(wgt_split_dyn, 1, MPI_REAL8, score, MPI_COMM_WORLD, ierr) 
		
		!av_wgtp = sum(dynamic_bank(M+1:size_bank)%wgt)/real(size_bank-M,8)
		!av_wgtd = sum(dynamic_bank(1:M)%wgt)/real(M,8)
		!w_c = av_wgtp / av_wgtd
		w_c = sum(prec_bank(:)%wgt)/size(prec_bank)
		psize = M 
		
		!if (icore==score) print *, size(dynamic_bank), size(prec_bank), w_c, sum(dynamic_bank(:)%wgt)
		

		
	end subroutine
	
	
	
	
	! ================================================== !
	!    collision_dynamic_MG()  
	! ================================================== !
    subroutine collision_dynamic_MG(p)
        type(Particle), intent(inout) :: p
        real(8) :: sig_tot, temp, rnum, wgt_s, uvw_temp(3)
        integer :: i,j, i_group, idx_group, n, bsize, imat
		real(8) :: beta
		real(8) :: t, t0, val, lambda_b
		integer :: isize 
		
        p % n_collision = p % n_collision + 1
        p % n_coord = 1
        
        imat = p%material
		beta = sum(MGD(imat)%beta(:))
        sig_tot = sum(XS_MG(imat)%sig_scat(p%g,:)) + XS_MG(imat)%sig_abs(p%g)
        
		
		temp = 0 
		do i = 1, npg
			temp = temp + MGD(p%material)%beta(i) / MGD(p%material)%lambda(i)
		enddo 
		lambda_b = beta / temp
		
        !> Precursor bank add
		temp = (1.0/w_c)*p%wgt*beta*XS_MG(imat)%sig_fis(p%g)*XS_MG(imat)%nu(p%g)/sig_tot
        !temp = p%wgt*(beta/lambda_b)*(XS_MG(imat)%nu(p%g)*XS_MG(imat)%sig_fis(p%g)/sig_tot)
        if (rang() < temp) then
			!print *, 'precursor sampled ', prec_idx
			!if (curr_timestep > 1) then 
			!	print '(a,3f10.3)', 'precursor sampled', p%wgt, temp, wgt_split_dyn
			!endif
			prec_idx = prec_idx + 1
			n = size(MGD(imat)%beta(:))
			prec_thread(prec_idx)%wgt 			= w_c
			prec_thread(prec_idx)%xyz 			= p%coord(1)%xyz
			prec_thread(prec_idx)%time 			= p%time
			prec_thread(prec_idx)%G 			= p%g
			prec_thread(prec_idx)%idx 			= p%material
			prec_thread(prec_idx)%beta(1:n)		= MGD(imat)%beta(:)
			prec_thread(prec_idx)%lambda(1:n)	= MGD(imat)%lambda(:)
			

			! ==========================================================
			! A precursor produced from the present time step 
			! must produce a delayed neutron within the timestep 
			! ==========================================================
			! Add to split bank 
			!t0 = p%time
			!t  = t0 + rang() * (curr_time + del_t - t0)
			!
			!split_idx = split_idx + 1
			!split_thread(split_idx)%xyz 	= p%coord(1)%xyz
			!split_thread(split_idx)%uvw 	= rand_vec()
			!split_thread(split_idx)%E		= real(p%material,8)
			!split_thread(split_idx)%G		= fission_G(imat,.true.)
			!split_thread(split_idx)%delayed = .true.
			!split_thread(split_idx)%time 	= t
			!
			!
			!temp = 0;
			!do j = 1, npg
			!	val = (MGD(imat)%beta(j)/beta) * MGD(imat)%lambda(j) &
			!			*exp(-MGD(imat)%lambda(j)*(t-t0))
			!	temp = temp + val 
			!enddo
			!
			!split_thread(split_idx)%wgt = (curr_time + del_t - p%time) * temp
			
			
			!if buffer is almost full -> add to the split_bank_temp
			if ( split_idx > 14000 ) then 
			  !$omp critical
				isize = size(split_bank_temp)
				if(allocated(temp_bank)) deallocate(temp_bank)
				allocate(temp_bank(1:isize+split_idx)) 
				if (isize>0) temp_bank(1:isize) = split_bank_temp(:)
				deallocate(split_bank_temp)
				temp_bank(isize+1:isize+split_idx) = split_thread(1:split_idx)
				call move_alloc(temp_bank, split_bank_temp)
			  !$omp end critical
				split_idx = 0
			endif
			
			
        endif
        
		
        
		temp = (1.0-beta)*XS_MG(imat)%sig_fis(p%g)*XS_MG(imat)%nu(p%g) + sum(XS_MG(imat)%sig_scat(p%g,:))
		!temp = XS_MG(i)%sig_fis(p%g)*XS_MG(i)%nu(p%g) + sum(XS_MG(i)%sig_scat(p%g,:))

		!> Collision estimator 
		!$omp atomic 
		DMC_loss = DMC_loss + p%wgt * XS_MG(imat)%sig_abs(p%g) / sig_tot
		
		!$omp atomic 
		DMC_prod = DMC_prod + p%wgt * XS_MG(imat)%sig_fis(p%g)*XS_MG(imat)%nu(p%g) / sig_tot
		
		
		
		!> Determine collision type (improved branchless method)
		if (rang() < sum(XS_MG(imat)%sig_scat(p%g,:)/temp)) then ! scattering
			rnum = rang()
			do i_group = 1, size(XS_MG(p%material)%sig_scat(p%g,:))
				if (rnum < sum(XS_MG(p%material)%sig_scat(p%g,1:i_group))/sum(XS_MG(p%material)%sig_scat(p%g,:))) then 
					idx_group = i_group
					exit
				endif 
			enddo 
			
		else ! prompt fission 
            rnum = rang()
            do i_group = 1, size(XS_MG(p%material)%chi(:))
                if (rnum < sum(XS_MG(p%material)%chi(1:i_group))/sum(XS_MG(p%material)%chi(:))) then 
                    idx_group = i_group
                    exit
                endif 
            enddo
		endif 
		
		p % g   = idx_group
		p % last_uvw(:) = p % coord(1)% uvw(:)
		p % coord(1)% uvw(:) = rand_vec()
		! improved branch-less
        p % wgt = p % wgt * temp / sig_tot
		
		!call Russian_Roulette(p)
        if (p%wgt < wgt_min_dyn) THEN 
            wgt_s = 2*wgt_min_dyn
            if ((p%wgt/wgt_s).ge.rang()) then
                p%wgt = wgt_s
            else
                p%alive = .false.
            endif 
        endif 
        
    end subroutine collision_dynamic_MG	
	
	
	
	

	! ====================================================================== !
	!	collision_dynamic_CE : Collision simulation for dynamic MC simulation 
	! ====================================================================== !
	subroutine collision_dynamic_CE (p) 
		type(particle), intent(inout) :: p
		integer :: iso, i, i_iso, xn, isab
		real(8) :: rn, el, inel, r, sigt_sum, temp, sum1, sum2
		real(8) :: micro_xs(6), macro_xs(5)
		! * microscopic cross section
		! 1 : total
		! 2 : elastic
		! 3 : absorption
		! 4 : fission
		! 5 : nufission
		! 6 : thermal elastic
		real(8) :: ipfac, pdf
		integer :: ierg
		integer :: nsplit, isize, n_iso 
		
		real(8) :: sig_sum
		real(8), allocatable ::  sig_arr(:)
		type (CrossSectionDataForm), pointer :: sigmt
		type (EnergyDist),  pointer :: eg
		integer :: iMT, MT
		integer :: pt1, pt2, pt3
		integer :: law, ilaw        ! collision law
		real(8) :: F         ! collision probability
		real(8) :: erg_out, mu
		real(8) :: val
		real(8) :: beta, lambda_b
		real(8) :: beta_g(8), lambda(8)
		integer :: ng
		integer :: NE
		real(8) :: nu_del, nu
		real(8) :: wgt_s

		
		p%n_collision = p%n_collision + 1
		p % n_coord = 1
		xn = 1
		!===============================================
		! Sample a target isotope in the mixture
		macro_xs = getMacroXS(materials(p%material), p%E, p%kT, p%urn)
		rn = rang(); temp = 0 
		n_iso = materials(p%material)%n_iso
		iso = materials(p%material)%ace_idx(n_iso)
        isab = materials(p%material) % sablist(i)
		do i = 1, n_iso
			micro_xs = getMicroXS( materials(p%material)%ace_idx(i), p%E)
			! S(a,b)
			call GET_SAB_MIC(materials(p%material),i,p%E,micro_xs,p%kT)
			temp = temp + micro_xs(1)*materials(p%material)%numden(i)*barn
			if ( rn < temp/macro_xs(1) ) then
				iso = materials(p%material)%ace_idx(i)
				i_iso = i
				exit
			endif
		enddo
		
		p%iso = iso 
		
		if ( materials(p%material)%sab .and. ace(iso)%sab_iso /= 0 &
			.and. p%E < 4D-6 ) then 
			p%yes_sab = .true.
		else 
			p%yes_sab = .false.
		endif
		
		! Precursor bank add ==============================
		
		! Calculate beta_g & lambda 
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
		
		
		
		beta = sum(beta_g(1:ng))
		temp = 0 
		do i = 1, ng 
			temp = temp + beta_g(i) / lambda(i)
		enddo 
		lambda_b = beta / temp	
		
		 !> Precursor bank add
		temp = (1.0/w_c)*p%wgt*beta*micro_xs(5)/micro_xs(1)

		if (rang() < temp) then
		
			prec_idx = prec_idx + 1
			prec_thread(prec_idx)%wgt 			= w_c
			prec_thread(prec_idx)%xyz 			= p%coord(1)%xyz
			prec_thread(prec_idx)%time 			= p%time
			prec_thread(prec_idx)%E 			= p%E
			prec_thread(prec_idx)%idx 			= iso
			prec_thread(prec_idx)%beta(1:ng)	= beta_g(1:ng)
			prec_thread(prec_idx)%lambda(1:ng)	= lambda(1:ng)
			
			
			!if buffer is almost full -> add to the split_bank_temp
			if ( split_idx > 14000 ) then 
			  !$omp critical
				isize = size(split_bank_temp)
				if(allocated(temp_bank)) deallocate(temp_bank)
				allocate(temp_bank(1:isize+split_idx)) 
				if (isize>0) temp_bank(1:isize) = split_bank_temp(:)
				deallocate(split_bank_temp)
				temp_bank(isize+1:isize+split_idx) = split_thread(1:split_idx)
				call move_alloc(temp_bank, split_bank_temp)
			  !$omp end critical
				split_idx = 0
			endif
			
		endif
		
		
		!!===============================================
		!Sampling reaction: elastic vs.non-elastic
		el   = micro_xs(2)
		inel = 0 
		call getierg(iso,ierg,p%E)
		ipfac = max(0.d0, min(1.d0,(p%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
		do i = 1, ace(iso)%NXS(5) !> through the reaction types...
			if (abs(ace(iso)%TY(i)) == 19 .or. ace(iso)%TY(i)==0) cycle 
			inel = inel + ace(iso)%sig_MT(i)%cx(ierg) & 
						+ ipfac*(ace(iso)%sig_MT(i)%cx(ierg+1) - ace(iso)%sig_MT(i)%cx(ierg))
		enddo 

		
		
		!$omp atomic 
		DMC_loss = DMC_loss + p%wgt * macro_xs(2) / macro_xs(1)
		
		!$omp atomic 
		DMC_prod = DMC_prod + p%wgt * macro_xs(4) / macro_xs(1)
		
		
		val = (1.0-beta)*micro_xs(5) + el + inel 
		r = rang()
		if( r < (el)/val ) then ! elastic scattering
            if ( p%yes_sab .and. isab > 0 ) then
                call SAB_CE(p,iso,isab,micro_xs(2),micro_xs(6))
            elseif( p % yes_sab .and. isab < 0) then
                call SAB_THERM_CE(p, iso, abs(isab), micro_xs(2), micro_xs(6))
            else
                call elastic_CE (p, iso)
            end if
		elseif (r < (el+inel)/val) then ! inelastic scattering
			call inElastic_CE (p,iso,xn)
			p%wgt = p%wgt * dble(xn)
		else ! prompt fission 
			p%E = fission_E (p%E, iso,.false.)
			do i = 1,p%n_coord
				p%coord(i)%uvw(:) = rand_vec()
			enddo
			!p%n_coord = 1
		endif
		
		
		
		p%wgt = p%wgt * (val/micro_xs(1))
		
		
		if (p%wgt < wgt_min_dyn) THEN !call Russian_Roulette(p)
			wgt_s = 2*wgt_min_dyn
			if ((p%wgt/wgt_s).ge.rang()) then
				p%wgt = wgt_s
			else
				p%alive = .false.
			endif
		endif
		
		
	end subroutine


end module 
