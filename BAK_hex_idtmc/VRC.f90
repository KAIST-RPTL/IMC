module VRC 
	use constants, 			only: wgt_min, TINY_BIT, INFINITY, K_B, barn
    use particle_header,    only: particle
	use material_header,	only: materials 
	use XS_header, 			only: XS_MG
	use variables, 			only: E_mode, k_vrc, fiss_vrc, loss_vrc
	use ace_xs, 			only: getMacroXS, getMicroXS, getnudel 
	use geometry, 			only: distance_to_boundary, find_cell, cross_surface
	use randoms, 			only: rang, rand_vec
    use omp_lib
	use bank_header, 		only: vrc_thread, vrc_idx
	use tally, 				only: mesh_distance
	use ace_xs,				only: getMicroXS, getierg, get_otf_db_mic, get_sab_mic
	use ace_header, 		only: ace 
	use ace_reactions, 		only: elastic_CE, inElastic_CE, WHAT_TEMPERATURE, SAB_CE, fission_E
	
	implicit none 
	
	type RayType
		real(8) :: wgt
		real(8) :: xyz(3)
		real(8) :: uvw(3) 
		real(8) :: E 
		integer :: G 
		real(8), allocatable :: Sigt(:)
		real(8), allocatable :: nuSigf(:) 
        real(8) :: kT
		
		contains
        procedure :: reset => reset_ray
	endtype 
	type(RayType) :: RayBuffer(1000) 
	
	integer :: m_pseudo 
	
	
	contains 
	
	!==============================================================
	! TRACE_PSUDORAY 
	!==============================================================
	subroutine trace_psudoray(p)
        type(Particle), intent(inout) :: p
		real(8) :: val
        integer :: j                      ! coordinate level
        integer :: surface_crossed        ! surface which particle is on
        real(8) :: d_boundary             ! distance to nearest boundary
        logical :: found_cell             ! found cell which particle is in?
        real(8) :: macro_xs(5)
        integer :: i_cell
		real(8) :: wgt0, wgt_s, wgt
		integer :: iter 
		
		wgt0 = p%wgt
		val = 0 
		!print *, 'psudo-ray start'
		!found_cell = .false.
		call find_cell(p, found_cell, i_cell)
        do while (p%alive == .true.)
			! calculate reduced weight 
			!if (.not. found_cell) call find_cell(p, found_cell, i_cell)
			!print *, p%coord(1)%xyz(1), XS_MG(p%material)%mat_id
			if (E_mode == 0) then 
				macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
				macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
				macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
			elseif (E_mode == 1) then 
				macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
			endif 
			
			! Sample a distance to boundary
			call distance_to_boundary(p, d_boundary, surface_crossed)
			p%wgt = p%wgt*exp(-val)
			val = d_boundary * macro_xs(1)
			!> Volumetric-ray-casting estimator
			!$omp atomic
			fiss_vrc = fiss_vrc + p%wgt*macro_xs(4)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
			!$omp atomic
			loss_vrc = loss_vrc + p%wgt*macro_xs(2)*(1-exp(-macro_xs(1)*d_boundary))/macro_xs(1)
			!> Advance particle
			do j = 1, p % n_coord
				p % coord(j) % xyz = p % coord(j) % xyz + d_boundary * p % coord(j) % uvw
			enddo
			
			call cross_surface(p, surface_crossed)

			!if (p%wgt < 0.000000001) then 
			!	p%alive = .false.
			!endif 
			! Russian Roulette for psudo-ray
			if (p%wgt < 0.0001) THEN !call Russian_Roulette(p)
				wgt_s = 2*0.0001
				if ((p%wgt/wgt_s).ge.rang()) then
					p%wgt = wgt_s
				else
					p%alive = .false.
				endif 
			endif 
			
		enddo 
		
	end subroutine	
	
	
	
	subroutine trace_psudoray_tet_vrc(p, dist)
        type(Particle), intent(inout) :: p
		real(8), intent(in) :: dist
		real(8) :: val
        integer :: j                      ! coordinate level
        integer :: surface_crossed        ! surface which particle is on
        real(8) :: d_boundary             ! distance to nearest boundary
        logical :: found_cell             ! found cell which particle is in?
        real(8) :: macro_xs(5)
        integer :: i_cell
		real(8) :: wgt0, wgt_s, wgt
		integer :: iter 
		real(8) :: dist_accum, d_mesh, distance, ddiff
		logical :: inside_mesh
		integer :: income_mesh, i_surf, i_xyz(3)
		
		dist_accum = 0 
		wgt0 = p%wgt
		val = 0 
		call find_cell(p, found_cell, i_cell)
        do while (p%alive == .true.)
			! calculate reduced weight 
			if (E_mode == 0) then 
				macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
				macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
				macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
			elseif (E_mode == 1) then 
				macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
			endif 
			
			! Sample a distance to boundary
			call distance_to_boundary(p, d_boundary, surface_crossed)
			p%wgt = p%wgt*exp(-val)
			val = d_boundary * macro_xs(1)
			
			! Sample a distance to mesh
			d_mesh = INFINITY
			inside_mesh = .false. 
			call MESH_DISTANCE(p,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)
			
			!> minimum distance
			ddiff = abs(d_boundary-d_mesh)/d_boundary
			if ( ddiff < TINY_BIT ) then
				d_mesh = d_boundary
			else if ( d_boundary < 5E-5 .and. ddiff < 1E-8 ) then
				d_mesh = d_boundary
			end if
			
			distance = min(d_boundary, d_mesh) 
			
			!> Advance particle
			do j = 1, p % n_coord
				p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
			enddo
			
			dist_accum = dist_accum + distance
			
			call cross_surface(p, surface_crossed)
			
			if (abs(dist_accum-dist) < 1.0d-5)  then 
				vrc_idx = vrc_idx + 1
				vrc_thread(vrc_idx)%wgt = p%wgt
				vrc_thread(vrc_idx)%xyz = p%coord(1)%xyz
				vrc_thread(vrc_idx)%uvw = p%coord(1)%uvw
				vrc_thread(vrc_idx)%E   = p%E
				vrc_thread(vrc_idx)%G   = p%material 
				vrc_thread(vrc_idx)%delayed = .True.
				exit
				
				
			elseif (dist_accum > dist ) then 
				print *, "psudo ray passed the length"
				print *, dist_accum, dist, dist_accum - dist
				print *, p%coord(1)%xyz
				stop 
			endif 
			
			
			! Russian Roulette for psudo-ray
			if (p%wgt < 1.0d-8) THEN !call Russian_Roulette(p)
				wgt_s = 2*1.0d-8
				if ((p%wgt/wgt_s).ge.rang()) then
					p%wgt = wgt_s
				else
					p%alive = .false.
				endif 
			endif 
			
		enddo 
		
	end subroutine		
	
	!==============================================================
	! CREATE_RAY_DYNAMIC
	!==============================================================
	subroutine create_ray_dynamic (p, p_psudo) 
        type(Particle), intent(in) :: p
        type(Particle) :: p_psudo
		integer :: i, n_mat
		integer :: xn, ierg
		real(8) :: rn, el, noel, r, sigt_sum, temp, sum1, sum2
		real(8) :: micro_xs(6), macro_xs(5)
		real(8) :: ipfac
		integer :: n_iso
		real(8) :: dtemp
		integer :: iso
		real(8) :: val, beta, pdf, nu_del 
		integer :: pt1, pt2, pt3
		integer :: ng, NE
		
		
		call p_psudo%initialize()
		p_psudo = p 
		
		if (E_mode == 0) then 
			
		else 
			! Sample a target isotope in the mixture
			call WHAT_TEMPERATURE(p_psudo)
			macro_xs = getMacroXS(materials(p_psudo%material), p_psudo%E,p_psudo%kT, p_psudo%urn)
			rn = rang(); temp = 0
			do i = 1, materials(p_psudo%material)%n_iso
				dtemp = abs(p_psudo%kT-ace(materials(p_psudo%material)%ace_idx(i))%temp)
				if ( materials(p_psudo%material)%db .and. dtemp > K_B .and. p_psudo%E < 1D0 ) then
				! On-the-fly Doppler broadening
				call GET_OTF_DB_MIC(p_psudo%kT,materials(p_psudo%material)%ace_idx(i),p_psudo%E,micro_xs)
				!if (  micro_xs(1) > 1E+30 ) stop
				else
				! point-wise data at the given temperature
				micro_xs = getMicroXS( materials(p_psudo%material)%ace_idx(i), p_psudo%E)
				end if
				! S(a,b)
				call GET_SAB_MIC(materials(p_psudo%material),i,p_psudo%E,micro_xs)
				temp = temp + micro_xs(1)*materials(p_psudo%material)%numden(i)*barn
				if ( rn < temp/macro_xs(1) ) then
					iso = materials(p_psudo%material)%ace_idx(i)
					if ( materials(p_psudo%material)%sab .and. ace(iso)%sab_iso /= 0 &
						.and. p_psudo%E < 4D-6 ) then
						p_psudo%yes_sab = .true.
					else
						p_psudo%yes_sab = .false.
					endif
					exit
				endif
			enddo		
		
			xn = 1
			micro_xs = getMicroXS( iso, p_psudo%E)
			el   = micro_xs(2)
			noel = 0 
			call getierg(iso,ierg,p_psudo%E)
			ipfac = max(0.d0, min(1.d0,(p_psudo%E-ace(iso)%E(ierg))/(ace(iso)%E(ierg+1)-ace(iso)%E(ierg))))
			do i = 1, ace(iso)%NXS(5) !> through the reaction types...
				if (abs(ace(iso)%TY(i)) == 19) cycle 
				noel = noel + ace(iso)%sig_MT(i)%cx(ierg) & 
							+ ipfac*(ace(iso)%sig_MT(i)%cx(ierg+1) - ace(iso)%sig_MT(i)%cx(ierg))
			enddo 
			
			
			! =========== steady state =================
			!r = rang()*(noel+el)-el
			!if( ace(iso)%nxs(5) == 0 .or. r <= 0.0d0 ) then 
			!	call elastic_CE (p_psudo, iso)
			!else
			!	call notElastic_CE (p_psudo, iso, xn)
			!end if
			!
			!p_psudo%wgt = p_psudo%wgt * ((el+noel)/micro_xs(1)) !(1 - micro_xs(3)/micro_xs(1))
			! ===============================================
			
			
			
			
			
			
			! Calculate beta
			nu_del = getnudel(iso,p_psudo%E)
			ng = ace(iso)%NXS(8)
			beta = 0 
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
				beta = beta + pdf * nu_del
			enddo 		
			
			
			val = (1.0-beta)*micro_xs(5) + el + noel 
			r = rang()
			if( ace(iso)%nxs(5) == 0 .or. r < (el)/val ) then ! elastic scattering
				if ( p_psudo%yes_sab ) then
					call SAB_CE(p_psudo,iso,micro_xs(2),micro_xs(6))
				else
					call elastic_CE (p_psudo, iso)
				end if
				
			elseif (r < (el+noel)/val) then ! inelastic scattering
				call inElastic_CE (p_psudo,iso,xn)
				p_psudo%wgt = p_psudo%wgt * dble(xn)
				
			else ! prompt fission 
				p_psudo%E = fission_E (p_psudo%E, iso,.false.)
				do i = 1,p_psudo%n_coord
					p_psudo%coord(i)%uvw(:) = rand_vec()
				enddo
				!p%n_coord = 1
			endif
			p_psudo%wgt = p_psudo%wgt * (val/micro_xs(1))
			
		endif
		p_psudo%n_coord = 1
		
	end subroutine
	
	
	!==============================================================
	! RESET_RAY clears data of a ray
	!==============================================================
    elemental subroutine reset_ray(this)
        class(RayType), intent(inout) :: this
		
        this % xyz = 0 
        this % uvw = 0 
        this % E   = 0 
        this % G   = 0 
		if (allocated(this%sigt)) deallocate(this%sigt)
		if (allocated(this%nuSigf)) deallocate(this%nuSigf)
		
    end subroutine reset_ray
	
end module 
