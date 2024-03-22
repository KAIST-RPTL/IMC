module tracking
    use omp_lib
    use MPI
    use variables
    use constants
    use surface_header
    use geometry_header,    only: cell, lattices, cells, cell_distance, universes
    use geometry,           only: cross_surface, distance_to_boundary, find_cell
    use particle_header,    only: particle
    use randoms,            only: rang
    use physics
    use XS_header 
    use tally,              only: TallyCoord, TallyFlux, FindTallyBin,&
                                  TallyPower, tallyon, MESH_DISTANCE, meshon, &
                                  MC_TRK, TALLY_SURF, MC_TRK_S, meshon_tet_vrc, mesh_power
    use ace_xs,             only: getMacroXS, getMacroXS_UEG
    use material_header,    only: materials, n_materials
    use ace_reactions,      only: collision_CE, WHAT_TEMPERATURE, WHAT_DENSITY
    use FMFD,               only: FMFD_TRK, FMFD_COL, FMFD_SURF, fmfdon, &
                                  DTMCBU, DTMC_BU_COL, DTMC_BU_TRK
    use FMFD_header,        only: acc_skip
    use DEPLETION_MODULE,   only: tally_burnup
    use VRC,                only: m_pseudo, trace_psudoray, trace_psudoray_tet_vrc, create_ray_dynamic
    use transient, 			only: del_t, curr_time, curr_timestep, n_col, n_cross
	use DMC, 				only: collision_dynamic_MG, collision_dynamic_CE
	use PCQS, 				only: collision_pcqs_MG, collision_pcqs_MG_init, collision_pcqs_CE, collision_pcqs_CE_init
	use bank_header
	use tetrahedral, 		only: tet, distance_tet, distance_tet_from_outside, in_the_tet, find_tet, node, find_tet_old, RayIntersectsTriangle
	use TH_HEADER,          only: th_on
    use TEMPERATURE,        only: TH_INSIDE, TH_COL
    use ace_header,         only: uresiso, n_unr
    implicit none

contains

!===============================================================================
! TRANSPORT - the main logic for moving a particle through geometry.
!===============================================================================

subroutine transport(p)
    implicit none

    type(Particle), intent(inout) :: p
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_mesh                 ! distance to FMFD grid
    real(8) :: d_gmsh                 ! distance to Gmsh grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)!, tmparr(3)
    real(8) :: xyz(3)
    integer :: i_cell, i_bin(4), i_lat, i_surf
    integer :: i_xyz(3), idx_xyz, j_xyz(3)
    logical :: inside_mesh
    integer :: income_mesh
    logical :: inside_th
    real(8) :: ddiff
    logical :: fm_crossed
	real(8) :: d_s, val 
	integer :: idx_surf
    integer :: next_tet
	integer :: bc
	integer :: tet_prev, tet_face
    integer :: pt1, pt2, pt3
    real(8) :: tmp_power
    real(8) :: speedn

	
	xyz = p%coord(1)%xyz
	
	
	
    found_cell = .false.
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    !print *, 'PTCL', p % kT/K_B, p % coord(1) % xyz(3)
	if (p%material < 1 .or. p%material > n_materials) then 
        ! print *, 'WTF happen? mat:', p%material
        p%alive = .false.
        return
	endif 
	
	if (p%in_tet) then
		if(.not. in_the_tet(p%coord(p%n_coord)%xyz, p%tet)) then 
			p%tet = find_tet(p%coord(p%n_coord)%xyz)
			p%tet_face = 0
			p%in_tet = .true.
			if (p%tet .le. 0) p%in_tet = .false.
		endif 
	endif 
	
	!> Tetrahedron mesh temperature treatment
	if (p%in_tet) then
		p%sqrtkT = sqrt(K_B * Tet(p%tet)%temperature ) 
		p%kT = K_B * Tet(p%tet)%temperature 
	elseif (E_mode == 1) then 
        call WHAT_TEMPERATURE(p)
        !call WHAT_DENSITY(p)
	endif

	
	
    !> Surface distance(boundary)
    call distance_to_boundary(p, d_boundary, surface_crossed)
    !> Sample distances from special boundaries in univ 0
	
	
	val = 1.0d0
    !> Sample a distance to collision
    if (E_mode == 0) then 
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) &
                    + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
        macro_xs(5) = XS_MG(p%material)%sig_fis(p%g)
        d_collision = -log(rang())/macro_xs(1)
    elseif (E_mode == 1) then
        if(do_ueg) then
            macro_xs = getMacroXS_UEG(materials(p%material), p%E,p%kT,p%urn)
        else
            macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
        endif
        macro_xs = macro_xs * p % dens
	    d_collision = -log(rang())/macro_xs(1)
        speedn = sqrt(2.d0*p%E*mevj/(m_u*m_n))*1.0d2 !> meter/sec
    endif 
    
    ! ===================================================
    !> CMFD distance 
    d_mesh = INFINITY
    if ( meshon ) &
        call MESH_DISTANCE(p,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)

    ! =========================================================================
    !> TH distance
!    d_TH = INFINITY
!    if ( th_on ) call TH_DISTANCE(p,j_xyz,d_TH)
    
    !> minimum distance
    ddiff = abs(d_boundary-d_mesh)/d_boundary
    if ( ddiff < TINY_BIT ) then
        d_mesh = d_boundary
        fm_crossed = .true.
    else if ( d_boundary < 5E-5 .and. ddiff < 1E-8 ) then
        d_mesh = d_boundary
        fm_crossed = .true.
    else if ( d_mesh < 1E-8 .and. d_boundary < d_mesh ) then
        d_mesh = d_boundary
        fm_crossed = .true.
    else
        fm_crossed = .false.
    end if
	
    ! ==================================================
	!> Tetrahedron mesh distance 
	d_gmsh = INFINITY
	if (do_gmsh .and. curr_cyc > n_inact) then 
		i_bin = FindTallyBin(p)
		if ( i_bin(1) > 0 ) then
			if (p%in_tet) then 
				call distance_tet (p, d_gmsh, next_tet, bc)
			else 
				call distance_tet_from_outside (p, d_gmsh, next_tet, tet_face) 
			endif
		endif
	endif
    distance = min(d_boundary, d_collision, d_mesh, d_gmsh)
    ! print *, 'TST', icore, p%n_cross,p%coord(1)%xyz, distance, d_boundary, d_collision, d_mesh
    p % trvltime = p % trvltime + distance / speedn 
    if(distance>TOOLONG) then
        print *, 'ESCAPED',distance,p%coord(1)%xyz(1:2),p%coord(1)%uvw(1:2)
        p%alive = .false.
    endif
	
    !> Track-length estimator
    !       >> Removed 'ATOMIC' due to Reduction: needs to be verified
    !           >> NO... While loop doesn't work for Reduction
    !$OMP ATOMIC 
    k_tl = k_tl + distance*p%wgt*macro_xs(4) 

    !> Flux Tally ===========================================================================
    if ( tally_switch > 0 .and. do_transient == .false.) then 
        i_bin = FindTallyBin(p)
		if (do_gmsh .and. i_bin(1) > 0) then 
			if (p%in_tet) then 
				i_bin(1) = p%tet
			else
				i_bin(1) = -1
			endif 
		endif 
		
		
        if ( i_bin(1) > 0 ) then 
            !$omp atomic
            TallyFlux(i_bin(1)) = TallyFlux(i_bin(1)) + distance*p%wgt
    !> ==== Power Tally =====================================================================
            !$omp atomic
            TallyPower(i_bin(1)) = TallyPower(i_bin(1)) + distance*p%wgt*macro_xs(5)
			
        endif
    endif 
    !> Cycle-power Tally ===================================================================  
    if(curr_cyc > n_inact .and. materials(p%material)%fissionable) then 
        !$omp atomic
        cyc_power = cyc_power + distance * p % wgt * macro_xs(5)

        if(do_mgtally) then
            ! 1. Find MG
            if(p%E > Ebin(1) .or. p%E < Ebin(n_mg+1)) then            
            else
            pt1 = 1; pt2 = n_mg+1
            do
            if(pt2-pt1 == 1 ) exit
            pt3 = (pt1 + pt2) / 2
            if(p % E < Ebin(pt3)) then
                pt1 = pt3
            else
                pt2 = pt3
            endif
            enddo
            !$omp atomic
            micro_flux(pt1) = micro_flux(pt1) + distance * p % wgt
            endif
        endif
        if(do_ifp)then
            ! ADJOINT related
            ! Effective beta calc.
            !$omp atomic
            denom= denom + distance*p%wgt*macro_xs(4)
            !$omp atomic
            gen_numer = gen_numer + distance*p%wgt*p%nlifearr(1)*macro_xs(4)
            if(p%delayedarr(1)>0) then
                !$omp atomic
                beta_numer(p%delayedarr(1)) = beta_numer(p%delayedarr(1)) + distance*p%wgt*macro_xs(4)
                !$omp atomic
                lam_denom(p%delayedarr(1))  = lam_denom(p%delayedarr(1)) + distance*p%wgt*macro_xs(4)/p%delayedlam(1)
            endif
        endif
    endif 

    !> MC Tally
    if ( tallyon .and. .not. fmfdon .and. inside_mesh ) &
        call MC_TRK(p%E,p%wgt,distance,macro_xs,i_xyz)
    
    !> Burn-up Tally ========================================================================
    call tally_burnup (p%material, distance, p%wgt, p%E)


    
    !> FMFD Tally (track length) 
    if ( fmfdon .and. inside_mesh ) then
        call FMFD_TRK(p%wgt,distance,macro_xs,i_xyz)
        if(DTMCBU .and. curr_cyc > acc_skip) then    
            i_cell = p%coord(p%n_coord)%cell
            call DTMC_BU_TRK(i_xyz, cells(i_cell)%dtmc, p%wgt*distance, macro_xs(5))
        endif
    endif
	
	
    !> Advance particle
    do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
    enddo
	
    if ( distance == d_collision ) then ! collision 
		if (do_PCQS .and. curr_cyc > n_inact) n_col = n_col + 1 
        if ( fmfdon .and. inside_mesh ) then
        call FMFD_COL(p%wgt,macro_xs,i_xyz)
        if ( DTMCBU .and. curr_cyc > acc_skip ) &
        call DTMC_BU_COL(i_xyz,cells(i_cell)%dtmc,p%wgt/macro_xs(1),macro_xs(5))
        end if
		p%tet_face = 0 
        if (E_mode == 0) then 
            call collision_MG(p)
        else !(E_mode == 1) 
            call collision_CE(p)
        endif
        if ( th_on .and. .not. fmfdon ) then
            call TH_INSIDE(p%coord(1)%xyz(:),j_xyz(:),inside_th)
            if ( inside_th ) call TH_COL(p%wgt,macro_xs(1),macro_xs(4),j_xyz(:))
        end if

    elseif  ( distance == d_mesh ) then 
        call FMFD_SURF(inside_mesh, income_mesh,i_surf, i_xyz, &
                        p%coord(1)%uvw, p%wgt, surfaces(surface_crossed)%bc)
!         if(surfaces(surface_crossed) % bc == 2) then
!             if((abs(p%coord(1)%xyz(1)-surfaces(surface_crossed)%parmtrs(3))>TINY_BIT) .and. (abs(p%coord(1)%xyz(2)-surfaces(surface_crossed)%parmtrs(3))>TINY_BIT) ) then
!                 print *, 'WTF22?', p % coord(1) % xyz(:)
!             endif
!         endif
!         call cross_surface(p, surface_crossed)
        if ( fm_crossed ) then
            call cross_surface(p, surface_crossed)
        else
            p%coord(1)%xyz = p%coord(1)%xyz + TINY_BIT * p%coord(1)%uvw
        end if
		
    elseif (abs(distance - d_boundary) < TINY_BIT) then
        p%n_cross = p%n_cross + 1 
!         if(surfaces(surface_crossed) % bc == 2) then
!             if((abs(p%coord(1)%xyz(1)-surfaces(surface_crossed)%parmtrs(3))>TINY_BIT) .and. (abs(p%coord(1)%xyz(2)-surfaces(surface_crossed)%parmtrs(3))>TINY_BIT) ) then
!                 print *, 'WTF?', p % coord(1) % xyz(:)
!             endif
!         endif
! 

        if (surface_crossed > 0) call cross_surface(p, surface_crossed)
        if(p%alive == .false.) then 
            !$omp atomic
            loss_vrc = loss_vrc + p%wgt! * exp(-macro_xs(1)*d_boundary)
        endif 
		
	elseif (abs(distance - d_gmsh) < TINY_BIT) then 
		tet_prev = p%tet
		p%tet_prev = p%tet
		if (.not. p%in_tet) then   ! outside -> tet 
			p%in_tet = .true. 
			p%tet = next_tet
			p%tet_face = tet_face
		elseif (next_tet .le. 0) then  ! tet -> outside 
			p%in_tet = .false.
			p%tet = -1
			p%tet_face = -1
			
		else                       ! tet1 -> tet2
			p%tet = next_tet
			do i = 1, 4
				if (tet(next_tet)%neighbor(i) == tet_prev) then 
					p%tet_face = i
					exit
				endif
			enddo
		endif 
		do j = 1, p % n_coord
			p % coord(j) % xyz = p % coord(j) % xyz + 1.0d-15 * p % coord(j) % uvw
		enddo
    endif

end subroutine transport

!===============================================================================
! transport_DT handles the Woodcock Delta-tracking algorithm 
!===============================================================================
subroutine transport_DT(p) 

    type(Particle), intent(inout) :: p
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_mesh                  ! distance to CMFD grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    
    real(8) :: macro_major              ! the global majorant cross-section
    real(8), allocatable :: macro_tot(:)
    real(8) :: macro_xs(5)
    real(8) :: xyz(3)
    integer :: i_cell, idx_surf
    integer :: i_xyz(3), idx_xyz, idx
    integer :: bc
    integer :: n_mat 
    
    
    found_cell = .false.
    i_xyz(:) = -1; idx_xyz = -1
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
    !> Determine macro_major 
    if (E_mode == 0) then 
        n_mat = size( XS_MG ) 
        allocate(macro_tot(n_mat))
        do i = 1, n_mat 
            macro_tot(i) = (sum(XS_MG(i)%sig_scat(p%g,:)) + XS_MG(i)%sig_abs(p%g))
        enddo 
    else 
        n_mat = size( materials ) 
        allocate(macro_tot(n_mat))
        do i = 1, n_mat 
            macro_xs = getMacroXS(materials(i), p%E,p%kT,p%urn)
            macro_tot(i) = macro_xs(1) 
        enddo 
    endif 
    macro_major = maxval(macro_tot, n_mat)
    deallocate(macro_tot)
    
    !> Sample a distance to collision
    d_collision = -log(rang())/macro_major
    
    !> Sample distances from special boundaries in univ 0
    idx = p % coord(1) % cell
    call cell_distance(cells(idx), p%coord(1)%xyz, p%coord(1)%uvw, surfaces, d_boundary, idx_surf)
    
    distance = min(d_boundary, d_collision)
    !print *, d_boundary, d_collision
    
    !> Determine Virtual / Real collsion OR Cross-surface
    if (d_collision < d_boundary) then ! collision 
        ! Advance particle
        p % n_coord = 1
        p % coord(1) % xyz = p % coord(1) % xyz + distance * p % coord(1) % uvw
        call find_cell(p, found_cell, i_cell)
        
        if (E_mode == 0) then 
            macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        elseif (E_mode == 1) then 
            macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
        endif 
        
        ! Reject? 
        if (rang() < macro_xs(1)/macro_major ) then !> real collision
            if (E_mode == 0) then 
                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
                call collision_MG_DT(p, macro_major)
            else 
                !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
                call collision_CE(p)
            endif
        endif 
    else
        call cross_surface(p, idx_surf)
    endif
    
    p%n_coord = 1 
    
end subroutine


!===============================================================================
! transport_VRC handles the Volumetric-Ray-Casting Method 
! in combination with the Delta-tracking 
!===============================================================================
subroutine transport_VRC(p)

    type(Particle), intent(inout) :: p
    type(Particle) :: p_psudo
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: d_CMFD                  ! distance to CMFD grid
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3)
    integer :: i_cell, i_bin, i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    logical :: inside_CMFD
    integer :: bc
    real(8) :: val
    
    
    found_cell = .false.
    i_xyz(:) = -1; idx_xyz = -1
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
    call distance_to_boundary(p, d_boundary, surface_crossed)
    
    ! Sample a distance to collision
    if (E_mode == 0) then 
        d_collision = -log(rang())/(sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
        
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E,p%kT,p%urn)
        d_collision = -log(rang())/macro_xs(1)
    endif 
            
    distance = min(d_boundary, d_collision)
    !> Track-length estimator
    !!$omp atomic 
    !k_tl = k_tl + distance*p%wgt*macro_xs(4)

    
    
    if (p % vrc_traced == .false.) then 
        p_psudo = p
        call trace_psudoray(p_psudo)
        call p_psudo%clear()
        p % vrc_traced = .true.
    endif
    
    
    !> Advance particle
    do j = 1, p % n_coord
        p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
    enddo
            
    if (distance == d_collision) then ! collision     
        if (E_mode == 0) then 
            !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
            call collision_MG(p)
        else !(E_mode == 1) 
            !call CMFD_tally_col(p%wgt, macro_xs, idx_xyz, inside_CMFD)
            call collision_CE(p)
        endif
        p_psudo = p
        call trace_psudoray(p_psudo)
        !call p_psudo%clear()
        
    else
        call cross_surface(p, surface_crossed)
        if(p%alive == .false.) then 
            !$omp atomic
            loss_vrc = loss_vrc + p%wgt! * exp(-macro_xs(1)*d_boundary)
        endif 
    endif
    
end subroutine transport_VRC



!===============================================================================
! TRANSPORT_DYNAMIC - moving a particle through geometry in DMC algorithm
!===============================================================================

subroutine transport_dynamic(p)
	
    type(Particle), intent(inout) :: p
    type(Particle) :: p_psudo
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3), uvw(3), E_prev, wgt
    integer :: i_cell, i_bin(4), i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    real(8) :: speedn, tfly
	integer :: isize, nsplit
	logical :: store
	integer :: idx_surf 
	real(8) :: d_s, val 
    real(8) :: d_gmsh                 ! distance to Gmsh grid
    integer :: next_tet
	integer :: bc
	integer :: tet_prev, tet_face
	real(8) :: d_mesh 
    logical :: inside_mesh
    integer :: income_mesh
	real(8) :: ddiff 
	logical :: fm_crossed 
	
	
	
	
	store = .false.
	
    found_cell = .false.
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
	
	
	!> Distance tet 
	if (p%in_tet) then
		if(.not. in_the_tet(p%coord(p%n_coord)%xyz, p%tet)) then 
			p%tet = find_tet(p%coord(p%n_coord)%xyz)
			p%tet_face = 0
			p%in_tet = .true.
			if (p%tet .le. 0) p%in_tet = .false.
		endif 
	endif 
	
	!> Tetrahedron mesh temperature treatment
	if (p%in_tet) then
		p%sqrtkT = sqrt(K_B * Tet(p%tet)%temperature ) 
		p%kT = K_B * Tet(p%tet)%temperature 
	else
		p%kT = materials(p%material)%temp
	endif
	
	
	
    !> Surface distance(boundary)
    call distance_to_boundary(p, d_boundary, surface_crossed)
    	
	
    !> Sample a distance to collision
	val = 1.0d0
    if (E_mode == 0) then 
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) &
                    + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
		speedn = MGD(p%material)%vel(p%g)*1.0d-2
		
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E, p%kT,p%urn)
		speedn = sqrt(2.0d0*p%E*mevj/(m_u*m_n))   ! m/s
    endif
	
	d_collision = -log(rang())/macro_xs(1)


	
	
	!> Tetrahedron mesh distance 
	d_gmsh = INFINITY
	if (do_gmsh .and. curr_cyc > n_inact) then 
		i_bin = FindTallyBin(p)
		if ( i_bin(1) > 0 ) then
			if (p%in_tet) then 
				call distance_tet (p, d_gmsh, next_tet, bc)
			else 
				call distance_tet_from_outside (p, d_gmsh, next_tet, tet_face) 
			endif
		endif
	endif
	
	
    !> minimum distance
    distance = min(d_boundary, d_collision, d_gmsh)
	tfly = (distance*1.0d-2) / speedn
	
	
	
	if ((p%time+tfly) .ge. (curr_time+del_t)) then    ! particle stop at time boundary
		tfly = curr_time+del_t-p%time
		distance = tfly*speedn*1.0d2
		store = .true. 
	endif 
	p%time = p%time + tfly
    
	
	
    !> Flux Tally ===========================================================================
    if ( tally_switch > 0 ) then 
        i_bin = FindTallyBin(p)
		if (do_gmsh .and. i_bin(1) > 0) then 
			if (p%in_tet) then 
				i_bin(1) = p%tet
			else
				i_bin(1) = -1
			endif 
		endif 
		
        if ( i_bin(1) > 0 ) then 
            !$omp atomic
            TallyFlux(i_bin(1)) = TallyFlux(i_bin(1)) + distance*p%wgt
		!> Power Tally =====================================================================
            !$omp atomic
            TallyPower(i_bin(1)) = TallyPower(i_bin(1)) + distance*p%wgt*macro_xs(5)
			
        endif
    endif 	
	
	
	

	!> Cycle-power Tally ====================================================================
	if(curr_cyc > n_inact) then
		if (E_mode == 1) then ! CE
			!$omp atomic
			cyc_power = cyc_power + distance*p%wgt*macro_xs(5)
			
			if (inside_mesh) then 
				!$omp atomic
				mesh_power = mesh_power + distance*p%wgt*macro_xs(5)
			endif 
		else 
			!$omp atomic
			cyc_power = cyc_power + distance*p%wgt*macro_xs(4)
			
			if (inside_mesh) then 
				!$omp atomic
				mesh_power = mesh_power + distance*p%wgt*macro_xs(4)
			endif 
		endif 
	endif
	
    !> Advance particle
	do j = 1, p % n_coord
		p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
	enddo
	
	
	
	!if (p%n_cross > 1e5) then 
	!	print *, p%n_cross, icore
	!	print *, distance, p%time, p%wgt 
	!	print *, p%coord(1)%xyz
	!	print *, p%coord(1)%uvw 
	!	print *, materials(p%material)%mat_name 
	!	print *, '' 
	!endif 
	
	
	
	!> TODO 
	
    if ( distance == d_collision ) then ! collision
	
		!$omp atomic
		n_col = n_col + 1
	
		!> Main particle
		if (E_mode == 0) then 
			call collision_dynamic_MG(p) 
        else !(E_mode == 1) 
			call collision_dynamic_CE(p)
        endif
		
		
    elseif (distance == d_boundary) then 
		if (surface_crossed > 0)  call cross_surface(p, surface_crossed)
		!write(prt_dynamic,*) "cross surface "
		if (.not. p%alive) then 
			!$omp atomic 
			DMC_loss = DMC_loss + p%wgt
		endif 
		p%n_cross = p%n_cross + 1 
		
		!$omp atomic 
		n_cross = n_cross + 1
		
		
		
	elseif (abs(distance - d_gmsh) < TINY_BIT) then 
		p%n_cross = p%n_cross + 1 
		
		tet_prev = p%tet
		p%tet_prev = p%tet
		if (.not. p%in_tet) then   ! outside -> tet 
			p%in_tet = .true. 
			p%tet = next_tet
			p%tet_face = tet_face

		elseif (next_tet .le. 0) then  ! tet -> outside 
			p%in_tet = .false.
			p%tet = -1
			p%tet_face = -1
			
		else                       ! tet1 -> tet2
			p%tet = next_tet
			do i = 1, 4
				if (tet(next_tet)%neighbor(i) == tet_prev) then 
					p%tet_face = i
					exit
				endif
			enddo
		endif 
		
		do j = 1, p % n_coord
			p % coord(j) % xyz = p % coord(j) % xyz + 1.0d-15 * p % coord(j) % uvw
		enddo

		
	elseif (store) then  ! stop the particle and store in the bank for the next time step 
		! Split time source 
		nsplit = int(p%wgt/wgt_split_dyn) + 1 
		do i = 1, nsplit
			bank_idx = bank_idx + 1
			thread_bank(bank_idx)%wgt 		= p%wgt / real(nsplit,8)
			thread_bank(bank_idx)%xyz 		= p%coord(1)%xyz
			thread_bank(bank_idx)%uvw 		= p%coord(1)%uvw
			thread_bank(bank_idx)%E			= p%E
			thread_bank(bank_idx)%G			= p%G
			thread_bank(bank_idx)%delayed 	= .false.
			thread_bank(bank_idx)%time 		= p%time
		enddo
		p%alive = .false.
	endif

	
	if (p%wgt > wgt_split_dyn .and. p%alive) then 
		! Split the particle
		nsplit = int(p%wgt/wgt_split_dyn) + 1 
		do i = 1, nsplit-1
			split_idx = split_idx + 1
			split_thread(split_idx)%wgt 	= p%wgt / real(nsplit,8)
			split_thread(split_idx)%xyz 	= p%coord(1)%xyz
			split_thread(split_idx)%uvw 	= p%coord(1)%uvw
			split_thread(split_idx)%E		= p%E
			split_thread(split_idx)%G		= p%G
			split_thread(split_idx)%delayed = .false.
			split_thread(split_idx)%time 	= p%time
			
			
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
		enddo
		p%wgt = p%wgt / real(nsplit,8)
		
	endif
	
	
	
	!if buffer is almost full -> add to the fission bank
	if ( bank_idx > 7500 ) then 
	!$omp critical
		isize = size(dynamic_bank)
		if(allocated(temp_bank)) deallocate(temp_bank)
		allocate(temp_bank(1:isize+bank_idx)) 
		if (isize>0) temp_bank(1:isize) = dynamic_bank(:)
		deallocate(dynamic_bank)
		temp_bank(isize+1:isize+bank_idx) = thread_bank(1:bank_idx)
		call move_alloc(temp_bank, dynamic_bank)
	!$omp end critical
		bank_idx = 0
	endif
	
end subroutine transport_dynamic








!===============================================================================
! TRANSPORT_PCQS - moving a particle through for PCQS MC 
!===============================================================================

subroutine transport_pcqs(p)
	use PCQS, only: PCQS_prod, PCQS_abs, PCQS_leak, PKE_gamma, n_pcqs_inact, PKE_tally 
	use transient, only: npg 
	
    type(Particle), intent(inout) :: p
    type(Particle) :: p_psudo
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3), uvw(3), E_prev, wgt
    integer :: i_cell, i_bin(4), i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    real(8) :: speedn, tfly
	integer :: isize, nsplit
	logical :: store
	integer :: idx_surf 
	real(8) :: d_s, val 
    real(8) :: d_gmsh                 ! distance to Gmsh grid
    integer :: next_tet
	integer :: bc
	integer :: tet_prev, tet_face
	real(8) :: d_mesh 
    logical :: inside_mesh
    integer :: income_mesh
	real(8) :: ddiff 
	logical :: fm_crossed 
	real(8) :: sigt_pcqs
	real(8) :: beta 
	
    found_cell = .false.
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
	
	if (p%material == 0 ) then 
		print *, 'killed', p%coord(1)%xyz 
		p%alive = .false. 
		return 
	endif 
	
	
    !> Surface distance(boundary)
    call distance_to_boundary(p, d_boundary, surface_crossed)
    
	
    !> Sample a distance to collision
	val = 1.0d0
    if (E_mode == 0) then 
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) &
                    + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
		speedn = MGD(p%material)%vel(p%g)
		
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E, p%kT,p%urn)
		speedn = sqrt(2.0d0*p%E*mevj/(m_u*m_n))*1.0d2   ! cm/s
    endif
	
	!> Defind sigt_pcqs 
	sigt_pcqs = macro_xs(1) + 1.0d0 / (speedn*del_t) + (PKE_gamma / speedn)
	
	d_collision = -log(rang())/sigt_pcqs

	
    !> Mesh distance 
    d_mesh = INFINITY
	inside_mesh = .false. 
    if ( meshon_tet_vrc ) &
        call MESH_DISTANCE(p,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)
    
	
    !> minimum distance
    ddiff = abs(d_boundary-d_mesh)/d_boundary
    if ( ddiff < TINY_BIT ) then
        d_mesh = d_boundary
        fm_crossed = .true.
    else if ( d_boundary < 5E-5 .and. ddiff < 1E-8 ) then
        d_mesh = d_boundary
        fm_crossed = .true.
    else
        fm_crossed = .false.
    end if
	
	
    !> minimum distance
    distance = min(d_boundary, d_collision, d_mesh)

	
	!$omp atomic
	PCQS_prod = PCQS_prod + p%wgt*macro_xs(4) * distance
	!$omp atomic
	PCQS_abs = PCQS_abs + p%wgt*macro_xs(2) * distance
	
	
	wgt = p%wgt * macro_xs(1) / sigt_pcqs
	beta = sum(MGD(p%material)%beta(:))
	if (curr_cyc > n_pcqs_inact .and. beta > 0) then 
		!> PKE parameter tally 
		!$omp critical 
		do i = 1, npg 
			PKE_tally(i) = PKE_tally(i)+ distance*wgt*MGD(p%material)%beta(i)*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)
			PKE_tally(i+npg) = PKE_tally(i+npg) + distance*wgt*MGD(p%material)%beta(i)*XS_MG(p%material)%nu(p%g) &
												*XS_MG(p%material)%sig_fis(p%g) / MGD(p%material)%lambda(i)
		enddo 
		PKE_tally(npg*2+1)     = PKE_tally(npg*2+1) + distance*wgt*XS_MG(p%material)%nu(p%g)*XS_MG(p%material)%sig_fis(p%g)
		PKE_tally(npg*2+2)     = PKE_tally(npg*2+2) + distance*wgt / speedn
		!$omp end critical 
	endif
	
	
	
	
	
	
	
	
	
    !> Flux Tally ===========================================================================
    if ( tally_switch > 0 ) then 
        i_bin = FindTallyBin(p)
		if (do_gmsh .and. i_bin(1) > 0) then 
			if (p%in_tet) then 
				i_bin(1) = p%tet
			else
				i_bin(1) = -1
			endif 
		endif 
		
        if ( i_bin(1) > 0 ) then 
            !$omp atomic
            TallyFlux(i_bin(1)) = TallyFlux(i_bin(1)) + distance*p%wgt
		!> Power Tally =====================================================================
            !$omp atomic
            TallyPower(i_bin(1)) = TallyPower(i_bin(1)) + distance*p%wgt*macro_xs(5)
			
        endif
    endif 
	

	!> Cycle-power Tally ====================================================================
	!if(curr_cyc > n_pcqs_inact) then
		if (E_mode == 1) then ! CE
			!$omp atomic
			cyc_power = cyc_power + distance*p%wgt*macro_xs(5)
			
			if (inside_mesh) then 
				!$omp atomic
				mesh_power = mesh_power + distance*p%wgt*macro_xs(5)
			endif 
			
		else 
			!$omp atomic
			cyc_power = cyc_power + distance*p%wgt*macro_xs(4)
			
			if (inside_mesh) then 
				!$omp atomic
				mesh_power = mesh_power + distance*p%wgt*macro_xs(4)
			endif 
			
		endif 
		
	!endif
	
    !> Advance particle
	do j = 1, p % n_coord
		p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
	enddo
	
	
    if ( distance == d_collision ) then ! collision
	
		!!> Create Psudo-ray 
		!if (.not.inside_mesh .and. meshon_tet_vrc) then 
		!	do i = 1, m_pseudo
		!		call create_ray_dynamic (p, p_psudo) 
        !        p_psudo % wgt = p_psudo % wgt / real(m_pseudo, 8)
		!		!> Mesh distance 
		!		d_mesh = INFINITY
		!		inside_mesh = .false. 
		!		call MESH_DISTANCE(p_psudo,i_xyz,d_mesh,inside_mesh,income_mesh,i_surf)
		!		
		!		if (d_mesh < toolong) then 
		!			call trace_psudoray_tet_vrc(p_psudo, d_mesh)
		!		endif 
		!		call p_psudo%clear()
		!	enddo 
		!endif 
	    
	    
		
		!> Main particle
		xyz = p%coord(1)%xyz
		uvw = p%coord(1)%uvw
		E_prev = p%E
		wgt = p%wgt
		if (E_mode == 0) then 
			call collision_pcqs_MG(p)
        else !(E_mode == 1) 
			call collision_pcqs_CE(p)
        endif
		
		!! material is saved 
		!if (inside_mesh) then 
		!	vrc_idx = vrc_idx + 1
		!	vrc_thread(vrc_idx)%wgt = wgt
		!	vrc_thread(vrc_idx)%xyz = xyz
		!	vrc_thread(vrc_idx)%uvw = uvw
		!	vrc_thread(vrc_idx)%E   = E_prev
		!	vrc_thread(vrc_idx)%G   = p%material
		!	vrc_thread(vrc_idx)%delayed = .false. ! collision source
		!endif 
		
		!if (curr_cyc > n_pcqs_inact) n_col = n_col + 1
		!$omp atomic
		n_col = n_col + 1
		
		
		
    !elseif  ( distance == d_mesh ) then 
	!	
    !    if ( fm_crossed ) then
    !        call cross_surface(p, surface_crossed)
    !    else
	!		do j = 1, p % n_coord
	!			p%coord(j)%xyz = p%coord(j)%xyz + TINY_BIT * p%coord(j)%uvw
	!		enddo 
    !    end if
		
    elseif (distance == d_boundary) then 
		if (surface_crossed > 0)  call cross_surface(p, surface_crossed)
		p%n_cross = p%n_cross + 1 
		if (.not. p%alive) then 
			!$omp atomic
			PCQS_leak = PCQS_leak + p%wgt
		else 
			!if (curr_cyc > n_pcqs_inact) n_cross = n_cross + 1
			!$omp atomic
			n_cross = n_cross + 1
		endif 
		
	endif

	
	
end subroutine transport_pcqs


!===============================================================================
! TRANSPORT_PCQS_INIT - Initialize banks for PCQS MC 
!===============================================================================

subroutine transport_pcqs_init(p)
	
    type(Particle), intent(inout) :: p
    type(Particle) :: p_psudo
    
    integer :: i 
    integer :: j                      ! coordinate level
    integer :: next_level             ! next coordinate level to check
    integer :: surface_crossed        ! surface which particle is on
    integer :: lattice_translation(3) ! in-lattice translation vector
    integer :: n_event                ! number of collisions/crossings
    real(8) :: d_boundary             ! distance to nearest boundary
    real(8) :: d_collision            ! distance to collision
    real(8) :: distance               ! distance particle travels
    logical :: found_cell             ! found cell which particle is in?
    real(8) :: macro_xs(5)
    real(8) :: xyz(3), uvw(3), E_prev, wgt
    integer :: i_cell, i_bin(4), i_lat, i_surf
    integer :: i_xyz(3), idx_xyz
    real(8) :: speedn
	integer :: idx_surf 
	real(8) :: d_s, val 
	integer :: bc
	real(8) :: sigt_pcqs
	
	
    found_cell = .false.
    if (p%n_coord == 1) call find_cell(p, found_cell, i_cell)
    
	
    !> Surface distance
    call distance_to_boundary(p, d_boundary, surface_crossed)
    	
	
    !> Sample a distance to collision
    if (E_mode == 0) then 
        macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) &
                    + XS_MG(p%material)%sig_abs(p%g))
        macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
        macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
        macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
		speedn = MGD(p%material)%vel(p%g)
    elseif (E_mode == 1) then 
        macro_xs = getMacroXS(materials(p%material), p%E, p%kT,p%urn)
		speedn = sqrt(2.0d0*p%E*mevj/(m_u*m_n))*1.0d2   ! cm/s
    endif
	
	d_collision = -log(rang())/macro_xs(1)
	
    !> minimum distance
    distance = min(d_boundary, d_collision)
	
    !> Advance particle
	do j = 1, p % n_coord
		p % coord(j) % xyz = p % coord(j) % xyz + distance * p % coord(j) % uvw
	enddo
	
    if ( distance == d_collision ) then ! collision
		if (E_mode == 0) then 
			call collision_pcqs_MG_init(p) 
        else !(E_mode == 1) 
			call collision_pcqs_CE_init(p)
        endif
		
    elseif (distance == d_boundary) then 
		if (surface_crossed > 0)  call cross_surface(p, surface_crossed)
		p%n_cross = p%n_cross + 1 
		
	endif

	
	
end subroutine transport_pcqs_init






end module tracking
