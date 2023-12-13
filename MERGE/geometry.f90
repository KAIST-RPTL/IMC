module geometry
    use constants,             only: TINY_BIT, MAX_COORD
    use surface_header 
    use geometry_header
    use particle_header
    use omp_lib
    implicit none
    
    
    
    contains 
    
    function cell_contains(c, p) result(in_cell)
        type(Cell), intent(in) :: c
        type(Particle), intent(in) :: p
        logical :: in_cell
        integer :: i,j, n 
        
        j = p % n_coord
        if (c%operand_flag >= 0) then   !> and 
            in_cell = .true.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),p%coord(j)%xyz) == .false.) in_cell = .false.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), p%coord(j)%xyz) == .true.) in_cell = .false.
            enddo
        else     !> or
            in_cell = .false.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),p%coord(j)%xyz) == .true.) in_cell = .true.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), p%coord(j)%xyz) == .false.) in_cell = .true.
            enddo
        endif
        
    end function 
    
    function cell_xyz(c, xyz) result(in_cell)
        type(Cell), intent(in) :: c
        real(8),    intent(in) :: xyz(3)
        logical :: in_cell
        integer :: i, j, n 
        
        ! print *, 'CXYZ: ', c%cell_id, xyz(1:3)
        if (c%operand_flag >= 0) then   !> and 
            in_cell = .true.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),xyz) == .false.) in_cell = .false.
                !print *, 'NEG: ', i, surfaces(c%neg_surf_idx(i))%surf_id, in_cell
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), xyz) == .true.) in_cell = .false.
                !print *, 'POS: ', i, surfaces(c%pos_surf_idx(i))%surf_id, in_cell
            enddo
        else     !> or
            in_cell = .false.
            n = size(c%neg_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%neg_surf_idx(i)),xyz) == .true.) in_cell = .true.
            enddo
            n = size(c%pos_surf_idx)
            do i = 1, n
                if (surf_neg_or_pos(surfaces(c%pos_surf_idx(i)), xyz) == .false.) in_cell = .true.
            enddo
        endif
        
    end function 
    
    
    !===============================================================================
    ! FIND_CELL determines what cell a source particle is in within a particular
    ! universe. If the base universe is passed, the particle should be found as long
    ! as it's within the geometry
    !===============================================================================
    !recursive subroutine find_cell(p, found, search_cells, cell_idx)
    recursive subroutine find_cell(p, found, cell_idx)
        type(Particle), intent(inout) :: p
        logical,        intent(inout) :: found
        !integer,        optional      :: search_cells(:)
        integer,         optional      :: cell_idx
        integer :: i                    ! index over cells
        integer :: j, k, idx            ! coordinate level index
        integer :: offset               ! instance # of a distributed cell
        integer :: distribcell_index
        integer :: i_xyz(3)             ! indices in lattice
        integer :: n                    ! number of cells to search
        integer :: i_cell               ! index in cells array
        integer :: i_univ_cell            ! index of cell in universe cell list
        integer :: i_universe           ! index in universes array

        
        do j = p % n_coord + 1, MAX_COORD
          call p % coord(j) % reset()
        enddo
        j = p % n_coord

        i_universe = p % coord(j) % universe
        
        
        if ( j == 1 ) then
            !p%univ => universes(0)
            p%coord(1)% universe = 0 
            idx = 0
        else
            idx = p % coord(j) % universe
            !p%univ => universes(idx)
            !p%coord(p%n_coord)% universe = idx
            !> coordinate translation for the new universe
            !> call last_univ_coord_tranlation(universes,p,j)
            !p%coord(j)%xyz = p%coord(j-1)%xyz에다가 연산.
            
        endif
        
        !n = p%univ%ncell; 
        n = universes(idx)%ncell; found = .false.
        !print *, 'num of cell', n
        CELL_LOOP: do i = 1, n
            ! select cells based on whether we are searching a universe or a provided
            ! list of cells (this would be for lists of neighbor cells)
            i_cell = universes(idx)%cell(i)
            ! Move on to the next cell if the particle is not inside this cell
            !print *, i, cells(i_cell)%cell_id
            !print *, p%coord(1)%xyz
            if(cell_contains(cells(i_cell), p)) then
                ! Set cell on this level
                p % coord(j) % cell = i_cell
                found = .true.
                !print *, 'in univ ',universes(idx)%univ_id , 'found in cell:', cells(i_cell)%cell_id
                !print *, i_cell
                !print *, p%coord(j)%xyz
                exit
            endif
            
        end do CELL_LOOP
        
        
                
        if ( found ) then
            associate(c => cells(i_cell))
                CELL_TYPE: if (c % filltype == FILL_MATERIAL) then
                    ! ======================================================================
                    ! AT LOWEST UNIVERSE, TERMINATE SEARCH
                    
                    p % coord(j) % cell = i_cell      ! index in cells(:) array
                    p % material = c%mat_idx
                    
                    if ( p%material == 0 ) then 
                        p%alive = .false. 
					else 
						p%alive = .true. 
                    endif
                    
                    if ( present(cell_idx) ) cell_idx = i_cell
                    
                    
                elseif (c % filltype == FILL_UNIVERSE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
                    !print *, 'FILL_UNIVERSE : ', universes(find_univ_idx(universes, c % fill))%univ_id
                    !print *, 'cell translation', c % translation
                    ! Store lower level coordinates
                    
                    p % coord(j + 1) % xyz = p % coord(j) % xyz
                    p % coord(j + 1) % uvw = p % coord(j) % uvw
                    
                    ! Move particle to next level and set universe
                    j = j + 1
                    p % n_coord = j
                                        
                    !p % coord(j) % universe = universes%find_idx(cells(i_cell) % fill)
                    p % coord(j) % universe = find_univ_idx(universes, c % fill)
                    !print *, 'universe', c % fill
                    
                    ! Apply translation
                    if (allocated(c % translation)) then
                        p % coord(j) % xyz = p % coord(j) % xyz - c % translation
                    end if
                    
                    call find_cell(p, found, cell_idx)
                    !j = p % n_coord
                
                
                elseif (c % filltype == FILL_LATTICE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
                    !print *, 'FILL_LATTICE : ', lattices(find_lat_idx(lattices, c % fill))%lat_id 
                    
                    associate (latptr => lattices(find_lat_idx(lattices, c % fill)))
                        ! Determine lattice indices
                        !i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)
                        i_xyz = lattice_coord(latptr,p % coord(j) % xyz)  !latptr%lat_pos(p % coord(j) % xyz)
                        !print *, p % coord(j) % xyz
                        !print *, 'lattice coordinate :', i_xyz(1:2)
                        ! Store lower level coordinates
                        p % coord(j + 1) % xyz = get_local_xyz(latptr, p % coord(j) % xyz, i_xyz)
                        p % coord(j + 1) % uvw = p % coord(j) % uvw
                        
                        ! set particle lattice indices
                        p % coord(j + 1) % lattice   = find_lat_idx(lattices, c % fill)
                        !print *, 'lattice   ', c%cell_id, c % fill, p % coord(j + 1) % lattice
                        p % coord(j + 1) % lattice_x = i_xyz(1)
                        p % coord(j + 1) % lattice_y = i_xyz(2)
                        p % coord(j + 1) % lattice_z = i_xyz(3)
                        
                        !print *, 'find cell', lattices%find_idx(c % fill), i_xyz(:)
                        
                        ! Set the next lowest coordinate level.
                        p % coord(j + 1) % universe = latptr % lat(i_xyz(1),i_xyz(2),i_xyz(3))
                        !print *, 'lattices universe', universes(p % coord(j + 1) % universe)%univ_id
                        
                    end associate
                    
                    ! Move particle to next level and search for the lower cells.
                    j = j + 1
                    p % n_coord = j
                    
                    call find_cell(p, found, cell_idx)
                    !j = p % n_coord
                    
                    
                end if CELL_TYPE
            end associate
        end if        
        
        
    end subroutine 
    
    !===============================================================================
    ! TRANSFORM_COORD transforms the xyz coordinate from univ1 to univ2 
    !===============================================================================
    subroutine transform_coord ()
        
        
    end subroutine
    
    
    
    
    !===============================================================================
    ! DISTANCE_TO_BOUNDARY calculates the distance to the nearest boundary for a
    ! particle 'p' traveling in a certain direction. For a cell in a subuniverse
    ! that has a parent cell, also include the surfaces of the edge of the universe.
    !===============================================================================
    subroutine distance_to_boundary(p, dist, surface_crossed)!, level)!, lattice_translation, next_level)
        type(Particle), intent(inout) :: p
        real(8),        intent(out)   :: dist
        !type(LocalCoord), intent(out) :: level(MAX_COORD) 
        integer,        intent(out)   :: surface_crossed
        !integer,        intent(out)   :: lattice_translation(3)
        !integer,        intent(out)   :: next_level
        
        integer :: idx
        integer :: i, j, j_idx
        integer :: i_xyz(3)           ! lattice indices
        integer :: level_surf_cross   ! surface crossed on current level
        integer :: level_lat_trans(3) ! lattice translation on current level
        real(8) :: xyz_t(3)           ! local particle coordinates
        real(8) :: d_lat              ! distance to lattice boundary
        real(8) :: d_surf             ! distance to surface
        real(8) :: xyz_cross(3)       ! coordinates at projected surface crossing
        real(8) :: surf_uvw(3)        ! surface normal direction

        integer :: idx_surf
        integer :: i_xyz_out(3)
        real(8) :: dist_temp
        integer :: univ_idx

        !do j = 1, MAX_COORD
        !  call level(j) % reset()
        !enddo
        
        ! initialize distance to infinity (huge)
        dist = INFINITY
        d_lat = INFINITY
        d_surf = INFINITY
        dist_temp = INFINITY
        i_xyz_out(:) = 0 
        
        !lattice_translation(:) = [0, 0, 0]

        !next_level = 0
        surface_crossed = -1 
        
        !> surface 탐색 순서
        ! 1. universe 내부 cell boundary
        ! 2. 만약 위 거리가 toolong 보다 길면 -> 상위 universe 탐색 (lattice or universe)
        p%coord(:)%dist = INFINITY            ! reset level distance
        univ_idx = 0 
        j = 0 
        LEVEL_LOOP: do 
            j = j+1

            ! get pointer to cell on this level
            idx = p % coord(j) % cell
            
            ! =======================================================================
            ! FIND MINIMUM DISTANCE TO SURFACE IN THIS CELL
            call cell_distance(cells(idx), p%coord(j)%xyz, p%coord(j)%uvw, surfaces, p%coord(j)%dist, idx_surf)
            
            ! =======================================================================
            ! FIND MINIMUM DISTANCE TO LATTICE SURFACES        
            !if ((p%coord(j) % dist > TOOLONG).and.(p % coord(j) % lattice /= NONE)) then 
            if (p % coord(j) % lattice /= NONE) then
                i_xyz(1) = p % coord(j) % lattice_x
                i_xyz(2) = p % coord(j) % lattice_y
                i_xyz(3) = p % coord(j) % lattice_z
                call lat_distance(lattices(p % coord(j) % lattice), surfaces, p % coord(j) % xyz, &
                                    p % coord(j) % uvw, i_xyz, dist_temp, idx_surf)
                if (dist_temp < p % coord(j) % dist) p % coord(j) % dist = dist_temp
            endif 
            
            
            if (j == 1) then 
                j_idx = 1
                surface_crossed = idx_surf
            elseif ((p%coord(j)%dist < p%coord(j_idx)%dist+TINY_BIT).and.&
                    (p%coord(j)%dist > p%coord(j_idx)%dist-TINY_BIT)) then  !> similar 
                if (surfaces(j)%bc == 2)  then 
                    surface_crossed = idx_surf
                    j_idx = j
                endif 
            elseif ((p%coord(j)%dist < p%coord(j_idx)%dist)) then 
                surface_crossed = idx_surf
                j_idx = j
            endif 
            
            if (j >= p%n_coord ) exit
        enddo LEVEL_LOOP
        
        dist = p%coord(j_idx)%dist
        
        !dist = p%coord(1)%dist; level(1) = p%coord(1)
        !do i = 1, p%n_coord 
        !    if (dist > p%coord(i)%dist) then 
        !        dist = p%coord(i)%dist
        !        level(1:i) =p%coord(1:i) 
        !    enddo
        !enddo
        
        !dist = minval(p%coord(:)%dist)
        
        
    end subroutine
    
    
!===============================================================================
! CROSS_SURFACE handles all surface crossings, whether the particle leaks out of
! the geometry, is reflected, or crosses into a new lattice or cell
!===============================================================================

    subroutine cross_surface(p, surface_crossed)
        type(Particle), intent(inout) :: p
        
        integer :: surface_crossed
        real(8) :: xyz(3)     ! Saved global coordinate
        real(8) :: uvw(3)     ! Saved global coordinate
        integer :: i_surface  ! index in surfaces
        logical :: rotational ! if rotational periodic BC applied
        logical :: found      ! particle found in universe?
        class(Surface), pointer :: surf
        class(Surface), pointer :: surf2 ! periodic partner surface
        integer :: i, i_cell, i_cell_prev
        !p % n_coord = 1
        !call find_cell(p, found)
        if (surfaces(surface_crossed)%bc == 1) then     !> Vacuum BC
            p % alive = .false.
            return 
        elseif (surfaces(surface_crossed)%bc == 2) then !> Reflective BC 
            
            uvw = p%coord(1)%uvw
            call reflective_bc(p%coord(1)%uvw, p%coord(1)%xyz, surface_crossed)
            p % n_coord = 1
            p % coord(1) % xyz(:) = p % coord(1) % xyz(:) - 100*TINY_BIT * uvw
            !p%last_material = p%material
            !p % coord(1) % xyz = p % coord(1) % xyz + TINY_BIT * p % coord(1) % uvw
        else
            p % n_coord = 1
            !p % coord(1) % xyz = p % coord(1) % xyz + 0.1*TINY_BIT * p % coord(1) % uvw
            p % coord(1) % xyz = p % coord(1) % xyz + 10*TINY_BIT * p % coord(1) % uvw
        endif
        call find_cell(p, found)
    end subroutine cross_surface
    
    subroutine reflective_bc (uvw, xyz, surface_crossed)
        real(8), intent(inout) :: uvw(3)
        real(8), intent(in)       :: xyz(3)
        real(8) :: xyz_(3), r, a, tmp, d(4), mind, dd(6)
        integer :: surface_crossed
        integer :: surf_type
        integer :: flag
        
        surf_type = surfaces(surface_crossed)%surf_type
        
        select case(surf_type) 
        
        case (1) !> px
            uvw(1) = -uvw(1) 
        case (2) !> py
            uvw(2) = -uvw(2) 
        case (3) !> pz
            uvw(3) = -uvw(3) 
            
        case (6) !> sqcz  (TO BE EDITTED)
            xyz_(3)   = xyz(3) 
            xyz_(1:2) = xyz(1:2) - surfaces(surface_crossed)%parmtrs(1:2) 
            r = surfaces(surface_crossed)%parmtrs(3) 
            
            if ((xyz_(2) >= -r-1.1*TINY_BIT).and.(xyz_(2) <= -r + 1.1*TINY_BIT)) then 
                uvw(2) = -uvw(2)
            elseif ((xyz_(1) >= -r-1.1*TINY_BIT).and.(xyz_(1) <= -r + 1.1*TINY_BIT)) then 
                uvw(1) = -uvw(1)
            elseif ((xyz_(2) >= r-1.1*TINY_BIT).and.(xyz_(2) <= r + 1.1*TINY_BIT)) then 
                uvw(2) = -uvw(2) 
            elseif ((xyz_(1) >= r-1.1*TINY_BIT).and.(xyz_(1) <= r + 1.1*TINY_BIT)) then 
                uvw(1) = -uvw(1) 
            else 
                print *, 'particle is not on the surface'
                print *, 'xyz', xyz
				print *, 'Surface ID : ', surfaces(surface_crossed)%surf_id
                stop 
            end if
            
        case (9) !> cylz
            xyz_(3)   = xyz(3) 
            xyz_(1:2) = xyz(1:2) - surfaces(surface_crossed)%parmtrs(1:2) 
            
            uvw(1) = uvw(1) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2))*xyz_(1)&
                        /(surfaces(surface_crossed)%parmtrs(3))**2
            uvw(2) = uvw(2) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2))*xyz_(2)&
                        /(surfaces(surface_crossed)%parmtrs(3))**2
            uvw(3) = uvw(3) 
            
        case (10) !> sph
            xyz_(1:3) = xyz(1:3) - surfaces(surface_crossed)%parmtrs(1:3) 
            
            uvw(1) = uvw(1) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(1) &
                        /(surfaces(surface_crossed)%parmtrs(4))**2
            uvw(2) = uvw(2) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(2) &
                        /(surfaces(surface_crossed)%parmtrs(4))**2
            uvw(3) = uvw(3) - 2*(xyz_(1)*uvw(1) + xyz_(2)*uvw(2) + xyz_(3)*uvw(3))*xyz_(3) &
                        /(surfaces(surface_crossed)%parmtrs(4))**2

        case (11) !> hexxc
            ! NOTICE: xyz_ for hexagonal surface /= xyz_ for others
            tmp = sqrt(3.d0)*0.5d0
            xyz_(1) = abs(xyz(1)-surfaces(surface_crossed)%parmtrs(1))
            xyz_(2) = abs((xyz(1)-surfaces(surface_crossed)%parmtrs(1))*0.5d0 &
                      - (xyz(2)-surfaces(surface_crossed)%parmtrs(2))*tmp)
            xyz_(3) = abs((xyz(1)-surfaces(surface_crossed)%parmtrs(1))*0.5d0 &
                      + (xyz(2)-surfaces(surface_crossed)%parmtrs(2))*tmp)
            r = surfaces(surface_crossed)%parmtrs(3)

            if ((xyz_(1) >= r-200*TINY_BIT).and.(xyz_(1) <= r+200*TINY_BIT))then
                uvw(1) = -uvw(1)
            elseif((xyz_(2)>=r-200*TINY_BIT).and.(xyz_(2)<=r+200*TINY_BIT)) then
                ! REFLECTION to <tmp,0.5> -> n = <-0.5,tmp>
                a = -uvw(1)+uvw(2)*tmp*2.d0 !2*DOT PRODUCT with n
                uvw(1) = uvw(1)+0.5d0*a
                uvw(2) = uvw(2)-tmp*a
            elseif((xyz_(3)>=r-200*TINY_BIT).and.(xyz_(3)<=r+200*TINY_BIT)) then
                ! REFLECTION to <tmp,-0.5> -> n = <0.5,tmp>
                a = uvw(1)+uvw(2)*tmp*2.d0
                uvw(1) = uvw(1)-0.5d0*a
                uvw(2) = uvw(2)-tmp*a
            else ! EXCEPTION MSG
                print *, 'particle not on surface'
                print *, 'ixyz', xyz_(1)-r,xyz_(2)-r,xyz_(3)-r
                print *, 'tolerance',400*TINY_BIT
                stop
            endif
        
        case(12) !> hexyc
            ! NOTICE: xyz_ for hexagonal surface /= xyz_ for others
            tmp = sqrt(3.d0)*0.5d0
            xyz_(1) = abs((xyz(1)-surfaces(surface_crossed)%parmtrs(1))*tmp &
                      -(xyz(2)-surfaces(surface_crossed)%parmtrs(2))*0.5d0)
            xyz_(2) = abs((xyz(1)-surfaces(surface_crossed)%parmtrs(1))*tmp &
                      +(xyz(2)-surfaces(surface_crossed)%parmtrs(2))*0.5d0)
            xyz_(3) = abs(xyz(2)-surfaces(surface_crossed)%parmtrs(2))
            r = surfaces(surface_crossed)%parmtrs(3)
            !print *, 'STUCK?', xyz(1:2)
            if ((xyz_(1) >= r-200*TINY_BIT).and.(xyz_(1) <= r+200*TINY_BIT))then
                ! REFLECTION to <tmp,-0.5>
                a = uvw(1)*tmp*2.d0-uvw(2)
                uvw(1) = uvw(1)-tmp*a
                uvw(2) = uvw(2)+0.5d0*a
            elseif((xyz_(2)>=r-200*TINY_BIT).and.(xyz_(2)<=r+200*TINY_BIT)) then
                ! REFLECTION to <tmp,0.5>
                a = uvw(1)*tmp*2.d0+uvw(2)
                uvw(1) = uvw(1)-tmp*a
                uvw(2) = uvw(2)-0.5d0*a
            elseif((xyz_(3)>=r-200*TINY_BIT).and.(xyz_(3)<=r+200*TINY_BIT)) then
                uvw(2) = -uvw(2)
            else ! EXCEPTION MSG
                print *, 'particle not on surface'
                print *, 'ixyz', xyz_(1)-r,xyz_(2)-r,xyz_(3)-r
                print *, 'tolerance',400*TINY_BIT
                stop
            endif
        case(13) !> RECT
            d(1) = abs(xyz(1)-surfaces(surface_crossed)%parmtrs(1))
            d(2) = abs(xyz(1)-surfaces(surface_crossed)%parmtrs(2))
            d(3) = abs(xyz(2)-surfaces(surface_crossed)%parmtrs(3))
            d(4) = abs(xyz(2)-surfaces(surface_crossed)%parmtrs(4))
            mind = minval(d(:))
            if(d(1)==mind .or. d(2)==mind) then
                uvw(1) = -uvw(1)
            else
                uvw(2) = -uvw(2)
            endif
        case(16) !> cuboid
            dd(1) = abs(xyz(1)-surfaces(surface_crossed)%parmtrs(1))
            dd(2) = abs(xyz(1)-surfaces(surface_crossed)%parmtrs(2))
            dd(3) = abs(xyz(2)-surfaces(surface_crossed)%parmtrs(3))
            dd(4) = abs(xyz(2)-surfaces(surface_crossed)%parmtrs(4))
            dd(5) = abs(xyz(3)-surfaces(surface_crossed)%parmtrs(5))
            dd(6) = abs(xyz(3)-surfaces(surface_crossed)%parmtrs(6))
            mind = minval(dd(:))
            if(dd(1)==mind .or. dd(2)==mind) then
                uvw(1) = -uvw(1)
            elseif(dd(3)==mind .or. dd(4)==mind) then
                uvw(2) = -uvw(2)
            else
                uvw(3) = -uvw(3)
            endif
        end select 
        
        
    end subroutine reflective_bc
    
    !===============================================================================
    ! NEIGHBOR_LISTS builds a list of neighboring cells to each surface to speed up
    ! searches when a cell boundary is crossed.
    !===============================================================================
    subroutine neighbor_lists()
        
    end subroutine 
    
    
    
    recursive subroutine find_cell_xyz(xyz, idx_univ, cell_idx)
        real(8),         intent(inout) :: xyz(3)
        logical                          :: found
        integer,         optional      :: cell_idx
        integer :: i                    ! index over cells
        integer :: j, k, idx            ! coordinate level index
        integer :: offset               ! instance # of a distributed cell
        integer :: distribcell_index
        integer :: i_xyz(3)             ! indices in lattice
        integer :: n                    ! number of cells to search
        integer :: i_cell               ! index in cells array
        integer :: i_univ_cell            ! index of cell in universe cell list
        integer :: i_universe           ! index in universes array
        integer,         intent(inout) :: idx_univ
        n = universes(idx_univ)%ncell; found = .false.
        CELL_LOOP: do i = 1, n

            ! select cells based on whether we are searching a universe or a provided
            ! list of cells (this would be for lists of neighbor cells)
            i_cell = universes(idx_univ)%cell(i)
            ! Move on to the next cell if the particle is not inside this cell
            if(cell_xyz(cells(i_cell), xyz)) then
                found = .true.
                exit
            endif
            !stop
        end do CELL_LOOP
        
                
        if ( found ) then
            associate(c => cells(i_cell))
                CELL_TYPE: if (c % filltype == FILL_MATERIAL) then
                    ! ======================================================================
                    ! AT LOWEST UNIVERSE, TERMINATE SEARCH
                    !print *, 'FILL_MATERIAL : ', trim(c%cell_id), xyz(:)
                    
                    if (present(cell_idx)) then 
                        cell_idx = i_cell
                    endif
                    
                elseif (c % filltype == FILL_UNIVERSE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LOWER UNIVERSE, RECURSIVELY FIND CELL
                    !print *, 'cell translation', c % translation
                    ! Store lower level coordinates
                    !p % coord(j + 1) % xyz = p % coord(j) % xyz
                    !p % coord(j + 1) % uvw = p % coord(j) % uvw

                    ! Move particle to next level and set universe
                    idx_univ = find_univ_idx(universes, c % fill)!universes%find_idx(c % fill)

                    ! Apply translation
                    if (allocated(c % translation)) then
                        xyz = xyz - c % translation
                    end if
                    !print *, 'FILL_UNIVERSE : ', universes(idx_univ)%univ_id, xyz(:)

                    call find_cell_xyz(xyz, idx_univ, cell_idx)
            
            
                elseif (c % filltype == FILL_LATTICE) then CELL_TYPE
                    ! ======================================================================
                    ! CELL CONTAINS LATTICE, RECURSIVELY FIND CELL
                    
                    associate (latptr => lattices(find_lat_idx(lattices,c % fill)))
                        ! Determine lattice indices
                        !i_xyz = lat % get_indices(p % coord(j) % xyz + TINY_BIT * p % coord(j) % uvw)
                        i_xyz = lattice_coord(latptr, xyz)!latptr%lat_pos(xyz)
                        !print *, 'lattice coordinate :', i_xyz(1:2)
                        ! Store lower level coordinates
                        xyz = get_local_xyz(latptr , xyz, i_xyz)
                        ! Set the next lowest coordinate level.
                        !print *, latptr % lat(i_xyz(1),i_xyz(2))
                        !print *, idx_univ
                        idx_univ = latptr % lat(i_xyz(1),i_xyz(2),i_xyz(3))
                         
                        !print *, 'FILL_LATTICE : ', latptr%lat_id, xyz(:)
                    end associate
                    
                    ! Move particle to next level and search for the lower cells.                    
                    call find_cell_xyz(xyz, idx_univ, cell_idx)
                    
                    
                end if CELL_TYPE
            end associate
            
        else 
            print *, 'cell not found from find_cell_xyz'
			print *, xyz
            stop
        end if        
        
        
    end subroutine     
    
    
    function getXYZ(index, a, b, c) result(xyz)
        integer, intent(in) :: index, a,b,c
        integer, dimension(3) :: xyz

        if ((index.le.0).or.(index.gt.(a*b*c))) then
            print *, 'function getXYZ() :: INDEX OUT OF RANGE'
            stop
        endif
        
        xyz(3) = ceiling(real(index, 8)/real(a*b,8))
        xyz(2) = ceiling(real(index - (xyz(3)-1)*a*b, 8)/a)
        xyz(1) = index - (xyz(3)-1)*a*b - (xyz(2)-1)*a
    
    end function getXYZ
    
end module
