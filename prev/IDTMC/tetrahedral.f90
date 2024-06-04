module tetrahedral 
    use omp_lib
    use variables
    use constants
    use particle_header,    only: particle
    use randoms,            only: rang
    use physics,            only: collision_MG
    use XS_header 
    use tally,              only: TallyFlux, FindTallyBin, TallyPower
    use ace_xs,             only: getMacroXS
    use material_header,    only: materials
    use ace_reactions,      only: collision_CE
    use depletion_module,   only: tally_burnup
    use geometry_header,    only: cells, surfaces
    use geometry,           only: cell_distance, find_cell
    
    implicit none

    type :: type_node 
        real(8) :: xyz(1:3)
    end type 
    type(type_node), allocatable :: node(:) 
    
    type :: type_tet 
        integer :: node(4)
        integer :: neighbor(4)
        real(8) :: abcd(4,4) 
        integer :: bc(4)
        real(8) :: vol
        integer :: mat
        integer :: part 
        real(8) :: temperature ! in Kelvin [K]
    end type
    type(type_tet), allocatable :: tet(:) 
    
    type type_tri
        integer :: node(3) 
        integer :: tet 
        integer :: idx
    endtype
    type(type_tri), allocatable :: tri_outside(:)
    
    
    
    integer :: num_node, num_tet, num_outsidetri
    integer :: tet_bc
    integer, parameter :: rd_tet = 99
    real(8) :: tet_xyz(3) 
    
contains

    
    
    subroutine transport_tet(p)
        type(Particle), intent(inout) :: p
        
        integer :: i 
        integer :: j                      ! coordinate level
        real(8) :: d_boundary             ! distance to nearest boundary
        real(8) :: d_collision            ! distance to collision
        real(8) :: distance               ! distance particle travels
        logical :: found_cell             ! found cell which particle is in?
        real(8) :: macro_xs(5)
        real(8) :: xyz(3), uvw(3)
        integer :: i_cell
        integer :: next_tet
        integer :: bc
        real(8) :: d_s, val 
        integer :: idx_surf
        
        ! Find distance to boundary
        call distance_tet(p, d_boundary, next_tet, bc)
        !p%last_material = p%material
        p%material = tet(p%coord(1)%cell)%mat
        
        
        val = 1.0d0
        ! Sample a distance to collision
        if (E_mode == 0) then 
            d_collision = -log(rang())/(sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
            macro_xs(1) = (sum(XS_MG(p%material)%sig_scat(p%g,:)) + XS_MG(p%material)%sig_abs(p%g))
            macro_xs(2) = XS_MG(p%material)%sig_abs(p%g)
            macro_xs(3) = XS_MG(p%material)%sig_fis(p%g)
            macro_xs(4) = XS_MG(p%material)%sig_fis(p%g)*XS_MG(p%material)%nu(p%g)
            
        elseif (E_mode == 1) then
            !if (p%last_material /= p%material)  p%macro_xs = getMacroXS(materials(p%material), p%E)
            !p%macro_xs = getMacroXS(materials(p%material), p%E)
            !macro_xs = p%macro_xs
            macro_xs = getMacroXS(materials(p%material), p%E, p%kT)
            
            ! Leakage free collision distance estimator 
            if (bc == 1) then 
                val = 1.0d0    - exp(-macro_xs(1)*d_boundary)
                d_collision = -log(1.0d0-rang()*val)/macro_xs(1)
                d_boundary = INFINITY 
                p%wgt = p%wgt * val
            else 
                d_collision = -log(rang())/macro_xs(1)
            endif 
            
        endif 
        
        distance = min(d_boundary, d_collision)
        
        !> Track-length estimator
        !$omp atomic 
        k_tl = k_tl + distance*p%wgt*macro_xs(4)
        
        !!> Tally ===========================================================================
        !if (tally_switch > 0 .and. curr_cyc > n_inact) then 
        !    i_cell = p%coord(1)%cell
        !    !$omp atomic
        !    TallyFlux(i_cell) = TallyFlux(i_cell) + p%wgt*distance; 
        !    !$omp atomic
        !    TallyPower(i_cell) = TallyPower(i_cell) + p%wgt*distance*macro_xs(5); 
        !endif 

        !if (bc == 1) p%wgt = p%wgt * val
        
        !> Advance particle =====================================================================
        p % coord(1) % xyz = p % coord(1) % xyz + distance * p % coord(1) % uvw

        !> Collision or Surface-cross ===========================================================
        if (distance == d_collision) then ! collision     
            !print *, 'col'
            
            !> Tally ===========================================================================
            if (tally_switch > 0 .and. curr_cyc > n_inact) then 
                i_cell = p%coord(1)%cell
                !$omp atomic
                TallyFlux(i_cell) = TallyFlux(i_cell) + p%wgt; 
                !$omp atomic
                TallyPower(i_cell) = TallyPower(i_cell) + p%wgt*macro_xs(5); 
            endif 
            if (E_mode == 0) then 
                call collision_MG(p)
            else !(E_mode == 1) 
                call collision_CE(p)
            endif
            !p%coord(1)%universe = 0 ! initialize next tetrahedron's face index
            p%tet_face = 0 ! initialize next tetrahedron's face index
            
        else
            p%n_cross = p%n_cross + 1 
            ! Reassign particle location 
            if (bc == 2) then  ! reflective bc
                !print *, 'ref'
                p % coord(1) % xyz = p % coord(1) % xyz - TINY_BIT * p % coord(1) % uvw
                call reflective_bc_tet(p, next_tet)
                if (.not. in_the_tet(p % coord(1) % xyz, p%coord(1)%cell)) &
                    p%coord(1)%cell = find_tet(p%coord(1)%xyz)
                
                p%tet_face = 0 ! initialize next tetrahedron's face index
                
                if (p%coord(1)%cell < 0) then 
                    print *, "KILLED B/C THE PARTICLE ESCAPE :: ", icore
                    print *, p%coord(1)%xyz(:) 
                    p%alive = .false. 
                endif 
                
                
            elseif (next_tet > 0) then 
                !print *, 'cross'
                p%coord(1)%cell = next_tet
                p%material = tet(p%coord(1)%cell)%mat
            elseif (bc == 1) then 
                p%alive = .false. 
            endif
            if (p%n_cross > 1e6) then 
                print *, "KILLED B/C TOO MANY CROSSING :: ",p%n_cross, icore 
                p%alive = .false. 
            endif 
        endif
        
        
    end subroutine transport_tet    

    
    subroutine read_msh() 
        integer :: i, j, i_tet
        real(8) :: xyz(3), T_tet
        real(8), allocatable :: T_node(:)

        open(rd_tet, file="./inputfile/msh.inp",action="read", status="old")
        if(icore==score) print *, "   Start reading pre-processed mesh data..."
        read(rd_tet,*) num_node 
        allocate(node(1:num_node))
        
        do i = 1, num_node 
            read(rd_tet,*) node(i)%xyz(1), node(i)%xyz(2), node(i)%xyz(3)
        enddo 
        
        read(rd_tet,*) num_tet
        allocate(tet(1:num_tet))
        do i = 1, num_tet 
        read(rd_tet,*) tet(i)%node(1), tet(i)%node(2), tet(i)%node(3), tet(i)%node(4), &
                        tet(i)%neighbor(1), tet(i)%neighbor(2), tet(i)%neighbor(3), tet(i)%neighbor(4), &
                        tet(i)%bc(1), tet(i)%bc(2), tet(i)%bc(3), tet(i)%bc(4), tet(i)%mat, tet(i)%part
                        
            do j = 1, 4
                if (tet(i)%bc(j) == 0) tet(i)%bc(j) = tet_bc
            enddo 
        enddo 
        
        
        read(rd_tet,*) num_outsidetri
        allocate(tri_outside(1:num_outsidetri))
        do i = 1, num_outsidetri 
            read(rd_tet,*) tri_outside(i)%node(1), tri_outside(i)%node(2), tri_outside(i)%node(3) &
                            , tri_outside(i)%tet, tri_outside(i)%idx
        enddo 
        
        close(rd_tet)
        
        
        ! Fix abcd and volume 
        do i_tet = 1, num_tet 
            Tet(i_tet)%vol = volume_tet(node(tet(i_tet)%node(1))%xyz, node(tet(i_tet)%node(2))%xyz, &
                                        node(tet(i_tet)%node(3))%xyz, node(tet(i_tet)%node(4))%xyz)
            
            !call get_bcd (node(tet(i_tet)%node(1))%xyz, node(tet(i_tet)%node(2))%xyz, &
            !              node(tet(i_tet)%node(3))%xyz, node(tet(i_tet)%node(4))%xyz, &
            !              Tet(i_tet)%abcd)
            
            Tet(i_tet)%vol = abs(Tet(i_tet)%vol)
        enddo 
        
        
        
        
        
        
        
        
        !===================== Setting tet temperature ==============================
        !do i = 1, num_tet 
        !    if (Tet(i)%mat == 1 .or. Tet(i)%mat == 2 .or. Tet(i)%mat == 4) then 
        !        Tet(i)%temperature = 600
        !    else 
        !        Tet(i)%temperature = 600
        !    endif 
        !enddo 
        !do i = 1, num_tet 
        !    Tet(i)%temperature = 600
        !enddo 

        
        allocate(T_node(num_node)) 
        open(rd_tet, file="./Temperature.out",action="read", status="old")
        if(icore==score) print *, "   Start reading Temperature data..."
        do i = 1, num_node 
            read(rd_tet,*) xyz(1:3), T_node(i)
            if (isnan(T_node(i)) .and. icore==score) then 
                print *, "NaN detected in Temperature.out line ", i 
                T_node(i) = 800d0   ! fix 필요 
            endif 
        enddo 
        close(rd_tet) 
        
        
        do i_tet = 1, num_tet 
            T_tet = 0 
            do i = 1, 4
                j = Tet(i_tet)%node(i)
                T_tet = T_tet + T_node(j)
            enddo 
            Tet(i_tet)%temperature = (T_tet / 4.0) + 273.0
        enddo 
        
        deallocate(T_node) 
        
        
        !========================================================================



        end subroutine
    
    
    
    
    subroutine distance_tet(p, d_boundary, next_tet, bc) 
        type(particle), intent(inout) :: p
        real(8), intent(inout) :: d_boundary
        integer, intent(inout) :: next_tet
        integer, intent(inout) :: bc
        integer :: i , j 
        real(8) :: a,b,c,d, a1,b1,c1,d1, a2,b2,c2,d2
        real(8) :: xyz(3), xyz1(3), xyz2(3), xyz3(3), uvw(3) 
        real(8) :: dist_temp
        integer :: idx(3), temp
        integer :: idx_tet
        integer :: vertex_idx
        real(8) :: test_dist(4)
        logical :: found_cell
        real(8) :: dr
        
        
        vertex_idx = -1
        test_dist(:) = -1 
        
        d_boundary = INFINITY
        idx_tet = p%tet
        
        
        xyz = p%coord(p%n_coord)%xyz
        uvw = p%coord(p%n_coord)%uvw
        
        do i = 1, 4 
            if (i == p%tet_face) cycle
            temp = 0
            loop: do j = 1, 4 
                if (i==j) cycle loop
                temp = temp+1
                idx(temp) = j
            enddo loop
            
            !print *, 'test', Tet(idx_tet)%node(idx(1)), idx_tet
            xyz1 = Node(Tet(idx_tet)%node(idx(1)))%xyz 
            xyz2 = Node(Tet(idx_tet)%node(idx(2)))%xyz 
            xyz3 = Node(Tet(idx_tet)%node(idx(3)))%xyz 
            
            a1 = xyz2(1) - xyz1(1)
            b1 = xyz2(2) - xyz1(2)
            c1 = xyz2(3) - xyz1(3)
            a2 = xyz3(1) - xyz1(1)
            b2 = xyz3(2) - xyz1(2)
            c2 = xyz3(3) - xyz1(3)
            a  = b1 * c2 - b2 * c1
            b  = a2 * c1 - a1 * c2
            c  = a1 * b2 - b1 * a2
            d  = a * xyz1(1) + b * xyz1(2) + c * xyz1(3)
            dist_temp = (d-a*xyz(1)-b*xyz(2)-c*xyz(3))/(a*uvw(1)+b*uvw(2)+c*uvw(3))
            
            !call RayIntersectsTriangle(xyz, uvw, xyz1, xyz2, xyz3, dist_temp) 
            test_dist(i) = dist_temp
            
            if (dist_temp >= 0 .and. d_boundary > dist_temp) then
                d_boundary = dist_temp
                vertex_idx = i
            endif
            
        enddo 

        if (vertex_idx < 0) then 
             
            !print *, in_the_tet(p%coord(p%n_coord)%xyz,  p%tet), p%tet
            
            
            do j = 1, p%n_coord
                p%coord(j)%xyz = p%coord(j)%xyz + TINY_BIT*p%coord(j)%uvw
            enddo 
            p%tet = find_tet(p%coord(p%n_coord)%xyz)
            p%tet_face = 0
            p%in_tet = .true.
            if (p%tet .le. 0) p%in_tet = .false.
            found_cell = .false.
            call find_cell(p, found_cell)
            next_tet = p%tet
            return 
        endif 
        
        
        next_tet = tet(idx_tet)%neighbor(vertex_idx)
        
                 
        
        !p%tet_face = -1
        !do i = 1, 4
        !    if (tet(next_tet)%neighbor(i) == idx_tet) then 
        !        p%tet_face = i
        !        exit
        !    endif
        !enddo
        !bc = 0 
        
    end subroutine 
    
    
    
    
    
    
    
    function in_the_tet(xyz, idx_tet) result (inside)
        real(8), intent(in) :: xyz(3)
        integer, intent(in) :: idx_tet
        real(8) :: vol0, vol(4)
        integer :: a,b,c,d
        integer :: i
        logical :: inside
        
        inside = .true. 
        a = tet(idx_tet)%node(1); b = tet(idx_tet)%node(2)
        c = tet(idx_tet)%node(3); d = tet(idx_tet)%node(4)
        
        vol0   = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
        vol(1) = volume_tet (xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
        vol(2) = volume_tet (Node(a)%xyz, xyz, Node(c)%xyz, Node(d)%xyz);
        vol(3) = volume_tet (Node(a)%xyz, Node(b)%xyz, xyz, Node(d)%xyz);
        vol(4) = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, xyz);
        
        do i = 1, 4 
            if (vol(i)/vol0 < -1.0d-10) then 
                inside = .false.
                !print *, i, vol(i)/vol0
                !print '(a,5e14.3)', 'inthetet ',vol(:), vol0 
            endif 
        enddo 
        
    end function 
    
    
    integer function find_tet(xyz) 
        real(8), intent(in) :: xyz(3)
        integer :: a,b,c,d
        real(8) :: vol0, vol(4) 
        real(8) :: dist, dist_temp
        integer :: i, j, k
        integer :: idx
        logical :: exist
        
        ! 1. 가장 가까운 점 찾기 
        ! 2. Tet(i)%node에 이게 없으면 cycle
        ! 2.1 여기서 못찾으면 전체에서 찾기
        
        idx = 0 
        dist = huge(0.0_8)
        do_node: do i = 1, num_node 
            dist_temp = sqrt((xyz(1)-Node(i)%xyz(1))**2+(xyz(2)-Node(i)%xyz(2))**2+(xyz(3)-Node(i)%xyz(3))**2)
            if (dist_temp < dist) then 
                idx = i
                dist = dist_temp 
            endif
        enddo do_node
        
        
        do_tet1: do i = 1, num_tet 
            exist = .false. 
            do j = 1, 4 
                if (Tet(i)%node(j) == idx) then 
                    exist = .true.
                    exit
                endif 
            enddo 
            if (exist) then 
                a = tet(i)%node(1); b = tet(i)%node(2)
                c = tet(i)%node(3); d = tet(i)%node(4)
                
                vol0   = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
                vol(1) = volume_tet (xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
                vol(2) = volume_tet (Node(a)%xyz, xyz, Node(c)%xyz, Node(d)%xyz);
                vol(3) = volume_tet (Node(a)%xyz, Node(b)%xyz, xyz, Node(d)%xyz);
                vol(4) = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, xyz);
                do j = 1, 4 
                    if (vol(j)/vol0 < 0) cycle do_tet1
                enddo 
                find_tet = i 
                return ! 대부분 여기서 return 
            endif 
        enddo do_tet1
        
        ! 못 찾은 경우 전체에서 찾기 
        do_tet2: do i = 1, num_tet 
            a = tet(i)%node(1); b = tet(i)%node(2)
            c = tet(i)%node(3); d = tet(i)%node(4)
            vol0   = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
            vol(1) = volume_tet (xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
            vol(2) = volume_tet (Node(a)%xyz, xyz, Node(c)%xyz, Node(d)%xyz);
            vol(3) = volume_tet (Node(a)%xyz, Node(b)%xyz, xyz, Node(d)%xyz);
            vol(4) = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, xyz);
            do j = 1, 4 
                if (vol(j)/vol0 < 0) cycle do_tet2
            enddo 
            find_tet = i 
            return 
        enddo do_tet2
        
        find_tet = -1 
        return        
    end function
    
    
    integer function find_tet_old(xyz) 
        real(8), intent(in) :: xyz(3)
        
        real(8) :: vol0, vol(4) 
        integer :: a,b,c,d
        integer :: i, j
        
        do_tet: do i = 1, num_tet 
            a = tet(i)%node(1); b = tet(i)%node(2)
            c = tet(i)%node(3); d = tet(i)%node(4)
            
            vol0   = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
            vol(1) = volume_tet (xyz, Node(b)%xyz, Node(c)%xyz, Node(d)%xyz);
            vol(2) = volume_tet (Node(a)%xyz, xyz, Node(c)%xyz, Node(d)%xyz);
            vol(3) = volume_tet (Node(a)%xyz, Node(b)%xyz, xyz, Node(d)%xyz);
            vol(4) = volume_tet (Node(a)%xyz, Node(b)%xyz, Node(c)%xyz, xyz);
            
            do j = 1, 4 
                if (vol(j)/vol0 < 0) cycle do_tet
            enddo 
            
            find_tet_old = i 
            return 
            
        enddo do_tet
        
        find_tet_old = -1 
        return
        
    end function
    
    real(8) function volume_tet (xyz1, xyz2, xyz3, xyz4 ) 
        real(8), dimension(3), intent(in) :: xyz1, xyz2, xyz3, xyz4 
        real(8) :: mat(4,4) 
        
        integer :: i
        
        do i = 1,3
            mat(i,1) = xyz1(i)
            mat(i,2) = xyz2(i)
            mat(i,3) = xyz3(i)
            mat(i,4) = xyz4(i)
        enddo 
        mat(4,:) = 1
        
        volume_tet = (1.0/6.0)*det(mat, 4)
    end function 
    
    real(8) function area_tri (xyz1, xyz2, xyz3)
        real(8), dimension(3), intent(in) :: xyz1, xyz2, xyz3
        real(8) :: a, b, c, s
        a = sqrt((xyz1(1)-xyz2(1))**2 + (xyz1(2)-xyz2(2))**2 + (xyz1(3)-xyz2(3))**2)
        b = sqrt((xyz2(1)-xyz3(1))**2 + (xyz2(2)-xyz3(2))**2 + (xyz2(3)-xyz3(3))**2)
        c = sqrt((xyz3(1)-xyz1(1))**2 + (xyz3(2)-xyz1(2))**2 + (xyz3(3)-xyz1(3))**2)
        s = (a + b + c) / 2.0
        area_tri = sqrt(s * (s - a) * (s - b) * (s - c))
        return 
    end function 



    real(8) function det(aa, n)
        real(8) :: aa(4,4)
        real(8) :: tmp,c(4,4)
        real(8) :: max
        integer i,j,k,l,m,n,num(4)
        
        det=1.    
        do k=1,n
            max=aa(k,k);num(k)=k;
            do i=k+1,n 
                if(abs(max)<abs(aa(i,k))) then
                    max=aa(i,k)
                    num(k)=i
                endif
            enddo
            if (num(k)/=k) then
                do l=k,n 
                    tmp=aa(k,l)
                    aa(k,l)=aa(num(k),l)
                    aa(num(k),l)=tmp
                enddo
                det=-1.*det
            endif
            do m=k+1,n
                c(m,k)=aa(m,k)/aa(k,k)
                do l=k,n 
                    aa(m,l)=aa(m,l)-c(m,k)*aa(k,l)
                enddo
            enddo !There we made matrix triangular!    
        enddo

        do i=1,n
        det=det*aa(i,i)
        enddo
        return
    end function    



    subroutine get_bcd(xyz1, xyz2, xyz3, xyz4, abcd)
        real(8), intent(in) :: xyz1(3), xyz2(3), xyz3(3), xyz4(3)
        real(8), intent(inout) :: abcd(4,4) 
        integer :: i 
        
        do i = 1, 4
            abcd(1,i) = 0 
        enddo 
        
        abcd(2,1) = xyz4(2)*(xyz3(3)-xyz2(3))+xyz3(2)*(xyz2(3)-xyz4(3))+xyz2(2)*(xyz4(3)-xyz3(3))
        abcd(2,2) = xyz4(2)*(xyz1(3)-xyz3(3))+xyz1(2)*(xyz3(3)-xyz4(3))+xyz3(2)*(xyz4(3)-xyz1(3))
        abcd(2,3) = xyz4(2)*(xyz2(3)-xyz1(3))+xyz2(2)*(xyz1(3)-xyz4(3))+xyz1(2)*(xyz4(3)-xyz2(3))
        abcd(2,4) = xyz3(2)*(xyz1(3)-xyz2(3))+xyz1(2)*(xyz2(3)-xyz3(3))+xyz2(2)*(xyz3(3)-xyz1(3))
        
        abcd(3,1) = xyz4(1)*(xyz2(3)-xyz3(3))+xyz2(1)*(xyz3(3)-xyz4(3))+xyz3(1)*(xyz4(3)-xyz2(3))
        abcd(3,2) = xyz4(1)*(xyz3(3)-xyz1(3))+xyz3(1)*(xyz1(3)-xyz4(3))+xyz1(1)*(xyz4(3)-xyz3(3))
        abcd(3,3) = xyz4(1)*(xyz1(3)-xyz2(3))+xyz1(1)*(xyz2(3)-xyz4(3))+xyz2(1)*(xyz4(3)-xyz1(3))
        abcd(3,4) = xyz3(1)*(xyz2(3)-xyz1(3))+xyz2(1)*(xyz1(3)-xyz3(3))+xyz1(1)*(xyz3(3)-xyz2(3))
        
        abcd(4,1) = xyz4(1)*(xyz3(2)-xyz2(2))+xyz3(1)*(xyz2(2)-xyz4(2))+xyz2(1)*(xyz4(2)-xyz3(2))
        abcd(4,2) = xyz4(1)*(xyz1(2)-xyz3(2))+xyz1(1)*(xyz3(2)-xyz4(2))+xyz3(1)*(xyz4(2)-xyz1(2))
        abcd(4,3) = xyz4(1)*(xyz2(2)-xyz1(2))+xyz2(1)*(xyz1(2)-xyz4(2))+xyz1(1)*(xyz4(2)-xyz2(2))
        abcd(4,4) = xyz3(1)*(xyz1(2)-xyz2(2))+xyz1(1)*(xyz2(2)-xyz3(2))+xyz2(1)*(xyz3(2)-xyz1(2))
        
    end subroutine 
    
    
    
    
    
    subroutine reflective_bc_tet(p, idx_surf)  
        type(particle), intent(inout) :: p 
        integer, intent(in) :: idx_surf
        integer :: idx_tet, i, idx(3), temp
        real(8) :: a,b,c, a1,b1,c1,d1, a2,b2,c2,d2
        real(8) :: xyz1(3), xyz2(3), xyz3(3), uvw(3) 
        real(8) :: val, norm 
    
        !idx_tet = p%coord(1)%cell 
        idx_tet = p%tet 
        
        temp = 0
        loop: do i = 1, 4 
            if (idx_surf==i) cycle loop
            temp = temp+1
            idx(temp) = i
        enddo loop
            
        xyz1 = Node(Tet(idx_tet)%node(idx(1)))%xyz 
        xyz2 = Node(Tet(idx_tet)%node(idx(2)))%xyz 
        xyz3 = Node(Tet(idx_tet)%node(idx(3)))%xyz 
        
        a1 = xyz2(1) - xyz1(1)
        b1 = xyz2(2) - xyz1(2)
        c1 = xyz2(3) - xyz1(3)
        a2 = xyz3(1) - xyz1(1)
        b2 = xyz3(2) - xyz1(2)
        c2 = xyz3(3) - xyz1(3)
        a  = b1 * c2 - b2 * c1
        b  = a2 * c1 - a1 * c2
        c  = a1 * b2 - b1 * a2
    
        uvw = p%coord(1)%uvw 
        val = 2.0*(a*uvw(1)+b*uvw(2)+c*uvw(3))/(a**2 + b**2 + c**2)
        
        p%coord(1)%uvw(1) = uvw(1) - val*a 
        p%coord(1)%uvw(2) = uvw(2) - val*b 
        p%coord(1)%uvw(3) = uvw(3) - val*c 
        
        norm = sqrt(p%coord(1)%uvw(1)**2 + p%coord(1)%uvw(2)**2 + p%coord(1)%uvw(3)**2)
        
        p%coord(1)%uvw = p%coord(1)%uvw / norm
            
    end subroutine 

    subroutine distance_tet_from_outside (p, d_boundary, next_tet, tet_face) 
        type(particle), intent(inout) :: p
        real(8), intent(inout) :: d_boundary
        integer, intent(inout) :: next_tet
        integer, intent(inout) :: tet_face
        integer :: i_tet, i, j
        integer :: face_out(4) 
        real(8) :: xyz(3), uvw(3)
        real(8) :: a,b,c,d, a1,b1,c1,d1, a2,b2,c2,d2
        real(8) :: xyz1(3), xyz2(3), xyz3(3) 
        real(8) :: dist_temp
        integer :: idx(3), temp
        real(8) :: xyz_new(3)
        real(8) :: area0, area(3), val 
        real(8) :: area_save, area_temp
        integer :: n_found
        
        xyz = p%coord(p%n_coord)%xyz
        uvw = p%coord(p%n_coord)%uvw
        
        d_boundary = INFINITY
        tet_face = -1
        
        area_save = INFINITY
        
        !n_found = 0 
        !do_tet: do i_tet = 1, num_tet 
        !    face_out(:) = 0 
        !    do i = 1, 4 
        !        if (Tet(i_tet)%neighbor(i) .le. 0) face_out(i) = 1 ! no neighboring tet
        !    enddo 
        !    if (sum(face_out(:)) == 0) cycle do_tet
        !    
        !    do_face: do i = 1, 4 
        !        if (face_out(i) == 0) cycle do_face
        !        
        !        ! 1. node index 3개 구하기 
        !        temp = 0
        !        loop: do j = 1, 4
        !            if (i==j) cycle loop
        !            temp = temp+1
        !            idx(temp) = j
        !        enddo loop
        !        
        !        ! 2. 3개의 xyz를 통해 plane equation 구하기 -> 거리 계산
        !        xyz1 = Node(Tet(i_tet)%node(idx(1)))%xyz 
        !        xyz2 = Node(Tet(i_tet)%node(idx(2)))%xyz 
        !        xyz3 = Node(Tet(i_tet)%node(idx(3)))%xyz 
        !        
        !        a1 = xyz2(1) - xyz1(1)
        !        b1 = xyz2(2) - xyz1(2)
        !        c1 = xyz2(3) - xyz1(3)
        !        a2 = xyz3(1) - xyz1(1)
        !        b2 = xyz3(2) - xyz1(2)
        !        c2 = xyz3(3) - xyz1(3)
        !        a  = b1 * c2 - b2 * c1
        !        b  = a2 * c1 - a1 * c2
        !        c  = a1 * b2 - b1 * a2
        !        d  = a * xyz1(1) + b * xyz1(2) + c * xyz1(3)
        !        dist_temp = (d-a*xyz(1)-b*xyz(2)-c*xyz(3))/(a*uvw(1)+b*uvw(2)+c*uvw(3))                
        !        !call RayIntersectsTriangle(xyz, uvw, xyz1, xyz2, xyz3, dist_temp)
        !
        !        
        !        ! 3. plane과 ray의 intersection 구하기 
        !        xyz_new(:) = xyz(:) + dist_temp * uvw(:) 
        !        
        !        ! 4. intersection이 삼각형 내부에 있는지 판단
        !        ! 5. 있다면 dist 저장 
        !        if (inside_tri(xyz1,xyz2,xyz3,xyz_new) .and. dist_temp > 0 .and. d_boundary > dist_temp) then 
        !            !n_found = n_found + 1
        !            !area_save = area_temp
        !            d_boundary = dist_temp! + TINY_BIT
        !            next_tet = i_tet
        !            tet_face = i
        !            
        !            !if (n_found == 2) then 
        !            !    !print *, 'exit'
        !            !    return
        !            !endif 
        !        endif
        !        
        !    enddo do_face
        !    
        !enddo do_tet
        
        
        do_face: do i = 1, num_outsidetri
            
            i_tet = tri_outside(i)%tet
            
            
            ! 2. 3개의 xyz를 통해 plane equation 구하기 -> 거리 계산
            xyz1 = Node(tri_outside(i)%node(1))%xyz 
            xyz2 = Node(tri_outside(i)%node(2))%xyz 
            xyz3 = Node(tri_outside(i)%node(3))%xyz 
            
            a1 = xyz2(1) - xyz1(1)
            b1 = xyz2(2) - xyz1(2)
            c1 = xyz2(3) - xyz1(3)
            a2 = xyz3(1) - xyz1(1)
            b2 = xyz3(2) - xyz1(2)
            c2 = xyz3(3) - xyz1(3)
            a  = b1 * c2 - b2 * c1
            b  = a2 * c1 - a1 * c2
            c  = a1 * b2 - b1 * a2
            d  = a * xyz1(1) + b * xyz1(2) + c * xyz1(3)
            dist_temp = (d-a*xyz(1)-b*xyz(2)-c*xyz(3))/(a*uvw(1)+b*uvw(2)+c*uvw(3))                
            !call RayIntersectsTriangle(xyz, uvw, xyz1, xyz2, xyz3, dist_temp)
        
            
            ! 3. plane과 ray의 intersection 구하기 
            xyz_new(:) = xyz(:) + dist_temp * uvw(:) 
            
            ! 4. intersection이 삼각형 내부에 있는지 판단
            ! 5. 있다면 dist 저장 
            if (inside_tri(xyz1,xyz2,xyz3,xyz_new) .and. dist_temp > 0 .and. d_boundary > dist_temp) then 
                !area_save = area_temp
                d_boundary = dist_temp! + TINY_BIT
                next_tet = i_tet
                tet_face = tri_outside(i)%idx
            endif
            
        enddo do_face
            
    end subroutine
    
    
    
    function inside_tri (xyz1,xyz2,xyz3,xyzp)
        real(8),dimension(3), intent(in) :: xyz1,xyz2,xyz3,xyzp
        logical :: inside_tri 
        real(8),dimension(3) :: v0, v1, v2
        real(8) :: dot00,dot01,dot02,dot11,dot12, invDenom, u, v, eps 
        
        eps = 1.0d-15
        
        ! Compute vectors
        v0 = xyz3 - xyz1
        v1 = xyz2 - xyz1
        v2 = xyzp - xyz1

        ! Compute dot products
        dot00 = dot_product(v0, v0)
        dot01 = dot_product(v0, v1)
        dot02 = dot_product(v0, v2)
        dot11 = dot_product(v1, v1)
        dot12 = dot_product(v1, v2)

        ! Compute barycentric coordinates
        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom

        ! Check if point is in triangle
        if ((u >= -eps) .and. (v >= -eps) .and. (u + v < 1.0d0+eps)) then 
            inside_tri = .true. 
        else 
            inside_tri = .false. 
        endif 
        
    end function
    
    
    
    subroutine RayIntersectsTriangle(xyz0, uvw, xyz1, xyz2, xyz3, distance) 
        real(8), dimension(3), intent(in) :: xyz0, uvw, xyz1, xyz2, xyz3
        real(8), intent(inout) :: distance
        real(8), dimension(3) :: edge1, edge2, h, s, q 
        real(8) :: a, f, u, v
        
        edge1 = xyz2 - xyz1
        edge2 = xyz3 - xyz1 
        
        h = cross_product(uvw, edge2) 
        a = dot_product(edge1, h) 
        
        if (abs(a) < 1.0d-10) then 
            distance = -INFINITY
            return 
        endif 
        
        f = 1.0d0 / a 
        s = xyz0 - xyz1
        u = f * dot_product(s,h) 
        if (u < 0.0d0 .or. u > 1.0d0) then 
            distance = -INFINITY
            return 
        endif 
        
        q = cross_product(s, edge1) 
        v = f * dot_product(uvw, q) 
        if (v < 0.0d0 .or. (u+v) > 1.0d0) then 
            distance = -INFINITY
            return 
        endif 
        distance = f * dot_product(edge2, q)
        
    end subroutine
    
    function cross_product(a,b) result(axb)
        implicit none
        integer,parameter :: wp=selected_real_kind(15, 307) !double precision
        real(wp),dimension(3) :: axb
        real(wp),dimension(3),intent(in) :: a
        real(wp),dimension(3),intent(in) :: b 

        axb(1) = a(2)*b(3) - a(3)*b(2)
        axb(2) = a(3)*b(1) - a(1)*b(3)
        axb(3) = a(1)*b(2) - a(2)*b(1)
    end function cross_product
    
    
    
    
    
    

end module
